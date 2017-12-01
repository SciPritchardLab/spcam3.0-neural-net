#include <misc.h>
#include <params.h>
! MSP on Kraken Cray XT5. For mysterious reasons, history files do not write correctly
! Unless this is defined:
#define HDEBUG 

module history
!----------------------------------------------------------------------- 
! 
! Purpose: History module.  Contains data and functions for writing history files.
!
! Public functions/subroutines:
!   addfld, add_default
!   intht
!   write_restart_history
!   read_restart_history
!   outfld
!   wshist
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
! $Id: history.F90,v 1.26.2.38 2003/12/15 18:52:35 hender Exp $
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use ppgrid,    only: pcols
   use constituents, only: pcnst, pnats, cnst_name, cnst_longname, &
                           dcconnam, sflxnam, tendnam, qphystendnam
   use filenames, only: mss_wpass, mss_irt, interpret_filename_spec, get_archivedir
   use cloudsimulatorparms, only: ntau, npres, prlim, taulim, doisccp
#if ( defined STAGGERED )
   use pmgrid,    only: masterproc, beglat, endlat, plat, plon, plev, plevp, dyngrid_set, splon, beglev, endlev, endlevp
#else
   use pmgrid,    only: masterproc, beglat, endlat, plat, plon, plev, plevp, dyngrid_set
#endif

   implicit none

PRIVATE

   include 'netcdf.inc'

   integer, parameter :: pflds = 1000           ! max number of fields 
   integer, parameter :: ptapes = 6             ! max number of tapes
   integer, parameter :: max_chars = 128        ! max chars for char variables

   real(r8), parameter :: fillvalue = 1.e36     ! fill value for reduced grid and others

   type field_info
      character*8 :: name                       ! field name
      character*(max_chars) :: long_name        ! long name
      character*(max_chars) :: units            ! units

      logical :: flag_xyfill                    ! non-applicable xy points flagged with fillvalue
      logical :: flag_isccplev                  ! levels dimension is ISCCP not CAM

      integer :: coldimin                       ! column dimension of model array
      integer :: numlev                         ! vertical dimension (.nc file and internal arr)
      integer :: begver                         ! on-node vert start index
      integer :: endver                         ! on-node vert end index
      integer :: begdim3                        ! on-node chunk or lat start index
      integer :: enddim3                        ! on-node chunk or lat end index
      integer :: decomp_type                    ! type of decomposition (physics or dynamics)
      integer, pointer :: colperdim3(:)         ! number of valid elements per chunk or lat
   end type field_info
!
! master_entry: elements of an entry in the master field list
!
   type master_entry
      type (field_info)     :: field            ! field information
      character*1           :: avgflag(ptapes)  ! averaging flag
      character*(max_chars) :: time_op(ptapes)  ! time operator (e.g. max, min, avg)
      logical               :: actflag(ptapes)  ! active/inactive flag
   end type master_entry

   type (master_entry) :: masterlist(pflds)     ! master field list
!
! hbuffer_2d, hbuffer_3d: 2-D and 3-D history buffer pointers.
!     Select either r4 or r8 kind buffer depending on hbuf_prec.
!
   type hbuffer_2d
      real(r8), pointer :: buf8(:,:)            ! 2-D history buffer for r8
      real(r4), pointer :: buf4(:,:)            ! 2-D history buffer for r4
   end type hbuffer_2d

   type hbuffer_3d
      real(r8), pointer :: buf8(:,:,:)          ! 3-D history buffer for r8
      real(r4), pointer :: buf4(:,:,:)          ! 3-D history buffer for r4
   end type hbuffer_3d
!
! arrays served as targets for history pointers
!
   integer,  target :: nothing_int(1,1)         ! 2-D integer target
   real(r8), target :: nothing_r8(1,1,1)        ! 3-D r8 target
   real(r4), target :: nothing_r4(1,1,1)        ! 3-D r4 target

   type column_info
      character*(max_chars) :: lat_name ! latitude name for this column or columns
      character*(max_chars) :: lon_name ! latitude name for this column or columns
      integer :: num_lats            ! number of lats in a group of contiguous columns
      integer :: num_lons            ! number of lons in a group of contiguous columns
      integer :: columnlat(2)       ! beginning and ending latitude (range) dimensioned by groups
      integer :: columnlon(2)       ! beginning and ending longitude (range) dimensioned by groups
   end type column_info
!
! hentry: elements of an entry in the list of active fields on a single history file
!
   type hentry
      type (field_info)     :: field            ! field information
      character*1           :: avgflag          ! averaging flag
      character*(max_chars) :: time_op          ! time operator (e.g. max, min, avg)
      character*(max_chars),pointer :: field_column_name(:) ! names of column groups written to tape

      integer :: hbuf_prec                      ! history buffer precision
      integer :: hwrt_prec                      ! history output precision

      type (hbuffer_3d)   :: hbuf               ! history buffer
      integer, pointer :: nacs(:,:)             ! accumulation counter
   end type hentry
!
! active_entry: vehicle for producing a ragged array
!
   type active_entry
      type (hentry) :: hlist(pflds)             ! array of history tape entries
      type (column_info),pointer :: column(:)             ! array of history tape entries
   end type active_entry

   type (active_entry) :: tape(ptapes)          ! history tapes
!
! dim_index_2d, dim_index_3d: 2-D & 3-D dimension index lower & upper bounds
!
   type dim_index_2d                   ! 2-D dimension index
      integer :: beg1, end1            ! lower & upper bounds of 1st dimension
      integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
   end type dim_index_2d

   type dim_index_3d                   ! 3-D dimension index
      integer :: beg1, end1            ! lower & upper bounds of 1st dimension
      integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
      integer :: beg3, end3            ! lower & upper bounds of 3rd dimension
   end type dim_index_3d

   integer :: ndm(12)                           ! number of days in each month (jan-dec)
   save ndm
   data ndm/31,28,31,30,31,30,31,31,30,31,30,31/

   integer :: nfmaster = 0             ! number of fields in master field list
   integer :: nflds(ptapes)            ! number of fields per tape

! per tape sampling frequency (0=monthly avg)

   integer :: i                        ! index for nhtfrq initialization
   integer :: nhtfrq(ptapes) = (/0, (-24, i=2,ptapes)/)  ! history write frequency (0 = monthly)
   integer :: mfilt(ptapes) = 30       ! number of time samples per tape
   integer :: nfils(ptapes)            ! Array of no. of files on current h-file
   integer :: ngroup(ptapes)           ! Array of no. of contiguous columns on current h-file
   integer :: mtapes = 0               ! index of max history file requested 
   integer :: nexcl(ptapes)            ! Actual number of excluded fields
   integer :: nincl(ptapes)            ! Actual number of included primary file fields
   integer :: nhstpr(ptapes) = 8       ! history buffer precision (8 or 4 bytes)
   integer :: ndens(ptapes) = 2        ! packing density (nf_float vs nf_double)
   integer :: ncprec(ptapes) = -999    ! netcdf packing parameter based on ndens
   real(r8) :: beg_time(ptapes)        ! time at beginning of an averaging interval
!
! Netcdf ids
!
   integer :: nfid(ptapes)             ! file id
   integer :: varid(pflds,ptapes)      ! variable ids
   integer :: mdtid(ptapes)            ! var id for timestep
   integer :: ndbaseid(ptapes)         ! var id for base day
   integer :: nsbaseid(ptapes)         ! var id for base seconds of base day
   integer :: nbdateid(ptapes)         ! var id for base date
   integer :: nbsecid(ptapes)          ! var id for base seconds of base date
   integer :: ndcurid(ptapes)          ! var id for current day
   integer :: nscurid(ptapes)          ! var id for current seconds of current day
   integer :: dateid(ptapes)           ! var id for current date
   integer :: datesecid(ptapes)        ! var id for curent seconds of current date
   integer :: nstephid(ptapes)         ! var id for current timestep
   integer :: timeid(ptapes)           ! var id for time
   integer :: tbndid(ptapes)           ! var id for time_bnds
   integer :: gwid(ptapes)             ! var id for gaussian weights
   integer :: date_writtenid(ptapes)   ! var id for date time sample written
   integer :: time_writtenid(ptapes)   ! var id for time time sample written
   integer :: nlonid(ptapes)           ! var id for number of longitudes
   integer :: wnummaxid(ptapes)        ! var id for cutoff fourier wavenumber (reduced grid)
   integer :: nscurf(ptapes)           ! First "current" second of day for each h-file
   integer :: ncsecf(ptapes)           ! First "current" second of date for each h-file

   logical :: rgnht(ptapes) = .false.  ! flag array indicating regeneration volumes
   logical :: hstwr(ptapes) = .false.  ! Flag for history writes
   logical :: empty_htapes  = .false.  ! Namelist flag indicates no default history fields
   logical :: htapes_defined = .false. ! flag indicates history contents have been defined

   integer, parameter :: nlen = 256    ! Length of strings
   character(len=nlen) :: hrestpath(ptapes) = (/(' ',i=1,ptapes)/) ! Full history restart pathnames
   character(len=nlen) :: nfpath(ptapes) = (/(' ',i=1,ptapes)/) ! Array of first pathnames, for header
   character(len=nlen) :: cpath(ptapes)                   ! Array of current pathnames
   character(len=nlen) :: nhfil(ptapes)                   ! Array of current file names
   character(len=1)  :: avgflag_pertape(ptapes) = (/(' ',i=1,ptapes)/) ! per tape averaging flag
   character(len=8)  :: logname             ! user name
   character(len=16) :: host                ! host name
   character(len=nlen) :: ctitle = ' '      ! Case title
   character(len=8)  :: inithist = 'YEARLY' ! If set to '6-HOURLY, 'DAILY', 'MONTHLY' or
                                            ! 'YEARLY' then write IC file
   character(len=10) :: fincl(pflds,ptapes) ! List of fields to add to primary h-file
   character(len=max_chars) :: fincllonlat(pflds,ptapes) ! List of fields to add to primary h-file
   character(len=8)  :: fexcl(pflds,ptapes) ! List of fields to rm from primary h-file
   character(len=10) :: fhstpr(pflds,ptapes) ! List of fields to change default hbuf size
   character(len=10) :: fwrtpr(pflds,ptapes) ! List of fields to change default history output prec
!
! Equivalence to please namelist on a wide variety of platforms
! NOTE: It is *ASSUMED* that ptapes is 6
!
#ifdef INTEL64
#else
   character*10 fincl1(pflds)
   character*10 fincl2(pflds)
   character*10 fincl3(pflds)
   character*10 fincl4(pflds)
   character*10 fincl5(pflds)
   character*10 fincl6(pflds)

   equivalence (fincl1,fincl(1,1))
   equivalence (fincl2,fincl(1,2))
   equivalence (fincl3,fincl(1,3))
   equivalence (fincl4,fincl(1,4))
   equivalence (fincl5,fincl(1,5))
   equivalence (fincl6,fincl(1,6))

   character(len=max_chars) fincl1lonlat(pflds)
   character(len=max_chars) fincl2lonlat(pflds)
   character(len=max_chars) fincl3lonlat(pflds)
   character(len=max_chars) fincl4lonlat(pflds)
   character(len=max_chars) fincl5lonlat(pflds)
   character(len=max_chars) fincl6lonlat(pflds)

   equivalence (fincl1lonlat,fincllonlat(1,1))
   equivalence (fincl2lonlat,fincllonlat(1,2))
   equivalence (fincl3lonlat,fincllonlat(1,3))
   equivalence (fincl4lonlat,fincllonlat(1,4))
   equivalence (fincl5lonlat,fincllonlat(1,5))
   equivalence (fincl6lonlat,fincllonlat(1,6))

   character*8 fexcl1(pflds)
   character*8 fexcl2(pflds)
   character*8 fexcl3(pflds)
   character*8 fexcl4(pflds)
   character*8 fexcl5(pflds)
   character*8 fexcl6(pflds)

   equivalence (fexcl1,fexcl(1,1))
   equivalence (fexcl2,fexcl(1,2))
   equivalence (fexcl3,fexcl(1,3))
   equivalence (fexcl4,fexcl(1,4))
   equivalence (fexcl5,fexcl(1,5))
   equivalence (fexcl6,fexcl(1,6))

   character*10 fhstpr1(pflds)
   character*10 fhstpr2(pflds)
   character*10 fhstpr3(pflds)
   character*10 fhstpr4(pflds)
   character*10 fhstpr5(pflds)
   character*10 fhstpr6(pflds)

   equivalence (fhstpr1,fhstpr(1,1))
   equivalence (fhstpr2,fhstpr(1,2))
   equivalence (fhstpr3,fhstpr(1,3))
   equivalence (fhstpr4,fhstpr(1,4))
   equivalence (fhstpr5,fhstpr(1,5))
   equivalence (fhstpr6,fhstpr(1,6))

   character*10 fwrtpr1(pflds)
   character*10 fwrtpr2(pflds)
   character*10 fwrtpr3(pflds)
   character*10 fwrtpr4(pflds)
   character*10 fwrtpr5(pflds)
   character*10 fwrtpr6(pflds)

   equivalence (fwrtpr1,fwrtpr(1,1))
   equivalence (fwrtpr2,fwrtpr(1,2))
   equivalence (fwrtpr3,fwrtpr(1,3))
   equivalence (fwrtpr4,fwrtpr(1,4))
   equivalence (fwrtpr5,fwrtpr(1,5))
   equivalence (fwrtpr6,fwrtpr(1,6))
#endif
!
! Overloading assignment operator
!
   interface assignment (=)
      module procedure hbuf_assigned_to_hbuf
      module procedure hbuf_assigned_to_real8
   end interface
!
! Generic procedures
!
   interface allocate_hbuf
      module procedure allocate_hbuf2d
      module procedure allocate_hbuf3d
   end interface

   interface deallocate_hbuf
      module procedure deallocate_hbuf2d
      module procedure deallocate_hbuf3d
   end interface

   interface nullify_hbuf
      module procedure nullify_hbuf2d
      module procedure nullify_hbuf3d
   end interface
!
! Public entities
!

!
! Filename specifiers for history, initial files and restart history files
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
!
   character(len=256) :: ifilename_spec = '%c.cam2.i.%y-%m-%d-%s.nc'  ! Initial files
   character(len=256) :: rhfilename_spec = '%c.cam2.rh%t.%y-%m-%d-%s' ! history restart
   character(len=256), public :: hfilename_spec(ptapes) = (/ (' ', i=1, ptapes) /) ! filename specifyer

! Needed by anyone calling addfld

   integer, parameter, public :: phys_decomp = 1     ! flag indicates physics decomposition
   integer, parameter, public :: dyn_decomp  = 2     ! flag indicates dynamics decomposition

! To allow parameterizations to initialize arrays to the fillvalue
! THIS NEEDS TO BE FIXED.  No parameterization should be allowed access to fillvalue

   public :: fillvalue

! Needed by cam

   public :: bldfld

! Needed by initext

   public :: nhtfrq, mfilt, inithist, ctitle

! Needed by read_namelist

#ifdef INTEL64
   public :: fincl
   public :: fincllonlat
   public :: fexcl
   public :: fhstpr
   public :: fwrtpr
#else
   public :: fincl, fincl1, fincl2, fincl3, fincl4, fincl5, fincl6
   public :: fincllonlat, fincl1lonlat, fincl2lonlat, fincl3lonlat, fincl4lonlat, fincl5lonlat, fincl6lonlat
   public :: fexcl, fexcl1, fexcl2, fexcl3, fexcl4, fexcl5, fexcl6
   public :: fhstpr, fhstpr1, fhstpr2, fhstpr3, fhstpr4, fhstpr5, fhstpr6
   public :: fwrtpr, fwrtpr1, fwrtpr2, fwrtpr3, fwrtpr4, fwrtpr5, fwrtpr6
#endif

   public :: pflds, ptapes, empty_htapes, nhstpr, ndens
   public :: avgflag_pertape

#ifdef CRM
   public :: ngroup
   public :: tape
#endif

! Needed by stepon

   public :: hstwr
   public :: nfils

! Functions

   public :: write_restart_history     ! Write restart history data
   public :: read_restart_history      ! Read restart history data
   public :: wshist                    ! Write files out
   public :: write_inithist            ! Write the initial file
   public :: outfldcol                 ! Output a single column of a field
   public :: outfld                    ! Output a field
   public :: intht                     ! Initialization
   public :: wrapup                    ! Archive history files at end of run
   public :: addfld                    ! Add a field to history file
   public :: add_default               ! Add the default fields
   public :: get_hfilepath             ! Return history filename
   public :: get_mtapes                ! Return the number of tapes being used
   public :: get_hist_restart_filepath ! Return the full filepath to the history restart file
#ifdef CRM
   public :: sec2hms
   public :: date2yyyymmdd
#endif

CONTAINS

   subroutine intht ()
!----------------------------------------------------------------------- 
! 
! Purpose: Initialize history file handler for initial or continuation run.
!          For example, on an initial run, this routine initializes "mtapes"
!          history files.  On a restart or regeneration  run, this routine 
!          only initializes history files declared beyond what existed on the 
!          previous run.  Files which already existed on the previous run have 
!          already been initialized (i.e. named and opened) in routine RESTRT.
! 
! Method: Loop over tapes and fields per tape setting appropriate variables and
!         calling appropriate routines
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ioFileMod
      use time_manager, only: get_curr_date, get_curr_time
#if ( defined SPMD )
      use mpishorthand
#endif
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer :: t, f              ! tape, field indices
      integer :: coldimin          ! column dimension of model array
      integer :: begver            ! on-node vert start index
      integer :: endver            ! on-node vert end index
      integer :: begdim3           ! on-node chunk or lat start index
      integer :: enddim3           ! on-node chunk or lat end index
      type (dim_index_3d) :: dimind  ! 3-D dimension index
      integer :: day, sec          ! day and seconds from base date
      integer :: i                 ! index
!
! Get users logname and machine hostname
!
      if ( masterproc )then
         logname = ' '
         call getenv ('LOGNAME',logname)
         if (logname=='        ') then
            write(6,*)'READ_NAMELIST: Cannot find LOGNAME environment variable'
            call endrun
         end if
         host = ' '
         call getenv ('HOST',host)
      end if
!
! Set default history history contents
!
      call h_default ()
!
! Override averaging flag for all fields on a particular tape if namelist input so specifies
!
      do t=1,ptapes
         if (avgflag_pertape(t) /= ' ') then
            call h_override (t)
         end if
      end do
!
! Define field list information for all history files.  
! Update mtapes to reflect *current* number of history files (note, 
! restart and regen runs can have additional auxiliary history files
! declared).
!
      call fldlst ()
!
! Loop over max. no. of history files permitted  
!
      call get_curr_time(day, sec)  ! elapased time since reference date
      do t=1,mtapes
         nfils(t) = 0            ! no. of time samples in hist. file no. t

! Time a beginning of current averaging interval.

         beg_time(t) = day + sec/86400._r8

      end do
!
! Check that the number of history files declared does not exceed
! the maximum allowed.
!
      if (mtapes > ptapes) then
         write(6,*) 'INTHT: Too many history files declared, max=',ptapes
         write(6,*)'To increase, change parameter ptapes.'
         call endrun
      end if
!
! Initialize history variables
!
      do t=1,mtapes
         do f=1,nflds(t)
            coldimin = tape(t)%hlist(f)%field%coldimin
            begver   = tape(t)%hlist(f)%field%begver
            endver   = tape(t)%hlist(f)%field%endver
            begdim3  = tape(t)%hlist(f)%field%begdim3
            enddim3  = tape(t)%hlist(f)%field%enddim3

            dimind = dim_index_3d (1,coldimin,begver,endver,begdim3,enddim3)
            call allocate_hbuf (tape(t)%hlist(f)%hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
            tape(t)%hlist(f)%hbuf = 0._r8

            allocate (tape(t)%hlist(f)%nacs(coldimin,begdim3:enddim3))
            do i = begdim3, enddim3
               tape(t)%hlist(f)%nacs(:coldimin,i) = 0
            end do
         end do
      end do
      
      return
   end subroutine intht

!#######################################################################

   subroutine write_restart_history (nrg, luhrest)

      use binary_io
      use ioFileMod
#ifdef SPMD
# if ( defined STAGGERED )
      use spmd_dyn, only: npes, compute_gsfactors, comm_y
      use pmgrid, only: myid_z, strip3dxzy, strip3dxzyp, strip2d, twod_decomp
      use parutilitiesmodule, only: commglobal, pargatherreal, pargatherint, pargatherreal4
# else
      use spmd_dyn, only: npes, compute_gsfactors
# endif
      use mpishorthand
#endif
      use phys_grid, only: gather_chunk_to_field_int
      use dycore, only: dycore_is

!--------------------------------------------------------------------------------------------------

!
! Arguments
!
      integer, intent(in) :: nrg        ! unit number
      integer, intent(in) :: luhrest    ! unit number
!
! Local workspace
!
      integer t,f                       ! Tape, field indices
      integer numlev                    ! number of vertical levels (dimension and loop)
      integer ioerr                     ! write error status
      integer coldimin                  ! column dimension of model array
      character(len=256) :: fname       ! History restart filename

#ifdef SPMD
      integer :: numsend                ! number of items to be sent
      integer :: numrecv(0:npes-1)      ! number of items to be received
      integer :: displs(0:npes-1)       ! displacement array
      integer :: numowned               ! number of items owned by this MPI task
      integer :: mpireal                ! MPI real data type
#endif

      type (hbuffer_3d) :: hbuf         ! full-field history buffer written by master
      integer, pointer :: fullnacs(:,:) ! full-field accumulation counter written by master
      type (dim_index_3d) :: dimind     ! 3-D dimension index
!
!-----------------------------------------------------------------------
! Write the history restart data if necessary
!-----------------------------------------------------------------------
!
      where (hstwr(:))
         rgnht(:) = .false.
      elsewhere
         rgnht(:) = .true.
      end where

      if (masterproc) then
         write (nrg) rgnht, mtapes, varid, fincl, fexcl
      end if
!
! If a history tape is not currently being disposed then write a history buffer restart file.
! History restart info: some f90 compilers (e.g. SGI) complain about I/O of derived types which 
! have pointer components, so explicitly write each one.
!
      do t=1,mtapes
         if (masterproc) then
            write (nrg, iostat=ioerr) nhtfrq(t), nflds(t), nfils(t), mfilt(t), &
                                      nfpath(t), cpath(t), nhfil(t),  &
                                      nhstpr(t), ndens(t), ncprec(t), beg_time(t)
            if (ioerr /= 0 ) then
               write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
               call endrun
            end if

            write (nrg,iostat=ioerr) ngroup(t)
            if (ngroup(t) .gt. 0) then 
               do f=1,ngroup(t)
                  write (nrg,iostat=ioerr) tape(t)%column(f)
               end do
            else
               write (nrg,iostat=ioerr) tape(t)%column(1)
            end if
            do f=1,nflds(t)
               write (nrg,iostat=ioerr) tape(t)%hlist(f)%field%name,        &
                                        tape(t)%hlist(f)%field%long_name,   &
                                        tape(t)%hlist(f)%field%units,       &
                                        tape(t)%hlist(f)%field%numlev,      &
                                        tape(t)%hlist(f)%field%decomp_type, &
                                        tape(t)%hlist(f)%avgflag,           &
                                        tape(t)%hlist(f)%time_op,           &
                                        tape(t)%hlist(f)%hbuf_prec,         &
                                        tape(t)%hlist(f)%hwrt_prec,         &
                                        tape(t)%hlist(f)%field%flag_xyfill, &
                                        tape(t)%hlist(f)%field%flag_isccplev 
               if (ioerr /= 0 ) then
                  write(6,*)'WRITE_RESTART_HISTORY: End or error condition writing ', &
                       'history restart field ',f,' from tape ',t
                  call endrun
               end if
            end do

         end if

         if (rgnht(t)) then
            if (masterproc) then
               fname = interpret_filename_spec( rhfilename_spec, number=(t-1) )
               hrestpath(t) = trim(get_archivedir('rest'))//fname
               write (nrg,iostat=ioerr) hrestpath(t)
               if (ioerr /= 0 ) then
                  write(6,*)'WRITE_RESTART_HISTORY: End or error condition writing ', &
                       'history restart filename to master restart file'
                  call endrun
               end if
               call opnfil (fname, luhrest, 'u')
            endif

            do f=1,nflds(t)
               coldimin =  tape(t)%hlist(f)%field%coldimin
               numlev   =  tape(t)%hlist(f)%field%numlev
#ifdef SPMD
               if (masterproc) then
                  dimind = dim_index_3d (1,plon,1,numlev,1,plat)
                  call allocate_hbuf (hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
                  allocate (fullnacs(plon,plat))
               else
                  call assoc_hbuf_with_nothing (hbuf,tape(t)%hlist(f)%hbuf_prec)
                  fullnacs => nothing_int
               end if

               if (tape(t)%hlist(f)%hbuf_prec == 8) then
                  mpireal = mpir8
               else
                  mpireal = mpir4
               end if

               select case (tape(t)%hlist(f)%field%decomp_type)
               case (phys_decomp)
                  call gather_chunk_to_field_hbuf (1, numlev, 1, plon, tape(t)%hlist(f)%hbuf, hbuf)
                  call gather_chunk_to_field_int (1, 1, 1, plon, tape(t)%hlist(f)%nacs, fullnacs)
               case (dyn_decomp)
                  if ( dycore_is('LR') )then
# if ( defined STAGGERED )
! NEW LR CODING
                     if (tape(t)%hlist(f)%hbuf_prec == 8) then
                        select case (numlev)
                        case (1)
                           if (myid_z .eq. 0) call pargatherreal(comm_y, 0,  &
                                  tape(t)%hlist(f)%hbuf%buf8, strip2d, hbuf%buf8)
                        case (plev)
                           call pargatherreal(commglobal, 0, tape(t)%hlist(f)%hbuf%buf8, &
                                       strip3dxzy, hbuf%buf8)
                        case (plevp)
                           call pargatherreal(commglobal, 0, tape(t)%hlist(f)%hbuf%buf8, &
                                       strip3dxzyp, hbuf%buf8)
                        case default
                           write(6,*)'WRITE_RESTART_HISTORY: invalid number of levels=', numlev
                           call endrun ()
                        end select
                        if (myid_z .eq. 0) call pargatherint(comm_y, 0, &
                           tape(t)%hlist(f)%nacs, strip2d, fullnacs)
                     else
                        select case (numlev)
                        case (1)
                           if (myid_z .eq. 0) call pargatherreal4(comm_y, 0,  &
                                  tape(t)%hlist(f)%hbuf%buf4, strip2d, hbuf%buf4)
                        case (plev)
                           call pargatherreal4(commglobal, 0, tape(t)%hlist(f)%hbuf%buf4, &
                                       strip3dxzy, hbuf%buf4)
                        case (plevp)
                           call pargatherreal4(commglobal, 0, tape(t)%hlist(f)%hbuf%buf4, &
                                       strip3dxzyp, hbuf%buf4)
                        case default
                           write(6,*)'WRITE_RESTART_HISTORY: invalid number of levels=', numlev
                           call endrun ()
                        end select
                        if (myid_z .eq. 0) call pargatherint(comm_y, 0, &
                           tape(t)%hlist(f)%nacs, strip2d, fullnacs)
                     endif
# endif
                  else
                     numowned = coldimin*numlev
                     call compute_gsfactors (numowned, numsend, numrecv, displs)
                     call mpigatherv_hbuf (tape(t)%hlist(f)%hbuf, numsend, mpireal, hbuf, numrecv, &
                                           displs, mpireal, 0, mpicom)

                     numowned = coldimin
                     call compute_gsfactors (numowned, numsend, numrecv, displs)
                     call mpigatherv (tape(t)%hlist(f)%nacs, numsend, mpiint, fullnacs, numrecv, &
                                      displs, mpiint, 0, mpicom)
                  endif
               end select

               if (masterproc) then
                  call write_hbuf (hbuf, luhrest, ioerr)
                  call deallocate_hbuf (hbuf)
                  write (luhrest,iostat=ioerr) fullnacs
                  deallocate (fullnacs)
               else
                  call nullify_hbuf (hbuf)
                  nullify (fullnacs)
               end if
#else
               select case (tape(t)%hlist(f)%field%decomp_type)
               case (phys_decomp)
                  dimind = dim_index_3d (1,plon,1,numlev,1,plat)
                  call allocate_hbuf (hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
                  call gather_chunk_to_field_hbuf (1, numlev, 1, plon, tape(t)%hlist(f)%hbuf, hbuf)
                  call write_hbuf (hbuf, luhrest, ioerr)
                  call deallocate_hbuf (hbuf)
                  allocate (fullnacs(plon,plat))
                  call gather_chunk_to_field_int (1, 1, 1, plon, tape(t)%hlist(f)%nacs, fullnacs)
                  write (luhrest,iostat=ioerr) fullnacs
                  deallocate (fullnacs)
               case (dyn_decomp)
                  call write_hbuf (tape(t)%hlist(f)%hbuf, luhrest, ioerr)
                  write (luhrest,iostat=ioerr) tape(t)%hlist(f)%nacs
               end select
#endif
            end do
            if (masterproc) then
               close (unit=luhrest)
               call putfil (fname, hrestpath(t), mss_wpass, mss_irt, .true.)
            end if
         end if
      end do

      return
   end subroutine write_restart_history

!#######################################################################

   subroutine read_restart_history (nrg, luhrest)

      use ppgrid,    only: begchunk, endchunk
      use phys_grid, only: get_ncols_p, scatter_field_to_chunk_int
      use rgrid,     only: nlon
      use dycore,    only: dycore_is
      use binary_io
      use ioFileMod
#ifdef SPMD
# if ( defined STAGGERED )
      use pmgrid, only: strip3dxzy, strip3dxzyp, strip2d, twod_decomp
      use restart_dynamics, only: lrreadin, lrreadini, lrreadin4
# endif
      use mpishorthand
#endif
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nrg            ! unit number
      integer, intent(in) :: luhrest        ! unit number
!
! Local workspace
!
      integer t, f                     ! tape, field indices
      integer c                        ! chunk or lat index
      integer lenc                     ! length of useful character data
      integer numlev                   ! number of vertical levels (dimension and loop)
      integer begver                   ! on-node vert start index
      integer endver                   ! on-node vert end index
      integer ioerr                    ! error code from read()
      integer coldimin                 ! column dimension of model array
      integer begdim3                  ! on-node chunk or lat start index
      integer enddim3                  ! on-node chunk or lat end index
      integer ncol                     ! number of active columns per chunk
      integer lenarr                   ! global size of array to be read
      integer g                        ! group index

      character(len=80)  :: locfn       ! Local filename

      integer, pointer :: nacs(:,:)    ! accumulation counter

      type (hbuffer_3d) :: hbuf             ! history buffer
      integer,  pointer :: fullnacs(:,:)    ! accumulation buffer
      type (dim_index_3d) :: dimind    ! 3-D dimension index

      if (masterproc) then

         read (nrg,iostat=ioerr) rgnht, mtapes, varid, fincl, fexcl
         if (ioerr /= 0 ) then
            write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if

         do t=1,mtapes
            read (nrg,iostat=ioerr) nhtfrq(t), nflds(t), nfils(t), mfilt(t), &
                                    nfpath(t), cpath(t), nhfil(t), &
                                    nhstpr(t), ndens(t), ncprec(t), beg_time(t)
            if (ioerr /= 0) then
               write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
               call endrun
            end if

            read (nrg,iostat=ioerr) ngroup(t)
            if (ioerr /= 0) then
               write (6,*) 'READ ioerror for group number',ioerr,' on i/o unit = ',nrg
               call endrun
            end if

            if (ngroup(t) .gt. 0) then
               allocate(tape(t)%column(ngroup(t)))
               do f=1,ngroup(t)
                  read (nrg,iostat=ioerr) tape(t)%column(f)
               end do
               if (ioerr /= 0) then
                  write (6,*) 'READ ioerror for column data',ioerr,' on i/o unit = ',nrg
                  call endrun
               end if
            else
               allocate(tape(t)%column(1))
               read (nrg,iostat=ioerr) tape(t)%column(1)
            end if

            do f=1,nflds(t)
               read (nrg,iostat=ioerr) tape(t)%hlist(f)%field%name,        &
                                       tape(t)%hlist(f)%field%long_name,   &
                                       tape(t)%hlist(f)%field%units,       &
                                       tape(t)%hlist(f)%field%numlev,      &
                                       tape(t)%hlist(f)%field%decomp_type, &
                                       tape(t)%hlist(f)%avgflag,           &
                                       tape(t)%hlist(f)%time_op,           &
                                       tape(t)%hlist(f)%hbuf_prec,         &
                                       tape(t)%hlist(f)%hwrt_prec,         &
                                       tape(t)%hlist(f)%field%flag_xyfill, &
                                       tape(t)%hlist(f)%field%flag_isccplev 

               if (ioerr /= 0) then
                  write(6,*)'READ_RESTART_HISTORY: ', &
                     'End or error condition reading history restart field ', &
                     f,' from tape ',t
                  write(6,*)'ioerr=',ioerr
                  call endrun
               end if
            end do
            if (rgnht(t)) then
               read (nrg,iostat=ioerr) hrestpath(t)
               if (ioerr /= 0) then
                  write (6,*) 'READ ioerror on read of filename ',ioerr,' on i/o unit = ',nrg
                  call endrun
               end if
            end if
         end do

      end if

#if ( defined SPMD )
      call mpibcast (rgnht  ,ptapes   ,mpilog ,0,mpicom)

      call mpibcast (mtapes ,1        ,mpiint ,0,mpicom)
      call mpibcast (nhtfrq ,ptapes   ,mpiint ,0,mpicom)
      call mpibcast (nflds  ,ptapes   ,mpiint ,0,mpicom)
      call mpibcast (nfils  ,ptapes   ,mpiint ,0,mpicom)
      call mpibcast (mfilt  ,ptapes   ,mpiint ,0,mpicom)
#ifdef CRM
      call mpibcast (ngroup ,ptapes   ,mpiint ,0,mpicom)
      ! We need to allocate our tape columns for each
      ! node; we've already done it for masterproc, and the mpibcast
      ! should include a barrier (right?) so the nodes don't run
      ! ahead of the master
      if (.not. masterproc) then
         do t=1,mtapes
            if (ngroup(t) .gt. 0) then
               allocate(tape(t)%column(ngroup(t)))
            else
               allocate(tape(t)%column(1))
            end if
         end do
      endif
#endif

      do t=1,mtapes
#ifdef CRM
         do g=1,ngroup(t)
            call mpibcast (tape(t)%column(g)%columnlon,2   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%columnlat,2   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%num_lats, 1   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%num_lons, 1   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%lon_name,max_chars   ,mpichar ,0,mpicom)
            call mpibcast (tape(t)%column(g)%lat_name,max_chars   ,mpichar ,0,mpicom)
         end do
#endif
         do f=1,nflds(t)
            call mpibcast (tape(t)%hlist(f)%field%numlev,       1,        mpiint, 0,mpicom)
            call mpibcast (tape(t)%hlist(f)%field%decomp_type,  1,        mpiint, 0,mpicom)
            call mpibcast (tape(t)%hlist(f)%field%name,         8,        mpichar,0,mpicom)
            call mpibcast (tape(t)%hlist(f)%field%units,        max_chars,mpichar,0,mpicom)
            call mpibcast (tape(t)%hlist(f)%avgflag,            1,        mpichar,0,mpicom)
            call mpibcast (tape(t)%hlist(f)%hbuf_prec,          1,        mpiint, 0,mpicom)
            call mpibcast (tape(t)%hlist(f)%hwrt_prec,          1,        mpiint, 0,mpicom)
            call mpibcast (tape(t)%hlist(f)%field%flag_xyfill,  1,        mpilog, 0,mpicom)
            call mpibcast (tape(t)%hlist(f)%field%flag_isccplev,1,        mpilog, 0,mpicom)
         end do
      end do
#endif
!
! Allocate space for history buffers and initialize
!
      do t=1,mtapes
         do f=1,nflds(t)
            numlev   = tape(t)%hlist(f)%field%numlev
            select case (tape(t)%hlist(f)%field%decomp_type)
            case (phys_decomp)
               tape(t)%hlist(f)%field%begdim3 = begchunk
               tape(t)%hlist(f)%field%enddim3 = endchunk
               tape(t)%hlist(f)%field%begver = 1
               tape(t)%hlist(f)%field%endver = numlev
               allocate (tape(t)%hlist(f)%field%colperdim3(begchunk:endchunk))
               do c=begchunk,endchunk
                  ncol = get_ncols_p(c)
                  tape(t)%hlist(f)%field%colperdim3(c) = ncol
               end do
               tape(t)%hlist(f)%field%coldimin = pcols
            case (dyn_decomp)
               tape(t)%hlist(f)%field%begdim3 = beglat
               tape(t)%hlist(f)%field%enddim3 = endlat
               if ( dycore_is('LR') )then
# if ( defined STAGGERED )
                  select case (numlev)
                  case (1)
                     tape(t)%hlist(f)%field%begver = 1
                     tape(t)%hlist(f)%field%endver = 1
                  case (plev)
                     tape(t)%hlist(f)%field%begver = beglev
                     tape(t)%hlist(f)%field%endver = endlev
                  case (plevp)
                     tape(t)%hlist(f)%field%begver = beglev
                     tape(t)%hlist(f)%field%endver = endlevp
                  case default
                     write(6,*)'READ_RESTART_HISTORY: invalid number of levels=', numlev
                     call endrun ()
                  end select
# endif
               else
                  tape(t)%hlist(f)%field%begver = 1
                  tape(t)%hlist(f)%field%endver = numlev
               endif
               allocate (tape(t)%hlist(f)%field%colperdim3(beglat:endlat))
               do c=beglat,endlat
                  tape(t)%hlist(f)%field%colperdim3(c) = nlon(c)
               end do
               tape(t)%hlist(f)%field%coldimin = plon
            case default
               write(6,*)'READ_RESTART_HISTORY: bad decomp_type=',tape(t)%hlist(f)%field%decomp_type
               call endrun ()
            end select

            coldimin = tape(t)%hlist(f)%field%coldimin
            begdim3  = tape(t)%hlist(f)%field%begdim3
            enddim3  = tape(t)%hlist(f)%field%enddim3
            begver   = tape(t)%hlist(f)%field%begver
            endver   = tape(t)%hlist(f)%field%endver

            dimind = dim_index_3d (1,coldimin,begver,endver,begdim3,enddim3)
            call allocate_hbuf (tape(t)%hlist(f)%hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
            tape(t)%hlist(f)%hbuf = 0._r8
            allocate (tape(t)%hlist(f)%nacs(coldimin,begdim3:enddim3))
            tape(t)%hlist(f)%nacs(:coldimin,begdim3:enddim3) = 0
         end do
      end do
!
! must allocate space and initialize field column names (this field is dimensioned by number of column groups)
!
      do t=1,mtapes
         do f=1,nflds(t)
            if (ngroup(t) .gt. 0) then
               allocate(tape(t)%hlist(f)%field_column_name(ngroup(t)))
               do i=1,ngroup(t)
                  tape(t)%hlist(f)%field_column_name(i) = trim(tape(t)%hlist(f)%field%name) // "_" // &
                       trim(tape(t)%column(i)%lon_name) // "_" // trim(tape(t)%column(i)%lat_name)
               end do
            else
               allocate(tape(t)%hlist(f)%field_column_name(1))
               tape(t)%hlist(f)%field_column_name(1) = ' '
            end if
         end do
      end do
!
!-----------------------------------------------------------------------
! Read history restart files
!-----------------------------------------------------------------------
!
! Loop over the total number of history files declared and
! read the pathname for any history restart files
! that are present (if any). Test to see if the run is a restart run
! AND if any history buffer regen files exist (rgnht=.T.). Note, rgnht 
! is preset to false, reset to true in routine WSDS if hbuf restart files
! are written and saved in the master restart file. Each history buffer
! restart file is then obtained.
! Note: some f90 compilers (e.g. SGI) complain about I/O of 
! derived types which have pointer components, so explicitly read each one.
! 
      do t=1,mtapes
         if (rgnht(t)) then
!
! Open history restart file
!
            if (masterproc) then
               call getfil (hrestpath(t), locfn)
               call opnfil (locfn, luhrest, 'u')
            end if
!
! Read history restart file
!
            do f=1,nflds(t)
               coldimin   =  tape(t)%hlist(f)%field%coldimin
               begdim3    =  tape(t)%hlist(f)%field%begdim3
               enddim3    =  tape(t)%hlist(f)%field%enddim3
               numlev     =  tape(t)%hlist(f)%field%numlev
               nacs       => tape(t)%hlist(f)%nacs(:coldimin,begdim3:enddim3)
#ifdef SPMD
               if (masterproc) then
                  dimind = dim_index_3d (1,plon,1,numlev,1,plat)
                  call allocate_hbuf (hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
                  allocate (fullnacs(plon,plat))
               else
                  call assoc_hbuf_with_nothing (hbuf,tape(t)%hlist(f)%hbuf_prec)
                  fullnacs => nothing_int
               end if

               select case (tape(t)%hlist(f)%field%decomp_type)
               case (phys_decomp)
                  if (masterproc) then
                     call read_hbuf (hbuf,luhrest,ioerr)
                     read (luhrest) fullnacs
                  end if
                  call scatter_field_to_chunk_hbuf (1, numlev, 1, plon, hbuf, tape(t)%hlist(f)%hbuf)
                  call scatter_field_to_chunk_int (1, 1, 1, plon, fullnacs, tape(t)%hlist(f)%nacs)
               case (dyn_decomp)
                  if ( dycore_is('LR') )then
# if ( defined STAGGERED )
! NEW LR CODING
                     if (tape(t)%hlist(f)%hbuf_prec == 8) then
                        lenarr = plon*numlev*plat
                        select case (numlev)
                        case (1)
                           call lrreadin(luhrest, strip2d, tape(t)%hlist(f)%hbuf%buf8,  &
                                         lenarr, 2)
                        case (plev)
                           call lrreadin(luhrest, strip3dxzy, tape(t)%hlist(f)%hbuf%buf8,  &
                                         lenarr, 3)
                        case (plevp)
                           call lrreadin(luhrest, strip3dxzyp, tape(t)%hlist(f)%hbuf%buf8,  &
                                         lenarr, 3)
                        case default
                           write(6,*)'READ_RESTART_HISTORY: invalid number of levels=', numlev
                           call endrun ()
                        end select
                        lenarr = plon*plat
                        call lrreadini(luhrest, strip2d, nacs,  &
                                      lenarr, 2)
                     else
                        lenarr = plon*numlev*plat
                        select case (numlev)
                        case (1)
                           call lrreadin4(luhrest, strip2d, tape(t)%hlist(f)%hbuf%buf4,  &
                                         lenarr, 2)
                        case (plev)
                           call lrreadin4(luhrest, strip3dxzy, tape(t)%hlist(f)%hbuf%buf4,  &
                                         lenarr, 3)
                        case (plevp)
                           call lrreadin4(luhrest, strip3dxzyp, tape(t)%hlist(f)%hbuf%buf4,  &
                                         lenarr, 3)
                        case default
                           write(6,*)'READ_RESTART_HISTORY: invalid number of levels=', numlev
                           call endrun ()
                        end select
                        lenarr = plon*plat
                        call lrreadini(luhrest, strip2d, nacs,  &
                                      lenarr, 2)
                     endif
# endif
                  else
                     call readin_hbuf(luhrest, tape(t)%hlist(f)%hbuf, coldimin*numlev)
                     call readin_int (luhrest, nacs, coldimin)
                  endif
               case default
                  write(6,*)'READ_RESTART_HISTORY: bad decomp_type=',tape(t)%hlist(f)%field%decomp_type
                  call endrun ()
               end select

               if (masterproc) then
                  call deallocate_hbuf (hbuf)
                  deallocate (fullnacs)
               else
                  call nullify_hbuf (hbuf)
                  nullify (fullnacs)
               end if
#else
               select case (tape(t)%hlist(f)%field%decomp_type)
               case (phys_decomp)
                  dimind = dim_index_3d (1,plon,1,numlev,1,plat)
                  call allocate_hbuf (hbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
                  call read_hbuf (hbuf,luhrest,ioerr)
                  call scatter_field_to_chunk_hbuf (1, numlev, 1, plon, hbuf, tape(t)%hlist(f)%hbuf)
                  call deallocate_hbuf (hbuf)
                  allocate (fullnacs(plon,plat))
                  read (luhrest) fullnacs
                  call scatter_field_to_chunk_int (1, 1, 1, plon, fullnacs, tape(t)%hlist(f)%nacs)
                  deallocate (fullnacs)
               case (dyn_decomp)
                  call readin_hbuf(luhrest, tape(t)%hlist(f)%hbuf, coldimin*numlev)
                  call readin_int (luhrest, nacs, coldimin)
               case default
                  write(6,*)'READ_RESTART_HISTORY: bad decomp_type=',tape(t)%hlist(f)%field%decomp_type
                  call endrun ()
               end select
#endif
            end do
!          
! Done reading this history restart file
!
            if (masterproc) close (luhrest)

         end if  ! rgnht(t)
      end do     ! end of do mtapes loop
!
! If the history files are partially complete (contain less than
! mfilt(t) time samples, then get the files and open them.
!
      do t=1,mtapes
         if (masterproc .and. nfils(t) > 0) then
            call getfil (cpath(t), locfn)
            call wrap_open (locfn, NF_WRITE, nfid(t))
            call h_inquire (t)
         end if
!
! If the history file is full, close the current unit
!
         if (nfils(t) >= mfilt(t)) then
            if (masterproc) then
               write(6,*)'READ_RESTART_HISTORY: nf_close(',nfid(t),')=',nhfil(t)
               call wrap_close (nfid(t))
            end if
            nfils(t) = 0
         end if
      end do
!
! set flag indicating h-tape contents are now defined (needed by addfld)
!
      htapes_defined = .true.      

      return
   end subroutine read_restart_history

!#######################################################################

   character(len=nlen) function get_hfilepath( tape )
!----------------------------------------------------------------------- 
! 
! Purpose: Return full filepath of history file for given tape number
! This allows public read access to the filenames without making
! the filenames public data.
!
!----------------------------------------------------------------------- 
  integer, intent(in) :: tape  ! Tape number

  get_hfilepath = cpath( tape )
  end function get_hfilepath

!#######################################################################

   character(len=nlen) function get_hist_restart_filepath( tape )
!----------------------------------------------------------------------- 
! 
! Purpose: Return full filepath of restart file for given tape number
! This allows public read access to the filenames without making
! the filenames public data.
!
!----------------------------------------------------------------------- 
  integer, intent(in) :: tape  ! Tape number

  get_hist_restart_filepath = hrestpath( tape )
  end function get_hist_restart_filepath

!#######################################################################

  integer function get_mtapes( )
!----------------------------------------------------------------------- 
! 
! Purpose: Return the number of tapes being used.
! This allows public read access to the number of tapes without making
! mtapes public data.
!
!----------------------------------------------------------------------- 
  get_mtapes = mtapes
  end function get_mtapes


!#######################################################################

   subroutine fldlst ()
!----------------------------------------------------------------------- 
! 
! Purpose: Define the contents of each history file based on namelist input for initial or branch
! run, and restart data if a restart run.
!          
! Method: Use arrays fincl and fexcl to modify default history tape contents.
!         Then sort the result alphanumerically for later use by OUTFLD to
!         allow an n log n search time.
! 
#ifdef SPMD
      use mpishorthand
#endif
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
      integer t, f                   ! tape, field indices
      integer ff                     ! index into include, exclude and fprec list
      character(len=8) :: name       ! field name portion of fincl (i.e. no avgflag separator)
      character(len=8) :: mastername ! name from masterlist field
      character(len=max_chars) :: tmpcolumn_name,latlonname,latlonnamep1 ! tmp char fields
      character(len=1) :: avgflag    ! averaging flag
      character(len=1) :: prec_acc   ! history buffer precision flag
      character(len=1) :: prec_wrt   ! history buffer write precision flag
      integer g  ! group index

      type (hentry) :: tmp           ! temporary used for swapping
      type (column_info) :: tmpcolumn ! temporary used for swapping
!
! First ensure contents of fincl, fexcl, fhstpr and fwrtpr are all valid names
!
      do t=1,ptapes
         f = 1
         do while (f < pflds .and. fincl(f,t) /= ' ')
            name = getname (fincl(f,t))
            do ff=1,nfmaster
               mastername = masterlist(ff)%field%name
               if (name == mastername) exit
            end do
            if (name /= mastername) then
               write(6,*)'FLDLST: ', trim(name), ' in fincl(', f, ') not found'
               call endrun
            end if
            f = f + 1
         end do

         f = 1
         do while (f < pflds .and. fexcl(f,t) /= ' ')
            do ff=1,nfmaster
               mastername = masterlist(ff)%field%name
               if (fexcl(f,t) == mastername) exit
            end do
            if (fexcl(f,t) /= mastername) then
               write(6,*)'FLDLST: ', fexcl(f,t), ' in fexcl(', f, ') not found'
               call endrun
            end if
            f = f + 1
         end do

         f = 1
         do while (f < pflds .and. fhstpr(f,t) /= ' ')
            name = getname (fhstpr(f,t))
            do ff=1,nfmaster
               mastername = masterlist(ff)%field%name
               if (name == mastername) exit
            end do
            if (name /= mastername) then
               write(6,*)'FLDLST: ', trim(name), ' in fhstpr(', f, ') not found'
               call endrun
            end if
            do ff=1,f-1                 ! If duplicate entry is found, stop
               if (trim(name) == trim(getname(fhstpr(ff,t)))) then
                  write(6,*)'FLDLST: Duplicate field ', name, ' in fhstpr'
                  call endrun
               end if
            end do
            f = f + 1
         end do

         f = 1
         do while (f < pflds .and. fwrtpr(f,t) /= ' ')
            name = getname (fwrtpr(f,t))
            do ff=1,nfmaster
               mastername = masterlist(ff)%field%name
               if (name == mastername) exit
            end do
            if (name /= mastername) then
               write(6,*)'FLDLST: ', trim(name), ' in fwrtpr(', f, ') not found'
               call endrun
            end if
            do ff=1,f-1                 ! If duplicate entry is found, stop
               if (trim(name) == trim(getname(fwrtpr(ff,t)))) then
                  write(6,*)'FLDLST: Duplicate field ', name, ' in fwrtpr'
                  call endrun
               end if
            end do
            f = f + 1
         end do
      end do
!
! If kind values r8 and r4 are identical, set accumulation precision to 8 bytes
!
      if (r4 == r8 .and. any(nhstpr == 4)) then
         nhstpr(:) = 8
         if (masterproc) then
            write(6,*) 'FLDLST: Set nhstpr to 8 because kind values r8 and r4 are identical'
         end if
      end if

      nflds(:) = 0

      do t=1,ptapes
!
! Add the field to the tape if specified via namelist (FINCL[1-ptapes]), or if
! it is on by default and was not excluded via namelist (FEXCL[1-ptapes]).
! Also set history buffer accumulation and output precision values according
! to the values specified via namelist (FHSTPR[1-ptapes] and FWRTPR[1-ptapes])
! or, if not on the list, to the default values given by ndens(t) and
! nhstpr(t), respectively.
!
         do f=1,nfmaster
            mastername = masterlist(f)%field%name
            call list_index (fhstpr(1,t), mastername, ff)
            if (ff > 0) then
               prec_acc = getflag(fhstpr(ff,t))
            else
               prec_acc = ' '
            end if

            call list_index (fwrtpr(1,t), mastername, ff)
            if (ff > 0) then
               prec_wrt = getflag(fwrtpr(ff,t))
            else
               prec_wrt = ' '
            end if

            call list_index (fincl(1,t), mastername, ff)
            if (ff > 0) then
               avgflag = getflag (fincl(ff,t))
               call inifld (t, f, avgflag, prec_acc, prec_wrt)
            else if (.not. empty_htapes) then
               call list_index (fexcl(1,t), mastername, ff)
               if (ff == 0 .and. masterlist(f)%actflag(t)) then
                  call inifld (t, f, ' ', prec_acc, prec_wrt)
               end if
            end if
         end do
!
! If column output is specified make sure there are some fields defined
! for that tape
!
         if (nflds(t) .eq. 0 .and. fincllonlat(1,t) .ne. ' ') then
            write(6,*) 'FLDLST: Column output is specified for tape ',t,' but no fields defined for that tape.'
            call endrun()
         end if


!
! Specification of tape contents now complete.  Sort each list of active 
! entries for efficiency in OUTFLD.  Simple bubble sort.
!
         do f=nflds(t)-1,1,-1
            do ff=1,f

               if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then
            
                  tmp = tape(t)%hlist(ff)
                  tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                  tape(t)%hlist(ff+1) = tmp

               else if (tape(t)%hlist(ff  )%field%name == tape(t)%hlist(ff+1)%field%name) then
                  
                  write(6,*)'FLDLST: Duplicate field ', tape(t)%hlist(ff  )%field%name
                  write(6,*)'t,ff,name=',t,ff,tape(t)%hlist(ff  )%field%name
                  call endrun
            
               end if

            end do
         end do
!
! Bubble sort columns check for duplicates and rebuild field_column_names in newly sorted order
!
         if (ngroup(t) .gt. 0) then
            do i=ngroup(t)-1,1,-1
               do ff=1,i
                  latlonname=trim(tape(t)%column(ff)%lon_name) // "_" // trim(tape(t)%column(ff)%lat_name)
                  latlonnamep1=trim(tape(t)%column(ff+1)%lon_name) // "_" // trim(tape(t)%column(ff+1)%lat_name)
                  if (trim(latlonname) > trim(latlonnamep1)) then
                     tmpcolumn = tape(t)%column(ff)
                     tape(t)%column(ff) = tape(t)%column(ff+1)
                     tape(t)%column(ff+1) = tmpcolumn
                  else if (trim(latlonname) == trim(latlonnamep1)) then
                     write(6,*)'FLDLST: Duplicate column entry for tape ',t,'duplicate column name is ',trim(latlonname)
                     call endrun
                  end if
               end do
            end do
            do f=1,nflds(t)
               do i=1,ngroup(t)
                  if (ngroup(t) .gt. 0) then
                     tape(t)%hlist(f)%field_column_name(i) = trim(tape(t)%hlist(f)%field%name) // "_" // &
                          trim(tape(t)%column(i)%lon_name) // "_" // trim(tape(t)%column(i)%lat_name)
                  else 
                     tape(t)%hlist(f)%field_column_name(1) = ' '
                  end if
               end do
            end do
         end if

         if (masterproc) then
            if (nflds(t) > 0) then
               write(6,*)'FLDLST: Included fields tape ',t,'=',nflds(t)
            end if
   
            do f=1,nflds(t)
               if (ngroup(t) .gt. 0) then
                  if (f.eq.1) write(6,*)'Fields on this tape will be output as column data (FIELD_LON_LAT)'
                  do i=1,ngroup(t)
                     ff=(f-1)*ngroup(t)+i
                     write(6,*) ff,' ',trim(tape(t)%hlist(f)%field_column_name(i)), tape(t)%hlist(f)%field%numlev, &
                       ' ',tape(t)%hlist(f)%avgflag
                  end do
               else
                  write(6,*) f,' ',tape(t)%hlist(f)%field%name, tape(t)%hlist(f)%field%numlev, &
                            ' ',tape(t)%hlist(f)%avgflag
               end if
            end do
         end if
      end do
!
! Determine total number of active history tapes
!
      mtapes = 0
      do t=ptapes,1,-1
         if (nflds(t) > 0) then
            mtapes = t
            exit
         end if
      end do
!
! Ensure there are no "holes" in tape specification, i.e. empty tapes. Enabling
! holes should not be difficult if necessary.
!
      do t=1,mtapes
         if (nflds(t)  ==  0) then
            write(6,*)'FLDLST: Tape ',t,' is empty'
            call endrun
         end if
      end do
      
!
! Packing density, ndens: With netcdf, only 1 (nf_double) and 2 (nf_float)
! are allowed
! Accumulation precision, nhstpr, must be either 8 (real*8) or 4 (real*4)
!
      do t=1,mtapes
         if (ndens(t) == 1) then
            ncprec(t) = nf_double
         else if (ndens(t) == 2) then
            ncprec(t) = nf_float
         else
            write(6,*)'FLDLST: ndens must be 1 or 2'
            call endrun
         end if

         if (nhstpr(t) /= 8 .and. nhstpr(t) /= 4) then
            write(6,*)'FLDLST: nhstpr must be 8 or 4'
            call endrun
         end if
      end do

      if (masterproc) then
         if (nhtfrq(1) == 0) then
            write(6,*)'History File 1 write frequency MONTHLY'
         else
            write(6,*)'History File 1 write frequency ',nhtfrq(1)
         end if
         
         do t=2,mtapes
            write(6,*)'History File ',t,' write frequency ',nhtfrq(t)
         end do

         do t=1,mtapes
            write(6,*)'Accumulation precision history file ', t, '=', nhstpr(t)
            write(6,*)'Packing density history file ', t, '=', ndens(t)
            write(6,*)'Number of time samples per file (MFILT) for history file ',t,' is ',mfilt(t)
         end do
      end if
!
! set flag indicating h-tape contents are now defined (needed by addfld)
!
      htapes_defined = .true.      
#ifdef CRM
#ifdef SPMD
      call mpibcast (ngroup ,ptapes   ,mpiint ,0,mpicom)
      ! We need to allocate our tape columns for each
      ! node; we've already done it for masterproc, and the mpibcast
      ! should include a barrier (right?) so the nodes don't run
      ! ahead of the master
      if (.not. masterproc) then
         do t=1,mtapes
            if (ngroup(t) .gt. 0) then
               allocate(tape(t)%column(ngroup(t)))
            else
               allocate(tape(t)%column(1))
            end if
         end do
      endif

      do t=1,mtapes
         do g=1,ngroup(t)
            call mpibcast (tape(t)%column(g)%columnlon,2   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%columnlat,2   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%num_lats, 1   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%num_lons, 1   ,mpiint ,0,mpicom)
            call mpibcast (tape(t)%column(g)%lon_name,max_chars   ,mpichar ,0,mpicom)
            call mpibcast (tape(t)%column(g)%lat_name,max_chars   ,mpichar ,0,mpicom)
         end do
      end do
#endif
#endif
      return
   end subroutine fldlst

!#######################################################################

   subroutine inifld (t, f, avgflag, prec_acc, prec_wrt)
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the active list for a history tape
! 
! Method: Copy the data from the master field list to the active list for the tape
!         Also: define mapping arrays from (col,chunk) -> (lon,lat)
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
     use commap, only: latdeg,londeg
!
! Arguments
!
      integer, intent(in) :: t            ! history tape index
      integer, intent(in) :: f            ! field index from master field list

      character*1, intent(in) :: avgflag  ! averaging flag
      character*1, intent(in) :: prec_acc ! history buffer precision flag
      character*1, intent(in) :: prec_wrt ! history output precision flag
!
! Local workspace
!
      integer :: n                        ! field index on defined tape
      integer :: ff            ! column index
      integer :: lonind(pflds,2) !   beginning and ending longitude range
      integer :: latind(pflds,2) !   beginning and ending latitude range
      integer :: group !   number of group
      character(len=max_chars) :: lonlatname(pflds), lonname(pflds), latname(pflds) ! variable name
!
! Ensure that it is not to late to add a field to the history tape
!
      if (htapes_defined) then
         write(6,*)'INIFLD: Attempt to add field ',masterlist(f)%field%name,' after history files set'
         call endrun ()
      end if

      nflds(t) = nflds(t) + 1
      n = nflds(t)
!
! Copy field info.
!
      tape(t)%hlist(n)%field = masterlist(f)%field
!
! Set history buffer size and its output data type flags. Set them to
! the default values given by, respective, nhstpr(t) and ndens(t)
! if the input flags prec_acc and prec_wrt are blank; otherwise set to
! the specified values.
!
#ifdef CRAY
      tape(t)%hlist(n)%hbuf_prec = 8
      tape(t)%hlist(n)%hwrt_prec = 8

      if (prec_acc /= ' ' .and. prec_acc /= '8') then
         if (prec_acc == '4') then
            if (masterproc) then
               write(6,*) 'INIFLD: Requested change in history buffer size for ', &
                          tape(t)%hlist(n)%field%name, ' ignored'
            end if
         else
            write(6,*)'INIFLD: unknown prec_acc=',prec_acc
            call endrun
         end if
      end if

      if (prec_wrt /= ' ' .and. prec_wrt /= '8') then
         if (prec_wrt == '4') then
            if (masterproc) then
               write(6,*) 'INIFLD: Requested change in history output size for ', &
                          tape(t)%hlist(n)%field%name, ' ignored'
            end if
         else
            write(6,*)'INIFLD: unknown prec_wrt=',prec_wrt
            call endrun
         end if
      end if
#else
      select case (prec_acc)
      case (' ')
         tape(t)%hlist(n)%hbuf_prec = nhstpr(t)
      case ('4')
         if (r4 /= r8) then
            tape(t)%hlist(n)%hbuf_prec = 4
            if (masterproc) then
               write(6,*) 'INIFLD: History buffer for ', tape(t)%hlist(n)%field%name, &
                          ' is real*4'
            end if
         else       ! if kind values r4 and r8 are identical, ignore the request
            tape(t)%hlist(n)%hbuf_prec = 8
            if (masterproc) then
               write(6,*) 'INIFLD: Requested change in history output size for ', &
                           tape(t)%hlist(n)%field%name, ' ignored'
               write(6,*) '        because kind values r8 and r4 are identical'
            end if
         end if
      case ('8')
         tape(t)%hlist(n)%hbuf_prec = 8
         if (masterproc) then
            write(6,*) 'INIFLD: History buffer for ', tape(t)%hlist(n)%field%name, &
                       ' is real*8'
         end if
      case default
         write(6,*)'INIFLD: unknown prec_acc=',prec_acc
         call endrun
      end select

      select case (prec_wrt)
      case (' ')
         if (ndens(t) == 1) then
            tape(t)%hlist(n)%hwrt_prec = 8
         else
            tape(t)%hlist(n)%hwrt_prec = 4
         end if
      case ('4')
         tape(t)%hlist(n)%hwrt_prec = 4
         if (masterproc) then
            write(6,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
                       ' is real*4'
         end if
      case ('8')
         tape(t)%hlist(n)%hwrt_prec = 8
         if (masterproc) then
            write(6,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
                       ' is real*8'
         end if
      case default
         write(6,*)'INIFLD: unknown prec_wrt=',prec_wrt
         call endrun
      end select
#endif
!
! Override the default averaging (masterlist) averaging flag if non-blank
!
      if (avgflag == ' ') then
         tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
         tape(t)%hlist(n)%time_op = masterlist(f)%time_op(t)
      else
         tape(t)%hlist(n)%avgflag = avgflag
         select case (avgflag)
         case ('A')
            tape(t)%hlist(n)%time_op = 'mean'
         case ('I')
            tape(t)%hlist(n)%time_op = ' '
         case ('X')
            tape(t)%hlist(n)%time_op = 'maximum'
         case ('M')
            tape(t)%hlist(n)%time_op = 'minimum'
         case default
            write(6,*)'INIFLD: unknown avgflag=', avgflag
            call endrun
         end select
      end if
!
! Setup column information if this field will be written as group
! First verify the column information in the namelist
!
      ff = 1
      do while (fincllonlat(ff,t) /= ' ')
         lonlatname(ff) = trim(fincllonlat(ff,t))
         call getlatind(lonlatname(ff),latind(ff,1),latind(ff,2),latname(ff))
         call getlonind(lonlatname(ff),lonind(ff,1),lonind(ff,2),lonname(ff))
#ifdef HDEBUG
         write(6,*)'closest lon lat range for group ',lonlatname(ff),' is ',lonind(ff,:),latind(ff,:)
#endif
         ff = ff + 1
      end do

      group = ff-1

      ngroup(t)=group

      if (group .gt. 0) then 

         allocate (tape(t)%column(group))
         allocate (tape(t)%hlist(n)%field_column_name(group))

         do ff = 1, group
            tape(t)%column(ff)%lat_name=latname(ff)
            tape(t)%column(ff)%lon_name=lonname(ff)
            tape(t)%column(ff)%columnlon(:)=lonind(ff,:)
            tape(t)%column(ff)%columnlat(:)=latind(ff,:)
            tape(t)%column(ff)%num_lats = &
                 tape(t)%column(ff)%columnlat(2)-tape(t)%column(ff)%columnlat(1)+1
            tape(t)%column(ff)%num_lons = &
                 tape(t)%column(ff)%columnlon(2)-tape(t)%column(ff)%columnlon(1)+1
            tape(t)%hlist(n)%field_column_name(ff) = &
                 trim(tape(t)%hlist(n)%field%name) // "_" // trim(lonname(ff)) // "_" // trim(latname(ff))
        end do


      else

         allocate (tape(t)%column(1))
         allocate (tape(t)%hlist(n)%field_column_name(1))

         tape(t)%hlist(n)%field_column_name(1) = ' '
         tape(t)%column(1)%lat_name=' '
         tape(t)%column(1)%lon_name=' '
         tape(t)%column(1)%columnlon(:)=0
         tape(t)%column(1)%columnlat(:)=0
         tape(t)%column(1)%num_lats=0
         tape(t)%column(1)%num_lons=0

      end if

#ifdef HDEBUG
      write(6,*)'INIFLD: field ', tape(t)%hlist(n)%field%name, ' added as ', 'field number ', n,' on tape ', t
      write(6,*)'units=',tape(t)%hlist(n)%field%units
      write(6,*)'numlev=',tape(t)%hlist(n)%field%numlev
      write(6,*)'avgflag=',tape(t)%hlist(n)%avgflag
      write(6,*)'time_op=',tape(t)%hlist(n)%time_op
      write(6,*)'hbuf_prec=',tape(t)%hlist(n)%hbuf_prec
      write(6,*)'hwrt_prec=',tape(t)%hlist(n)%hwrt_prec
#endif

      return
   end subroutine inifld

!#######################################################################

   character(len=8) function getname (inname)
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve name portion of inname
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          lop it off
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*) inname
!
! Local workspace
!
      integer :: length
      integer :: i
   
      length = len (inname)
      
      if (length < 8 .or. length > 10) then
         write(6,*) 'getname: bad length=',length
         call endrun
      end if
   
      getname = ' '
      do i=1,8
         if (inname(i:i) == ':') exit
         getname(i:i) = inname(i:i)
      end do
      
      return
   end function getname

!#######################################################################

   subroutine getlatind(inname,beglatind,endlatind,latname)
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve closes model lat index given north south latitude string
!          
! Author: John Truesdale
! 
!-------------------------------------------------------------------------------
     use commap, only: latdeg
!
! Arguments
!
      character(len=*) inname,latname
      integer :: beglatind,endlatind
!
! Local workspace
!
      character(len=max_chars) str1,str2,tmpstr(4)
      integer :: i,j,marker,latind(2),tmplen(4)
      real(r8) :: latdegree(2)

!
!-------------------------------------------------------------------------------
!
      str1=' '
      str2=' '
      tmpstr(:)= ' '

!
! make sure _ separator is present
!      
      if (scan(inname,'_').eq.0) then
         write(6,*)'HISTORY: Improperly formatted column string.  Missing underscore character (xxxE_yyyS) ',inname,scan(inname,'_')
         call endrun
      end if
!
! split inname string into lat lon strings
!      
      marker=scan(inname,'_')
      str1=inname(:marker-1)
      str2=inname(marker+1:)

!
! split ranges of lats( or lons) into seperate substrings. Substrings 1,2 will contain lats(lons) from portion of string before underscore
! if a single column and not a range is specified substrings 1 and 2 will be the same
! and substrings 3,4 will contain lats(lons) from portion of main string after underscore character
! if a single column and not a range is specified substrings 3 and 4 will be the same
!
      if (scan(str1,':') .ne. 0) then
         marker=scan(str1,':')
         tmpstr(1)=str1(:marker-1)
         tmpstr(2)=str1(marker+1:)
      else
         tmpstr(1)=str1
         tmpstr(2)=str1
      end if
      if (scan(str2,':') .ne. 0) then
         marker=scan(str2,':')
         tmpstr(3)=str2(:marker-1)
         tmpstr(4)=str2(marker+1:)
      else
         tmpstr(3)=str2
         tmpstr(4)=str2
      end if

! check format of substrings - (number followed by single character north/south east/west designation)

      do i = 1,4
         tmplen(i)=len_trim(tmpstr(i))
         if (verify(tmpstr(i),"0123456789.").ne.tmplen(i) .or. verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS") .ne. 0) then
            write(6,*)'HISTORY1: Improperly formatted column string. ',&
                 inname,verify(tmpstr(i),"0123456789."),tmplen(i),&
                 verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS")
            call endrun
         end if
      end do

! find latitude substrings and put into temporary work space tmpstr(1:2)

      if (verify(tmpstr(1)(tmplen(1):tmplen(1)),"nNsSs").eq.0  .and. verify(tmpstr(2)(tmplen(2):tmplen(2)),"nNsSs").eq.0 ) then

      else if (verify(tmpstr(3)(tmplen(3):tmplen(3)),"nNsSs").eq.0 .and. verify(tmpstr(3)(tmplen(3):tmplen(3)),"nNsSs").eq.0) then
         tmpstr(1) = tmpstr(3)
         tmplen(1)=tmplen(3)
         tmpstr(2) = tmpstr(4)
         tmplen(2)=tmplen(4)
      else
         write(6,*)'HISTORY3: Improperly formatted column string. ',inname
         call endrun
      end if
!
! convert lat substrings to real and if in southern hemisphere make negative
!
      do i = 1,2
         read(tmpstr(i)(1:tmplen(i)-1),*) latdegree(i)
         if (verify(tmpstr(i)(tmplen(i):tmplen(i)),'sS').eq.0) then
            latdegree(i)=(-1.)*latdegree(i)
         end if
!
! Make sure specified latitudes is in bounds
!
         if (latdegree(i) .lt. -90. .or. latdegree(i) .gt. 90.) then
            write(6,*)'HISTORY: latitude for column namelist is out of range (-90 .. 90) value=',latdegree(i)
            call endrun
         endif
!
! Find closest lat index for each substring
!
         if (latdegree(i).ge.latdeg(plat)) then
            latind(i)=plat
         else
            do j = 1, plat-1
               if ( abs(latdeg(j)-latdegree(i)) .lt.  &
                    abs(latdeg(j+1)-latdegree(i))) then
                  latind(i)=j
                  exit
               endif
            enddo
         end if
      end do
!
! output begining and ending latitude indicies.   If just a column is specified then beginning latitude index will be the same as the
! ending latitude index  

!      
      if (latind(1) .le. latind(2) ) then
         beglatind =latind(1)
         endlatind =latind(2)
      else 
         beglatind = latind(2)
         endlatind = latind(1)
      end if

      if (beglatind .eq. endlatind) then
         latname = trim(tmpstr(1))
      else
         if (latind(1) .le. latind(2) ) then
            latname = trim(tmpstr(1)) // "_to_" // trim(tmpstr(2))
         else 
            latname = trim(tmpstr(2)) // "_to_" // trim(tmpstr(1))
         end if
      end if
         return
   end subroutine getlatind
!#######################################################################

   subroutine getlonind (inname,beglonind,endlonind,lonname)
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve closes model lat index given north south latitude string
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          lop it off
! 
!-------------------------------------------------------------------------------
     use commap, only: londeg
!
! Arguments
!
      character(len=*) inname,lonname
      integer :: beglonind,endlonind
!
! Local workspace
!
      character(len=max_chars) str1,str2,tmpstr(4)
      integer :: i,j,marker,lonind(2),tmplen(4)
      real(r8) :: londegree(2)

!
!-------------------------------------------------------------------------------
!
      str1=' '
      str2=' '
      tmpstr(:)= ' '
!
! make sure _ separator is present
!      
      if (scan(inname,'_').eq.0) then
         write(6,*)'HISTORY4: Improperly formatted column string.  Missing underscore character (xxxE_yyyS) ',&
              inname,scan(inname,'_')
         call endrun
      end if
!
! split string in to lat lon strings
!      
      marker=scan(inname,'_')
      str1=inname(:marker-1)
      str2=inname(marker+1:)

!
! split ranges of lats( or lons) into seperate substrings. Substrings 1,2 will contain lats(lons) from portion of string before underscore
! if a single column and not a range is specified substrings 1 and 2 will be the same
! and substrings 3,4 will contain lats(lons) from portion of main string after underscore character
! if a single column and not a range is specified substrings 3 and 4 will be the same
!
      if (scan(str1,':') .ne. 0) then
         marker=scan(str1,':')
         tmpstr(1)=str1(:marker-1)
         tmpstr(2)=str1(marker+1:)
      else
         tmpstr(1)=str1
         tmpstr(2)=str1
      end if

      if (scan(str2,':') .ne. 0) then
         marker=scan(str2,':')
         tmpstr(3)=str2(:marker-1)
         tmpstr(4)=str2(marker+1:)
      else
         tmpstr(3)=str2
         tmpstr(4)=str2
      end if

! check format of substrings - (number followed by single character north/south east/west designation)

      do i = 1,4
         tmplen(i)=len_trim(tmpstr(i))
         if (verify(tmpstr(i),"0123456789.").ne.tmplen(i) .or. verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS") .ne. 0) then
            write(6,*)'HISTORY1: Improperly formatted column string. ', &
                 inname,verify(tmpstr(i),"0123456789."),tmplen(i), &
                 verify(tmpstr(i)(tmplen(i):tmplen(i)),"eEwWnNsS")
            call endrun
         end if
      end do

! find longitude substrings and put into temporary work space tmpstr(1:2)

      if (verify(tmpstr(1)(tmplen(1):tmplen(1)),"eEwW").eq.0  .and. verify(tmpstr(2)(tmplen(2):tmplen(2)),"eEwW").eq.0 ) then

      else if (verify(tmpstr(3)(tmplen(3):tmplen(3)),"eEwW").eq.0 .and. verify(tmpstr(3)(tmplen(3):tmplen(3)),"eEwW").eq.0) then
         tmpstr(1) = tmpstr(3)
         tmplen(1)=tmplen(3)
         tmpstr(2) = tmpstr(4)
         tmplen(2)=tmplen(4)
      else
         write(6,*)'HISTORY3: Improperly formatted column string. ',inname
         call endrun
      end if


!
! convert lon substrings to real and make sure its degrees east
!
      do i = 1,2
         read(tmpstr(i)(1:tmplen(i)-1),*) londegree(i)
         if (verify(tmpstr(i)(tmplen(i):tmplen(i)),'wW').eq.0) then
            londegree(i) = 360. - londegree(i)
         end if
!
! Make sure specified longitudes are in bounds
!
         if (londegree(i) .lt. 0 .or. londegree(i) .gt. 360.) then
            write(6,*)'HISTORY: longitude for column namelist is out of range (0 .. 360) value=',londegree(i)
            call endrun
         endif

!
! Find closest lon index for each substring.  If just a column is specified then beginning longitude index will be the same as the
! ending longitude index  
!
         if (londegree(i).ge.londeg(plon,1)) then
            lonind(i)=plon
         else
            do j = 1, plon-1
               if ( abs(londeg(j,1)-londegree(i)) .le. abs(londeg(j+1,1)-londegree(i))) then
                  lonind(i)=j
                  exit
               endif
            end do
         end if
      end do
!
! output begining and ending longitude indicies.   If just a column is specified then beginning longitude index will be the same as the
! ending longitude index  

!      
      if (lonind(1) .le. lonind(2) ) then
         beglonind =lonind(1)
         endlonind =lonind(2)
      else 
         beglonind = lonind(2)
         endlonind = lonind(1)
      end if


      if (beglonind .eq. endlonind) then
         lonname = trim(tmpstr(1))
      else
         if (lonind(1) .le. lonind(2) ) then
            lonname = trim(tmpstr(1)) // "_to_" // trim(tmpstr(2))
         else 
            lonname = trim(tmpstr(2)) // "_to_" // trim(tmpstr(1))
         end if
      end if


      return
    end subroutine getlonind

!#######################################################################

   character(len=1) function getflag (inname)
!----------------------------------------------------------------------- 
! 
! Purpose: retrieve flag portion of inname
!          
! Method:  If an averaging flag separater character is present (":") in inname, 
!          return the character after it as the flag
! 
!-------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*) inname   ! character string
!
! Local workspace
!
      integer :: length         ! length of inname
      integer :: i              ! loop index

      length = len (inname)

      if (length /= 10) then
         write(6,*) 'getflag: bad length=',length
         call endrun
      end if

      getflag = ' '
      do i=1,9
         if (inname(i:i) == ':') then
            getflag = inname(i+1:i+1)
            exit
         end if
      end do

      return
   end function getflag

!#######################################################################

   subroutine list_index (list, name, index)
!
! Input arguments
!
      character(len=*) :: list(pflds) ! input list of names, possibly ":" delimited
      character(len=8) :: name        ! name to be searched for
!
! Output arguments
!
      integer :: index                ! index of "name" in "list"
!
! Local workspace
!
      character(len=8) :: listname    ! input name with ":" stripped off.
      integer f                       ! field index

      index = 0
      do f=1,pflds
!
! Only list items
!
         listname = getname (list(f))
         if (listname == ' ') exit
         if (listname == name) then
            index = f
            exit
         end if
      end do
      
      return
   end subroutine list_index

!#######################################################################

   subroutine outfld (fname, field, idim, c)
!-----------------------------------------------------------------------
!
! Purpose: Accumulate (or take min, max, etc. as appropriate) input field
!          into its history buffer for appropriate tapes
!
! Method: Search for fname among fields on history tapes.  If found, do the
!         accumulation.  If not found, return silently.
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

      integer, intent(in) :: idim           ! Longitude dimension of field array
      integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

      real(r8), intent(in) :: field(idim,*) ! Array containing field values
!
! Local variables
!
      integer :: t, f                ! tape, field indices
      integer coldimin               ! column dimension of model array
      integer :: fl, fu              ! upper, lower indices used in binary search thru sorted list
      integer :: begver              ! on-node vert start index
      integer :: endver              ! on-node vert end index
      integer :: lexdiff             ! lexical difference (analagous to output from strcmp in C)
      integer :: endi                ! ending longitude index (reduced grid)

      character*8 :: fname8          ! 8-char equivalent of fname
      character*1 :: avgflag         ! averaging flag

      type (hbuffer_2d) :: hbuf      ! history buffer
      integer, pointer :: nacs(:)    ! accumulation counter
      type (dim_index_2d) :: dimind  ! 2-D dimension index
      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
!-----------------------------------------------------------------------
!      call t_startf ('outfld')
      fname8 = fname
!
! Note, the field may be on any or all of the history files (primary
! and auxiliary).
!
!      write(6,*)'fname8=',fname8
      do 40 t=1,mtapes
!
! Search the sorted list of fields.  Algorithm taken from Numerical
! Recipes in C
!
         fl = 0
         fu = nflds(t) + 1
         
         lexdiff = -1
         do while (fu - fl .gt. 1)
            f = (fl + fu)/2
            lexdiff = strcmpf (fname8, tape(t)%hlist(f)%field%name)
            if (lexdiff < 0) then
               fu = f
            else if (lexdiff > 0) then
               fl = f
            else
               exit
            end if
         end do

         if (lexdiff /= 0) then
!          write(6,*)'OUTFLD: field ',fname,' not found on tape ', t
            goto 40
         end if
!
! Field name found.  Update history buffer
!
         begver  = tape(t)%hlist(f)%field%begver
         endver  = tape(t)%hlist(f)%field%endver
         endi    = tape(t)%hlist(f)%field%colperdim3(c)
         avgflag = tape(t)%hlist(f)%avgflag
         flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill
         coldimin= tape(t)%hlist(f)%field%coldimin
         nacs   => tape(t)%hlist(f)%nacs(:coldimin,c)
         call assoc_hbuf2d_with_hbuf3d (hbuf, tape(t)%hlist(f)%hbuf, c)

         dimind = dim_index_2d (1,endi,begver,endver)

         select case (avgflag)

         case ('I') ! Instantaneous

            call hbuf_accum_inst (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('A') ! Time average

            call hbuf_accum_add (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('X') ! Maximum over time

            call hbuf_accum_max (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('M') ! Minimum over time

            call hbuf_accum_min (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case default

            write(6,*)'OUTFLD: invalid avgflag=', avgflag
            call endrun ()

         end select
40    continue
!      call t_stopf ('outfld')
      return
   end subroutine outfld

   subroutine outfldcol (fname, field, idim, icol, c)
!-----------------------------------------------------------------------
!
! Purpose: Accumulate (or take min, max, etc. as appropriate) single column input field
!          into its history buffer for appropriate tapes
!    doesn't work right
!
! Method: Search for fname among fields on history tapes.  If found, do the
!         accumulation.  If not found, return silently.
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

      integer, intent(in) :: idim           ! Longitude dimension of field array
      integer, intent(in) :: icol           ! column index of field array
      integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

      real(r8), intent(in) :: field(idim,*) ! Array containing field values
!
! Local variables
!
      integer :: t, f                ! tape, field indices
      integer coldimin               ! column dimension of model array
      integer :: fl, fu              ! upper, lower indices used in binary search thru sorted list
      integer :: begver              ! on-node vert start index
      integer :: endver              ! on-node vert end index
      integer :: lexdiff             ! lexical difference (analagous to output from strcmp in C)
      integer :: begi                ! ending longitude index (reduced grid)
      integer :: endi                ! ending longitude index (reduced grid)
      integer k

      character*8 :: fname8          ! 8-char equivalent of fname
      character*1 :: avgflag         ! averaging flag

      type (hbuffer_2d) :: hbuf      ! history buffer
      integer, pointer :: nacs(:)    ! accumulation counter
      type (dim_index_2d) :: dimind  ! 2-D dimension index
      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
!-----------------------------------------------------------------------
!      call t_startf ('outfld')
      fname8 = fname
!
! Note, the field may be on any or all of the history files (primary
! and auxiliary).
!
!      write(6,*)'fname8=',fname8
      do 40 t=1,mtapes
!
! Search the sorted list of fields.  Algorithm taken from Numerical
! Recipes in C
!
         fl = 0
         fu = nflds(t) + 1
         
         lexdiff = -1
         do while (fu - fl .gt. 1)
            f = (fl + fu)/2
            lexdiff = strcmpf (fname8, tape(t)%hlist(f)%field%name)
            if (lexdiff < 0) then
               fu = f
            else if (lexdiff > 0) then
               fl = f
            else
               exit
            end if
         end do

         if (lexdiff /= 0) then
!          write(6,*)'OUTFLD: field ',fname,' not found on tape ', t
            goto 40
         end if
!
! Field name found.  Update history buffer
!
         begver  = tape(t)%hlist(f)%field%begver
         endver  = tape(t)%hlist(f)%field%endver
	 begi=icol
         endi    = icol
         avgflag = tape(t)%hlist(f)%avgflag
         flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill
         coldimin= tape(t)%hlist(f)%field%coldimin
         nacs   => tape(t)%hlist(f)%nacs(:coldimin,c)
         call assoc_hbuf2d_with_hbuf3d (hbuf, tape(t)%hlist(f)%hbuf, c)

         dimind = dim_index_2d (begi,endi,begver,endver)

         select case (avgflag)

         case ('I') ! Instantaneous

            call hbuf_accum_inst (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('A') ! Time average

            call hbuf_accum_add (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('X') ! Maximum over time

            call hbuf_accum_max (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case ('M') ! Minimum over time

            call hbuf_accum_min (hbuf%buf4, hbuf%buf8, field, nacs, dimind, idim, flag_xyfill)

         case default

            write(6,*)'OUTFLD: invalid avgflag=', avgflag
            call endrun ()

         end select
40    continue
!      call t_stopf ('outfld')
      return
   end subroutine outfldcol

!#######################################################################

   integer function strcmpf (name1, name2)
!----------------------------------------------------------------------- 
! 
! Purpose: Return the lexical difference between two strings
! 
! Method: Use ichar() intrinsic as we loop through the names
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=8), intent(in) :: name1, name2  ! strings to compare
      integer n                                     ! loop index
   
      do n=1,8
         strcmpf = ichar(name1(n:n)) - ichar(name2(n:n))
         if (strcmpf /= 0) exit
      end do

      return
   end function strcmpf

!#######################################################################

   subroutine h_inquire (t)
!----------------------------------------------------------------------- 
! 
! Purpose: Ensure that the proper variables are on a history file
! 
! Method: Issue the appropriate netcdf wrapper calls
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: t   ! tape index
!
! Local workspace
!
      integer f,ff               ! field index
      integer ret                ! return value from function call
      integer londim             ! longitude dimension id
      integer latdim             ! latitude dimension id
      integer levdim             ! level dimension id
      integer ilevdim            ! intfc dimension id
      integer tbnddim            ! time_bnds dimension id
      integer old_mode           ! returned from nf_set_fill
      integer marker             ! index for string marker
!
! Setup netcdf file - create the dimensions of lat,lon,time,level
!     
      ret = nf_set_fill (nfid(t), nf_nofill, old_mode)
!
! Dimension id's
!
      call wrap_inq_dimid (nfid(t), 'lat', latdim)
      call wrap_inq_dimid (nfid(t), 'lon', londim)
      call wrap_inq_dimid (nfid(t), 'lev', levdim)
      call wrap_inq_dimid (nfid(t), 'ilev', ilevdim)
      call wrap_inq_dimid (nfid(t), 'tbnd', tbnddim)
!
! Create variables for model timing and header information 
!
      call wrap_inq_varid (nfid(t),'ndcur   ',    ndcurid(t))
      call wrap_inq_varid (nfid(t),'nscur   ',    nscurid(t))
      call wrap_inq_varid (nfid(t),'date    ',    dateid(t))
      call wrap_inq_varid (nfid(t),'datesec ',    datesecid(t))
      call wrap_inq_varid (nfid(t),'nsteph  ',    nstephid(t))
      call wrap_inq_varid (nfid(t),'time    ',    timeid(t))
      call wrap_inq_varid (nfid(t),'time_bnds',   tbndid(t))
      call wrap_inq_varid (nfid(t),'date_written',date_writtenid(t))
      call wrap_inq_varid (nfid(t),'time_written',time_writtenid(t))
!
! Obtain variable name from ID which was read from restart file
!
      do f=1,nflds(t)
!
! If this field will be put out as columns then get column names for field
!
         if (ngroup(t) .gt. 0) then
            do i=1,ngroup(t)
               ff=(f-1)*ngroup(t)+i
               call wrap_inq_varname (nfid(t), varid(ff,t), tape(t)%hlist(f)%field_column_name(i))
               marker=scan(tape(t)%hlist(f)%field_column_name(i),'_')
               tape(t)%hlist(f)%field%name=tape(t)%hlist(f)%field_column_name(i)(:marker-1)
            end do
         else	
            call wrap_inq_varname (nfid(t), varid(f,t), tape(t)%hlist(f)%field%name)
         end if
      end do
!
      ret = nf_enddef(nfid(t))
      write(6,*)'H_INQUIRE: Successfully opened netcdf file '

      return
   end subroutine h_inquire

!#######################################################################

   subroutine h_default ()
!----------------------------------------------------------------------- 
! 
! Purpose: Define default contents of history files
! 
! Method: Call add_default for each field.  Arguments are field name, tape index, and
!         modification to averaging flag (blank means no modification)
! 
!-----------------------------------------------------------------------
      use dycore, only: dycore_is

#include <comctl.h>

!
! Local workspace
!
      integer m      ! tracer index
!
! First tape (monthly output by default)
!
      call add_default ('PHIS    ', 1, ' ')
      call add_default ('PS      ', 1, ' ')
      call add_default ('T       ', 1, ' ')
      call add_default ('Q       ', 1, ' ')
      call add_default ('U       ', 1, ' ')
      call add_default ('V       ', 1, ' ')

#if ( defined STAGGERED )
      call add_default ('US      ', 1, ' ')
      call add_default ('VS      ', 1, ' ')
#endif

#if ( defined COUP_CSM )
      call add_default ('CPLRAINC', 1, ' ')
      call add_default ('CPLRAINL', 1, ' ')
      call add_default ('CPLSNOWC', 1, ' ')
      call add_default ('CPLSNOWL', 1, ' ')
      call add_default ('CPLPRCER', 1, ' ')
#endif

      call add_default (dcconnam(1), 1, ' ')

      if ( .not. dycore_is('LR') )then
         call add_default ('DTH     ', 1, ' ')
      end if
      call add_default ('SNOWHICE', 1, ' ')
      call add_default ('SNOWHLND', 1, ' ')
      call add_default ('PRECL   ', 1, ' ')
      call add_default ('PRECC   ', 1, ' ')
      call add_default ('PRECSH  ', 1, ' ')
      call add_default ('PRECSL  ', 1, ' ')
      call add_default ('PRECSC  ', 1, ' ')
      call add_default ('SHFLX   ', 1, ' ')
      call add_default ('LHFLX   ', 1, ' ')
      call add_default ('QFLX    ', 1, ' ')
      call add_default ('PBLH    ', 1, ' ')
      call add_default ('TREFHT  ', 1, ' ')
#if ( defined COUP_CSM )
      call add_default ('QREFHT  ', 1, ' ')
#endif

      call add_default ('DTV     ', 1, ' ')
      call add_default ('FSNS    ', 1, ' ')
      call add_default ('FLNS    ', 1, ' ')
      call add_default ('FLNT    ', 1, ' ')
      call add_default ('FLUT    ', 1, ' ')
      call add_default ('FLUTC   ', 1, ' ')
      call add_default ('FSNT    ', 1, ' ')
      call add_default ('FSNTOA  ', 1, ' ')
      call add_default ('FSNTOAC ', 1, ' ')
      call add_default ('FLWDS   ', 1, ' ')
      call add_default ('FLWDSC  ', 1, ' ')
      call add_default ('FSUTOA  ', 1, ' ')
      call add_default ('FSUTOAC ', 1, ' ')
      call add_default ('LWCF    ', 1, ' ')
      call add_default ('SWCF    ', 1, ' ')
      call add_default ('CLOUD   ', 1, ' ')
      call add_default ('CLDTOT  ', 1, ' ')
      call add_default ('CLDLOW  ', 1, ' ')
      call add_default ('CLDMED  ', 1, ' ')
      call add_default ('CLDHGH  ', 1, ' ')
      call add_default ('FLNTC   ', 1, ' ')
      call add_default ('FSNTC   ', 1, ' ')
      call add_default ('FLNSC   ', 1, ' ')
      call add_default ('FSNSC   ', 1, ' ')
      call add_default ('FSDSC   ', 1, ' ')
      call add_default ('OMEGA   ', 1, ' ')
      call add_default ('TAUX    ', 1, ' ')
      call add_default ('TAUY    ', 1, ' ')
      call add_default ('SRFRAD  ', 1, ' ')
      call add_default ('QRS     ', 1, ' ')
      call add_default ('QRL     ', 1, ' ')
#ifdef QRLDAMP
      call add_default ('DQRL     ', 1, ' ')
#endif
#ifdef FLUXDAMP
      call add_default ('DUBOT    ', 1, ' ')
      call add_default ('DVBOT    ', 1, ' ')
      call add_default ('DWINDBOT ', 1, ' ')
      call add_default ('DLHFLX   ', 1, ' ')
      call add_default ('DSHFLX   ', 1, ' ')
#endif
      call add_default ('TS      ', 1, ' ')
      call add_default ('TSMN    ', 1, ' ')
      call add_default ('TSMX    ', 1, ' ')
      call add_default ('SOLIN   ', 1, ' ')
      call add_default ('DTCOND  ', 1, ' ')
      call add_default ('TTEND   ', 1, ' ')
!      call add_default ('STEND   ', 1, ' ')
      call add_default ('FSDS    ', 1, ' ')
      call add_default ('RELHUM  ', 1, ' ')
      call add_default ('TPHYSTND', 1, ' ')
      call add_default ('UPHYSTND', 1, ' ')
      call add_default ('VPHYSTND', 1, ' ')
#ifndef CRM
      do m=1,pcnst+pnats
         call add_default (cnst_name(m), 1, ' ')
         call add_default (tendnam(m), 1, ' ')
         call add_default (qphystendnam(m), 1, ' ')
      end do
      call add_default ('CMFDT   ', 1, ' ')
      call add_default ('CMFDQ   ', 1, ' ')
      call add_default ('CMFDQR  ', 1, ' ')
      call add_default ('QC      ', 1, ' ')
      call add_default ('CMFMC   ', 1, ' ')
      call add_default ('CONCLD  ', 1, ' ')
#endif
#ifdef CRM
      call add_default ('PRES    ', 1, ' ')
      call add_default ('DPRES   ', 1, ' ')
      call add_default ('HEIGHT  ', 1, ' ')
      call add_default ('_QRS    ', 1, ' ')
      call add_default ('_QRL    ', 1, ' ')
      call add_default ('_FSDS   ', 1, ' ')
      call add_default ('_FSNS   ', 1, ' ')
      call add_default ('_FLNS   ', 1, ' ')
      call add_default ('_FLUT   ', 1, ' ')
      call add_default ('_FSNTOA ', 1, ' ')
      call add_default ('_LWCF   ', 1, ' ')
      call add_default ('_SWCF   ', 1, ' ')
      call add_default ('SPDT    ', 1, ' ')
      call add_default ('SPDQ    ', 1, ' ')
      call add_default ('SPDQC   ', 1, ' ')
      call add_default ('SPDQI   ', 1, ' ')
      call add_default ('SPMC    ', 1, ' ')
      call add_default ('SPMCUP  ', 1, ' ')
      call add_default ('SPMCDN  ', 1, ' ')
      call add_default ('SPMCUUP ', 1, ' ')
      call add_default ('SPMCUDN ', 1, ' ')
      call add_default ('SPQC    ', 1, ' ')
      call add_default ('SPQI    ', 1, ' ')
      call add_default ('SPQS    ', 1, ' ')
      call add_default ('SPQG    ', 1, ' ')
      call add_default ('SPQR    ', 1, ' ')
      call add_default ('SPQTFLX ', 1, ' ')
      call add_default ('SPQTFLXS', 1, ' ')
      call add_default ('TKE     ', 1, ' ')
      call add_default ('TKES    ', 1, ' ')
      call add_default ('SPQPFLX ', 1, ' ')
      call add_default ('SPPFLX  ', 1, ' ')
      call add_default ('SPQTLS  ', 1, ' ')
      call add_default ('SPQTTR  ', 1, ' ')
      call add_default ('SPQPTR  ', 1, ' ')
      call add_default ('SPQPEVP ', 1, ' ')
      call add_default ('SPQPFALL', 1, ' ')
      call add_default ('SPQPSRC ', 1, ' ')
      call add_default ('SPTLS   ', 1, ' ')
      call add_default ('SPULS   ', 1, ' ')
      call add_default ('SPVLS   ', 1, ' ')
      call add_default ('CLOUDTOP', 1, ' ')
      call add_default ('CLOUDR  ', 1, ' ')
      call add_default ('_CONCLD ', 1, ' ')
      call add_default ('_CLOUD  ', 1, ' ')
      call add_default ('_CLDTOT ', 1, ' ')
      call add_default ('_CLDLOW ', 1, ' ')
      call add_default ('_CLDMED ', 1, ' ')
      call add_default ('_CLDHGH ', 1, ' ')
     if (doisccp) then
         call add_default ('_FISCCP1', 1, ' ')
         call add_default ('_TCLDAR ', 1, ' ')
         call add_default ('_MNPTOP ', 1, ' ')
         call add_default ('_MNTAU  ', 1, ' ')
         call add_default ('_MNTTOP ', 1, ' ')
         call add_default ('_CLOUDY ', 1, ' ')
      endif
      call add_default ('TIMINGF ', 1, ' ')
#endif

!++scyc  add fields associated with sulfur cycle

      if ( doramp_so4 ) then
         call add_default ('SULFBIO ', 1, ' ')
         call add_default ('SULFANT ', 1, ' ')
         call add_default ('SULFMMR ', 1, ' ')
      end if

      if ( indirect ) then
         call add_default ('MSO4    ', 1, ' ')
         call add_default ('LWC     ', 1, ' ')
         call add_default ('CLDFRQ  ', 1, ' ')
         call add_default ('WREL    ', 1, ' ')
         call add_default ('WLWC    ', 1, ' ')
      end if

      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('UU      ', 1, ' ')
      call add_default ('Z3      ', 1, ' ')
      call add_default ('TMQ     ', 1, ' ')
      call add_default ('PSL     ', 1, ' ')

      call add_default ('ICEFRAC ', 1, ' ')
      call add_default ('OCNFRAC ', 1, ' ')
      call add_default ('LANDFRAC ', 1, ' ')

#if ( defined COUP_SOM )
      call add_default ('MELTB   ', 1, ' ')
      call add_default ('MELTT   ', 1, ' ')
      call add_default ('MELTL   ', 1, ' ')
      call add_default ('GROWB   ', 1, ' ')
      call add_default ('FRAZIL  ', 1, ' ')
      call add_default ('FLOOD   ', 1, ' ')
      call add_default ('FRZMLT  ', 1, ' ')
      call add_default ('QFLUX   ', 1, ' ')
      call add_default ('QFLUX_FT', 1, ' ')
      call add_default ('QFLUX_TH', 1, ' ')
      call add_default ('QFLUX_A2', 1, ' ')
      call add_default ('FOCN    ', 1, ' ')
      call add_default ('EICEIN  ', 1, ' ')
      call add_default ('EICEOUT ', 1, ' ')
      call add_default ('F_ICE   ', 1, ' ')
      call add_default ('F_OCN   ', 1, ' ')
      call add_default ('FRZMLTMX', 1, ' ')
      call add_default ('DELTAICE', 1, ' ')
      call add_default ('IMBAL   ', 1, ' ')
      call add_default ('NRGERROR', 1, ' ')
      call add_default ('MLDANN  ', 1, ' ')
      call add_default ('ONF     ', 1, ' ')
      call add_default ('OIE     ', 1, ' ')
      call add_default ('OIERATE ', 1, ' ')
      call add_default ('NRGICE  ', 1, ' ')
      call add_default ('IIERATE ', 1, ' ')
#endif
      call add_default ('FSNSOI', 1, ' ')
      call add_default ('FLNSOI', 1, ' ')
      call add_default ('LHFLXOI', 1, ' ')
      call add_default ('SHFLXOI', 1, ' ')

! ISCCP fields

      if (doisccp) then
         call add_default ('FISCCP1 ', 1, ' ')
         call add_default ('TCLDAREA', 1, ' ')
         call add_default ('LCLDAREA', 1, ' ')
         call add_default ('MCLDAREA', 1, ' ')
         call add_default ('HCLDAREA', 1, ' ')
         call add_default ('MEANPTOP', 1, ' ')
         call add_default ('MEANTAU ', 1, ' ')
         call add_default ('MEANTTOP', 1, ' ')
         call add_default ('CLOUDY  ', 1, ' ')
#ifdef FILLDEBUG
         call add_default ('FAKEFLD ', 1, ' ')
         call add_default ('FAKEFLD1', 1, ' ')
         call add_default ('FAKEFLD2', 1, ' ')
#endif
      end if


! Mike Pritchard - for closing budgets:

! Snapshots for explicitly capturing the tendency term (to know residual)
#ifdef PVBUDGET
      call add_default ('PVAP    ',1,' ')
     call add_default ('PVPHYSTN',1,' ')
     call add_default ('PVTEND  ',1,' ')
#endif
      call add_default ('TAP     ',1,' ')
      call add_default ('QAP     ',1,' ')
      call add_default ('QCAP    ',1,' ')
      call add_default ('QIAP    ',1,' ')
      call add_default ('UAP     ',1,' ')
      call add_default ('VAP     ',1,' ')

! Information about tendency contirubtions from physics packages.
     call add_default ('TPHYSTND',1,' ')
     call add_default ('UPHYSTND',1,' ')
     call add_default ('VPHYSTND',1,' ')

     call add_default ('PHQ     ',1,' ')
     call add_default ('PHCLDLIQ',1,' ')
     call add_default ('PHCLDICE',1,' ')
     call add_default ('UTGWORO ',1,' ')
     call add_default ('VTGWORO ',1,' ')
     call add_default ('TTGWORO ',1,' ')
     call add_default ('DTVKE   ',1,' ')
     call add_default ('DTV     ',1,' ')
     call add_default ('DUV     ',1,' ')
     call add_default ('DVV     ',1,' ')
    ! Nb. SPDT, SPDQ, etc. are all already defaults   

    ! Pritch variables capturing dycore tendency:
     call add_default ('TTEND   ',1,' ')
     call add_default ('QTEND   ',1,' ')
     call add_default ('QCTEND  ',1,' ')
     call add_default ('QITEND  ',1,' ')
     call add_default ('UTEND   ',1,' ')
     call add_default ('VTEND   ',1,' ')
#ifdef QVORTDAMP
     call add_default ('XXXU3   ',1,' ')
     call add_default ('XXXU3AUX',1,' ')
     call add_default ('DQFCST',1,' ')
#endif
! These fields removed 10/30/2003 per CRB decision.
!      call add_default ('CLDST   ', 1, ' ')
!      call add_default ('CNVCLD  ', 1, ' ')
!      call add_default ('EVAPPCT ', 1, ' ')
!      call add_default ('MQ      ', 1, ' ')
!      call add_default ('O3VMR   ', 1, ' ')
!      call add_default ('OMEGAU  ', 1, ' ')
!      call add_default ('PRECTMX ', 1, ' ')
!      call add_default ('QPERT   ', 1, ' ')
!      call add_default ('SICTHK  ', 1, ' ')
!      call add_default ('TPERT   ', 1, ' ')
!      call add_default ('TREFMNAV', 1, ' ')
!      call add_default ('TREFMXAV', 1, ' ')
!      call add_default ('TSICE   ', 1, ' ')
!      call add_default ('TT      ', 1, ' ')
!      call add_default ('USTAR   ', 1, ' ')
!      call add_default ('VZ      ', 1, ' ')
!      call add_default ('WSPEED  ', 1, ' ')
!      call add_default ('ZMDQ    ', 1, ' ')
!      call add_default ('ZMDT    ', 1, ' ')
!      call add_default ('ZZ      ', 1, ' ')
!      call add_default ('FLNT    ', 2, ' ')
!      call add_default ('OMEGA500', 2, ' ')
!      call add_default ('OMEGA850', 2, ' ')
!      call add_default ('PRECT   ', 2, ' ')
!      call add_default ('PSL     ', 2, 'I')
!      call add_default ('QFLX    ', 2, ' ')
!      call add_default ('QREFHT  ', 2, ' ')
!      call add_default ('T300    ', 2, 'I')
!      call add_default ('T850    ', 2, 'I')
!      call add_default ('TMQ     ', 2, 'I')
!      call add_default ('TREFHT  ', 2, ' ')
!      call add_default ('TREFHTMN', 2, ' ')
!      call add_default ('TREFHTMX', 2, ' ')
!      call add_default ('TS      ', 2, ' ')
!      call add_default ('TSMN    ', 2, ' ')
!      call add_default ('TSMX    ', 2, ' ')
!      call add_default ('U200    ', 2, 'I')
!      call add_default ('U850    ', 2, 'I')
!      call add_default ('V200    ', 2, 'I')
!      call add_default ('V850    ', 2, 'I')
!      call add_default ('Z300    ', 2, 'I')
!      call add_default ('Z500    ', 2, 'I')
!      call add_default ('Z700    ', 2, 'I')

      return
   end subroutine h_default

!#######################################################################

   subroutine add_default (name, tindex, flag)
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the default "on" list for a given history file
! 
! Method: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name  ! field name
      character(len=1), intent(in) :: flag  ! averaging flag

      integer, intent(in) :: tindex         ! history tape index
!
! Local workspace
!
      integer :: f            ! field index
      logical :: found        ! flag indicates field found in masterlist
!
! Check validity of input arguments
!
      if (tindex > ptapes) then
         write(6,*)'ADD_DEFAULT: tape index=', tindex, ' is too big'
         call endrun
      end if

      if (flag /= ' ' .and. flag /= 'A' .and. flag /= 'I' .and. &
          flag /= 'X' .and. flag /= 'M') then

         write(6,*)'ADD_DEFAULT: unknown averaging flag=', flag
         call endrun
      end if
!
! Look through master list for input field name.  When found, set active
! flag for that tape to true.  Also set averaging flag if told to use other
! than default.
!
      found = .false.
      do f=1,nfmaster
         if (trim(name) == trim(masterlist(f)%field%name)) then
            masterlist(f)%actflag(tindex) = .true.
            if (flag /= ' ') then
               masterlist(f)%avgflag(tindex) = flag
               select case (flag)
               case ('A')
                  masterlist(f)%time_op(tindex) = 'mean'
               case ('I')
                  masterlist(f)%time_op(tindex) = ' '
               case ('X')
                  masterlist(f)%time_op(tindex) = 'maximum'
               case ('M')
                  masterlist(f)%time_op(tindex) = 'minimum'
               case default
                  write(6,*)'h_default: unknown avgflag=',flag
                  call endrun
               end select
            end if
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         write(6,*)'ADD_DEFAULT: field=', name, ' not found'
         call endrun
      end if

      return
   end subroutine add_default

!#######################################################################

   subroutine h_override (t)
!----------------------------------------------------------------------- 
! 
! Purpose: Override default history tape contents for a specific tape
!
! Method: Copy the flag into the master field list
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: t         ! history tape index
!
! Local workspace
!
      integer :: f                     ! index over fields
      character(len=1) :: avgflg       ! lcl equiv of avgflag_pertape(t) (to address xlf90 compiler bug)

      avgflg = avgflag_pertape(t)
      do f=1,nfmaster
         select case (avgflg)
         case ('A')
            masterlist(f)%avgflag(t) = avgflag_pertape(t)
            masterlist(f)%time_op(t) = 'mean'
         case ('I')
            masterlist(f)%avgflag(t) = avgflag_pertape(t)
            masterlist(f)%time_op(t) = ' '
         case ('X')
            masterlist(f)%avgflag(t) = avgflag_pertape(t)
            masterlist(f)%time_op(t) = 'maximum'
         case ('M')
            masterlist(f)%avgflag(t) = avgflag_pertape(t)
            masterlist(f)%time_op(t) = 'minimum'
         case default
            write(6,*)'h_override: unknown avgflag=',avgflag_pertape(t)
            call endrun ()
         end select
      end do
   end subroutine h_override
         
!#######################################################################

   subroutine h_define (t)
!----------------------------------------------------------------------- 
! 
! Purpose: Define contents of history file t
! 
! Method: Issue the required netcdf wrapper calls to define the history file contents
! 
!-----------------------------------------------------------------------
      use pspect
      use rgrid
      use commap
      use time_manager, only: get_step_size, get_ref_date
      use filenames,    only: caseid
#ifdef CRM
      use crmdims, only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy
#endif

#include <comctl.h>
#include <comhyb.h>

!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: t   ! tape index
!
! Local workspace
!
      integer :: i, j            ! longitude, latitude indices
      integer :: k               ! ISCCP vertical index
      integer :: l               ! ISCCP optical depth index
      integer :: kl              ! ISCCP merged k and l indices
      integer :: f               ! field index
      integer :: ff              ! varid index for fields output by column
      integer :: numlev          ! number of vertical levels (dimension and loop)
      integer :: ncreal          ! netCDF real data type
      integer :: dtime           ! timestep size
      integer :: ndbase = 0      ! days component of base time
      integer :: nsbase = 0      ! seconds component of base time
      integer :: nbdate          ! base date in yyyymmdd format
      integer :: nbsec           ! time of day component of base date [seconds]
      integer :: yr, mon, day    ! year, month, day components of a date

#ifdef STAGGERED
      integer :: slatdim         ! staggered latitude dimension
      integer :: slondim         ! staggered longitude dimension
      integer :: slatvar         ! variable id for staggered lat
      integer :: slonvar         ! variable id for staggered lon
      integer :: dimen4us(4)     ! dimension arrary for staggered U winds
      integer :: dimen4vs(4)     ! dimension arrary for staggered V winds
      integer :: wsid            ! Staggered latitude weight ID
      
      real(r8) slons( splon )    ! Staggered grid point array (lon)
      real(r8) slats( plat-1 )   ! Staggered grid point array (lat)
#endif

      real(r8) ailev(plevp)      ! interface level values
      real(r8) gausslat(plat)    ! gaussian latitudes
      real(r8) pie               ! 3.14159...
      
      character(len=max_chars) str ! character temporary 
!
! netcdf variables
!
      integer ret                ! function return value
      integer timdim             ! unlimited dimension id
      integer latvar             ! latitude variable id
      integer lonvar             ! longitude variable id
      integer glatvar(pflds)     ! column latitude variable id
      integer glonvar(pflds)     ! column longitude variable id
      integer rlonvar            ! reduced longitude variable id
      integer ps0var             ! variable id for PS0
      integer chardim            ! character dimension id
      
      real(r8) alon(plon)        ! longitude values (degrees)
      real(r8) alev(plev)        ! level values (pascals)
      real(r8) rlon(plon,plat)   ! reduced longitudes (degrees)
      real(r8) alat(plat)        ! latitude values (degrees)
      real(r8) prmid(npres)      ! pressure midpoints of ISCCP data
      real(r8) taumid(ntau)      ! optical depth midpoints of ISCCP data
      real(r8) prstau(npres*ntau)! prmid + taumid/1000

      integer dimenchar(2)       ! character dimension ids
      integer dimen1(1)          ! dimension ids (1d)
      integer dimen2(2)          ! dimension ids (2d)
      integer dimen2t(2)         ! dimension ids (2d time boundaries)
      integer dimen3(3)          ! dimension ids (3d)
      integer dimen4f(4)         ! dimension ids (4d at levels)
      integer dimen4g(4)         ! temp array holding dimension ids for groups of contiguous columns (4d)
      integer dimen4i(4)         ! dimension ids (4d at interfaces)
      integer dimen4n(4)         ! dimension ids (4d at isccp pressure levels)
      integer londim             ! longitude dimension id
      integer latdim             ! latitude dimension id
      integer grouplondim(pflds) ! longitude dimension id
      integer grouplatdim(pflds) ! latitude dimension id
      integer levdim             ! level dimension id
      integer ilevdim            ! interface dimension id
      integer isccp_prs_dim      ! dimension variable for ISCCP pressure levels
      integer isccp_tau_dim      ! dimension variable for ISCCP tau values
      integer isccp_prstau_dim   ! dimension variable for ISCCP pressure*tau levels
      integer tbnddim            ! time_bnds dimension id
      integer levvar             ! level variable id
      integer ilevvar            ! intfc variable id
      integer isccp_prs_var      ! ISCCP mean pressure variable id
      integer isccp_tau_var      ! ISCCP mean optical depth variable id
      integer isccp_prstau_var   ! ISCCP mixed variable id
      integer old_mode           ! returned mode from netcdf call
      integer hyaiid             ! hybrid A coef. intfc var id
      integer hybiid             ! hybrid B coef. intfc var id
      integer hyamid             ! hybrid A coef. level var id
      integer hybmid             ! hybrid B coef. level var id
      integer ntrmid             ! M truncation parameter var id
      integer ntrnid             ! N truncation parameter var id
      integer ntrkid             ! K truncation parameter var id
#ifdef CRM
      integer dimen5(5)          ! dimension ids (5d)
      integer dimen5g(5)         ! temp array holding dimension ids for groups of contiguous columns (4d)
      integer dimen6(6)          ! dimension ids (5d)
      integer dimen6g(6)         ! temp array holding dimension ids for groups of contiguous columns (4d)
      integer crm_x_var
      integer crm_y_var
      integer crm_z_var
      integer crm_x_dim
      integer crm_y_dim
      integer crm_z_dim
      real(r8) crm_x(crm_nx)
      real(r8) crm_y(crm_ny)
      real(r8) crm_z(crm_nz)
#endif
!
      write(6,*)'Opening netcdf history file ', trim(nhfil(t))
      call wrap_create (nhfil(t), nf_clobber, nfid(t))
!
! Setup netcdf file - create the dimensions of lat,lon,time,level
!
      ret = nf_set_fill (nfid(t), nf_nofill, old_mode)
      ret = nf_def_dim (nfid(t), 'lat', plat, latdim)
      ret = nf_def_dim (nfid(t), 'lon', plon, londim)

#ifdef STAGGERED
      ret = nf_def_dim (nfid(t), 'slat', plat-1, slatdim)
      ret = nf_def_dim (nfid(t), 'slon', splon, slondim)
#endif

      ret = nf_def_dim (nfid(t), 'lev', plev, levdim)
      ret = nf_def_dim (nfid(t), 'ilev', plevp, ilevdim)
      ret = nf_def_dim (nfid(t), 'isccp_prs', npres, isccp_prs_dim)
      ret = nf_def_dim (nfid(t), 'isccp_tau', ntau, isccp_tau_dim)
      ret = nf_def_dim (nfid(t), 'isccp_prstau', npres*ntau, isccp_prstau_dim)
      ret = nf_def_dim (nfid(t), 'time', nf_unlimited, timdim)
      ret = nf_def_dim (nfid(t), 'tbnd', 2, tbnddim)
      ret = nf_def_dim (nfid(t), 'chars', 8, chardim)
#ifdef CRM
      ret = nf_def_dim (nfid(t), 'crm_x', crm_nx, crm_x_dim)
      ret = nf_def_dim (nfid(t), 'crm_y', crm_ny, crm_y_dim)
      ret = nf_def_dim (nfid(t), 'crm_z', crm_nz, crm_z_dim)
#endif
!
! create dimensions for groups of contiguous columns
! variables with for single columns are created later in this routine
!
      if (ngroup(t).ne.0) then
         do i = 1, ngroup(t)
!	    write(6,*)tape(t)%column(i)%lat_name,tape(t)%column(i)%lon_name
            ret = nf_def_dim (nfid(t), tape(t)%column(i)%lat_name, tape(t)%column(i)%num_lats, grouplatdim(i))
            call wrap_def_var (nfid(t),tape(t)%column(i)%lat_name,nf_double,1,grouplatdim(i),glatvar(i))
            call wrap_put_att_text (nfid(t), glatvar(i),'long_name','latitude')
            call wrap_put_att_text (nfid(t), glatvar(i),'units','degrees_north')
            ret = nf_def_dim (nfid(t), tape(t)%column(i)%lon_name, tape(t)%column(i)%num_lons, grouplondim(i))
            call wrap_def_var (nfid(t),tape(t)%column(i)%lon_name,nf_double,1,grouplondim(i),glonvar(i))
            call wrap_put_att_text (nfid(t), glonvar(i),'long_name','longitude')
            call wrap_put_att_text (nfid(t), glonvar(i),'units','degrees_east')
         end do
      end if
!
! setup dimension arrays for 1,2,3,4d variables 
!     
 write(6,*)'define dimensions'
      dimen1(1) = timdim

      dimen2(1) = londim
      dimen2(2) = latdim
      
      dimen2t(1) = tbnddim
      dimen2t(2) = timdim

      dimen3(1) = londim
      dimen3(2) = latdim
      dimen3(3) = timdim
!            
      dimen4f(1) = londim
      dimen4f(4) = timdim

      dimen4i(1) = londim
      dimen4i(4) = timdim

      dimen4f(2) = latdim
      dimen4f(3) = levdim
         
      dimen4i(2) = latdim
      dimen4i(3) = ilevdim

      dimen4n(1) = londim
      dimen4n(2) = latdim
      dimen4n(3) = isccp_prstau_dim
      dimen4n(4) = timdim
#ifdef CRM
      dimen5(1) = londim
      dimen5(2) = latdim
      dimen5(3) = crm_x_dim
      dimen5(4) = crm_y_dim
      dimen5(5) = timdim
      dimen5g(3) = crm_x_dim
      dimen5g(4) = crm_y_dim
      dimen5g(5) = timdim
      dimen6(1) = londim
      dimen6(2) = latdim
      dimen6(3) = crm_x_dim
      dimen6(4) = crm_y_dim
      dimen6(5) = crm_z_dim
      dimen6(6) = timdim
      dimen6g(3) = crm_x_dim
      dimen6g(4) = crm_y_dim
      dimen6g(5) = crm_z_dim
      dimen6g(6) = timdim
#endif

#ifdef STAGGERED
      dimen4us(1) = londim
      dimen4us(2) = slatdim
      dimen4us(3) = levdim
      dimen4us(4) = timdim
      
      dimen4vs(1) = slondim
      dimen4vs(2) = latdim
      dimen4vs(3) = levdim
      dimen4vs(4) = timdim
#endif

! define variables to label the dimensions, use the same names as the
! dimensions

      call wrap_def_var (nfid(t),'P0',nf_double,0,0,ps0var)
      str = 'reference pressure'
      call wrap_put_att_text (nfid(t), ps0var, 'long_name', str)
      call wrap_put_att_text (nfid(t), ps0var, 'units', 'Pa')
      
      call wrap_def_var (nfid(t),'lat',nf_double,1,LATDIM,latvar)
      call wrap_put_att_text (nfid(t), latvar, 'long_name', 'latitude')
      call wrap_put_att_text (nfid(t), latvar, 'units', 'degrees_north')
      
      if (fullgrid) then

         call wrap_def_var (nfid(t),'lon',nf_double,1,LONDIM,lonvar)
         call wrap_put_att_text (nfid(t), lonvar,'long_name','longitude')
         call wrap_put_att_text (nfid(t), lonvar,'units','degrees_east')

      else

         call wrap_def_var (nfid(t),'rlon',nf_double,2,dimen2,rlonvar)
         call wrap_put_att_text (nfid(t), rlonvar, 'long_name', 'reduced_longitude')
         call wrap_put_att_text (nfid(t), rlonvar,'units','degrees_east')

      end if

#ifdef CRM
      call wrap_def_var (nfid(t),'crm_x',nf_double,1,crm_x_dim,crm_x_var)
      call wrap_put_att_text (nfid(t), crm_x_var, 'long_name', 'crm horiz. distance in x')
      call wrap_put_att_text (nfid(t), crm_x_var, 'units', 'km')

      call wrap_def_var (nfid(t),'crm_y',nf_double,1,crm_y_dim,crm_y_var)
      call wrap_put_att_text (nfid(t), crm_y_var, 'long_name', 'crm horiz. distance in y')
      call wrap_put_att_text (nfid(t), crm_y_var, 'units', 'km')

      call wrap_def_var (nfid(t),'crm_z',nf_double,1,crm_z_dim,crm_z_var)
      call wrap_put_att_text (nfid(t), crm_z_var, 'long_name', 'crm nominal pressure levels')
      call wrap_put_att_text (nfid(t), crm_z_var, 'units', 'mb')
#endif
! If a staggered grid is in use, output the lat and lon arrays variables.
! Note: the staggered grid is currently independent of the reduced grid,
! and the staggered grid is a full grid.

#ifdef STAGGERED
      call wrap_def_var (nfid(t), 'slat', nf_double,1,SLATDIM,slatvar)
      call wrap_put_att_text (nfid(t), slatvar, 'long_name', 'staggered latitude')
      call wrap_put_att_text (nfid(t), slatvar, 'units', 'degrees_north')

      call wrap_def_var (nfid(t),'slon',nf_double,1,SLONDIM,slonvar)
      call wrap_put_att_text (nfid(t), slonvar, 'long_name', 'staggered longitude')
      call wrap_put_att_text(nfid(t), slonvar, 'units', 'degrees_east')
      
      call wrap_def_var (nfid(t),'w_stag',nf_double,1,slatdim,wsid)
      str = 'staggered latitude weights'
      call wrap_put_att_text (nfid(t), wsid, 'long_name', str)
#endif

      call wrap_def_var (nfid(t),'lev',nf_double,1,LEVDIM,levvar)
      str = 'hybrid level at midpoints (1000*(A+B))'
      call wrap_put_att_text (nfid(t), levvar, 'long_name', str)
      str = 'level'
      call wrap_put_att_text (nfid(t), levvar, 'units', str)
      call wrap_put_att_text (nfid(t), levvar, 'positive', 'down')
      call wrap_put_att_text (nfid(t), levvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
      call wrap_put_att_text (nfid(t), levvar, 'formula_terms', 'a: hyam b: hybm p0: P0 ps: PS')

      call wrap_def_var (nfid(t),'ilev',nf_double,1,ILEVDIM,ilevvar)
      str = 'hybrid level at interfaces (1000*(A+B))'
      call wrap_put_att_text (nfid(t), ilevvar, 'long_name', str)
      str = 'level'
      call wrap_put_att_text (nfid(t), ilevvar, 'units', str)
      call wrap_put_att_text (nfid(t), ilevvar, 'positive', 'down')
      call wrap_put_att_text (nfid(t), ilevvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
      call wrap_put_att_text (nfid(t), ilevvar, 'formula_terms', 'a: hyai b: hybi p0: P0 ps: PS')

! ISCCP pressure, optical depth, and mixed dimension

      call wrap_def_var (nfid(t), 'isccp_prs', NF_DOUBLE, 1, isccp_prs_dim, isccp_prs_var)
      str = 'Mean ISCCP pressure'
      call wrap_put_att_text (nfid(t), isccp_prs_var, 'long_name', str)
      str = 'mb'
      call wrap_put_att_text (nfid(t), isccp_prs_var, 'units', str)
      call wrap_put_att_realx (nfid(t), isccp_prs_var, 'isccp_prs_bnds', NF_DOUBLE, npres+1, prlim)

      call wrap_def_var (nfid(t), 'isccp_tau', NF_DOUBLE, 1, isccp_tau_dim, isccp_tau_var)
      str = 'Mean ISCCP optical depth'
      call wrap_put_att_text (nfid(t), isccp_tau_var, 'long_name', str)
      str = 'unitless'
      call wrap_put_att_text (nfid(t), isccp_tau_var, 'units', str)
      call wrap_put_att_realx (nfid(t), isccp_tau_var, 'isccp_tau_bnds', NF_DOUBLE, ntau+1, taulim)

      call wrap_def_var (nfid(t), 'isccp_prstau', NF_DOUBLE, 1, isccp_prstau_dim, isccp_prstau_var)
      str = 'Mean pressure (mb).mean optical depth (unitless)/1000'
      call wrap_put_att_text (nfid(t), isccp_prstau_var, 'long_name', str)
      str = 'mixed'
      call wrap_put_att_text (nfid(t), isccp_prstau_var, 'units', str)

      call get_ref_date(yr, mon, day, nbsec)
      nbdate = yr*10000 + mon*100 + day
      call wrap_def_var (nfid(t),'time',nf_double,1,TIMDIM,timeid(t))
      call wrap_put_att_text (nfid(t), timeid(t), 'long_name', 'time')
      str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
      call wrap_put_att_text (nfid(t), timeid(t), 'units', str)
      call wrap_put_att_text (nfid(t), timeid(t), 'calendar', 'noleap')
      call wrap_put_att_text (nfid(t), timeid(t), 'bounds', 'time_bnds')

      call wrap_def_var (nfid(t),'time_bnds',nf_double,2,dimen2t,tbndid(t))
      call wrap_put_att_text (nfid(t), tbndid(t), 'long_name', 'time interval endpoints')
      str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
      call wrap_put_att_text (nfid(t), tbndid(t), 'units', str)
!
! Character
!
      dimenchar(1) = chardim
      dimenchar(2) = timdim
      call wrap_def_var (nfid(t),'date_written',NF_CHAR,2,dimenchar, date_writtenid(t))
      call wrap_def_var (nfid(t),'time_written',NF_CHAR,2,dimenchar, time_writtenid(t))
!
! Integer Header
!
      call wrap_def_var (nfid(t),'ntrm',NF_INT,0,0,ntrmid)
      str = 'spectral truncation parameter M'
      call wrap_put_att_text (nfid(t), ntrmid, 'long_name', str)

      call wrap_def_var (nfid(t),'ntrn',NF_INT,0,0,ntrnid)
      str = 'spectral truncation parameter N'
      call wrap_put_att_text (nfid(t), ntrnid, 'long_name', str)

      call wrap_def_var (nfid(t),'ntrk',NF_INT,0,0,ntrkid)
      str = 'spectral truncation parameter K'
      call wrap_put_att_text (nfid(t), ntrkid, 'long_name', str)

      call wrap_def_var (nfid(t),'ndbase',NF_INT,0,0,ndbaseid(t))
      str = 'base day'
      call wrap_put_att_text (nfid(t), ndbaseid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'nsbase',NF_INT,0,0,nsbaseid(t))
      str = 'seconds of base day'
      call wrap_put_att_text (nfid(t), nsbaseid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'nbdate',NF_INT,0,0,nbdateid(t))
      str = 'base date (YYYYMMDD)'
      call wrap_put_att_text (nfid(t), nbdateid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'nbsec',NF_INT,0,0,nbsecid(t))
      str = 'seconds of base date'
      call wrap_put_att_text (nfid(t), nbsecid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'mdt',NF_INT,0,0,mdtid)
      call wrap_put_att_text (nfid(t), mdtid, 'long_name', 'timestep')
      call wrap_put_att_text (nfid(t), mdtid, 'units', 's')

      call wrap_def_var (nfid(t),'nlon',NF_INT,1,latdim,nlonid(t))
      str = 'number of longitudes'
      call wrap_put_att_text (nfid(t), nlonid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'wnummax',NF_INT,1,latdim,wnummaxid(t))
      str = 'cutoff Fourier wavenumber'
      call wrap_put_att_text (nfid(t), wnummaxid(t), 'long_name', str)
!
! Floating point time-invariant
!
      call wrap_def_var (nfid(t),'hyai',NF_DOUBLE,1,ilevdim,hyaiid)
      str = 'hybrid A coefficient at layer interfaces'
      call wrap_put_att_text (nfid(t), hyaiid, 'long_name', str)

      call wrap_def_var (nfid(t),'hybi',NF_DOUBLE,1,ilevdim,hybiid)
      str = 'hybrid B coefficient at layer interfaces'
      call wrap_put_att_text (nfid(t), hybiid, 'long_name', str)

      call wrap_def_var (nfid(t),'hyam',NF_DOUBLE,1,levdim,hyamid)
      str = 'hybrid A coefficient at layer midpoints'
      call wrap_put_att_text (nfid(t), hyamid, 'long_name', str)

      call wrap_def_var (nfid(t),'hybm',NF_DOUBLE,1,levdim,hybmid)
      str = 'hybrid B coefficient at layer midpoints'
      call wrap_put_att_text (nfid(t), hybmid, 'long_name', str)

      call wrap_def_var (nfid(t),'gw',NF_DOUBLE,1,latdim,gwid)
      str = 'gauss weights'
      call wrap_put_att_text (nfid(t), gwid, 'long_name', str)
!
! Character header information
!
      str = 'CF-1.0'
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'Conventions', str)
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'source', 'CAM')
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'case',caseid)
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'title',ctitle)
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'logname',logname)
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'host', host)
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'Version', &
           '$Name:  $')
      call wrap_put_att_text (nfid(T), NF_GLOBAL, 'revision_Id', &
           '$Id: history.F90,v 1.26.2.38 2003/12/15 18:52:35 hender Exp $')
!
! Create variables for model timing and header information
!
      call wrap_def_var (nfid(t),'ndcur   ',nf_int,1,dimen1,ndcurid(t))
      str = 'current day (from base day)'
      call wrap_put_att_text (nfid(t), ndcurid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'nscur   ',nf_int,1,dimen1,nscurid(t))
      str = 'current seconds of current day'
      call wrap_put_att_text (nfid(t), nscurid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'date    ',nf_int,1,dimen1,dateid(t))
      str = 'current date (YYYYMMDD)'
      call wrap_put_att_text (nfid(t), dateid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'datesec ',nf_int,1,dimen1, datesecid(t))
      str = 'current seconds of current date'
      call wrap_put_att_text (nfid(t), datesecid(t), 'long_name', str)

      call wrap_def_var (nfid(t),'nsteph  ',nf_int,1,dimen1,nstephid(t))
      str = 'current timestep'
      call wrap_put_att_text (nfid(t), nstephid(t), 'long_name', str)
!
! Create variables and attributes for field list
!
      do f=1,nflds(t)
         numlev = tape(t)%hlist(f)%field%numlev
         if (tape(t)%hlist(f)%hwrt_prec == 8) then
            ncreal = nf_double
         else
            ncreal = nf_float
         end if
!
!  Create variables and atributes for fields written out as columns
!         
         if (ngroup(t).ne.0) then
            do i=1,ngroup(t)
               ff=(f-1)*ngroup(t)+i
               dimen4g(1)=grouplondim(i)
               dimen4g(2)=grouplatdim(i)
               dimen4g(4)=timdim
#ifdef CRM
               dimen5g(1)=grouplondim(i)
               dimen5g(2)=grouplatdim(i)
               dimen6g(1)=grouplondim(i)
               dimen6g(2)=grouplatdim(i)
#endif

               if (numlev == npres*ntau .and. tape(t)%hlist(f)%field%flag_isccplev) then
                  dimen4g(3)=isccp_prstau_dim
                  call wrap_def_var(nfid(t), tape(t)%hlist(f)%field_column_name(i),ncreal,4,dimen4g,varid(ff,t))
               else if (numlev == 1) then
                  dimen4g(3)=timdim
                  call wrap_def_var(nfid(t),trim(tape(t)%hlist(f)%field_column_name(i)),ncreal,3,dimen4g,varid(ff,t))

               else if (numlev == plev) then
                  dimen4g(3)=levdim
                  call wrap_def_var(nfid(t), trim(tape(t)%hlist(f)%field_column_name(i)), ncreal,4,dimen4g,varid(ff,t))
                  
               else if (numlev == plevp) then
                  dimen4g(3)=ilevdim
                  call wrap_def_var(nfid(t), trim(tape(t)%hlist(f)%field_column_name(i)), ncreal,4,dimen4g,varid(ff,t))
#ifdef CRM
               else if (numlev == crm_nx*crm_ny) then

                  call wrap_def_var(nfid(t), trim(tape(t)%hlist(f)%field_column_name(i)), ncreal,5,dimen5g,varid(ff,t))

	       else if (numlev == crm_nx*crm_ny*crm_nz) then
                  call wrap_def_var(nfid(t), trim(tape(t)%hlist(f)%field_column_name(i)), ncreal,6,dimen6g,varid(ff,t))
#endif
               else
                  write(6,*)'H_DEFINE: bad numlev=',numlev
                  call endrun
               end if

               if (.not. fullgrid .or. tape(t)%hlist(f)%field%flag_xyfill) then
                  call wrap_put_att_realx (nfid(t), varid(ff,t), '_FillValue', ncreal, 1, fillvalue)
               end if

!JR Add missing_value for nco operators

               if (tape(t)%hlist(f)%field%flag_xyfill) then
                  call wrap_put_att_realx (nfid(t), varid(ff,t), 'missing_value', ncreal, 1, fillvalue)
               end if
               
               str = tape(t)%hlist(f)%field%units
               if ( str(1:1) /= ' ' ) then
                  call wrap_put_att_text (nfid(t), varid(ff,t), 'units', str)
               end if

               str = tape(t)%hlist(f)%field%long_name
               call wrap_put_att_text (nfid(t), varid(ff,t), 'long_name', str)
!
! Assign field attributes defining valid levels and averaging info
!
               str = tape(t)%hlist(f)%time_op
               select case (str)
               case ('mean', 'maximum', 'minimum' )
                  call wrap_put_att_text (nfid(t), varid(ff,t),'cell_method', 'time: '//str)
               end select
            end do
!
!  else create variables and atributes for fields written out as a full model grid
!         
         else
         
            if (numlev == npres*ntau .and. tape(t)%hlist(f)%field%flag_isccplev) then

               call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncreal,4,dimen4n,varid(f,t))

            else if (numlev == 1) then

               call wrap_def_var(nfid(t),tape(t)%hlist(f)%field%name,ncreal,3,dimen3,varid(f,t))

            else if (numlev == plev) then	

#ifdef STAGGERED
               select case (tape(t)%hlist(f)%field%name)
               case ('US')
                  call wrap_def_var(nfid(t),tape(t)%hlist(f)%field%name,ncreal,4,dimen4us,varid(f,t))
               case ('VS')
                  call wrap_def_var(nfid(t),tape(t)%hlist(f)%field%name,ncreal,4,dimen4vs,varid(f,t))
               case default
                  call wrap_def_var(nfid(t),tape(t)%hlist(f)%field%name,ncreal,4,dimen4f,varid(f,t))
               end select
#else
               call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncreal,4,dimen4f,varid(f,t))
               
#endif
            else if (numlev == plevp) then	
               
               call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncreal,4,dimen4i,varid(f,t))
#ifdef CRM
            else if (numlev == crm_nx*crm_ny) then

               call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncreal,5,dimen5,varid(f,t))

            else if (numlev == crm_nx*crm_ny*crm_nz) then

               call wrap_def_var(nfid(t), tape(t)%hlist(f)%field%name, ncreal,6,dimen6,varid(f,t))
#endif
            else
               write(6,*)'H_DEFINE: bad numlev=',numlev
               call endrun
            end if
!	
            if (.not. fullgrid .or. tape(t)%hlist(f)%field%flag_xyfill) then
               call wrap_put_att_realx (nfid(t), varid(f,t), '_FillValue', ncreal, 1, fillvalue)
            end if

!JR Add missing_value for nco operators

            if (tape(t)%hlist(f)%field%flag_xyfill) then
               call wrap_put_att_realx (nfid(t), varid(f,t), 'missing_value', ncreal, 1, fillvalue)
            end if

            str = tape(t)%hlist(f)%field%units
            if ( str(1:1) /= ' ' ) then
               call wrap_put_att_text (nfid(t), varid(f,t), 'units', str)
            end if

            str = tape(t)%hlist(f)%field%units
            if ( str(1:1) /= ' ' ) then
               call wrap_put_att_text (nfid(t), varid(f,t), 'units', str)
            end if
            
            str = tape(t)%hlist(f)%field%long_name
            call wrap_put_att_text (nfid(t), varid(f,t), 'long_name', str)
!
! Assign field attributes defining valid levels and averaging info
!
            str = tape(t)%hlist(f)%time_op
            select case (str)
            case ('mean', 'maximum', 'minimum' )
               call wrap_put_att_text (nfid(t), varid(f,t),'cell_method', 'time: '//str)
            end select
         endif
      end do
!
      ret = nf_enddef(nfid(t))
      write(6,*)'H_DEFINE: Successfully opened netcdf file '
!
! Write time-invariant portion of history header
!
      call wrap_put_var_realx (nfid(t), ps0var, ps0)

      pie = 4.*atan(1.)
      do j=1,plat
         gausslat(j) = (180./pie)*clat(j)
      end do
      call wrap_put_var_realx (nfid(t), latvar, gausslat)
!
! Coordinate var defined only in full grid: otherwise define 2-d var.
!
      if (fullgrid) then
         do i=1,plon
            alon(i) = (i-1) * 360.0 / plon
         end do
         call wrap_put_var_realx (nfid(t), lonvar, alon)

      else

         do j=1,plat
            rlon(:,j) = fillvalue
            do i=1,nlon(j)
               rlon(i,j) = (i-1) * 360.0 / nlon(j)
            end do
         end do
         call wrap_put_var_realx (nfid(t), rlonvar, rlon)
      end if
#ifdef CRM
      do i=1,crm_nx
        crm_x(i) = (i-1)*crm_dx*0.001
      end do
      call wrap_put_var_realx (nfid(t), crm_x_var, crm_x)
      do i=1,crm_ny
        crm_y(i) = (i-1)*crm_dy*0.001
      end do
      call wrap_put_var_realx (nfid(t), crm_y_var, crm_y)
#endif

#ifdef STAGGERED

! Calculate latitudes for the staggered grid.

      do j = 1, plat-1

! Staggered latitudes fall in between the regular latitudes. Note that
! for Lin-Rood dynamics, gausslat latitudes are not gaussian!

         slats(j) = (180./pie) * clat_staggered(j)

      end do

! Sanity check.

      if (slats(plat-1) .gt. 90.0) then
         write(6,*) "H_DEFINE: WARNING - last staggered grid latitude"
         write(6,*) "was calculated to be: ", slats(plat-1)
         write(6,*) "Point has been redefined to be 90.0"
         slats(plat-1) = 90.0
      endif

      call wrap_put_var_realx (nfid(t), slatvar, slats)
      
      do i = 1, splon
         slons(i) = ((i-1) * 360.0 / splon) - 180.0/splon
      enddo

      call wrap_put_var_realx (nfid(t), slonvar, slons)
      call wrap_put_var_realx (nfid(t), wsid, w_staggered)
#endif
!
! write out coordinate variables for columns
!
!jt This needs to be fixed for reduced grid : note the bogus dimension 1 on londeg
!

      if (ngroup(t).ne.0) then
         do i = 1, ngroup(t)
            call wrap_put_var_realx (nfid(t), glonvar(i), londeg(tape(t)%column(i)%columnlon(1):tape(t)%column(i)%columnlon(2),1))
            call wrap_put_var_realx (nfid(t), glatvar(i), latdeg(tape(t)%column(i)%columnlat(1):tape(t)%column(i)%columnlat(2)))
         end do
      end if
!
! 0.01 converts Pascals to millibars
!
      alev(:plev) = 0.01*ps0*(hyam(:plev) + hybm(:plev))
      ailev(:plevp) = 0.01*ps0*(hyai(:plevp) + hybi(:plevp))

      do k=1,npres
         prmid(k) = 0.5*(prlim(k) + prlim(k+1))
      end do

      do k=1,ntau
         taumid(k) = 0.5*(taulim(k) + taulim(k+1))
      end do

!JR Kludgey way of combining pressure and optical depth into a single dimension:
!JR pressure in millibars will show up on the left side of the decimal point, 
!JR optical depth/1000 will show up on the right

      do k=1,npres
         do l=1,ntau
            kl = (k-1)*ntau + l
            prstau(kl) = prmid(k) + taumid(l)*0.001
         end do
      end do

      call wrap_put_var_realx (nfid(t), levvar, alev)
      call wrap_put_var_realx (nfid(t), ilevvar, ailev)
      call wrap_put_var_realx (nfid(t), hyaiid, hyai)
      call wrap_put_var_realx (nfid(t), hybiid, hybi)
      call wrap_put_var_realx (nfid(t), hyamid, hyam)
      call wrap_put_var_realx (nfid(t), hybmid, hybm)
      call wrap_put_var_realx (nfid(t), gwid, w)
      call wrap_put_var_realx (nfid(t), isccp_prs_var, prmid)
      call wrap_put_var_realx (nfid(t), isccp_tau_var, taumid)
      call wrap_put_var_realx (nfid(t), isccp_prstau_var, prstau)
#ifdef CRM
      do i=1,crm_nz
        crm_z(i) = alev(plev-i+1)
      end do
      call wrap_put_var_realx (nfid(t), crm_z_var, crm_z)
#endif

      call wrap_put_var_int (nfid(t), ntrmid, ptrm)
      call wrap_put_var_int (nfid(t), ntrnid, ptrn)
      call wrap_put_var_int (nfid(t), ntrkid, ptrk)
      dtime = get_step_size()
      call wrap_put_var_int (nfid(t), mdtid, dtime)
!
! Model date info
!
      call wrap_put_var_int (nfid(t), ndbaseid(t), ndbase)
      call wrap_put_var_int (nfid(t), nsbaseid(t), nsbase)

      call wrap_put_var_int (nfid(t), nbdateid(t), nbdate)
      call wrap_put_var_int (nfid(t), nbsecid(t), nbsec)
!
! Reduced grid info
!
      call wrap_put_var_int (nfid(t), nlonid(t), nlon)
      call wrap_put_var_int (nfid(t), wnummaxid(t), wnummax)

      return
   end subroutine h_define

!#######################################################################

character(len=10) function date2yyyymmdd (date)

   implicit none

! Input arguments

   integer, intent(in) :: date

! Local workspace

   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

   if (date < 0) then
      write(6,*)'DATE2YYYYMMDD: negative date not allowed'
      call endrun ()
   end if

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   if (month < 1 .or. month > 12) then
      write(6,*)'DATE2YYYYMMDD: invalid month encountered:', month
      call endrun ()
   end if

   if (day < 1 .or. day > ndm(month)) then
      write(6,*)'DATE2YYYYMMDD: invalid day encountered:', day
      call endrun ()
   end if

   write(date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)
   return
end function date2yyyymmdd

!#######################################################################

character(len=8) function sec2hms (seconds)

   implicit none
   
! Input arguments

   integer, intent(in) :: seconds

! Local workspace

   integer :: hours     ! hours of hh:mm:ss
   integer :: minutes   ! minutes of hh:mm:ss
   integer :: secs      ! seconds of hh:mm:ss

   if (seconds < 0 .or. seconds > 86400) then
      write(6,*)'SEC2HRS: bad input seconds:', seconds
      call endrun ()
   end if

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   if (minutes < 0 .or. minutes > 60) then
      write(6,*)'SEC2HRS: bad minutes = ',minutes
      call endrun ()
   end if

   if (secs < 0 .or. secs > 60) then
      write(6,*)'SEC2HRS: bad secs = ',secs
      call endrun ()
   end if

   write(sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)
   return
end function sec2hms

!#######################################################################

   subroutine h_normalize (f, t)
!----------------------------------------------------------------------- 
! 
! Purpose: Normalize fields on a history file by the number of accumulations
! 
! Method: Loop over fields on the tape.  Need averaging flag and number of
!         accumulations to perform normalization.
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: f       ! field index
      integer, intent(in) :: t       ! tape index
!
! Local workspace
!
      integer begver                 ! on-node vert start index
      integer endver                 ! on-node vert end index
      integer c                      ! chunk (or lat) index
      integer coldimin               ! column dimension of model array
      integer endi                   ! terminating column index
      integer begdim3                ! on-node chunk or lat start index
      integer enddim3                ! on-node chunk or lat end index  
      integer, pointer :: nacs(:)    ! accumulation counter
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer :: ieu, jeu  ! number of elements in each dimension
      integer k

      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
      character*1 avgflag            ! averaging flag

      type (hbuffer_2d) :: hbuf      ! Per field history buffer
      type (dim_index_2d) :: dimind  ! 2-D dimension index
      
      call t_startf ('h_normalize')

      begver  = tape(t)%hlist(f)%field%begver
      endver  = tape(t)%hlist(f)%field%endver
      avgflag = tape(t)%hlist(f)%avgflag
!
! normalize by number of accumulations for averaged case
!
      begdim3 = tape(t)%hlist(f)%field%begdim3
      enddim3 = tape(t)%hlist(f)%field%enddim3
      flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill

      do c=begdim3,enddim3
         call assoc_hbuf2d_with_hbuf3d (hbuf, tape(t)%hlist(f)%hbuf, c)
         coldimin = tape(t)%hlist(f)%field%coldimin
         nacs => tape(t)%hlist(f)%nacs(:coldimin,c)
         endi =  tape(t)%hlist(f)%field%colperdim3(c)

         dimind = dim_index_2d (1,endi,begver,endver)
         ib = dimind%beg1
         ie = dimind%end1
         jb = dimind%beg2
         je = dimind%end2
         
         ieu = ie-ib+1
         jeu = je-jb+1

         if (flag_xyfill) then
            if (associated(hbuf%buf8)) then
               do k=1,jeu
                  where (nacs(1:ieu) == 0)
                     hbuf%buf8(1:ieu,k) = fillvalue
                  endwhere
               end do
            else if (associated(hbuf%buf4)) then
               do k=1,jeu
                  where (nacs(1:ieu) == 0)
                     hbuf%buf4(1:ieu,k) = fillvalue
                  endwhere
               end do
            end if
         end if
         
         if (avgflag == 'A') then
            if (associated(hbuf%buf8)) then
               do k=1,jeu
                  where (nacs(1:ieu) /= 0)
                     hbuf%buf8(1:ieu,k) = hbuf%buf8(1:ieu,k) / nacs(1:ieu)
                  endwhere
               end do
            else if (associated(hbuf%buf4)) then
               do k=1,jeu
                  where (nacs(1:ieu) /= 0)
                     hbuf%buf4(1:ieu,k) = hbuf%buf4(1:ieu,k) / nacs(1:ieu)
                  endwhere
               end do
            end if
         end if
      end do

      call t_stopf ('h_normalize')
      
      return
   end subroutine h_normalize

!#######################################################################

   subroutine h_zero (f, t)
!----------------------------------------------------------------------- 
! 
! Purpose: Zero out accumulation buffers for a tape
! 
! Method: Loop through fields on the tape
! 
!-----------------------------------------------------------------------

      integer, intent(in) :: f     ! field index
      integer, intent(in) :: t     ! tape index
!
! Local workspace
!
      integer begver               ! on-node vert start index
      integer endver               ! on-node vert end index
      integer c                    ! chunk index
      integer endi                 ! terminating column index
      integer begdim3              ! on-node chunk or lat start index
      integer enddim3              ! on-node chunk or lat end index  
      type (dim_index_3d) :: dimind ! 3-D dimension index
      
      call t_startf ('h_zero')

      begdim3 = tape(t)%hlist(f)%field%begdim3
      enddim3 = tape(t)%hlist(f)%field%enddim3
      begver  = tape(t)%hlist(f)%field%begver
      endver  = tape(t)%hlist(f)%field%endver

      do c=begdim3,enddim3
         endi   = tape(t)%hlist(f)%field%colperdim3(c)

         dimind = dim_index_3d (1,endi,begver,endver,c,c)
         call set_hbuf_section_to_val (tape(t)%hlist(f)%hbuf,dimind,0._r8)
         tape(t)%hlist(f)%nacs(:endi,c) = 0
      end do

      call t_stopf ('h_zero')

      return
   end subroutine h_zero

!#######################################################################

   subroutine dump_field (f, t)
!----------------------------------------------------------------------- 
! 
! Purpose: Write a variable to a history tape
! 
! Method: If SPMD, first gather the data to the master processor.  
!         Next, transpose the data to COORDS order (the default).
!         Finally, issue the netcdf call to write the variable
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
# if ( defined STAGGERED )
      use spmd_dyn, only: npes, compute_gsfactors, comm_y
      use pmgrid, only: myid_z, strip3dxzy, strip3dxzyp, strip2d, twod_decomp
      use parutilitiesmodule, only: commglobal, pargatherreal, pargatherreal4
# else
      use spmd_dyn, only: npes, compute_gsfactors
# endif
      use mpishorthand
#endif
      use dycore, only: dycore_is
#ifdef CRM
      use crmdims, only: crm_nx, crm_ny, crm_nz
#endif
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: t, f ! tape, field indices
!
! Local workspace
!
      integer :: numlev           ! number of vertical levels (dimension and loop)
#ifdef CRM
      integer :: start(6)         ! array of starting field indices for the field in the netcdf file
      integer :: count(6)         ! array of count values for the field in the netcdf file
#else
      integer :: start(4)         ! array of starting field indices for the field in the netcdf file
      integer :: count(4)         ! array of count values for the field in the netcdf file
#endif
      integer :: columnlon1       ! tmp
      integer :: columnlon2       ! tmp
      integer :: columnlat1       ! tmp
      integer :: columnlat2       ! tmp

      type (hbuffer_3d) :: tmpxyzbuf ! temporary array to hold contiguous column data in COORDS orDer
      type (hbuffer_3d) :: xyzbuf ! temporary array to hold data in COORDS order
      type (hbuffer_3d) :: xzybuf ! temporary array holding full x, y, z dims

#ifdef SPMD
      integer :: coldimin         ! column dimension of model array
      integer :: numsend          ! number of items to be sent
      integer :: numrecv(0:npes-1)! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numowned         ! number of items owned by this MPI task
      integer :: mpireal          ! MPI real data type
#endif
      integer :: ff,i,lon,lat,k              ! varid index for fields with groups of columns
      type (dim_index_3d) :: dimind  ! 3-D dimension index
      logical :: flag_xyfill         ! non-applicable xy points flagged with fillvalue
      character*1 avgflag            ! averaging flag

      avgflag     = tape(t)%hlist(f)%avgflag
      flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill
      numlev      = tape(t)%hlist(f)%field%numlev
#ifdef HDEBUG
      write(6,*)'DUMP_FIELD: writing ',tape(t)%hlist(f)%field%name,masterproc
#endif
      start(:) = 0
      start(:) = 0

      start(1) = 1
      start(2) = 1
      
      count(1) = plon

      if (numlev == 1) then

         start(3) = nfils(t)
         count(2) = plat
         count(3) = 1

      else if (numlev > 1) then

         start(3) = 1
         start(4) = nfils(t)
         count(2) = plat
         count(3) = numlev
         count(4) = 1

      else

         write(6,*)'DUMP_FIELD: bad numlev=', numlev
         stop 999
         
      end if

#ifdef STAGGERED
      if (tape(t)%hlist(f)%field%name == 'US') then
         count(2) = plat - 1
      else if (tape(t)%hlist(f)%field%name == 'VS') then
         count(1) = splon
      endif
#endif
!------------------------------
! CRM is here:

#ifdef CRM

      if(numlev == crm_nx*crm_ny) then
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = 1
        start(5) = nfils(t)
        count(1) = plon
        count(2) = plat
        count(3) = crm_nx
        count(4) = crm_ny
        count(5) = 1
      end if
      if(numlev == crm_nx*crm_ny*crm_nz) then
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = 1
        start(5) = 1
        start(6) = nfils(t)
        count(1) = plon
        count(2) = plat
        count(3) = crm_nx
        count(4) = crm_ny
        count(5) = crm_nz
        count(6) = 1
      end if

#endif

#ifdef HDEBUG
      write(6,*)'DUMP_FIELD: writing time indx ', nfils(t), ' field id', varid(f,t), &
                ' name ', tape(t)%hlist(f)%field%name,masterproc
      write(6,*)'start=',start
      write(6,*)'count=',count
      call print_memusage ('dump_field')
#endif
!
! Transpose to COORDS order
!
! Note: this may or may not work when outputting the US variable
! in Lin-Rood because it has one less latitude (plat-1).
! xyzbuf is the memory for field to be written.  Initialize entire array to fillvalue since 
! history module does not know how to map unused portions of the physics domain
!
      if (masterproc) then
         dimind = dim_index_3d (1,plon,1,numlev,1,plat)
         call allocate_hbuf (xzybuf,dimind,tape(t)%hlist(f)%hbuf_prec)
         dimind = dim_index_3d (1,plon,1,plat,1,numlev)
         ! fvitt - for US field want xyzbuf to start at 2nd latitude
         if (tape(t)%hlist(f)%field%name == 'US') dimind = &
           dim_index_3d (1,plon,2,plat,1,numlev)
         call allocate_hbuf (xyzbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
      else
         call assoc_hbuf_with_nothing (xzybuf,tape(t)%hlist(f)%hbuf_prec)
         call nullify_hbuf (xyzbuf)
      end if
!
! physics decomposition: convert from col,lev,chunk -> lon,lev,lat -> lon,lat,lev
! dynamics decomposition: convert from lon,lev,lat -> lon,lat,lev
!
      select case (tape(t)%hlist(f)%field%decomp_type)
      case (phys_decomp)
         call gather_chunk_to_field_hbuf (1, numlev, 1, plon, tape(t)%hlist(f)%hbuf, xzybuf)
      case (dyn_decomp)
#ifdef SPMD
         if (tape(t)%hlist(f)%hbuf_prec == 8) then
            mpireal = mpir8
         else
            mpireal = mpir4
         end if

         if ( dycore_is('LR') )then
# if ( defined STAGGERED )
! NEW LR CODING
            if (tape(t)%hlist(f)%hbuf_prec == 8) then
               select case (numlev)
               case (1)
                  if (myid_z .eq. 0) call pargatherreal(comm_y, 0,  &
                         tape(t)%hlist(f)%hbuf%buf8, strip2d, xzybuf%buf8)
               case (plev)
                  call pargatherreal(commglobal, 0, tape(t)%hlist(f)%hbuf%buf8, &
                         strip3dxzy, xzybuf%buf8)
               case (plevp)
                  call pargatherreal(commglobal, 0, tape(t)%hlist(f)%hbuf%buf8, &
                         strip3dxzyp, xzybuf%buf8)
               case default
                  write(6,*)'DUMP_FIELD: invalid number of levels=', numlev
                  call endrun ()
               end select
            else
               select case (numlev)
               case (1)
                  if (myid_z .eq. 0) call pargatherreal4(comm_y, 0,  &
                         tape(t)%hlist(f)%hbuf%buf4, strip2d, xzybuf%buf4)
               case (plev)
                  call pargatherreal4(commglobal, 0, tape(t)%hlist(f)%hbuf%buf4, &
                         strip3dxzy, xzybuf%buf4)
               case (plevp)
                  call pargatherreal4(commglobal, 0, tape(t)%hlist(f)%hbuf%buf4, &
                         strip3dxzyp, xzybuf%buf4)
               case default
                  write(6,*)'DUMP_FIELD: invalid number of levels=', numlev
                  call endrun ()
               end select
            endif
# endif
         else
            coldimin = tape(t)%hlist(f)%field%coldimin
            numowned = coldimin*numlev
            call compute_gsfactors (numowned, numsend, numrecv, displs)
            call mpigatherv_hbuf (tape(t)%hlist(f)%hbuf, numsend, mpireal, xzybuf, numrecv, &
                                  displs, mpireal, 0, mpicom)
         endif
#else
!
! Need to copy instead of pointer associate to keep phys_decomp and dyn_decomp consistent
!
         xzybuf = tape(t)%hlist(f)%hbuf       ! Overloaded assignment
#endif
      case default
         write(6,*) 'DUMP_FIELD: bad decomp_type=',tape(t)%hlist(f)%field%decomp_type
         call endrun
      end select

      if (masterproc) then
         dimind = dim_index_3d (1,plon,1,numlev,1,plat) ! xzybuf order
         if (tape(t)%hlist(f)%field%name == 'US') dimind = dim_index_3d (1,plon,1,numlev,2,plat)
         call xzy_to_xyz (xyzbuf, xzybuf, dimind)
 !
! Max/min fields got initialized to -huge and huge, respectively.  If fillvalue was in play
! these might still exist and need to be changed to fillvalue before output
!
         if (flag_xyfill .and. (avgflag == 'X' .or. avgflag == 'M')) then
             dimind = dim_index_3d (1,plon,1,plat,1,numlev)
             if (tape(t)%hlist(f)%field%name == 'US') dimind = &
                   dim_index_3d (1,plon,2,plat,1,numlev)
              print *,'flag_xyfill=',flag_xyfill,'avgflag=',avgflag
	      write(6,*)'stopping because this seems to cause trouble for some unknown reason'
	      call endrun()
!             call fill_unset (xyzbuf%buf4, xyzbuf%buf8, dimind)
         end if
!
!  Check whether we are writing out an entire field or partial field
!  If writing partial field determine whether it is a single column or group of contiguous columns
!
#ifdef HDEBUG
         if (ngroup(t) .gt. 0) then
               columnlon1=tape(t)%column(1)%columnlon(1)
               columnlat1=tape(t)%column(1)%columnlat(1)
	       write(6,*)'writing xyzbuf again'
              write(6,'(a,g12.2)')'xyzbuf=',xyzbuf%buf8(columnlon1,columnlat1,1)
         endif
#endif

         if (ngroup(t) .gt. 0) then
            do i=1,ngroup(t)
               ff=(f-1)*ngroup(t)+i     ! this is the variable id for the columm (or groups of columns)
#ifdef HDEBUG
               write(6,*)'DUMP_FIELD: writing column time indx ', &
                    nfils(t), ' field id', varid(ff,t), &
                    ' name ', trim(tape(t)%hlist(f)%field_column_name(i)), &
                    'numlons=',tape(t)%column(i)%num_lons,'numlats=', &
                    tape(t)%column(i)%num_lats,' base name ', &
                    trim(tape(t)%hlist(f)%field%name),' columnlon=', &
                    tape(t)%column(i)%columnlon,' columnlat=', &
                    tape(t)%column(i)%columnlat
#endif
               dimind = dim_index_3d (1,tape(t)%column(i)%num_lons,1,tape(t)%column(i)%num_lats,1,numlev)
               ! fvitt - for US field want tmpxyzbuf to have one less latitude ??
               if (tape(t)%hlist(f)%field%name == 'US')                   &
                 dimind = dim_index_3d (1,tape(t)%column(i)%num_lons,1,   &
                                          tape(t)%column(i)%num_lats-1,1, &
                                          numlev)
               call allocate_hbuf (tmpxyzbuf,dimind,tape(t)%hlist(f)%hbuf_prec)
               columnlon1=tape(t)%column(i)%columnlon(1)
               columnlon2=tape(t)%column(i)%columnlon(2)
               columnlat1=tape(t)%column(i)%columnlat(1)
               columnlat2=tape(t)%column(i)%columnlat(2)
#ifdef HDEBUG
               write(6,*)'dimind=',dimind
	       write(6,'(a,4g12.2)')'columnlon,columnlat',columnlon1,columnlon2,columnlat1,columnlat2
               write(6,'(a,g12.2)')'xyzbuf=',xyzbuf%buf8(columnlon1,columnlat1,1)
#endif

              if (tape(t)%hlist(f)%hbuf_prec == 8) then
                  tmpxyzbuf%buf8(1:tape(t)%column(i)%num_lons,1:tape(t)%column(i)%num_lats,:)= &
                       xyzbuf%buf8(columnlon1:columnlon2,columnlat1:columnlat2,:)
               else if(tape(t)%hlist(f)%hbuf_prec == 4) then
                  tmpxyzbuf%buf4(1:tape(t)%column(i)%num_lons,1:tape(t)%column(i)%num_lats,:)= &
                       xyzbuf%buf4(columnlon1:columnlon2,columnlat1:columnlat2,:)
               end if
               count(1)=tape(t)%column(i)%num_lons
               count(2)=tape(t)%column(i)%num_lats
#ifdef HDEBUG
               write(6,*)'DUMP_FIELD: writing column variable ', &
                    varid(ff,t),' start=',start,' count=',count
#endif
               call wrap_put_vara_hbuf (nfid(t), varid(ff,t), start, count, tmpxyzbuf)
               call deallocate_hbuf (tmpxyzbuf)
            end do
!
! else if writing full fields
!
         else 

#ifdef HDEBUG
            write(6,*)'DUMP_FIELD: writing full field for time indx ', nfils(t), ' field id', varid(f,t), &
                 ' name ', tape(t)%hlist(f)%field%name
#endif
            call wrap_put_vara_hbuf (nfid(t), varid(f,t), start, count, xyzbuf)
         end if

         call deallocate_hbuf (xzybuf)
         call deallocate_hbuf (xyzbuf)
      else
         call nullify_hbuf (xzybuf)
      end if

      return
   end subroutine dump_field

!#######################################################################

   subroutine wshist ()
!----------------------------------------------------------------------- 
! 
! Purpose: Driver routine to write fields on history tape t
! 
! Method: For variables which do not need to be gathered (SPMD) just issue the netcdf call
!         For those that do need to be gathered, call "dump_field" to do the operation.
!         Finally, zero the history buffers for each field written.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

      use time_manager, only: get_nstep, get_curr_date, get_curr_time
!
! Local workspace
!
      character(len=8) :: cdate  ! system date
      character(len=8) :: ctime  ! system time

      integer t, f               ! tape, field indices
      integer ret                ! return value from netcdf call
      integer start              ! starting index required by nf_put_vara
      integer count1             ! count values required by nf_put_vara
      integer startc(2)          ! start values required by nf_put_vara (character)
      integer countc(2)          ! count values required by nf_put_vara (character)
#ifdef HDEBUG
      integer begdim3
      integer enddim3
#endif
      
      integer :: yr, mon, day      ! year, month, and day components of a date
      integer :: nstep             ! current timestep number
      integer :: ncdate            ! current date in integer format [yyyymmdd]
      integer :: ncsec             ! current time of day [seconds]
      integer :: ndcur             ! day component of current time
      integer :: nscur             ! seconds component of current time
      real(r8) :: time             ! current time
      real(r8) :: tdata(2)         ! time interval boundaries
      character(len=nlen) :: fname ! Filename
      logical :: prev              ! Label file with previous date rather than current
!-----------------------------------------------------------------------

      nstep = get_nstep()
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      call get_curr_time(ndcur, nscur)
!
! Write time-varying portion of history file header
!
      do t=1,mtapes
         if (nhtfrq(t) == 0) then
            hstwr(t) = nstep /= 0 .and. day == 1 .and. ncsec == 0
            prev = .true.
         else
            hstwr(t) = mod(nstep,nhtfrq(t)) == 0
            prev = .false.
         end if

         if (hstwr(t)) then
            if (masterproc) then
               write(6,100) nfils(t),t,yr,mon,day,ncsec
100            format('WSHIST: writing time sample ',i3,' to h-file ', &
                      i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
               write(6,*)
!
! Starting a new volume => define the metadata
!
               if (nfils(t)==0) then
                  fname    = interpret_filename_spec( hfilename_spec(t), number=(t-1), &
                             prev=prev )
!
! Check that this new filename isn't the same as a previous or current filename
!
                  do f = 1, mtapes
                    if ( trim(fname) == trim(nhfil(f)) )then
                       write(6,*)'WSHIST: New filename same as old file = ', trim(fname)
                       write(6,*)'Is there an error in your filename specifiers?'
                       write(6,*)'hfilename_spec(', t, ') = ', hfilename_spec(t)
                       if ( t /= f )then
                          write(6,*)'hfilename_spec(', f, ') = ', hfilename_spec(f)
                       end if
                       call endrun
                    end if
                  end do
                  nhfil(t) = fname
                  write(6,*)'WSHIST: nhfil(',t,')=',trim(nhfil(t))
                  cpath(t) = trim(get_archivedir('hist')) // nhfil(t)
                  if ( len_trim(nfpath(t)) == 0 ) nfpath(t) = cpath(t)
                  call h_define (t)
               end if
            end if
            
            nfils(t) = nfils(t) + 1
            start = nfils(t)
            count1 = 1
      
            if (masterproc) then
               call wrap_put_vara_int (nfid(t), ndcurid(t),start, count1,ndcur)
               call wrap_put_vara_int (nfid(t), nscurid(t),start, count1,nscur)
               call wrap_put_vara_int (nfid(t), dateid(t),start, count1,ncdate)
               call wrap_put_vara_int (nfid(t), datesecid(t),start,count1,ncsec)
               call wrap_put_vara_int (nfid(t), nstephid(t),start, count1,nstep)
               time = ndcur + nscur/86400._r8
               call wrap_put_vara_realx (nfid(t), timeid(t), start, count1,time)

               startc(1) = 1
               startc(2) = nfils(t)
               countc(1) = 2
               countc(2) = 1
               tdata(1) = beg_time(t)
               tdata(2) = time
               call wrap_put_vara_realx (nfid(t), tbndid(t), startc, countc, tdata)
               beg_time(t) = time  ! update beginning time of next interval

               startc(1) = 1
               startc(2) = nfils(t)
               countc(1) = 8
               countc(2) = 1
               call datetime (cdate, ctime)
               call wrap_put_vara_text (nfid(t), date_writtenid(t), startc, countc, cdate)
               call wrap_put_vara_text (nfid(t), time_writtenid(t), startc, countc, ctime)
#ifdef HDEBUG
               ret = nf_sync (nfid(t))
#endif
            end if

!$OMP PARALLEL DO PRIVATE (F)

            do f=1,nflds(t)
               call h_normalize (f, t)  ! Normalized averaged fields and/or put fillvalue
            end do
!
! Write field to history tape.  Note that this is NOT threaded due to netcdf limitations
!
            call t_startf ('dump_field')
            do f=1,nflds(t)
#ifdef HDEBUG
               begdim3 = tape(t)%hlist(f)%field%begdim3
               enddim3 = tape(t)%hlist(f)%field%enddim3
               if (tape(t)%hlist(f)%hbuf_prec == 8) then
                  write(6,*)'WSHIST:',tape(t)%hlist(f)%field%name,'(1,1,1)=',tape(t)%hlist(f)%hbuf%buf8(1,1,begdim3),masterproc
               else
                  write(6,*)'WSHIST:',tape(t)%hlist(f)%field%name,'(1,1,1)=',tape(t)%hlist(f)%hbuf%buf4(1,1,begdim3),masterproc
               end if
#endif
               call dump_field (f, t)
            end do
            call t_stopf ('dump_field')
!
! Zero history buffers and accumulators now that the fields have been written.
!
! JR: recoded away from array syntax because pgf90 generated bad code. Found
! by compiling with -Mbounds.  Appeared to be only multitasking that caused
! problems.
!
!$OMP PARALLEL DO PRIVATE (F)

            do f=1,nflds(t)
               call h_zero (f, t)
            end do
      
         end if
      end do

   end subroutine wshist

!#######################################################################
#ifdef CRM
   subroutine write_inithist (crminitsave,crmsavechunks)
#else 
   subroutine write_inithist ()
#endif
!----------------------------------------------------------------------- 
! 
! Purpose: Write an initial dataset
! 
! Method: Make appropriate netcdf calls to put required fields on tape
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use pspect
      use comsrf
      use buffer
      use prognostics
      use rgrid
      use ioFileMod
      use commap
      use phys_buffer, only: pbuf, pbuf_next_tim_idx, pbuf_old_tim_idx, pbuf_get_fld_idx
#if ( defined SPMD )
      use mpishorthand
# if ( defined STAGGERED )
      use pmgrid, only: strip3dxzy, strip2d, myid_z
      use spmd_dyn, only: comm_y
      use parutilitiesmodule, only: commglobal, pargatherreal
# endif
#endif
      use dycore, only: dycore_is
      use phys_grid, only: gather_chunk_to_field
      use time_manager, only: get_step_size, get_nstep, get_ref_date, &
                              get_curr_date, get_curr_time

      use filenames,    only: caseid
#if ( ! defined COUP_SOM ) && ( ! defined COUP_CSM )
      use sst_data, only: sst
#endif
#ifdef CRM
      use crmdims, only: crm_nx,crm_ny,crm_nz
#endif

#include <comctl.h>
#include <comhyb.h>
!
! Local variables
!

#ifdef CRM
      logical, intent(in) :: crminitsave
      integer, intent(in) :: crmsavechunks
#endif
      character*80 name              ! variable name
      character*80 filename          ! filename to write
      character*256 path             ! mass store path
      character str*(max_chars)      ! string to hold netcdf attribute
      character*8 cdate              ! date
      character*8 ctime              ! time
      
      integer :: dtime               ! timestep size
      integer :: nstep               ! current timestep number
      integer :: ndbase = 0          ! days component of base time
      integer :: nsbase = 0          ! seconds component of base time
      integer :: nbdate              ! base date in integer format [yyyymmdd]
      integer :: nbsec               ! time of day relative to base date [seconds]
      integer :: ncdate              ! current date in integer format [yyyymmdd]
      integer :: ncsec               ! time of day relative to current date [seconds]
      integer :: ndcur               ! days component of current time
      integer :: nscur               ! seconds component of current time
      integer :: yr, mon, day        ! year, month, day components of a date

      real(r8) alon(plon)            ! longitude values
      real(r8) rlon(plon,plat)       ! reduced longitude values
      real(r8) alev(plev)            ! vertical levels
      real(r8) ailev(plevp)          ! vertical interfaces
      real(r8) gausslat(plat)        ! gaussian latitudes
      real(r8) pie                   ! 3.14...
      real(r8) time                  ! current time value
      real(r8) qpert_tmp (pcols,begchunk:endchunk)

      real(r8), allocatable :: arr2d(:,:)   ! 2d array to output
      real(r8), allocatable :: arr3d(:,:,:) ! 3d array to output
      real(r8), allocatable :: arr2dm(:,:,:)! multicomponent 2d array to output
      real(r8), allocatable :: arr2dchnk(:,:)   ! 2d array to output
      real(r8), allocatable :: arr2dmchnk(:,:,:)! multicomponent 2d array to output
!
! pointer to physics buffer fields
!
      integer itim, ifld
      real(r8), pointer, dimension(:,:,:) :: arr3dchnk_ptr

      integer i, m, j                ! indices
      integer ret                    ! return from netcdf lib call

      integer id                     ! netcdf file id
      integer old_mode               ! returned from nf_set_fill
      integer londim                 ! longitude dimension id
      integer latdim                 ! latitude dimension id
      integer levdim                 ! level dimension id
      integer ilevdim                ! interface dimension id
      integer timdim                 ! time dimension id
      integer dimenchar(2)           ! character dimension ids
      integer dimen1                 ! 1d dimension id
      integer dimen2(2)              ! 2d dimension ids
      integer dimen3(3)              ! 3d dimension ids
      integer dimen4(4)              ! 4d dimension ids
      integer ps0var                 ! ps0 variable id
      integer chardim                ! character dimension id
      integer lonvar                 !----------------------------------------------------
      integer rlonvar                ! 
      integer latvar                 ! dimension variable ids
      integer levvar                 ! 
      integer ilevvar                !----------------------------------------------------

      integer timeid                 !------------------------------------------------------
      integer date_writtenid         ! 
      integer time_writtenid         ! 
      integer ntrmid                 ! 
      integer ntrnid                 ! 
      integer ntrkid                 ! 
      integer ndbaseid               ! 
      integer nsbaseid               ! 
      integer nbdateid               ! 
      integer nbsecid                ! 
      integer mdtid                  ! 
      integer hyaiid                 !  
      integer hybiid                 ! 
      integer hyamid                 ! 
      integer hybmid                 ! 
      integer tracid(pcnst+pnats)    ! netcdf variable ids for variables of the given name
      integer gwid                   ! with "id" tacked on to the end
      integer ndcurid                ! 
      integer nscurid                ! 
      integer dateid                 ! 
      integer datesecid              ! 
      integer nstephid               ! 
      integer uid                    ! 
      integer vid                    ! 
      integer tid                    ! 
      integer pblhtid                ! 
      integer tpertid                ! 
      integer qpertid                ! 
      integer cldid                  ! 
      integer qcwatid                ! 
      integer tcwatid                ! 
      integer lcwatid                ! 
      integer qid                    ! 
      integer psid                   ! 
      integer phisid                 ! 
      integer sghid                  ! 
      integer landmid                ! 
      integer tsid                   ! 
      integer tsiceid                ! 
      integer tsice_rad_id           ! 
      integer tbotid                 ! 
      integer tssubid(4)             ! 
      integer sicthkid               !
      integer snowhlndid             !
      integer snowhiceid             !
      integer landfracid             ! land-fraction id
      integer icefracid              ! ice-fraction id
      integer tsocnid
      integer nlonid, wnummaxid      !------------------------------------------------------

      integer start(1)               ! starting index array for nf_put_vara (1d vars)
      integer count(1)               ! count array for nf_put_vara (1d vars)
      data start/1/
      data count/1/

      integer start2d(3)             ! starting index array for nf_put_vara (2d vars)
      integer count2d(3)             ! count array for nf_put_vara (2d vars)
      data start2d/1,1,1/
      data count2d/plon,plat,1/

      integer start3d(4)             ! starting index array for nf_put_vara (3d vars)
      integer count3d(4)             ! count array for nf_put_vara (3d vars)
      data start3d/1,1,1,1/
      data count3d/plon,plev,plat,1/

      logical haslon                 ! flag indicates longitude dimension is present
      logical haslat                 ! flag indicates latitude dimension is present

      integer n, nn                  ! indices
      integer nvars                  ! number of variables returned from nf_inq_nvars
      integer ndims                  ! number of dimensions returned from nf_inq_var
      integer natts                  ! number of attributes
      integer xtype                  ! variable type flag (netcdf)
      integer dimids(nf_max_var_dims)! dimension ids

#ifdef STAGGERED
      real(r8) slons(splon)
      real(r8) slats(plat-1)

      integer slonvar                ! Staggered longitude variable
      integer slatvar                ! Staggered latitude variable
      integer usid                   ! Staggered latidude wind variable id
      integer vsid                   !     "     longitude "      "     "
      integer slatdim                ! Dimension of the staggered lat array
      integer slondim                !     "     "   "     "      lon   "
      integer dimen4us(4)            ! NetCDF array holding array dimensions
      integer dimen4vs(4)            !    "
      integer count3dus(4)           !    "
      integer count3dvs(4)           !    "
      integer wsid                   ! Staggered latitude weight ID

      data count3dvs/splon,plev,plat,1/

      real(r8), allocatable :: arr3dxyz_local(:,:,:) ! help array
      real(r8), allocatable :: arr3dxzy_local(:,:,:) ! help array
      real(r8), allocatable :: arr3dxyz(:,:,:) ! help array
      real(r8), allocatable :: arr3dxzy(:,:,:) ! help array
      real(r8), allocatable :: arr3dus(:,:,:)
      real(r8), allocatable :: arr2dxy_local(:,:) ! help array
#endif

      integer startc(2)              ! starting index array for nf_put_vara (character)
      integer countc(2)              ! count array for nf_put_vara (character)

      integer k                      ! level index
      integer endi                   ! ending longitude index
      integer jcen                   ! latitude index (extended array)
      integer numi                   ! number of longitude values
      integer numperlat              ! number of values per latitude band
#ifdef CRM
      character(1) numb
      integer :: crm3dchunkdim,crm2dchunkdim
      integer :: dimen4crmvar3d(4),dimen4crmvar2d(4)
      integer :: dimen4crmradbfr(4)
      integer :: ucrmid(crmsavechunks), vcrmid(crmsavechunks), wcrmid(crmsavechunks)
      integer :: tcrmid(crmsavechunks), qcrmid(crmsavechunks)
      integer :: qncrmid(crmsavechunks), qpcrmid(crmsavechunks)
      integer :: qrlcrmid(crmsavechunks),qrscrmid(crmsavechunks)

      integer :: fsdscrmid,fsnscrmid,fsdsccrmid
      integer :: fsntoacrmid,fsntoaccrmid
      integer :: fsutoacrmid,flwdscrmid,flnscrmid
      integer :: flutccrmid,flnsccrmid,flutcrmid

      integer :: crmradbfrid,qrs1id,qrl1id,radbufdim
      real(r8), allocatable :: arr3dcrmradbfr(:,:,:) ! 3d array to output
      integer start3dcrmradbfr(4)             
      integer count3dcrmradbfr(4)            
      data start3dcrmradbfr/1,1,1,1/
      data count3dcrmradbfr/plon,nrad_buffer,plat,1/


#endif
!-----------------------------------------------------------------------


#ifdef STAGGERED
!
! Set count3dus in executable code since "plat-1" cannot be computed in a data stmt
!
      count3dus(1) = plon
      count3dus(2) = plev
      count3dus(3) = plat - 1
      count3dus(4) = 1
#endif
      if (masterproc) then

         dtime = get_step_size()
         nstep = get_nstep()
         call get_ref_date(yr, mon, day, nbsec)
         nbdate = yr*10000 + mon*100 + day
         call get_curr_date(yr, mon, day, ncsec)
         if (dycore_is ('SLD') ) then
            call get_curr_date(yr, mon, day, ncsec,offset=dtime)
         endif
         ncdate = yr*10000 + mon*100 + day
         call get_curr_time(ndcur, nscur)

         write(6,*)'WRITE_INITHIST: ncdate=',ncdate
         filename = interpret_filename_spec( ifilename_spec )
         path = trim(get_archivedir( 'init' )) // trim(filename)
         
         ret = nf_create (filename, nf_clobber, id)
         if (ret /= nf_noerr) call handle_error (ret)
         
         ret = nf_set_fill (id, nf_nofill, old_mode)
         ret = nf_def_dim (id, 'lat', plat, latdim)
         ret = nf_def_dim (id, 'lon', plon, londim)
         
#ifdef STAGGERED
         ret = nf_def_dim (id, 'slat', plat-1, slatdim)
         ret = nf_def_dim (id, 'slon', splon, slondim)
#endif

         ret = nf_def_dim (id, 'lev', plev, levdim)
         ret = nf_def_dim (id, 'ilev', plevp, ilevdim)
         ret = nf_def_dim (id, 'time', nf_unlimited, timdim)
         ret = nf_def_dim (id, 'chars', 8, chardim)

#ifdef CRM
         if (crminitsave) then
           ret = nf_def_dim (id,'crm3dchunk', crm_nx*crm_ny*crm_nz/crmsavechunks,crm3dchunkdim)
           ret = nf_def_dim (id,'radbufvars',nrad_buffer,radbufdim)
           ret = nf_def_dim (id,'crm2dchunk', crm_nx*crm_ny,crm2dchunkdim)
!         ret = nf_def_dim (id,'crmy', crm_ny,crmydim)
!         ret = nf_def_dim (id,'crmz', crm_nz/crmsavechunks,crmzchunkdim)
         end if
#endif
! setup dimension arrays for 1,2,3,4d variables 
!     
         dimen1 = timdim

         dimen2(1) = londim
         dimen2(2) = latdim

         dimen3(1) = londim
         dimen3(2) = latdim
         dimen3(3) = timdim
!            
         dimen4(1) = londim
         dimen4(2) = levdim
         dimen4(3) = latdim
         dimen4(4) = timdim

#ifdef STAGGERED
         dimen4us(1) = londim
         dimen4us(2) = levdim
         dimen4us(3) = slatdim
         dimen4us(4) = timdim

         dimen4vs(1) = slondim
         dimen4vs(2) = levdim
         dimen4vs(3) = latdim
         dimen4vs(4) = timdim
#endif


#ifdef CRM
         if (crminitsave) then
           dimen4crmvar3d(1) = londim
           dimen4crmvar3d(2) = crm3dchunkdim
           dimen4crmvar3d(3) = latdim
           dimen4crmvar3d(4) = timdim

           dimen4crmvar2d(1) = londim
           dimen4crmvar2d(2) = crm2dchunkdim
           dimen4crmvar2d(3) = latdim
           dimen4crmvar2d(4) = timdim

           dimen4crmradbfr(1) = londim
           dimen4crmradbfr(2) = radbufdim
           dimen4crmradbfr(3) = latdim
           dimen4crmradbfr(4) = timdim
         end if

#endif
!           
! define variables to label the dimensions, use the same names as the
! dimensions
!
         call wrap_def_var (id,'P0',nf_double,0,0,ps0var)
         str = 'reference pressure'
         call wrap_put_att_text (id, ps0var, 'long_name', str)
         call wrap_put_att_text (id, ps0var, 'units', 'Pa')

         call wrap_def_var (id,'lat',nf_double,1,latdim,latvar)
         call wrap_put_att_text (id, latvar, 'long_name', 'latitude')
         call wrap_put_att_text (id, latvar, 'units', 'degrees_north')

         if (fullgrid) then

            call wrap_def_var (id,'lon',nf_double,1,londim,lonvar)
            call wrap_put_att_text (id, lonvar, 'long_name', 'longitude')
            call wrap_put_att_text (id, lonvar, 'units', 'degrees_east')

         else

            call wrap_def_var (id,'rlon',nf_double,2,dimen2,rlonvar)
            call wrap_put_att_text (id, rlonvar, 'long_name', 'reduced_longitude')
            call wrap_put_att_text (id, rlonvar, 'units', 'degrees_east')

         end if

! If a staggered grid is in use, output the lat and lon arrays variables.
! Note: the staggered grid is currently independent of any reduced grid
! information.

#ifdef STAGGERED
         call wrap_def_var (id, 'slat', nf_double, 1,slatdim, slatvar)
         call wrap_put_att_text (id, slatvar, 'long_name','staggered latitude')
         call wrap_put_att_text(id, slatvar, 'units', 'degrees_north')

         call wrap_def_var(id, 'slon', nf_double, 1, slondim, slonvar)
         call wrap_put_att_text (id, slonvar, 'long_name', 'staggered longitude')
         call wrap_put_att_text (id, slonvar, 'units', 'degrees_east')

         call wrap_def_var (id,'w_stag',nf_double,1,slatdim,wsid)
         str = 'staggered latitude weights'
         call wrap_put_att_text (id, wsid, 'long_name', str)
#endif

         call wrap_def_var (id,'lev',nf_double,1,levdim,levvar)
         str = 'hybrid level at midpoints (1000*(A+B))'
         call wrap_put_att_text (id, levvar, 'long_name', str)
         str = 'level'
         call wrap_put_att_text (id, levvar, 'units', str)
         call wrap_put_att_text (id, levvar, 'positive', 'down')
         call wrap_put_att_text (id, levvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
         call wrap_put_att_text (id, levvar, 'formula_terms', 'a: hyam b: hybm p0: P0 ps: PS')

         call wrap_def_var (id,'ilev',nf_double,1,ilevdim,ilevvar)
         str = 'hybrid level at interfaces (1000*(A+B))'
         call wrap_put_att_text (id, ilevvar, 'long_name', str)
         str = 'level'
         call wrap_put_att_text (id, ilevvar, 'units', str)
         call wrap_put_att_text (id, ilevvar, 'positive', 'down')
         call wrap_put_att_text (id, ilevvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
         call wrap_put_att_text (id, ilevvar, 'formula_terms', 'a: hyai b: hybi p0: P0 ps: PS')

         call wrap_def_var (id,'time',nf_double,1,timdim,timeid)
         call wrap_put_att_text (id, timeid, 'long_name', 'time')
         str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
         call wrap_put_att_text (id, timeid, 'units', str)
         call wrap_put_att_text (id, timeid, 'calendar', 'noleap')
!
! Character
!
         dimenchar(1) = chardim
         dimenchar(2) = timdim
         call wrap_def_var (id,'date_written',nf_char,2,dimenchar, date_writtenid)
         call wrap_def_var (id,'time_written',nf_char,2,dimenchar, time_writtenid)
!
! Integer Header
!
         call wrap_def_var (id,'ntrm',nf_int,0,0,ntrmid)
         str = 'spectral truncation parameter M'
         call wrap_put_att_text (id, ntrmid, 'long_name', str)

         call wrap_def_var (id,'ntrn',nf_int,0,0,ntrnid)
         str = 'spectral truncation parameter N'
         call wrap_put_att_text (id, ntrnid, 'long_name', str)

         call wrap_def_var (id,'ntrk',nf_int,0,0,ntrkid)
         str = 'spectral truncation parameter K'
         call wrap_put_att_text (id, ntrkid, 'long_name', str)

         call wrap_def_var (id,'ndbase',nf_int,0,0,ndbaseid)
         str = 'base day'
         call wrap_put_att_text (id, ndbaseid, 'long_name', str)

         call wrap_def_var (id,'nsbase',nf_int,0,0,nsbaseid)
         str = 'seconds of base day'
         call wrap_put_att_text (id, nsbaseid, 'long_name', str)

         call wrap_def_var (id,'nbdate',nf_int,0,0,nbdateid)
         str = 'base date (YYYYMMDD)'
         call wrap_put_att_text (id, nbdateid, 'long_name', str)

         call wrap_def_var (id,'nbsec',nf_int,0,0,nbsecid)
         str = 'seconds of base date'
         call wrap_put_att_text (id, nbsecid, 'long_name', str)

         call wrap_def_var (id,'mdt',nf_int,0,0,mdtid)
         call wrap_put_att_text (id, mdtid, 'long_name', 'timestep')
         call wrap_put_att_text (id, mdtid, 'units', 's')

         if (.not. fullgrid) then
            call wrap_def_var (id,'nlon',nf_int,1,latdim,nlonid)
            str = 'number of longitudes'
            call wrap_put_att_text (id, nlonid, 'long_name', str)

            call wrap_def_var (id,'wnummax',nf_int,1,latdim,wnummaxid)
            str = 'cutoff Fourier wavenumber'
            call wrap_put_att_text (id, wnummaxid, 'long_name', str)
         end if
!
! Floating point time-invariant
!
         call wrap_def_var (id,'hyai',nf_double,1,ilevdim,hyaiid)
         str = 'hybrid A coefficient at layer interfaces'
         call wrap_put_att_text (id, hyaiid, 'long_name', str)

         call wrap_def_var (id,'hybi',nf_double,1,ilevdim,hybiid)
         str = 'hybrid B coefficient at layer interfaces'
         call wrap_put_att_text (id, hybiid, 'long_name', str)

         call wrap_def_var (id,'hyam',nf_double,1,levdim,hyamid)
         str = 'hybrid A coefficient at layer midpoints'
         call wrap_put_att_text (id, hyamid, 'long_name', str)

         call wrap_def_var (id,'hybm',nf_double,1,levdim,hybmid)
         str = 'hybrid B coefficient at layer midpoints'
         call wrap_put_att_text (id, hybmid, 'long_name', str)

         call wrap_def_var (id,'gw',nf_double,1,latdim,gwid)
         str = 'gauss weights'
         call wrap_put_att_text (id, gwid, 'long_name', str)
!     
! Character header information 
!
         str = 'CF-1.0'
         call wrap_put_att_text (id, nf_global, 'Conventions', str)
         call wrap_put_att_text (id, nf_global, 'source', 'CAM')
         call wrap_put_att_text (id, nf_global, 'case',caseid)
         call wrap_put_att_text (id, nf_global, 'title',ctitle)
         call wrap_put_att_text (id, nf_global, 'logname',logname)
         call wrap_put_att_text (id, nf_global, 'host', host)
!
! Create variables for model timing and header information 
!
         call wrap_def_var (id,'ndcur   ',nf_int,1,dimen1,ndcurid)
         str = 'current day (from base day)'
         call wrap_put_att_text (id, ndcurid, 'long_name', str)

         call wrap_def_var (id,'nscur   ',nf_int,1,dimen1,nscurid)
         str = 'current seconds of current day'
         call wrap_put_att_text (id, nscurid, 'long_name', str)

         call wrap_def_var (id,'date    ',nf_int,1,dimen1,dateid)
         str = 'current date (YYYYMMDD)'
         call wrap_put_att_text (id, dateid, 'long_name', str)

         call wrap_def_var (id,'datesec ',nf_int,1,dimen1, datesecid)
         str = 'current seconds of current date'
         call wrap_put_att_text (id, datesecid, 'long_name', str)

         call wrap_def_var (id,'nsteph  ',nf_int,1,dimen1,nstephid)
         str = 'current timestep'
         call wrap_put_att_text (id, nstephid, 'long_name', str)
!
! Call wrapper routine to add attributes as well as define variable for 
! fields which are in master field list
!

#ifdef STAGGERED
         call addvar (id,'US      ',nf_double,4,dimen4us,usid)
         call addvar (id,'VS      ',nf_double,4,dimen4vs,vsid)
#else
         call addvar (id,'U       ',nf_double,4,dimen4,uid)
         call addvar (id,'V       ',nf_double,4,dimen4,vid)
#endif
         call addvar (id,'T       ',nf_double,4,dimen4,tid)
         call addvar (id,'Q       ',nf_double,4,dimen4,qid)
         call addvar (id,'PS      ',nf_double,3,dimen3,psid)
         call addvar (id,'PHIS    ',nf_double,3,dimen3,phisid)
         call addvar (id,'SGH     ',nf_double,3,dimen3,sghid)
         call addvar (id,'LANDM_COSLAT',nf_double,3,dimen3,landmid)
         call addvar (id,'PBLH    ',nf_double,3,dimen3,pblhtid)
         call addvar (id,'TPERT   ',nf_double,3,dimen3,tpertid)
         call addvar (id,'QPERT   ',nf_double,3,dimen3,qpertid)
         call addvar (id,'CLOUD   ',nf_double,4,dimen4,cldid)
         call addvar (id,'QCWAT   ',nf_double,4,dimen4,qcwatid)
         call addvar (id,'TCWAT   ',nf_double,4,dimen4,tcwatid)
         call addvar (id,'LCWAT   ',nf_double,4,dimen4,lcwatid)

#if ( ! defined COUP_CSM )
         call addvar (id,'TSICERAD',nf_double,3,dimen3,tsice_rad_id)
         call addvar (id,'TS      ',nf_double,3,dimen3,tsid)
         call addvar (id,'TSICE   ',nf_double,3,dimen3,tsiceid)
         call addvar (id,'TS1     ',nf_double,3,dimen3,tssubid(1))
         call addvar (id,'TS2     ',nf_double,3,dimen3,tssubid(2))
         call addvar (id,'TS3     ',nf_double,3,dimen3,tssubid(3))
         call addvar (id,'TS4     ',nf_double,3,dimen3,tssubid(4))
         call addvar (id,'SNOWHICE',nf_double,3,dimen3,snowhiceid)
         call addvar (id,'LANDFRAC',nf_double,3,dimen3,landfracid)
         call addvar (id,'TBOT'    ,nf_double,3,dimen3,tbotid)
         call addvar (id,'ICEFRAC ',nf_double,3,dimen3,icefracid)
         call addvar (id,'SICTHK  ',nf_double,3,dimen3,sicthkid)
         call addvar (id,'TSOCN   ',nf_double,3,dimen3,tsocnid)
#endif

#ifdef CRM
   if(crminitsave) then
    do m=1,crmsavechunks
     write(numb,'(i1)') m
     numb = adjustl(numb)
     call addvar (id,'U_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,ucrmid(m))
     call addvar (id,'V_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,vcrmid(m))
     call addvar (id,'W_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,wcrmid(m))
     call addvar (id,'T_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,tcrmid(m))
     call addvar (id,'Q_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,qcrmid(m))
     call addvar (id,'QN_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,qncrmid(m))
     call addvar (id,'QP_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,qpcrmid(m))
     call addvar (id,'QRS_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,qrscrmid(m))
     call addvar (id,'QRL_CRM'//trim(numb),nf_double,4,dimen4crmvar3d,qrlcrmid(m))
    ! Chunking happens over the CRM's z dimension so for 2D crm variables don't need  separate variables:
     call addvar (id,'FSDS_CRM',nf_double,4,dimen4crmvar2d,fsdscrmid)
     call addvar (id,'FSNS_CRM',nf_double,4,dimen4crmvar2d,fsnscrmid)
     call addvar (id,'FSDSC_CRM',nf_double,4,dimen4crmvar2d,fsdsccrmid)
     call addvar (id,'FSNTOA_CRM',nf_double,4,dimen4crmvar2d,fsntoacrmid)
     call addvar (id,'FSNTOAC_CRM',nf_double,4,dimen4crmvar2d,fsntoaccrmid)
     call addvar (id,'FSUTOA_CRM',nf_double,4,dimen4crmvar2d,fsutoacrmid)
     call addvar (id,'FLWDS_CRM',nf_double,4,dimen4crmvar2d,flwdscrmid)
     call addvar (id,'FLNS_CRM',nf_double,4,dimen4crmvar2d,flnscrmid)
     call addvar (id,'FLNSC_CRM',nf_double,4,dimen4crmvar2d,flnsccrmid)
     call addvar (id,'FLUT_CRM',nf_double,4,dimen4crmvar2d,flutcrmid)
     call addvar (id,'FLUTC_CRM',nf_double,4,dimen4crmvar2d,flutccrmid)

    end do
! Note that physics/read_restart_physics has CRM logic that suggests we also need these variables on the initial conditions tape:
     call addvar (id,'CRMRADBFR',nf_double,4,dimen4crmradbfr,crmradbfrid)
     call addvar (id,'QRS1',nf_double,4,dimen4,qrs1id) ! only 80% sure dimen4 works here.
     call addvar (id,'QRL1',nf_double,4,dimen4,qrl1id)

   end if

#endif



!
! Advected and non-advected constituents.
!
         do m=2,pcnst+pnats
            call addvar (id,cnst_name(m),nf_double,4,dimen4,tracid(m))
         end do
!
! Add fillvalue attribute to arrays which have lat and lon dimensions
!
         if (nf_inq_nvars (id, nvars) /= nf_noerr) then
            write(6,*)'WRITE_INITHIST: bad call to nf_inq_nvars'
            call endrun
         end if

         do n=1,nvars
            call wrap_inq_var (id, n, name, xtype, ndims, dimids, natts)

            haslon = .false.
            haslat = .false.

            do nn=1,ndims
               if (dimids(nn) == londim) then
                  haslon = .true.
               else if (dimids(nn) == latdim) then
                  haslat = .true.
               end if

#ifdef STAGGERED
               if (dimids(nn) == slondim) then
                  haslon = .true.
               else if (dimids(nn) == slatdim) then
                  haslat = .true.
               end if
#endif

            end do

            if (haslon .and. haslat .and. .not. fullgrid) then
               call wrap_put_att_realx (id, n, '_FillValue', nf_double, 1, fillvalue)
            end if
         end do

         ret = nf_enddef(id)
         write(6,*)'WRITE_INITHIST: Successfully opened netcdf file'
!
! Write time-invariant portion of history header
!
         call wrap_put_var_realx (id, ps0var, ps0)

         ret = nf_put_vara_int (id, ndcurid,start, count,ndcur)
         ret = nf_put_vara_int (id, nscurid,start, count,nscur)
         ret = nf_put_vara_int (id, dateid,start, count,ncdate)
         ret = nf_put_vara_int (id, datesecid,start, count,ncsec)
         ret = nf_put_vara_int (id, nstephid,start, count,nstep)
         time = ndcur + nscur/86400._r8
         call wrap_put_vara_realx (id, timeid, start, count, time)

         startc(1) = 1
         startc(2) = 1
         countc(1) = 8
         countc(2) = 1
         call datetime (cdate, ctime)
         ret = nf_put_vara_text (id, date_writtenid, startc, countc, cdate)
         ret = nf_put_vara_text (id, time_writtenid, startc, countc, ctime)

         pie = 4.*atan(1.)
         do j=1,plat
            gausslat(j) = (180./pie)*clat(j)
         end do
         call wrap_put_var_realx (id, latvar, gausslat)
!
! Coordinate var defined only in full grid: otherwise define 2-d var.
!
         if (fullgrid) then
            do i=1,plon
               alon(i) = (i-1) * 360.0 / plon
            end do
            call wrap_put_var_realx (id, lonvar, alon)

         else

            do j=1,plat
               rlon(:,j) = fillvalue
               do i=1,nlon(j)
                  rlon(i,j) = (i-1) * 360.0 / nlon(j)
               end do
            end do
            call wrap_put_var_realx (id, rlonvar, rlon)

            call wrap_put_var_int (id, nlonid, nlon)
            call wrap_put_var_int (id, wnummaxid, wnummax)
         end if

#ifdef STAGGERED
! Calculate the staggered grid latitudes.
! Staggered latitudes fall in between the regular latitudes. Note that
! for Lin-Rood, gausslat latitudes are not gaussian!

         do j = 1, plat-1
            slats(j) = (gausslat(j) + gausslat(j+1))/2.
         end do

! Sanity check.

         if (slats(plat-1) .gt. 90.0) then
            write(6,*)"WRITE_INITHIST: WARNING - last staggered grid latitude"
            write(6,*) "was calculated to be:", slats(plat-1)
            write(6,*) "Point has been redefined to be 90.0"
            slats(plat-1) = 90.0
         endif

         call wrap_put_var_realx (id, slatvar, slats)

         do i = 1, splon
            slons(i) = ((i-1) * 360.0 / splon) - 180.0/splon
         enddo

         call wrap_put_var_realx (id, slonvar, slons)
         call wrap_put_var_realx (id, wsid, w_staggered)
#endif

!
! 0.01 converts Pascals to millibars
!
         alev(:plev) = 0.01*ps0*(hyam(:plev) + hybm(:plev))
         ailev(:plevp) = 0.01*ps0*(hyai(:plevp) + hybi(:plevp))

         call wrap_put_var_realx (id, levvar, alev)
         call wrap_put_var_realx (id, ilevvar, ailev)
         call wrap_put_var_realx (id, hyaiid, hyai)
         call wrap_put_var_realx (id, hybiid, hybi)
         call wrap_put_var_realx (id, hyamid, hyam)
         call wrap_put_var_realx (id, hybmid, hybm)
         call wrap_put_var_realx (id, gwid, w)

         call wrap_put_var_int (id, ntrmid, ptrm)
         call wrap_put_var_int (id, ntrnid, ptrn)
         call wrap_put_var_int (id, ntrkid, ptrk)
         call wrap_put_var_int (id, mdtid, dtime)
!
! Model date info
!
         call wrap_put_var_int (id, ndbaseid, ndbase)
         call wrap_put_var_int (id, nsbaseid, nsbase)
         call wrap_put_var_int (id, nbdateid, nbdate)
         call wrap_put_var_int (id, nbsecid, nbsec)

      end if  ! end of if-masterproc
!
! Allocate space for arrays to be written to initial history file.  In SPMD mode,
! masterproc needs all plat latitudes for writing
!
! If the staggered arrays are a different size than the regular arrays,
! yer gonna need something like "begsplat" and "endsplat" to define the
! portions used for MPI output below. -gg

      if (masterproc) then
         allocate (arr2d(plon,plat))
         allocate (arr3d(plon,plev,plat))
         allocate (arr2dchnk(pcols,begchunk:endchunk))
         allocate (arr2dmchnk(pcols,plevmx,begchunk:endchunk))
         allocate (arr2dm(plon,plevmx,plat))
#ifdef CRM
         if (crminitsave) then
           allocate (arr3dcrmradbfr(plon,nrad_buffer,plat))
         endif
#endif
#ifdef STAGGERED
         allocate (arr3dxyz_local(plon,beglat:endlat,beglev:endlev))
         allocate (arr3dxzy_local(plon,beglev:endlev,beglat:endlat))
         allocate (arr3dxyz(plon,plat,plev))
         allocate (arr3dxzy(plon,plev,plat))
         allocate (arr3dus(plon,plev,plat-1))
         allocate (arr2dxy_local(plon,beglat:endlat))
#endif
      else
         allocate (arr2d(plon,beglat:endlat))
         allocate (arr3d(plon,plev,beglat:endlat))
         allocate (arr2dchnk(pcols,begchunk:endchunk))
         allocate (arr2dmchnk(pcols,plevmx,begchunk:endchunk))
         allocate (arr2dm(1,1,1))
#ifdef CRM
         if (crminitsave) then
           allocate (arr3dcrmradbfr(pcols,nrad_buffer,begchunk:endchunk))
         endif
#endif
#ifdef STAGGERED
         allocate (arr3dxyz_local(plon,beglat:endlat,beglev:endlev))
         allocate (arr3dxzy_local(plon,beglev:endlev,beglat:endlat))
         allocate (arr3dxyz(1,1,1))
         allocate (arr3dxzy(1,1,1))
         allocate (arr3dus(1,1,1))
         allocate (arr2dxy_local(plon,beglat:endlat))
#endif

      end if

!----------------------------------------------------------------------------
! output variables
! 3-d
!----------------------------------------------------------------------------

      numperlat = plnlv
      arr3d = fillvalue

#ifdef STAGGERED
!
! output u
!
!
! WS 2001.01.29 :  Fixed Glenn's "DANGER, WILL ROBINSON! DANGER!" problem
! was:  call gather_put (id, usid, arr3dus, num, start3d, count3dus)
! Temporary: requires better solution
!

#if defined( SPMD )
!
!  Copy U3S (ghosted) to arr3dxyz_local (hack)
!
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               arr3dxyz_local(i,j,k) = u3s(i,j,k)
            enddo
         enddo
      enddo
      call gather( arr3dxyz_local, strip3dxyz, arr3dxyz )
#else
      do k=1,plev
         do j=1,plat
            do i=1,plon
               arr3dxyz(i,j,k) = u3s(i,j,k)
            enddo
         enddo
      enddo
#endif

      if (masterproc) then
         do j = 1, plat-1         ! shift the indices (what a mess...)
            do k = 1, plev
               do i = 1, plon
                  arr3dus(i,k,j) = arr3dxyz(i,j+1,k)
               enddo
            enddo
         enddo
         call wrap_put_vara_realx (id, usid, start3d, count3dus, arr3dus)
      end if

!
! output v
!
      do j = beglat, endlat
         do k = beglev, endlev
            arr3dxzy_local(:,k,j) = v3s(:,j,k)
         enddo
      enddo
      
#if defined (SPMD)
      call pargatherreal(commglobal, 0, arr3dxzy_local, strip3dxzy, arr3dxzy)
#else
      arr3dxzy(:,:,:) = arr3dxzy_local(:,:,:)
#endif
      if (masterproc) then
         call wrap_put_vara_realx (id, vsid, start3d, count3dvs, arr3dxzy)
      end if

#else

      arr3d(:,:,:) = fillvalue
      arr2d(:,:) = fillvalue

!
! output u
!
      do j=beglat,endlat
         endi = i1+nlon(j)-1
         numi = nlon(j)
         jcen = j1 - 1 + j
         arr3d(:numi,:plev,j) = u3(i1:endi,:plev,jcen,n3m1)
      end do
      call gather_put (id, uid, arr3d, numperlat, start3d, count3d)

!
! output v
!
      do j=beglat,endlat
         endi = i1+nlon(j)-1
         numi = nlon(j)
         jcen = j1 - 1 + j
         arr3d(:numi,:plev,j) = v3(i1:endi,:plev,jcen,n3m1)
      end do
      call gather_put (id, vid, arr3d, numperlat, start3d, count3d)

#endif

!
! output t
!
#ifdef STAGGERED
      do j=beglat,endlat
         do k=beglev,endlev
            arr3dxzy_local(:plon,k,j) = t3(:plon,k,j)
         enddo
      end do

#if defined (SPMD)
      call pargatherreal(commglobal, 0, arr3dxzy_local, strip3dxzy, arr3dxzy)
#else
      arr3dxzy(:,:,:) = arr3dxzy_local(:,:,:)
#endif

      if (masterproc) then
         call wrap_put_vara_realx (id, tid, start3d, count3d, arr3dxzy)
      end if
#else

      do j=beglat,endlat
         endi = i1+nlon(j)-1
         numi = nlon(j)
         jcen = j1 - 1 + j
         arr3d(:numi,:plev,j) = t3(i1:endi,:plev,jcen,n3m1)
      end do
      call gather_put (id, tid, arr3d, numperlat, start3d, count3d)
#endif

!
! output q (first component)
!
#ifdef STAGGERED
      do j=beglat,endlat
         if ( dycore_is('LR') )then
            do k=beglev,endlev
               do i=1,plon
                  arr3dxzy_local(i,k,j) = q3(i,j,k,1)
               enddo
            enddo
         else
            arr3dxzy_local(:plon,:plev,j) = q3(:plon,:plev,1,j)
         end if
      end do

#if defined (SPMD)
      call pargatherreal(commglobal, 0, arr3dxzy_local, strip3dxzy, arr3dxzy)
#else
      arr3dxzy(:,:,:) = arr3dxzy_local(:,:,:)
#endif

      if (masterproc) then
         call wrap_put_vara_realx (id, qid, start3d, count3d, arr3dxzy)
      end if
#else

      do j=beglat,endlat
         endi = i1+nlon(j)-1
         numi = nlon(j)
         jcen = j1 - 1 + j
         arr3d(:numi,:plev,j) = q3(i1:endi,:plev,1,jcen,n3m1)
      end do
      call gather_put (id, qid, arr3d, numperlat, start3d, count3d)
#endif

!
! output q (other components)
!
#ifdef STAGGERED
      do m=2,pcnst+pnats
         do j=beglat,endlat
            if ( dycore_is('LR') )then
               do k=beglev,endlev
                  do i=1,plon
                     arr3dxzy_local(i,k,j) = q3(i,j,k,m)
                  enddo
               enddo
            else
               arr3dxzy_local(:plon,:plev,j) = q3(:plon,:plev,m,j)
            endif
         end do

#if defined (SPMD)
         call pargatherreal(commglobal, 0, arr3dxzy_local, strip3dxzy, arr3dxzy)
#else
         arr3dxzy(:,:,:) = arr3dxzy_local(:,:,:)
#endif

         if (masterproc) then
            call wrap_put_vara_realx (id, tracid(m), start3d, count3d, arr3dxzy)
         end if
      end do
#else
      do m=2,pcnst+pnats
         do j=beglat,endlat
            endi = i1+nlon(j)-1
            numi = nlon(j)
            jcen = j1 - 1 + j
            arr3d(:numi,:plev,j) = q3(i1:endi,:plev,m,jcen,n3m1)
         end do
         call gather_put (id, tracid(m), arr3d, numperlat, start3d, count3d)
      end do
#endif

!
! 2-d fields
!
      numperlat = plon
      arr2d = fillvalue

#ifdef STAGGERED
      do j=beglat,endlat
         arr2dxy_local(:plon,j) = ps(:plon,j)
      end do

#if defined (SPMD)
      if (myid_z .eq. 0) call pargatherreal(comm_y, 0,      &
                  arr2dxy_local, strip2d, arr2d)
#else
      arr2d(:,:) = arr2dxy_local(:,:)
#endif

      if (masterproc) then
         call wrap_put_vara_realx (id, psid, start2d, count2d, arr2d)
      end if
#else

      do j=beglat,endlat
         numi = nlon(j)
         arr2d(:numi,j) = ps(:numi,j,n3m1)
      end do
      call gather_put (id, psid, arr2d, numperlat, start2d, count2d)
#endif

      do j=beglat,endlat
         numi = nlon(j)
#ifdef STAGGERED
         arr2dxy_local(:numi,j) = phis(:numi,j)
#else
         arr2d(:numi,j) = phis(:numi,j)
#endif
      end do

#ifdef STAGGERED
#if defined (SPMD)
      if (myid_z .eq. 0) call pargatherreal(comm_y, 0,      &
                  arr2dxy_local, strip2d, arr2d)
#else
      arr2d(:,:) = arr2dxy_local(:,:)
#endif

      if (masterproc) then
         call wrap_put_vara_realx (id, phisid, start2d, count2d, arr2d)
      end if
#else

      call gather_put (id, phisid, arr2d, numperlat, start2d, count2d)
#endif

      call gather_chunk_to_field(1,1,1,plon,sgh,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, sghid, start2d, count2d, arr2d)
      end if

      call gather_chunk_to_field(1,1,1,plon,landm,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, landmid, start2d, count2d, arr2d)
      end if
!
! 2d physics arrays
!
      call gather_chunk_to_field(1,1,1,plon,pblht,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, pblhtid, start2d, count2d, arr2d)
      end if

      call gather_chunk_to_field(1,1,1,plon,tpert,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tpertid, start2d, count2d, arr2d)
      end if
         
      qpert_tmp(:,:) = qpert(:,1,:)
      call gather_chunk_to_field(1,1,1,plon,qpert_tmp,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, qpertid, start2d, count2d, arr2d)
      end if
!
! 3d physics arrays
!
      itim = pbuf_old_tim_idx()
      itim = pbuf_next_tim_idx(itim)
 
      ifld = pbuf_get_fld_idx('CLD')
      arr3dchnk_ptr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,itim)
      call gather_chunk_to_field(1,plev,1,plon,arr3dchnk_ptr(:,:,:),arr3d)
      if (masterproc) then
         call wrap_put_vara_realx (id, cldid  , start3d, count3d, arr3d)
      end if

      ifld = pbuf_get_fld_idx('QCWAT')
      arr3dchnk_ptr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,itim)
      call gather_chunk_to_field(1,plev,1,plon,arr3dchnk_ptr(:,:,:),arr3d)
      if (masterproc) then
         call wrap_put_vara_realx (id, qcwatid, start3d, count3d, arr3d)
      end if

      ifld = pbuf_get_fld_idx('TCWAT')
      arr3dchnk_ptr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,itim)
      call gather_chunk_to_field(1,plev,1,plon,arr3dchnk_ptr(:,:,:),arr3d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tcwatid, start3d, count3d, arr3d)
      end if

      ifld = pbuf_get_fld_idx('LCWAT')
      arr3dchnk_ptr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,itim)
      call gather_chunk_to_field(1,plev,1,plon,arr3dchnk_ptr(:,:,:),arr3d)
      if (masterproc) then
         call wrap_put_vara_realx (id, lcwatid, start3d, count3d, arr3d)
      end if

#ifdef CRM
      if (crminitsave) then
        call write_collapse_crmvar3d (id,ucrmid,u_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,vcrmid,v_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,tcrmid,t_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,wcrmid,w_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,qcrmid,q_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,qncrmid,qn_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,qpcrmid,qp_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,qrscrmid,qrs_crm,crmsavechunks)
        call write_collapse_crmvar3d (id,qrlcrmid,qrl_crm,crmsavechunks)

        call write_collapse_crmvar2d (id,fsdscrmid,fsds_crm)
        call write_collapse_crmvar2d (id,fsnscrmid,fsns_crm)
        call write_collapse_crmvar2d (id,fsdsccrmid,fsdsc_crm)
        call write_collapse_crmvar2d (id,fsntoacrmid,fsntoa_crm)
        call write_collapse_crmvar2d (id,fsntoaccrmid,fsntoac_crm)
        call write_collapse_crmvar2d (id,fsutoacrmid,fsutoa_crm)
        call write_collapse_crmvar2d (id,flwdscrmid,flwds_crm)
        call write_collapse_crmvar2d (id,flnscrmid,flns_crm)
        call write_collapse_crmvar2d (id,flutccrmid,flutc_crm)
        call write_collapse_crmvar2d (id,flutcrmid,flut_crm)
        call write_collapse_crmvar2d (id,flnsccrmid,flnsc_crm)

        call gather_chunk_to_field(1,plev,1,plon,qrl1(1:pcols,1:pver,begchunk:endchunk),arr3d)
        if (masterproc) then
          call wrap_put_vara_realx(id,qrl1id,start3d,count3d,arr3d) ! Feelin' a bit shaky about this.
        end if

        call gather_chunk_to_field(1,plev,1,plon,qrs1(1:pcols,1:pver,begchunk:endchunk),arr3d)
        if (masterproc) then
          call wrap_put_vara_realx(id,qrs1id,start3d,count3d,arr3d)
        end if

        call gather_chunk_to_field(1,nrad_buffer,1,plon,rad_buffer(1:pcols,1:nrad_buffer,begchunk:endchunk),arr3dcrmradbfr)
        if (masterproc) then
          call wrap_put_vara_realx(id,crmradbfrid,start3dcrmradbfr,count3dcrmradbfr,arr3dcrmradbfr)
        end if
!       call crmic_debug_dump ('[write_init] at write')
     end if ! if crminitsave
#endif

#if ( ! defined COUP_CSM )
      call gather_chunk_to_field(1,1,1,plon,tsice_rad,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tsice_rad_id, start2d, count2d, arr2d)
      end if

      call gather_chunk_to_field(1,1,1,plon,tsice,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tsiceid, start2d, count2d, arr2d)
      end if

      do j=begchunk,endchunk
         arr2dchnk(:,j) = srfflx_state2d(j)%ts(:)
      end do
      call gather_chunk_to_field(1,1,1,plon,arr2dchnk,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tsid, start2d, count2d, arr2d)
      end if


      do k=1,4
         do j=begchunk,endchunk
            arr2dmchnk(:,k,j) = surface_state2d(j)%tssub(:,k)
         end do
      end do
      call gather_chunk_to_field(1,plevmx,1,plon,arr2dmchnk,arr2dm)
      if (masterproc) then
         do k=1,plevmx
            do j=1,plat
               numi = nlon(j)
               arr2d(:numi,j) = arr2dm(:numi,k,j)
            end do
            call wrap_put_vara_realx (id, tssubid(k), start2d, count2d, arr2d)
         end do
      end if

      call gather_chunk_to_field (1,1,1,plon,sicthk,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, sicthkid, start2d, count2d, arr2d)
      end if
      call gather_chunk_to_field(1,1,1,plon,snowhice,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, snowhiceid, start2d, count2d, arr2d)
      end if
      call gather_chunk_to_field(1,1,1,plon,landfrac,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, landfracid, start2d, count2d, arr2d)
      end if
      do j=begchunk,endchunk
         arr2dchnk(:,j) = surface_state2d(j)%tbot(:)
      end do
      call gather_chunk_to_field(1,1,1,plon,arr2dchnk,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, tbotid, start2d, count2d, arr2d)
      end if
      call gather_chunk_to_field(1,1,1,plon,icefrac,arr2d)
      if (masterproc) then
         call wrap_put_vara_realx (id, icefracid, start2d, count2d, arr2d)
      end if
!
!JR Write out tsocn only in full-physics mode.  Problem is array sst has not been
!JR allocated in adiabatic and ideal_phys modes.
!
      if (.not. adiabatic .and. .not. ideal_phys) then
#if ( defined COUP_SOM )
         call gather_chunk_to_field(1,1,1,plon,tsocn,arr2d)
#else
         call gather_chunk_to_field(1,1,1,plon,sst,arr2d)
#endif
         if (masterproc) then
            call wrap_put_vara_realx (id, tsocnid, start2d, count2d, arr2d)
         end if
      end if
#endif

      if (masterproc) then
         ret = nf_close (id)
         call putfil (filename, path, mss_wpass, mss_irt, .true.)
      end if

      deallocate (arr2d)
      deallocate (arr3d)
      deallocate (arr2dm)
#ifdef STAGGERED
      deallocate (arr3dxyz_local)
      deallocate (arr3dxzy_local)
      deallocate (arr3dxyz)
      deallocate (arr3dxzy)
      deallocate (arr3dus)
      deallocate (arr2dxy_local)
#endif      
      return
   end subroutine write_inithist

#ifdef CRM
! The following utilities collapse interior CRM dimensions into a central
! dimension to define a fake "3D" variable that can be chunked from physics
! domain decomp into files using existing CAM wrapper functions gath_chunk*

! Slightly different logic is required if the CRM variable is crm_nx*crm_ny*crm_nz
! Versus crm_nx*crm_ny

subroutine write_collapse_crmvar2d (id,crmvar2did,crmvar2d)
        use phys_grid
        use crmdims, only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy

        implicit none

        real (r8), intent(in) :: crmvar2d(pcols,crm_nx,crm_ny,begchunk:endchunk)
        integer, intent(in) :: crmvar2did ! netcdf variable id
        integer, intent (in) :: id ! netcdf id

      integer :: k1,k2,cc,m
      integer start3dcrm(4)             ! starting index array for nf_put_vara
      integer count3dcrm(4)
      character (1) :: numb
      integer :: icol,icrmx,icrmy,kk
      data start3dcrm/1,1,1,1/
      real (r8), allocatable :: arr3dcrm(:,:,:), arr3dcrmchnk(:,:,:)

      allocate (arr3dcrm(plon,crm_nx*crm_ny,plat))
      allocate (arr3dcrmchnk(pcols,crm_nx*crm_ny,begchunk:endchunk))
      count3dcrm(1) = plon
      count3dcrm(2) = crm_nx*crm_ny
      count3dcrm(3) = plat
      count3dcrm(4) = 1
      do icol=1,pcols
        cc=0
        do icrmx=1,crm_nx
            do icrmy=1,crm_ny
                cc=cc+1
                arr3dcrmchnk (icol,cc,begchunk:endchunk) = crmvar2d(icol,icrmx,icrmy,begchunk:endchunk)
            end do
        end do
        call gather_chunk_to_field(1,crm_nx*crm_ny,1,plon,arr3dcrmchnk(:,:,:),arr3dcrm)
        if (masterproc) then
          call wrap_put_vara_realx(id,crmvar2did,start3dcrm,count3dcrm,arr3dcrm)
        end if
       end do ! Next chunk of each variable.
       deallocate (arr3dcrm)
       deallocate (arr3dcrmchnk)
end subroutine write_collapse_crmvar2d

subroutine write_collapse_crmvar3d (id,crmvarid,crmvar,crmsavechunks)
        use phys_grid
        use crmdims, only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy
 
        implicit none

        real (r8), intent(in) :: crmvar(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk)
        !integer, intent(in) :: crmvarid(crmsavechunks),crmsavechunks
        integer, intent(in) :: crmsavechunks
        integer, intent(in) :: crmvarid(crmsavechunks)
        integer, intent (in) :: id

      integer :: k1,k2,cc,m
      integer start3dcrm(4)             ! starting index array for nf_put_vara
      integer count3dcrm(4)
      character (1) :: numb
      integer :: icol,icrmx,icrmy,kk
      data start3dcrm/1,1,1,1/
      real (r8), allocatable :: arr3dcrm(:,:,:), arr3dcrmchnk(:,:,:)

      allocate (arr3dcrm(plon,crm_nx*crm_ny*crm_nz/crmsavechunks,plat))
      allocate (arr3dcrmchnk(pcols,crm_nx*crm_ny*crm_nz/crmsavechunks,begchunk:endchunk))
      count3dcrm(1) = plon
      count3dcrm(2) = crm_nx*crm_ny*crm_nz/crmsavechunks
      count3dcrm(3) = plat
      count3dcrm(4) = 1
        do m=1,crmsavechunks
          k1 = 1+(m-1)*crm_nz/crmsavechunks
          k2 = m*crm_nz/crmsavechunks
          do icol=1,pcols
            cc=0
            do icrmx=1,crm_nx
               do icrmy=1,crm_ny
                 do kk=k1,k2
                   cc=cc+1
                   arr3dcrmchnk (icol,cc,begchunk:endchunk) = crmvar(icol,icrmx,icrmy,kk,begchunk:endchunk)
                 end do
               end do
             end do
           end do
           call gather_chunk_to_field(1,crm_nx*crm_ny*crm_nz/crmsavechunks,1,plon,arr3dcrmchnk(:,:,:),arr3dcrm)
          if (masterproc) then
            call wrap_put_vara_realx(id,crmvarid(m),start3dcrm,count3dcrm,arr3dcrm)
          end if 
       end do ! Next chunk of each variable.
       deallocate (arr3dcrm)
       deallocate (arr3dcrmchnk)
   end subroutine write_collapse_crmvar3d
#endif
!#######################################################################

   subroutine gather_only (arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: gather an array on the master processor
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
      use mpishorthand
      use spmd_dyn, only: npes, compute_gsfactors
#endif
!
! Input Arguments
!
      integer, intent(in) :: numperlat ! number of array values to gather per latitude

      real(r8), intent(inout) :: arr(*)! array to be gathered then written

!-----------------------------------------------------------------------

#ifdef SPMD
      integer :: numrecv(0:npes-1)     ! number of items to be received
      integer :: displs(0:npes-1)      ! displacement array
      integer :: numsend               ! number of items to send
      
      call compute_gsfactors (numperlat, numsend, numrecv, displs)
      call mpigatherv (arr, numsend, mpir8, arr, numrecv, displs, mpir8, 0, mpicom)
#endif

      return
   end subroutine gather_only

!#######################################################################

   subroutine gather_put (fileid, fieldid, arr, numperlat, start, count)
!----------------------------------------------------------------------- 
! 
! Purpose: gather an array on the master processor and issue the proper netcdf call
!          to write it out
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
      use mpishorthand
      use spmd_dyn, only: npes, compute_gsfactors
#endif
!
! Input Arguments
!
      integer, intent(in) :: fileid    ! netcdf file id
      integer, intent(in) :: fieldid   ! netcdf field id
      integer, intent(in) :: numperlat ! number of array values to gather per latitude
      integer, intent(in) :: start(*)  ! array of start values for netcdf call
      integer, intent(in) :: count(*)  ! array of count values for netcdf call

      real(r8), intent(inout) :: arr(*)! array to be gathered then written

!-----------------------------------------------------------------------

#ifdef SPMD
      integer :: numrecv(0:npes-1)     ! number of items to be received
      integer :: displs(0:npes-1)      ! displacement array
      integer :: numsend               ! number of items to send
      
      call compute_gsfactors (numperlat, numsend, numrecv, displs)
      call mpigatherv (arr, numsend, mpir8, arr, numrecv, displs, mpir8, 0, mpicom)
#endif

      if (masterproc) then
         call wrap_put_vara_realx (fileid, fieldid, start, count, arr)
      end if
      
      return
   end subroutine gather_put

!#######################################################################

   subroutine addvar (ncid, name, xtype , ndims , dimids, vid)
!----------------------------------------------------------------------- 
! 
! Purpose: Issue the netcdf call to add a variable to the dataset
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: ncid       ! netcdf file id
      integer, intent(in) :: xtype      ! netcdf type flag
      integer, intent(in) :: ndims      ! number of dimensions
      integer, intent(in) :: dimids(:)  ! dimension ids
      integer, intent(in) :: vid        ! variable ids
      
      character(len=*), intent(in) :: name ! variable name
!
! Local workspace
!
      integer f                         ! field index

      call wrap_def_var (ncid, name, xtype , ndims , dimids, vid)
      
      do f=1,nfmaster
         if (masterlist(f)%field%name == name) then
            call wrap_put_att_text (ncid, vid, 'long_name', masterlist(f)%field%long_name)
            if (masterlist(f)%field%units(1:1) /= ' ') then
               call wrap_put_att_text (ncid, vid, 'units', masterlist(f)%field%units)
            end if
            return
         end if
      end do
!
! Field not found in masterlist: long_name and units attributes unknown so just
! return
!
      return
   end subroutine addvar

   subroutine bldfld ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Build Master Field List of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute.
! 
! Method: Call a subroutine to add each field
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ppgrid, only: pver, pverp
      use comsrf, only: plevmx, tsnam
#ifdef CRM
      use crmdims, only: crm_nx, crm_ny, crm_nz
#endif
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!
! Local variables
!
      integer m,j                   ! Indices
      character(len=20) :: string   ! Character string to use
!
! Call addfld to add each field to the Master Field List.
!
      call addfld ('PHIS    ','m2/s2   ',1,    'I','Surface geopotential',phys_decomp)
      call addfld ('PS      ','Pa      ',1,    'A','Surface pressure',phys_decomp)
      call addfld ('T       ','K       ',plev, 'A','Temperature',phys_decomp)
      call addfld ('U       ','m/s     ',plev, 'A','Zonal wind',phys_decomp)
      call addfld ('V       ','m/s     ',plev, 'A','Meridional wind',phys_decomp)
      call addfld ('ETADOT  ','1/s     ',plevp,'A','Vertical (eta) velocity',dyn_decomp)
      call addfld ('SGH     ','m       ',1,    'I','Standard deviation of orography',phys_decomp)
      call addfld ('US      ','m/s     ',plev, 'A','Zonal wind, staggered',dyn_decomp)
      call addfld ('VS      ','m/s     ',plev, 'A','Meridional wind, staggered',dyn_decomp)

! MSP, track SLD dycore tendencies:
      call addfld ('TTEND','K/s',pver,'A','Temperature tendency due to SLD core',phys_decomp)
      call addfld ('STEND','W/kg',pver,'A','Dry static energy tendency due to SLD core',phys_decomp)
      call addfld ('PVTEND','(K m2/kg)/s2',pver,'A','PV tendency due to SLD core',phys_decomp)
      call addfld ('UTEND','m/s2',pver,'A','Zonal wind tendency due to SLD core',phys_decomp)
      call addfld ('VTEND','m/s2',pver,'A','Meridional wind tendency due to SLD core',phys_decomp)
      call addfld ('OMGATEND','Pa/s2',pver,'A','Pressure vertical velocity tendency due to SLD core',phys_decomp)
      call addfld ('QTEND','kg/kg/s',pver,'A','Water vapor tendency due to SLD core',phys_decomp) 
#ifdef QVORTDAMP
      call addfld ('XXXU3','kg/kg/s',pver,'A','Zonal velocity inside SLD dycore',phys_decomp) 
      call addfld ('XXXU3AUX','kg/kg/s',pver,'A','Auxilliary zonal velocity inside SLD dycore',phys_decomp) 
!      call addfld ('QFCST','kg/kg',pver,'A','SLT moisture forecast from scanslt',phys_decomp) 
      call addfld ('DQFCST','kg/kg/s',pver,'A','Difference in SLT moisture forecast tendency due to vortical inteference.',phys_decomp)
      call addfld ('QTENDQFCST','kg/kg',pver,'A','Component of vapor adiabatic tendency from SLT forecast', dyn_decomp)
      call addfld ('QTENDQFIX','kg/kg',pver,'A','Component of vapor adiabatic tendency from obscure fixer', dyn_decomp)
#endif
      call addfld ('QCTEND','kg/kg/s',pver,'A','Liquid water tendency due to SLD core',phys_decomp);
      call addfld ('QITEND','kg/kg/s',pver,'A','Ice water tendency due to SLD core',phys_decomp); 


        call addfld ('QLEVER','kg/kg/s',pver,'A','Subregional moisture interference tendency',phys_decomp)
!
! Constituent tracers
!
      do m=1,pcnst+pnats
         call addfld (cnst_name(m),'kg/kg   ',pver, 'A',cnst_longname(m),phys_decomp)
         call addfld (tendnam(m),'kg/kg/s   ',pver, 'A',trim(cnst_longname(m))//' tendency',dyn_decomp)
         call addfld (qphystendnam(m),'kg/kg/s   ',pver, 'A',trim(cnst_longname(m))//' physics tendency',phys_decomp)
         call addfld (sflxnam(m),  'kg/m2/s ',1,    'A',trim(cnst_name(m))//' surface flux',phys_decomp)
         call addfld (dcconnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' tendency due to moist processes',phys_decomp)
      end do

      call addfld ('DUH     ','K/s     ',plev, 'A','U horizontal diffusive heating',dyn_decomp)
      call addfld ('DVH     ','K/s     ',plev, 'A','V horizontal diffusive heating',dyn_decomp)
      call addfld ('DTH     ','K/s     ',plev, 'A','T horizontal diffusive heating',dyn_decomp)
      call addfld ('SNOWHICE','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
      call addfld ('SNOWHLND','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
      call addfld ('SICTHK  ','m       ',1,    'A','Sea ice thickness',phys_decomp)
      call addfld ('PRECL   ','m/s     ',1,    'A','Large-scale (stable) precipitation rate',phys_decomp)
      call addfld ('PRECC   ','m/s     ',1,    'A','Convective precipitation rate',phys_decomp)
      call addfld ('PRECSH  ','m/s     ',1,    'A','Shallow Convection precipitation rate',phys_decomp)
      write(string,'(F4.2)') precc_thresh*3600.*1000.0
      call addfld ('PRECCFRQ','fraction',1,    'A',&
      'Convective precipitation frequency (fraction of time where rate is > '//trim(string)//'mm/hr)' &
      ,phys_decomp)
      call addfld ('PRECCINT','mm/hr   ',1,    'A',&
      'Convective precipitation rate (less than '//trim(string)// &
      'mm/hr is set to zero -- to get intensity divide by PRECCFRQ)',phys_decomp)
      write(string,'(F4.2)') precl_thresh*3600.*1000.0
      call addfld ('PRECLFRQ','fraction',1,    'A',&
      'Large-scale (stable) precipitation frequency (fraction of time where rate is > '//trim(string)//'mm/hr)' &
      ,phys_decomp)
      call addfld ('PRECLINT','mm/hr   ',1,    'A',&
      'Large-scale (stable) precipitation rate (less than '//trim(string)// &
      'mm/hr is set to zero -- to get intensity divide by PRECLFRQ)',phys_decomp)
      call addfld ('PRECT   ','m/s     ',1,    'A','Total (convective and large-scale) precipitation rate',phys_decomp)
      call addfld ('BRAINRAIN','mm/day',1,'A','Neural net estimated rain rate (mm/day, I think',phys_decomp) 
      call addfld ('BRAINOLR','W/m2',1,'A','Neural net estimated OLR (W/m2, I think',phys_decomp) 

      call addfld ('EVAPPCT ','percent ',1,    'A','Percentage of Zhang-McFarlane precipitation going into evaporation',phys_decomp)
      call addfld ('PRECTMX ','m/s     ',1,    'X','Maximum (convective and large-scale) precipitation rate',phys_decomp)
      call addfld ('PRECSL  ','m/s     ',1,    'A','Large-scale (stable) snow rate (water equivalent)',phys_decomp)
      call addfld ('PRECSC  ','m/s     ',1,    'A','Convective snow rate (water equivalent)',phys_decomp)
      call addfld ('SHFLX   ','W/m2    ',1,    'A','Surface sensible heat flux',phys_decomp)
      call addfld ('LHFLX   ','W/m2    ',1,    'A','Surface latent heat flux',phys_decomp)
      call addfld ('QFLX    ','kg/m2/s ',1,    'A','Surface water flux',phys_decomp)
      call addfld ('PBLH    ','m       ',1,    'A','PBL height',phys_decomp)
      call addfld ('USTAR   ','m/s     ',1,    'A','Surface friction velocity',phys_decomp)
      call addfld ('TREFHT  ','K       ',1,    'A','Reference height temperature',phys_decomp)
#if ( defined COUP_CSM )
      call addfld ('QREFHT  ','kg/kg   ',1,    'A','Reference height humidity',phys_decomp)
#endif
      call addfld ('CGH     ','K/m     ',pverp,'A','Counter-gradient term for heat in PBL',phys_decomp)
      call addfld ('CGQ     ','1/m     ',pverp,'A','Counter-gradient term for moisture in PBL',phys_decomp)
      call addfld ('CGS     ','s/m2    ',pverp,'A','Counter-gradient coeff on surface kinematic fluxes',phys_decomp)
      call addfld ('TPERT   ','K       ',1,    'A','Perturbation temperature (eddies in PBL)',phys_decomp)
      call addfld ('QPERT   ','kg/kg   ',1,    'A','Perturbation specific humidity (eddies in PBL)',phys_decomp)
      call addfld ('KVH     ','m2/s    ',pverp,'A','Vertical diffusion diffusivities (heat/moisture)',phys_decomp)
      call addfld ('KVM     ','m2/s    ',pverp,'A','Vertical diffusion diffusivities (momentum)',phys_decomp)
      call addfld ('DUV     ','m/s2    ',pver, 'A','U vertical diffusion',phys_decomp)
      call addfld ('DVV     ','m/s2    ',pver, 'A','V vertical diffusion',phys_decomp)
      call addfld ('DTV     ','K/s     ',pver, 'A','T vertical diffusion',phys_decomp)
      call addfld ('UCONVMOM','m/s2    ',pver, 'A','U Tendency due to convective transport',phys_decomp)
      call addfld ('VCONVMOM','m/s2    ',pver, 'A','V Tendency due to convective transport',phys_decomp)
      call addfld ('DTVKE   ','K/s     ',pver, 'A','dT/dt vertical diffusion KE dissipation',phys_decomp)
      call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp)
      call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp)
      call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp)
      call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp)
      call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp)
      call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp)
      call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp)
      call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp)
      call addfld ('FLWDS   ','W/m2    ',1,    'A','Downwelling longwave flux at the surface',phys_decomp)
      call addfld ('FLWDSC  ','W/m2    ',1,    'A','Clearsky downwelling longwave flux at the surface',phys_decomp)
      call addfld ('FSUTOA  ','W/m2    ',1,    'A','Upwelling shortwave flux at the TOA',phys_decomp)
      call addfld ('FSUTOAC ','W/m2    ',1,    'A','Clearsky upwelling shortwave flux at the TOA',phys_decomp)
      call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp)
      call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp)
      call addfld ('CLOUD   ','fraction',pver, 'A','Cloud fraction',phys_decomp)
      call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp)
      call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp)
      call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp)
      call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp)
      call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp)
      call addfld ('OMEGA   ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)
      call addfld ('DQP     ','kg/kg/s ',pver, 'A','Specific humidity tendency due to precipitation',phys_decomp)
      call addfld ('TAUX    ','N/m2    ',1,    'A','Zonal surface stress',phys_decomp)
      call addfld ('TAUY    ','N/m2    ',1,    'A','Meridional surface stress',phys_decomp)
      call addfld ('SRFRAD  ','W/m2    ',1,    'A','Net radiative flux at surface',phys_decomp)
      call addfld ('QRS     ','K/s     ',pver, 'A','Solar heating rate',phys_decomp)
      call addfld ('QRL     ','K/s     ',pver, 'A','Longwave heating rate',phys_decomp)
#ifdef QRLDAMP
      call addfld ('DQRL    ','K/s     ',pver, 'A','Interference in longwave heating rate',phys_decomp)
#endif
#ifdef FLUXDAMP
      call addfld ('DUBOT   ','m/s2    ',1, 'A','Interference in surface zonal winds used for surface fluxes',phys_decomp)
      call addfld ('DVBOT   ','m/s2    ',1, 'A','Interference in suface meridional winds used for surface fluxes',phys_decomp)
      call addfld ('DWINDBOT','m/s2    ',1, 'A','Interference in suface wind speed used for surface fluxes',phys_decomp)
      call addfld ('DLHFLX  ','W/m2    ',1, 'A','Interference in latent heat flux resulting from surface wind interference', phys_decomp)
      call addfld ('DSHFLX  ','W/m2    ',1, 'A','Interference in sensible heat flux resulting from surface wind interference', phys_decomp)
#endif
      call addfld ('CLDTOT  ','fraction',1,    'A','Vertically-integrated total cloud if lwp>0.02 kg/m2',phys_decomp)
      call addfld ('CLDLOW  ','fraction',1,    'A','Vertically-integrated low cloud ptop > 700 mb',phys_decomp)
      call addfld ('CLDMED  ','fraction',1,    'A','Vertically-integrated mid-level cloud 700 mb > ptop > 400 mb',phys_decomp)
      call addfld ('CLDHGH  ','fraction',1,    'A','Vertically-integrated high cloud ptop < 400 mb',phys_decomp)
      call addfld ('ENGYCORR','W/m2    ',plev, 'A','Energy correction for over-all conservation',dyn_decomp)
!
! Put TS1 on history file since it will be different than TS in a coupled
! run
!
      do m=1,plevmx
         call addfld (tsnam(m),'K       ',1,    'A',tsnam(m)//' subsoil temperature',phys_decomp)
      end do

      call addfld ('TS      ','K       ',1,    'A','Surface temperature (radiative)',phys_decomp)
      call addfld ('TSMN    ','K       ',1,    'M','Minimum surface temperature over output period',phys_decomp)
      call addfld ('TSMX    ','K       ',1,    'X','Maximum surface temperature over output period',phys_decomp)
      call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp)
      call addfld ('FU      ','m/s     ',plev, 'I','Zonal wind forcing term',dyn_decomp)
      call addfld ('FV      ','m/s     ',plev, 'I','Meridional wind forcing term',dyn_decomp)
!      call addfld ('UTEND   ','m/s2    ',plev, 'A','U tendency',dyn_decomp)
!      call addfld ('VTEND   ','m/s2    ',plev, 'A','V tendency',dyn_decomp)
!      call addfld ('TTEND   ','K/s     ',plev, 'A','T tendency',dyn_decomp)
     call addfld ('TPHYSTND','K/s     ',plev, 'A','T physics tendency',phys_decomp)
     call addfld ('UPHYSTND','m/s2     ',plev, 'A','U physics tendency',phys_decomp)
     call addfld ('VPHYSTND','m/s2     ',plev, 'A','V physics tendency',phys_decomp)
     call addfld ('PVPHYSTND','K m2/kg/s2',plev, 'A','PV physics tendency',phys_decomp)
     call addfld ('PV1PHYSTND','K m2/kg/s2',plev, 'A','Physics tendency of first component in PV = -G*(vr*dthdp - (dthdx*dvdp-dthdy*dudp)',phys_decomp)
     call addfld ('PV2PHYSTND','K m2/kg/s2',plev, 'A','Physics tendency of second component in PV = -G*(vr*dthdp - (dthdx*dvdp-dthdy*dudp)',phys_decomp)
     call addfld ('PV3PHYSTND','K m2/kg/s2',plev, 'A','Physics tendency of third component in PV = -G*(vr*dthdp - (dthdx*dvdp-dthdy*dudp)',phys_decomp)


! MSP: for redundantly tracking total tendency after physics package (so we know MSE budget residual magnitude exactly)
      call addfld ('PVAP','(K m2/kg)/s',plev,'I','PV after physics',phys_decomp)
      call addfld ('TAP','K     ',plev, 'I','T after physics',phys_decomp)
      call addfld ('QAP','kg/kg  ',plev, 'I','Q after physics',phys_decomp)
      call addfld ('QCAP','kg/kg  ',plev, 'I','QC after physics',phys_decomp)
      call addfld ('QIAP','kg/kg  ',plev, 'I','QI after physics',phys_decomp)
      call addfld ('UAP','m/s ',plev, 'I','U after physics',phys_decomp)
      call addfld ('VAP','kg/kg  ',plev, 'I','V after physics',phys_decomp)

      call addfld ('LPSTEN  ','Pa/s    ',1,    'A','Surface pressure tendency',dyn_decomp)
      call addfld ('DTCOND  ','K/s     ',pver, 'A','T tendency - moist processes',phys_decomp)
      call addfld ('CMFDT   ','K/s     ',pver, 'A','T tendency - Hack convection',phys_decomp)
      call addfld ('CMFDQ   ','kg/kg/s ',pver, 'A','Q tendency - Hack convection',phys_decomp)
      call addfld ('ZMDT    ','K/s     ',pver, 'A','T tendency - Zhang-McFarlane moist convection',phys_decomp)
      call addfld ('ZMDQ    ','kg/kg/s ',pver, 'A','Q tendency - Zhang-McFarlane moist convection',phys_decomp)
      call addfld ('BRAINDT ','K/s     ',pver, 'A','T tendency - Gentine-Pritchard neural net for moist convection',phys_decomp)
      call addfld ('BRAINDT2','K/s     ',pver, 'A','T tendency - neural net (after energy fix)',phys_decomp)
      call addfld ('BRAINDQ ','kg/kg/s ',pver, 'A','Q tendency - neural net for moist convection',phys_decomp)
      call addfld ('BRAINDQ2','kg/kg/s ',pver, 'A','Q tendency - neural net post energy fix',phys_decomp)
      call addfld ('CMFDQR  ','kg/kg/s ',pver, 'A','Q tendency - shallow convection rainout',phys_decomp)
      call addfld ('QC      ','kg/kg/s ',pver, 'A','Q tendency - shallow convection LW export',phys_decomp)
      call addfld ('CMFMC   ','kg/m2/s ',pverp,'A','Moist convection mass flux',phys_decomp)
      call addfld ('CMFSL   ','W/m2    ',pverp,'A','Moist convection liquid water static energy flux',phys_decomp)
      call addfld ('CMFLQ   ','W/m2    ',pverp,'A','Moist convection total water flux',phys_decomp)
      call addfld ('CNVCLD  ','fraction',1,    'A','Vertically integrated convective cloud amount',phys_decomp)
      call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp)
      call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp)
      call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp)
      call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp)
      call addfld ('O3VMR   ','m3/m3   ',pver, 'A','Ozone volume mixing ratio',phys_decomp)
      call addfld ('CLDST   ','fraction',pver, 'A','Stratus cloud fraction',phys_decomp)
      call addfld ('CONCLD  ','fraction',pver, 'A','Convective cloud cover',phys_decomp)
#ifdef CRM
      call addfld ('PRES    ','Pa      ',pver, 'A','Pressure',phys_decomp)
      call addfld ('DPRES   ','Pa      ',pver, 'A','Pressure thickness of layer',phys_decomp)
      call addfld ('HEIGHT  ','m       ',pver, 'A','Height above surface',phys_decomp)
      call addfld ('_QRS    ','K/s     ',pver, 'A','Solar heating rate',phys_decomp)
      call addfld ('_QRL    ','K/s     ',pver, 'A','Longwave heating rate',phys_decomp)
      call addfld ('_FSDS   ','W/m2    ',1,    'A','Flux Shortwave Downwelling Surface',phys_decomp)
      call addfld ('_FSNS   ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp)
      call addfld ('_FLNS   ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp)
      call addfld ('_FLUT   ','W/m2    ',1,    'A','Upward LW flux at TOM',phys_decomp)
      call addfld ('_FSNTOA ','W/m2    ',1,    'A','Net solar flx at top of Atmos',phys_decomp)
      call addfld ('_LWCF   ','W/m2    ',1,    'A','LW cloud forcing',phys_decomp)
      call addfld ('_SWCF   ','W/m2    ',1,    'A','SW cloud forcing',phys_decomp)
      call addfld ('_FLWDS  ','W/m2    ',1,    'A','Downwelling longwave flux at the surface',phys_decomp)
      call addfld ('_FLWDSC ','W/m2    ',1,    'A','Clearsky downwelling longwave flux at the surface',phys_decomp)
      call addfld ('SPDT    ','K/s     ',pver, 'A','T tendency due to CRM',phys_decomp)
      call addfld ('SPDQ    ','kg/kg/s ',pver, 'A','Q tendency due to CRM',phys_decomp)
      call addfld ('SPDQC   ','kg/kg/s ',pver, 'A','QC tendency due to CRM',phys_decomp)
      call addfld ('SPDQI   ','kg/kg/s ',pver, 'A','QI tendency due to CRM',phys_decomp)
      call addfld ('SPMC    ','kg/m2/s ',pver, 'A','Total mass flux from CRM',phys_decomp)
      call addfld ('SPMCUP  ','kg/m2/s ',pver, 'A','Updraft mass flux from CRM',phys_decomp)
      call addfld ('SPMCDN  ','kg/m2/s ',pver, 'A','Downdraft mass flux from CRM',phys_decomp)
      call addfld ('SPMCUUP ','kg/m2/s ',pver, 'A','Unsaturated updraft mass flux from CRM',phys_decomp)
      call addfld ('SPMCUDN ','kg/m2/s ',pver, 'A','Unsaturated downdraft mass flux from CRM',phys_decomp)
      call addfld ('SPQC    ','kg/kg   ',pver, 'A','Cloud water from CRM',phys_decomp)
      call addfld ('SPQI    ','kg/kg   ',pver, 'A','Cloud ice from CRM',phys_decomp)
      call addfld ('SPQS    ','kg/kg   ',pver, 'A','Snow from CRM',phys_decomp)
      call addfld ('SPQG    ','kg/kg   ',pver, 'A','Graupel from CRM',phys_decomp)
      call addfld ('SPQR    ','kg/kg   ',pver, 'A','Rain from CRM',phys_decomp)
      call addfld ('SPQTFLX ','kg/m2/s ',pver, 'A','Nonprecip. water flux from CRM',phys_decomp)
      call addfld ('SPUFLX  ','m2/s2   ',pver, 'A','x-momentum flux from CRM',phys_decomp)
      call addfld ('SPVFLX  ','m2/s2   ',pver, 'A','y-momentum flux from CRM',phys_decomp)
      call addfld ('SPQTFLXS','kg/m2/s ',pver, 'A','SGS Nonprecip. water flux from CRM',phys_decomp)
      call addfld ('TKE     ','kg/m/s2 ',pver, 'A','Total TKE in CRM',phys_decomp)
      call addfld ('TKES    ','kg/m/s2 ',pver, 'A','SGS TKE in CRM',phys_decomp)
      call addfld ('SPQPFLX ','kg/m2/s ',pver, 'A','Precip. water flux from CRM',phys_decomp)
      call addfld ('SPPFLX  ','m/s     ',pver, 'A','Precipitation flux from CRM',phys_decomp)
      call addfld ('SPQTLS  ','kg/kg/s ',pver, 'A','L.S. Vapor Tendency from CRM',phys_decomp)
      call addfld ('SPQTTR  ','kg/kg/s ',pver, 'A','Nonprec. water transport from CRM',phys_decomp)
      call addfld ('SPQPTR  ','kg/kg/s ',pver, 'A','Prec. water transport from CRM',phys_decomp)
      call addfld ('SPQPEVP ','kg/kg/s ',pver, 'A','Prec. water evaporation from CRM',phys_decomp)
      call addfld ('SPQPFALL','kg/kg/s ',pver, 'A','Prec. water fall-out from CRM',phys_decomp)
      call addfld ('SPQPSRC ','kg/kg/s ',pver, 'A','Prec. water source from CRM',phys_decomp)
      call addfld ('SPTLS   ','kg/kg/s ',pver, 'A','L.S. LIWSE Tendency from CRM',phys_decomp)
      call addfld ('SPULS   ','m/s2',pver, 'A','L.S. U Tendency from CRM',phys_decomp)
      call addfld ('SPVLS   ','m/s2 ',pver, 'A','L.S. V Tendency from CRM',phys_decomp)
      call addfld ('CLOUDTOP','        ',pver, 'A','Cloud Top PDF',phys_decomp)
      call addfld ('CLOUDR  ','        ',pver, 'A','Cloud fraction based on -30dBZ radar reflectivity',phys_decomp)
      call addfld ('_CLDST  ','fraction',pver, 'A','Stratus cloud fraction convenctional cloud param',phys_decomp)
      call addfld ('_CONCLD ','fraction',pver, 'A','Convective cloud cover convenctional cloud param',phys_decomp)
      call addfld ('_CNVCLD ','fraction',1,    'A','Vertically int convective cloud amount conventional cloud param',phys_decomp)
      call addfld ('_CLOUD  ','        ',pver, 'A','Cloud fraction from convenctional cloud param',phys_decomp)
      call addfld ('_CLDTOT ','        ',1,    'A','Vertically-integrated total cloud convenctional cloud param', &
                   phys_decomp)
      call addfld ('_CLDLOW ','        ',1,    'A','Vertically-integrated low cloud convenctional cloud param', &
                   phys_decomp)
      call addfld ('_CLDMED ','        ',1,    'A','Vertically-integrated mid-level cloud convenctional cloud param', &
                   phys_decomp)
      call addfld ('_CLDHGH ','        ',1,    'A','Vertically-integrated high cloud convenctional cloud param', &
                   phys_decomp)
      call addfld ('CRM_U   ','m/s     ',crm_nx*crm_ny*crm_nz, 'A','CRM x-wind',phys_decomp)
      call addfld ('CRM_V   ','m/s     ',crm_nx*crm_ny*crm_nz, 'A','CRM y-wind',phys_decomp)
      call addfld ('CRM_W   ','m/s     ',crm_nx*crm_ny*crm_nz, 'A','CRM z-wind',phys_decomp)
      call addfld ('CRM_TABS','K       ',crm_nx*crm_ny*crm_nz, 'A','CRM Temperature',phys_decomp)
      call addfld ('CRM_QV  ','g/kg    ',crm_nx*crm_ny*crm_nz, 'A','CRM Water Vapor',phys_decomp)
      call addfld ('CRM_QC  ','g/kg    ',crm_nx*crm_ny*crm_nz, 'A','CRM Cloud Water',phys_decomp)
      call addfld ('CRM_QI  ','g/kg    ',crm_nx*crm_ny*crm_nz, 'A','CRM Cloud Ice',phys_decomp)
      call addfld ('CRM_QPC ','g/kg    ',crm_nx*crm_ny*crm_nz, 'A','CRM Precipitating Water',phys_decomp)
      call addfld ('CRM_QPI ','g/kg    ',crm_nx*crm_ny*crm_nz, 'A','CRM Precipitating Ice',phys_decomp)
      call addfld ('CRM_PREC','m/s     ',crm_nx*crm_ny, 'A','CRM Precipitation Rate',phys_decomp)
      call addfld ('CRM_QRS ','K/s     ',crm_nx*crm_ny*crm_nz, 'A','CRM Shortwave radiative heating rate',phys_decomp)
      call addfld ('CRM_QRL ','K/s     ',crm_nx*crm_ny*crm_nz, 'A','CRM Longwave radiative heating rate',phys_decomp)
      call addfld ('CRM_FSNT','W/m2    ',crm_nx*crm_ny, 'A','CRM Net shortwave flux at top',phys_decomp)
      call addfld ('CRMFSNTC','W/m2    ',crm_nx*crm_ny, 'A','CRM clearsky net shortwave flux at top',phys_decomp)
      call addfld ('CRM_FSUT','W/m2    ',crm_nx*crm_ny, 'A','CRM Upward shortwave flux at top of atmosphere',phys_decomp)
      call addfld ('CRM_FSNS','W/m2    ',crm_nx*crm_ny, 'A','CRM Net shortwave flux at surface',phys_decomp)
      call addfld ('CRMFSDSC','W/m2    ',crm_nx*crm_ny, 'A','CRM clearsky downward shortwave flux at surface',phys_decomp)
      call addfld ('CRM_FSDS','W/m2    ',crm_nx*crm_ny, 'A','CRM Downward shortwave flux at surface',phys_decomp)
      call addfld ('CRM_FLUT','W/m2    ',crm_nx*crm_ny, 'A','CRM Upward longwave flux at top of model',phys_decomp)
      call addfld ('CRMFLUTC','W/m2    ',crm_nx*crm_ny, 'A','CRM clearsky upward longwave flux at top of model',phys_decomp)
      call addfld ('CRM_FLNS','W/m2    ',crm_nx*crm_ny, 'A','CRM Net upward longwave flux at surface',phys_decomp)
      call addfld ('CRMFLNSC','W/m2    ',crm_nx*crm_ny, 'A','CRM clearksy net upward longwave flux at surface',phys_decomp)
      call addfld ('CRM_FLDS','W/m2    ',crm_nx*crm_ny, 'A','CRM Downward longwave flux at surface',phys_decomp)
! ISCCP fields

      call addfld ('_FISCCP1','mixed   ',npres*ntau,'A', &
                   'grid box fraction covered by each ISCCP D level cloud types conventional cloud param', &
                   phys_decomp,flag_xyfill=.true., flag_isccplev=.true.)
      call addfld ('_TCLDAR ','fraction',1,'A','Total cloud area conventional cloud param', &
                   phys_decomp,flag_xyfill=.true.)
      call addfld ('_MNPTOP ','mb      ',1,'A','Mean cloud top pressure conventional cloud param', &
                   phys_decomp,flag_xyfill=.true.)
      call addfld ('_MNTAU  ','unitless',1,'A','Mean optical thickness conventional cloud param', &
                   phys_decomp,flag_xyfill=.true.)
      call addfld ('_MNTTOP ','K       ',1,'A','Mean cloud top temperature conventional cloud param', &
                    phys_decomp,flag_xyfill=.true.)
      call addfld ('_CLOUDY ','fraction',1,'A','Binary flag for cltot >=cldmin conventional cloud param', &
                    phys_decomp)

! ARM diagnostics (from one column of CRM):

      call addfld ('CLDTOTR ','        ',1, 'A','Total Cloud Fraction (CRM column)',phys_decomp)
      call addfld ('CLDHGHR ','        ',1, 'A','High-Cloud Fraction (CRM column)',phys_decomp)
      call addfld ('CLDMEDR ','        ',1, 'A','Mid-Cloud Fraction (CRM column)',phys_decomp)
      call addfld ('CLDLOWR ','        ',1, 'A','Low-Cloud Fraction (CRM column)',phys_decomp)
      call addfld ('CLDTOTRD','        ',1, 'A','Total Cloud Fraction (CRM domain)',phys_decomp)
      call addfld ('CLDHGHRD','        ',1, 'A','High-Cloud Fraction (CRM domain)',phys_decomp)
      call addfld ('CLDMEDRD','        ',1, 'A','Mid-Cloud Fraction (CRM domain)',phys_decomp)
      call addfld ('CLDLOWRD','        ',1, 'A','Low-Cloud Fraction (CRM domain)',phys_decomp)
      call addfld ('LWPR    ','kg/m2    ',1, 'A','Liquid Water Path (CRM column)',phys_decomp)
      call addfld ('IWPR    ','kg/m2    ',1, 'A','Ice Water Path (CRM column)',phys_decomp)
      call addfld ('LOBR    ','m       ',1, 'A','Lowest cloud boundary height (CRM column)',phys_decomp)
      call addfld ('HIBR    ','m       ',1, 'A','Highest cloud boundary height (CRM column)',phys_decomp)
      call addfld ('PRECR   ','m/s     ',1, 'A','Precipitation Rate (CRM column)',phys_decomp)
      call addfld ('FSDSR   ','W/m2    ',1, 'A','Downward shortwave flux at the surface (CRM column)',phys_decomp)
      call addfld ('FLDSR   ','W/m2    ',1, 'A','Downward longwave flux at the surface (CRM column)',phys_decomp)
      call addfld ('FSUTR   ','W/m2    ',1, 'A','Upward shortwave flux at the TOA (CRM column)',phys_decomp)
      call addfld ('FLUTR   ','W/m2    ',1, 'A','Upward longwave flux at the TOA (CRM column)',phys_decomp)
      call addfld ('FSDSCR  ','W/m2    ',1, 'A','Clear-sky downward shortwave flux at the surface (CRM column)',phys_decomp)
      call addfld ('FSUTCR  ','W/m2    ',1, 'A','Clear-sky upward shortwave flux at the TOA (CRM column)',phys_decomp)
      call addfld ('Z0M     ','m       ',1, 'A','Momentum roughness length for CRM',phys_decomp)
      call addfld ('BFLX    ','K*kg/m2/s',1, 'A','Buoyancy flux for CRM',phys_decomp)
      call addfld ('TAUX_CRM','N/m2    ',1, 'A','zonal CRM stress perturbation',phys_decomp)
      call addfld ('TAUY_CRM','N/m2    ',1, 'A','merid CRM stress perturbation',phys_decomp)
      call addfld ('TIMINGF ','        ',1, 'A','CRM CPU usage efficiency factor: 1 - ideal',phys_decomp)


#endif
!++pcw
      call addfld ('RELHUM  ','percent ',pver, 'A','Relative humidity',phys_decomp)
!--pcw

      call addfld ('SULFBIO ','kg/kg   ',pver, 'A','Biogenic sulfate mass mixing ratio' ,phys_decomp)
      call addfld ('SULFANT ','kg/kg   ',pver, 'A','Anthropogenic sulfate mass mixing ratio' ,phys_decomp)
      call addfld ('SULFMMR ','kg/kg   ',pver, 'A','Sulfate mass mixing ratio' ,phys_decomp)
!++tls
      call addfld ('AERMMR  ','kg/kg   ',pver, 'A','Sulfate + bg aerosol mass mixing ratio',phys_decomp)
      call addfld ('REL     ','um      ',pver, 'A','Effective radius of liquid water droplets',phys_decomp)
      call addfld ('TAUVIS  ','unitless',1,    'A','Total column aerosol extinction, vis band [aerosol optical depth]',phys_decomp)
      call addfld ('MSO4    ','gram/cm3',pver, 'A','Mass concentration of SO4',phys_decomp)
      call addfld ('LWC     ','kg/m3   ',pver, 'A','Liquid Water Content',phys_decomp)
      call addfld ('CLDFRQ  ','fraction',pver, 'A','Frequency of occurance of clouds (CLOUD > 0.01)',phys_decomp)
      call addfld ('WREL    ','um      ',pver, 'A','Weighted effective radius (by CLDFRQ)',phys_decomp)
      call addfld ('WLWC    ','kg/m3   ',pver, 'A','Weighted Liquid Water Content, prognostic (by CLDFRQ)',phys_decomp)
      call addfld ('TREFHTMN','K       ',1,    'M','Minimum reference height temperature over output period',phys_decomp)
      call addfld ('TREFHTMX','K       ',1,    'X','Maximum reference height temperature over output period',phys_decomp)
      call addfld ('TREFMNAV','K       ',1,    'A','Average of TREFHT daily minimum',phys_decomp)
      call addfld ('TREFMXAV','K       ',1,    'A','Average of TREFHT daily maximum',phys_decomp)
      call addfld ('LANDFRAC','fraction',1,    'A','Fraction of sfc area covered by land',phys_decomp)
      call addfld ('ICEFRAC ','fraction',1,    'A','Fraction of sfc area covered by sea-ice',phys_decomp)
      call addfld ('OCNFRAC ','fraction',1,    'A','Fraction of sfc area covered by ocean',phys_decomp)
      call addfld ('UBOT    ','m/s     ',1,    'A','Lowest model level zonal wind',phys_decomp)
      call addfld ('VBOT    ','m/s     ',1,    'A','Lowest model level meridional wind',phys_decomp)
      call addfld ('QBOT    ','kg/kg   ',1,    'A','Lowest model level water vapor mixing ratio',phys_decomp)
      call addfld ('TBOT    ','K       ',1,    'A','Lowest model level temperature', phys_decomp)
      call addfld ('PBOT    ','Pa      ',1,    'A','Lowest model level pressure', phys_decomp)
      call addfld ('ZBOT    ','m       ',1,    'A','Lowest model level height', phys_decomp)
#if ( defined COUP_CSM )
      call addfld ('CPLRAINC','kg/m2/s ',1,    'A','Convective rainfall sent to coupler' ,phys_decomp)
      call addfld ('CPLRAINL','kg/m2/s ',1,    'A','Large-scale rainfall sent to coupler' ,phys_decomp)
      call addfld ('CPLSNOWC','kg/m2/s ',1,    'A','Convective snowfall sent to coupler' ,phys_decomp)
      call addfld ('CPLSNOWL','kg/m2/s ',1,    'A','Large-scale snowfall sent to coupler' ,phys_decomp)
      call addfld ('CPLPRCER','kg/m2/s ',1,    'A','Error in precipitation state (rain or snow) sent to coupler' ,phys_decomp)
#endif
      call addfld ('PRECCav ','m/s     ',1,    'A','Average large-scale precipitation',phys_decomp)
      call addfld ('PRECLav ','m/s     ',1,    'A','Average convective precipitation',phys_decomp)


 call addfld ('UNTEND  ','m/s2    ',pver, 'A','U-analysis nudge tendency',phys_decomp)
  call addfld ('VNTEND  ','m/s2    ',pver, 'A','V-analysis nudge tendency',phys_decomp)
  call addfld ('SNTEND  ','W/kg     ',pver, 'A','S-analysis nudge tendency',phys_decomp)
  call addfld ('QNTEND  ','kg/kg/s ',pver, 'A','Q-analysis nudge tendency',phys_decomp)
  call addfld ('LPSNTEN ','Pa/s    ',1,    'A','Surface pressure analysis nudge tendency',phys_decomp)
  call addfld ('MSPSA','J/kg',pver,'A','S in the analysis at nudging',phys_decomp)
  call addfld ('MSPS','J/kg',pver,'A','state%s in the model at nudging',phys_decomp)
  call addfld ('MSPZ','m',pver,'A','state%zm in the model at nudging',phys_decomp)

!
! Variables output on specific pressure surfaces
!
      call addfld ('Z700    ','m       ',1,    'A','Geopotential Z at 700 mbar pressure surface',phys_decomp)
      call addfld ('Z500    ','m       ',1,    'A','Geopotential Z at 500 mbar pressure surface',phys_decomp)
      call addfld ('Z300    ','m       ',1,    'A','Geopotential Z at 300 mbar pressure surface',phys_decomp)
      call addfld ('Z050    ','m       ',1,    'A','Geopotential Z at 50 mbar pressure surface',phys_decomp)
      call addfld ('Q850    ','kg/kg   ',1,    'A','Specific Humidity at 850 mbar pressure surface',phys_decomp)
      call addfld ('Q200    ','kg/kg   ',1,    'A','Specific Humidity at 700 mbar pressure surface',phys_decomp)
      call addfld ('T850    ','K       ',1,    'A','Temperature at 850 mbar pressure surface',phys_decomp)
      call addfld ('T300    ','K       ',1,    'A','Temperature at 300 mbar pressure surface',phys_decomp)
      call addfld ('U850    ','m/s     ',1,    'A','Zonal wind at 850 mbar pressure surface',phys_decomp)
      call addfld ('U200    ','m/s     ',1,    'A','Zonal wind at 200 mbar pressure surface',phys_decomp)
      call addfld ('V850    ','m/s     ',1,    'A','Meridional wind at 850 mbar pressure surface',phys_decomp)
      call addfld ('V200    ','m/s     ',1,    'A','Meridional wind at 200 mbar pressure surface',phys_decomp)
      call addfld ('OMEGA850','Pa/s    ',1,    'A','Vertical velocity at 850 mbar pressure surface',phys_decomp)
      call addfld ('OMEGA500','Pa/s    ',1,    'A','Vertical velocity at 500 mbar pressure surface',phys_decomp)
!
! Fields needed such that all variables passed to flux coupler
! also appear on history tape
!
      call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp)
      call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp)
      call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp)
      call addfld ('SOLSD   ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp)
!
! Heating rate needed for d(theta)/dt computation
!
      call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)
!
! Fields which go out active monthly average history files
!
      call addfld ('VT      ','K m/s   ',plev, 'A','Meridional heat transport',phys_decomp)
      call addfld ('VU      ','m2/s2   ',plev, 'A','Meridional flux of zonal momentum' ,phys_decomp)
      call addfld ('VV      ','m2/s2   ',plev, 'A','Meridional velocity squared' ,phys_decomp)
      call addfld ('WSPEED  ','m/s     ',plev, 'X','Horizontal total wind speed' ,phys_decomp)
      call addfld ('UU      ','m2/s2   ',plev, 'A','Zonal velocity squared' ,phys_decomp)
      call addfld ('ZZ      ','m2      ',pver ,'A','Eddy height variance' ,phys_decomp)
      call addfld ('TT      ','K2      ',plev ,'A','Eddy temperature variance' ,phys_decomp)
      call addfld ('OMEGAT  ','K Pa/s  ',pver ,'A','Vertical heat flux' ,phys_decomp)
      call addfld ('OMEGAU  ','K Pa/s  ',pver ,'A','Vertical flux of zonal momentum' ,phys_decomp)
      call addfld ('VZ      ','m2/s    ',plev, 'A','Meridional transport of geopotential energy',phys_decomp)
      call addfld ('VQ      ','m/skg/kg',pver, 'A','Meridional water transport',phys_decomp)
      call addfld ('Z3      ','m       ',pver, 'A','Geopotential Height (above sea level)',phys_decomp)
      call addfld ('MQ      ','kg/m2   ',pver, 'A','Water vapor mass in layer',phys_decomp)
      call addfld ('TMQ     ','kg/m2   ',1,    'A','Total (vertically integrated) precipitatable water',phys_decomp)
      call addfld ('PSL     ','Pa      ',1,    'A','Sea level pressure',phys_decomp)
!
! Fields which are only needed for recording attributes in master field list
!
      call addfld ('LANDMCOS','unitless',1,    'I', &
                   'Land ocean transition mask: ocean (0), continent (1), transition (0-1)',phys_decomp)
#if ( defined COUP_SOM )
      call addfld ('MELTB   ','m/s   ',1,'A','Sea ice basal melt rate',phys_decomp)
      call addfld ('MELTT   ','m/s   ',1,'A','Sea ice top surface melt rate',phys_decomp)
      call addfld ('MELTL   ','m/s   ',1,'A','Sea ice lateral surface melt rate',phys_decomp)
      call addfld ('GROWB   ','m/s   ',1,'A','Sea ice basal growth rate',phys_decomp)
      call addfld ('FRAZIL  ','m/s   ',1,'A','Sea ice frazil growth rate',phys_decomp)
      call addfld ('FLOOD   ','units?',1,'A','Description of FLOOD goes here',phys_decomp)
      call addfld ('FRZMLT  ','W/m2  ',1,'A','Potential to grow/melt sea ice',phys_decomp)
      call addfld ('QFLUX   ','W/m2  ',1,'A','Ocean mixed layer heat flux',phys_decomp)
      call addfld ('QFLUX_FT','W/m2  ',1,'A','M.L. heat flux after ft adjustment',phys_decomp)
      call addfld ('QFLUX_TH','W/m2  ',1,'A','M.L. heat flux after thickness factor adjustment',phys_decomp)
      call addfld ('QFLUX_A2','W/m2  ',1,'A','Adjusted ocean mixed layer heat flux 2',phys_decomp)
      call addfld ('FOCN    ','W/m2  ',1,'A','Melting heat flux from ice model',phys_decomp)
      call addfld ('EICEIN  ','W/m2  ',1,'A','Ice internal energy init/dtime',phys_decomp)
      call addfld ('EICEOUT ','W/m2  ',1,'A','Ice internal energy final/dtime',phys_decomp)
      call addfld ('F_ICE   ','W/m2  ',1,'A','net ice sfc flx over water/ice',phys_decomp)
      call addfld ('F_OCN   ','W/m2  ',1,'A','ice melt flx',phys_decomp)
      call addfld ('FRZMLTMX','W/m2  ',1,'A','ice formation flx (positives zeroed)',phys_decomp)
      call addfld ('DELTAICE','W/m2  ',1,'A','change in ice area',phys_decomp)
      call addfld ('IMBAL   ','W/m2  ',1,'A','local ice imbalance',phys_decomp)
      call addfld ('NRGERROR','W/m2  ',1,'A','local ice imbalance',phys_decomp)
      call addfld ('MLDANN  ','m     ',1,'I','mixed layer depth',phys_decomp)
      call addfld ('ONF     ','W/m2  ',1,'A','net ocn sfc flx over water/ice',phys_decomp)
      call addfld ('OIE     ','J/m2  ',1,'A','ocean internal energy',phys_decomp)
      call addfld ('OIERATE ','W/m2  ',1,'A','ocean internal energy change rate',phys_decomp)
      call addfld ('NRGICE  ','J/m2  ',1,'A','ice internal energy',phys_decomp)
      call addfld ('IIERATE ','W/m2  ',1,'A','ice internal energy change rate',phys_decomp)
#endif
      call addfld ('TSICE   ','K     ',1,'A','Ice temperature',phys_decomp)
!
! Fields to calculate surface energy budget for each component model
!
      call addfld ('FSNSLND ','W/m2  ',1,'A','FSNS over land',phys_decomp)
      call addfld ('FLNSLND ','W/m2  ',1,'A','FLNS over land',phys_decomp)
      call addfld ('LHFLXLND','W/m2  ',1,'A','LHFLX over land',phys_decomp)
      call addfld ('SHFLXLND','W/m2  ',1,'A','SHFLX over land',phys_decomp)

      call addfld ('FSNSOCN ','W/m2  ',1,'A','FSNS over open ocn',phys_decomp)
      call addfld ('FLNSOCN ','W/m2  ',1,'A','FLNS over open ocn',phys_decomp)
      call addfld ('LHFLXOCN','W/m2  ',1,'A','LHFLX over open ocn',phys_decomp)
      call addfld ('SHFLXOCN','W/m2  ',1,'A','SHFLX over open ocn',phys_decomp)

      call addfld ('FSNSICE ','W/m2  ',1,'A','FSNS over sea ice',phys_decomp)
      call addfld ('FLNSICE ','W/m2  ',1,'A','FLNS over sea ice',phys_decomp)
      call addfld ('LHFLXICE','W/m2  ',1,'A','LHFLX over sea ice',phys_decomp)
      call addfld ('SHFLXICE','W/m2  ',1,'A','SHFLX over sea ice',phys_decomp)

! Fields to calculate ocean/ice forcing for som

      call addfld ('FSNSOI ','W/m2  ',1,'A','FSNS over open ocn and ice',phys_decomp)
      call addfld ('FLNSOI ','W/m2  ',1,'A','FLNS over open ocn and ice',phys_decomp)
      call addfld ('LHFLXOI','W/m2  ',1,'A','LHFLX over open ocn and ice',phys_decomp)
      call addfld ('SHFLXOI','W/m2  ',1,'A','SHFLX over open ocn and ice',phys_decomp)

! ISCCP fields

      call addfld ('FISCCP1 ','mixed   ',npres*ntau,'A', &
                   'grid box fraction covered by each ISCCP D level cloud types',phys_decomp,&
                   flag_xyfill=.true., flag_isccplev=.true.)
      call addfld ('TCLDAREA','fraction',1,'A','Total cloud area (ISCCP)',phys_decomp,flag_xyfill=.true.)
      call addfld ('HCLDAREA','fraction',1,'A','High cloud area (ISCCP)',phys_decomp,flag_xyfill=.true.)
      call addfld ('MCLDAREA','fraction',1,'A','Middle cloud area (ISCCP)',phys_decomp,flag_xyfill=.true.)
      call addfld ('LCLDAREA','fraction',1,'A','Low cloud area',phys_decomp,flag_xyfill=.true.)
      call addfld ('MEANPTOP','mb      ',1,'A','Mean cloud top pressure',phys_decomp,flag_xyfill=.true.)
      call addfld ('MEANTAU ','unitless',1,'A','Mean optical thickness',phys_decomp,flag_xyfill=.true.)
      call addfld ('MEANTTOP','K       ',1,'A','Mean cloud top temperature',phys_decomp,flag_xyfill=.true.)
      call addfld ('CLOUDY  ','fraction',1,'A','Binary flag for daytime cltot >=cldmin',phys_decomp)
!
! Add NSTEP to history tape for debugging
!
      call addfld ('NSTEP   ','timestep',1,'A','Model timestep',phys_decomp)
#ifdef FILLDEBUG
      call addfld ('FAKEFLD ','none'    ,plev,'A','1/nstep if nstep > x-index',phys_decomp,flag_xyfill=.true.)
      call addfld ('FAKEFLD1','none'    ,1,'A','1=day cloudy 0=day sunny',phys_decomp,flag_xyfill=.true.)
      call addfld ('FAKEFLD2','none'    ,1,'A','1=day cloudy',phys_decomp,flag_xyfill=.true.)
#endif
!
! Print master field list
!
      if (masterproc) then
         write(6,*)'BLDFLD:nfmaster=',nfmaster
         write(6,*)' ******* MASTER FIELD LIST *******'

         do j=1,nfmaster
            write(6,9000)j,masterlist(j)%field%name,masterlist(j)%field%units
9000        format (i5,1x,a8,1x,a8)
         end do
      end if

      return
   end subroutine bldfld

!#######################################################################

   subroutine addfld (fname, units, numlev, avgflag, long_name, &
                      decomp_type, flag_xyfill, flag_isccplev)
!----------------------------------------------------------------------- 
! 
! Purpose: Add a field to the master field list
! 
! Method: Put input arguments of field name, units, number of levels, averaging flag, and 
!         long name into a type entry in the global master field list (masterlist).
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ppgrid,    only: begchunk, endchunk
      use rgrid,     only: nlon
      use phys_grid, only: get_ncols_p, physgrid_set
      use dycore, only: dycore_is
!
! Arguments
!
      character(len=*), intent(in) :: fname      ! field name--should be 8 characters
      character(len=*), intent(in) :: units      ! units of fname--should be 8 chars
      character(len=1), intent(in) :: avgflag    ! averaging flag
      character(len=*), intent(in) :: long_name  ! long name of field
      
      integer, intent(in) :: numlev              ! number of vertical levels (dimension and loop)
      integer, intent(in) :: decomp_type         ! decomposition type

      logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue
      logical, intent(in), optional :: flag_isccplev ! levels are ISCCP levels not vertical
!
! Local workspace
!
      integer :: n     ! loop index
      integer :: c     ! chunk (physics ) or latitude (dynamics) index
      integer :: ncol  ! number of columns per chunk
!
! Ensure that required grid info is available now
!
      select case (decomp_type)
      case (phys_decomp)
         if (.not. physgrid_set) then
            write(6,*)'ADDFLD: Attempt to add field ',fname,' to masterlist before physics grid set'
            call endrun ()
         end if
      case (dyn_decomp)
         if (.not. dyngrid_set) then
            write(6,*)'ADDFLD: Attempt to add field ',fname,' to masterlist before dynamics grid set'
            call endrun ()
         end if
      end select
!
! Ensure that new field is not all blanks
!
      if (fname == ' ') then
         write(6,*)'ADDFLD: blank field name not allowed'
      end if
!
! Ensure that new field doesn't already exist
!
      do n=1,nfmaster
         if (masterlist(n)%field%name == fname) then
            write(6,*)'ADDFLD:', fname, ' already on list'
            call endrun ()
         end if
      end do
!
! Add field to Master Field List arrays fieldn and iflds
!
      nfmaster = nfmaster + 1  ! Increase number of fields on Master F.L.
!
! Check number of fields against max for master list.
!
      if (nfmaster > pflds) then
         write(6,*)'ADDFLD: Too many fields for primary history file -- pflds,nfmaster=', &
                   pflds,nfmaster
         stop 999
      end if
      masterlist(nfmaster)%field%name        = fname
      masterlist(nfmaster)%field%long_name   = long_name
      masterlist(nfmaster)%field%units       = units
      masterlist(nfmaster)%field%numlev      = numlev
      masterlist(nfmaster)%field%decomp_type = decomp_type
!
! Whether to apply xy fillvalue: default is false
!
      if (present(flag_xyfill)) then
         masterlist(nfmaster)%field%flag_xyfill = flag_xyfill
      else
         masterlist(nfmaster)%field%flag_xyfill = .false.
      end if
!
! Whether level dimension is ISCCP (currently 49) or CAM
!
      if (present(flag_isccplev)) then
         masterlist(nfmaster)%field%flag_isccplev = flag_isccplev
      else
         masterlist(nfmaster)%field%flag_isccplev = .false.
      end if

      if (masterlist(nfmaster)%field%flag_isccplev .and. numlev /= npres*ntau) then
         write(6,*)'ADDFLD: Number of ISCCP levels must be ',npres*ntau, ' got ', numlev
         call endrun ()
      end if
!
! Dimension history info based on decomposition type (dynamics or physics)
!
      select case (decomp_type)
      case (phys_decomp)
         masterlist(nfmaster)%field%begdim3  = begchunk
         masterlist(nfmaster)%field%enddim3  = endchunk
         masterlist(nfmaster)%field%begver   = 1
         masterlist(nfmaster)%field%endver   = numlev
         allocate (masterlist(nfmaster)%field%colperdim3(begchunk:endchunk))
         do c=begchunk,endchunk
            ncol = get_ncols_p(c)
            masterlist(nfmaster)%field%colperdim3(c) = ncol
         end do
         masterlist(nfmaster)%field%coldimin   = pcols
      case (dyn_decomp)
         masterlist(nfmaster)%field%begdim3  = beglat
         masterlist(nfmaster)%field%enddim3  = endlat

         if ( dycore_is('LR') )then
# if ( defined STAGGERED )
            select case (numlev)
            case (1)
               masterlist(nfmaster)%field%begver   = 1
               masterlist(nfmaster)%field%endver   = 1
            case (plev)
               masterlist(nfmaster)%field%begver   = beglev
               masterlist(nfmaster)%field%endver   = endlev
            case (plevp)
               masterlist(nfmaster)%field%begver   = beglev
               masterlist(nfmaster)%field%endver   = endlevp
            case default
               write(6,*)'ADDFLD: invalid number of levels=', numlev
               call endrun ()
            end select
# endif
         else
            masterlist(nfmaster)%field%begver   = 1
            masterlist(nfmaster)%field%endver   = numlev
         endif

         allocate (masterlist(nfmaster)%field%colperdim3(beglat:endlat))
         do c=beglat,endlat
            masterlist(nfmaster)%field%colperdim3(c) = nlon(c)
         end do
         masterlist(nfmaster)%field%coldimin = plon
      case default
         write(6,*)'ADDFLD: unknown decomp_type=', decomp_type
         call endrun ()
      end select
!
! These 2 fields are used only in master field list, not runtime field list
!
      masterlist(nfmaster)%avgflag(:) = avgflag
      masterlist(nfmaster)%actflag(:) = .false.

      select case (avgflag)
      case ('A')
         masterlist(nfmaster)%time_op(:) = 'mean'
      case ('I')
         masterlist(nfmaster)%time_op(:) = ' '
      case ('X')
         masterlist(nfmaster)%time_op(:) = 'maximum'
      case ('M')
         masterlist(nfmaster)%time_op(:) = 'minimum'
      case default
         write(6,*)'ADDFLD: unknown avgflag=',avgflag
         call endrun ()
      end select
      
      return
   end subroutine addfld

!#######################################################################

   subroutine wrapup ()
!-----------------------------------------------------------------------
!
! Purpose: 
! Close and dispose of history files.
! 
! Method: 
! This routine will close and dispose to Mass Store any full hist. files
! or any hist. file that has data on it when restart files are being 
! written.
! If a partially full history file was disposed (for restart 
! purposes), then wrapup will open that unit back up and position 
! it for appending new data. 
!
! Original version: CCM2
!
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pspect
      use ioFileMod
      use time_manager, only: get_nstep, get_curr_date, get_curr_time

#include <comctl.h>
#include <comlun.h>

!
! Local workspace
!
      integer :: nstep          ! current timestep number
      integer :: ncdate         ! current date in integer format [yyyymmdd]
      integer :: ncsec          ! time of day relative to current date [seconds]
      integer :: ndcur          ! days component of current time
      integer :: nscur          ! seconds component of current time
      integer :: yr, mon, day   ! year, month, day components of a date

      logical lfill   (ptapes)  ! Is history file ready to dispose?
      logical hdispose(ptapes)  ! Primary file disposed
      logical lremov            ! Rm file after dispose? 
      logical lhdisp            ! true => history file is disposed
      logical lhfill            ! true => history file is full

      integer t                 ! History file number

      real(r8) tday             ! Model day number for printout
!
! Externals
!
      logical, external :: rstwr ! whether it is time to write restarts
!-----------------------------------------------------------------------

      nstep = get_nstep()
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      call get_curr_time(ndcur, nscur)
!
!-----------------------------------------------------------------------
! Dispose history files.
!-----------------------------------------------------------------------
!
! Begin loop over mtapes (the no. of declared history files - primary
! and auxiliary).  This loop disposes a history file to Mass Store
! when appropriate.
!
      do t=1,mtapes

         hdispose(t) = .false.
         lfill(t) = .false.
!
! Find out if file is full
!
         if (hstwr(t) .and. nfils(t) >= mfilt(t)) then
            lfill(t) = .true.
         endif
!
! Dispose history file to mass store if 
!    1) file is filled or 
!    2) this is the end of run and file has data on it or
!    3) restarts are being put out and history file has data on it
!
         if (lfill(t) .or. (nlend .and. nfils(t) >= 1) .or. (rstwr() .and. nfils(t) >= 1)) then
!
! Dispose history file
!
            hdispose(t) = .true.
!
! Remove the disk history file unless this is the end of the run, or
! if the history file is not full.
!
            lremov = .true.
            if (nlend.or..not.lfill(t)) lremov = .false.
!
! Is this the 0 timestep data of a monthly run?
! If so, just close primary unit do not dispose. 
!
            if (masterproc) then
               write(6,*)'WRAPUP: nf_close(',nfid(t),')=',nhfil(t)
               call wrap_close (nfid(t))

               if (nhtfrq(t) /= 0 .or. nstep > 0) then

! Dispose history file.  

                  call putfil (nhfil(t), cpath(t), mss_wpass, mss_irt, lremov)
! 
! Print information concerning model output.
! Model day number = iteration number of history file data * delta-t / (seconds per day)
! 
                  tday = ndcur + nscur/86400._r8
                  if (t==1) then
                     write(6,*)'   Primary history file'
                  else
                     write(6,*)'   Auxiliary history file number ', t-1
                  end if
                  write(6,9003)nstep,nfils(t),tday
                  write(6,9004)
!                      
! Auxilary files may have been closed and saved off without being full. 
! We must reopen the files and position them for more data.
! Must position auxiliary files if not full
!              
                  if (.not.nlend .and. .not.lfill(t)) then
                     call wrap_open (nhfil(t), NF_WRITE, nfid(t))
                  end if
               endif                 ! if 0 timestep of montly run****
            end if                   ! if (masterproc)
         end if                      ! if time dispose history fiels***   
      end do                         ! do mtapes
!
! Reset number of files on each history tape
!
      do t=1,mtapes
         lhfill = hstwr(t) .and. nfils(t) >= mfilt(t)
         lhdisp = lhfill .or. (nlend .and. nfils(t) >= 1) .or. &
                              (rstwr() .and. nfils(t) >= 1)
         if (lhfill.and.lhdisp) then
            nfils(t) = 0
         endif
      end do
      
      return
9003  format('    Output at NSTEP     = ',i10,/, &
             '    Number of time samples on this file = ',i10,/, &
             '    Model Day           = ',f10.2)
9004  format('---------------------------------------')
   end subroutine wrapup

!#######################################################################

   subroutine allocate_hbuf2d (hbuf, dimind, hbuf_prec)
!-----------------------------------------------------------------------
!
! Purpose: Allocate memory for the 2-D history buffer of the right precision.
!          Nullify the other buffer.
!
!-----------------------------------------------------------------------
      type (hbuffer_2d) :: hbuf                   ! history buffer
      type (dim_index_2d) :: dimind               ! 2-D dimension index
      integer, intent(in) :: hbuf_prec            ! precision for history buffer

      if (hbuf_prec == 8) then
         nullify  (hbuf%buf4)
         allocate (hbuf%buf8(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2))
      else
         allocate (hbuf%buf4(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2))
         nullify  (hbuf%buf8)
      end if

   end subroutine allocate_hbuf2d

!#######################################################################

   subroutine allocate_hbuf3d (hbuf, dimind, hbuf_prec)
!-----------------------------------------------------------------------
!
! Purpose: Allocate memory for the 3-D history buffer of the right precision
!          Nullify the other buffer.
!
!-----------------------------------------------------------------------
      type (hbuffer_3d) :: hbuf                   ! history buffer
      type (dim_index_3d) :: dimind               ! 3-D dimension index
      integer, intent(in) :: hbuf_prec            ! precision for history buffer

      if (hbuf_prec == 8) then
         nullify  (hbuf%buf4)
         allocate (hbuf%buf8(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2, &
                             dimind%beg3:dimind%end3))
      else
         allocate (hbuf%buf4(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2, &
                             dimind%beg3:dimind%end3))
         nullify  (hbuf%buf8)
      end if

   end subroutine allocate_hbuf3d

!#######################################################################

   subroutine deallocate_hbuf2d (hbuf)
!-----------------------------------------------------------------------
!
! Purpose: Deallocate memory for 2-D history buffer
!
!-----------------------------------------------------------------------
      type (hbuffer_2d) :: hbuf                   ! history buffer

      if (associated(hbuf%buf8)) then
         deallocate (hbuf%buf8)
      else if (associated(hbuf%buf4)) then
         deallocate (hbuf%buf4)
      end if

   end subroutine deallocate_hbuf2d

!#######################################################################

   subroutine deallocate_hbuf3d (hbuf)
!-----------------------------------------------------------------------
!
! Purpose: Deallocate memory for 3-D history buffer
!
!-----------------------------------------------------------------------
      type (hbuffer_3d) :: hbuf                   ! history buffer

      if (associated(hbuf%buf8)) then
         deallocate (hbuf%buf8)
      else if (associated(hbuf%buf4)) then
         deallocate (hbuf%buf4)
      end if

   end subroutine deallocate_hbuf3d

!#######################################################################

   subroutine nullify_hbuf2d (hbuf)
!-----------------------------------------------------------------------
!
! Purpose: Nullify 2-D history buffer pointers
!
!-----------------------------------------------------------------------
      type (hbuffer_2d) :: hbuf                   ! history buffer

      nullify (hbuf%buf8, hbuf%buf4)

   end subroutine nullify_hbuf2d

!#######################################################################

   subroutine nullify_hbuf3d (hbuf)
!-----------------------------------------------------------------------
!
! Purpose: Nullify 3-D history buffer pointers
!
!-----------------------------------------------------------------------
      type (hbuffer_3d) :: hbuf                   ! history buffer

      nullify (hbuf%buf8, hbuf%buf4)

   end subroutine nullify_hbuf3d

!#######################################################################

   subroutine read_hbuf (hbuf, luhrest, ioerr)
!-----------------------------------------------------------------------
!
! Purpose: Read history information using the appropriate buffer
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(inout) :: hbuf    ! history buffer
      integer, intent(in ) :: luhrest
      integer, intent(out) :: ioerr

      if (associated(hbuf%buf8)) then
         read (luhrest,iostat=ioerr) hbuf%buf8
      else if (associated(hbuf%buf4)) then
         read (luhrest,iostat=ioerr) hbuf%buf4
      end if

   end subroutine read_hbuf

!#######################################################################

   subroutine write_hbuf (hbuf, luhrest, ioerr)
!-----------------------------------------------------------------------
!
! Purpose: Write (unformatted) history information to file unit luhrest
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(in) :: hbuf       ! history buffer
      integer, intent(in ) :: luhrest             ! file unit number
      integer, intent(out) :: ioerr               ! I/O error status

      if (associated(hbuf%buf8)) then
         write (luhrest,iostat=ioerr) hbuf%buf8
      else if (associated(hbuf%buf4)) then
         write (luhrest,iostat=ioerr) hbuf%buf4
      end if

   end subroutine write_hbuf

!#######################################################################

   subroutine hbuf_assigned_to_hbuf (hbuf1, hbuf2)
!-----------------------------------------------------------------------
!
! Purpose: Set hbuf1 to hbuf2 (copy the contents).
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(inout) :: hbuf1   ! history buffer
      type (hbuffer_3d), intent(in   ) :: hbuf2   ! history buffer

      if (associated(hbuf2%buf8)) then
         hbuf1%buf8 = hbuf2%buf8
      else if (associated(hbuf2%buf4)) then
         hbuf1%buf4 = hbuf2%buf4
      end if

   end subroutine hbuf_assigned_to_hbuf

!#######################################################################

   subroutine hbuf_assigned_to_real8 (hbuf, scalar)
!-----------------------------------------------------------------------
!
! Purpose: Set appropriate history buffer to the value, scalar
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(inout) :: hbuf    ! history buffer
      real(r8), intent(in) :: scalar              ! scalar real

      if (associated(hbuf%buf8)) then
         hbuf%buf8 = scalar
      else if (associated(hbuf%buf4)) then
         hbuf%buf4 = scalar
      end if

   end subroutine hbuf_assigned_to_real8

!#######################################################################

   subroutine set_hbuf_section_to_val (hbuf, dimind, scalar)
!-----------------------------------------------------------------------
!
! Purpose: Set section of appropriate history buffer to scalar
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(inout) :: hbuf    ! history buffer
      type (dim_index_3d), intent(in ) :: dimind  ! 3-D dimension index
      real(r8), intent(in) :: scalar              ! scalar real

      if (associated(hbuf%buf8)) then
         hbuf%buf8(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2, &
                   dimind%beg3:dimind%end3) = scalar
      else if (associated(hbuf%buf4)) then
         hbuf%buf4(dimind%beg1:dimind%end1, dimind%beg2:dimind%end2, &
                   dimind%beg3:dimind%end3) = scalar
      end if

   end subroutine set_hbuf_section_to_val

!#######################################################################

   subroutine assoc_hbuf2d_with_hbuf3d (hbuf2d, hbuf3d, c)
!-----------------------------------------------------------------------
!
! Purpose: Associate 2-D hbuf2d with 3-D hbuf3d of column c.
!          Nullify the other buffer of hbuf2d.
!
!-----------------------------------------------------------------------
      type (hbuffer_2d), intent(out) :: hbuf2d    ! 2-D history buffer
      type (hbuffer_3d), intent(in ) :: hbuf3d    ! 3-D history buffer
      integer, intent(in) :: c                    ! chunk (or lat) index

      if (associated(hbuf3d%buf8)) then
         hbuf2d%buf8 => hbuf3d%buf8(:,:,c)
         nullify (hbuf2d%buf4)
      else if (associated(hbuf3d%buf4)) then
         hbuf2d%buf4 => hbuf3d%buf4(:,:,c)
         nullify (hbuf2d%buf8)
      end if

   end subroutine assoc_hbuf2d_with_hbuf3d

!#######################################################################

   subroutine assoc_hbuf_with_nothing (hbuf, hbuf_prec)
!-----------------------------------------------------------------------
!
! Purpose: Associate the appropriate 3-D hbuf pointer with nothing.
!          Nullify the other pointer.
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(out) :: hbuf      ! 3-D history buffer
      integer, intent(in) :: hbuf_prec            ! hbuffer precision

      if (hbuf_prec == 8) then
         hbuf%buf8 => nothing_r8
         nullify (hbuf%buf4)
      else if (hbuf_prec == 4) then
         hbuf%buf4 => nothing_r4
         nullify (hbuf%buf8)
      end if
      return
   end subroutine assoc_hbuf_with_nothing

!#######################################################################

   subroutine xzy_to_xyz (xyzbuf, xzybuf, dimind)
!-----------------------------------------------------------------------
!
! Purpose: Set the appropriate buffer of xyzbuf to the transpose of
!          xzybuf with respect to the last two indices.
!
!-----------------------------------------------------------------------
   use rgrid,     only: nlon

   type (hbuffer_3d), intent(inout) :: xyzbuf      ! data in COORDS order
   type (hbuffer_3d), intent(in   ) :: xzybuf      ! data in xzy order
   type (dim_index_3d), intent(in ) :: dimind      ! dim index (xzybuf order)
   integer j, k                                    ! loop indices
   integer endi                                    ! ending longitude index

   if (associated(xzybuf%buf8)) then

      do k=dimind%beg2,dimind%end2
         do j=dimind%beg3,dimind%end3
            endi = nlon(j)
            xyzbuf%buf8(dimind%beg1:endi,j,k) = xzybuf%buf8(dimind%beg1:endi,k,j)
            if (endi < dimind%end1) then
               xyzbuf%buf8(endi+1:,j,k) = fillvalue
            end if
         end do
      end do

   else if (associated(xzybuf%buf4)) then

      do k=dimind%beg2,dimind%end2
         do j=dimind%beg3,dimind%end3
            endi = nlon(j)
            xyzbuf%buf4(dimind%beg1:endi,j,k) = xzybuf%buf4(dimind%beg1:endi,k,j)
            if (endi < dimind%end1) then
               xyzbuf%buf4(endi+1:,j,k) = fillvalue
            end if
         end do
      end do

   end if
   return
   end subroutine xzy_to_xyz

!#######################################################################

   subroutine fill_unset (buf4, buf8, dimind)
!-----------------------------------------------------------------------
!
! Purpose: For max/min calcs which initialized hbuf to +-huge, set those values to fillvalue
!
!-----------------------------------------------------------------------
   use rgrid,     only: nlon

   real(r4), pointer :: buf4(:,:,:)
   real(r8), pointer :: buf8(:,:,:)
   type (dim_index_3d), intent(in ) :: dimind      ! dim index (xzybuf order)
   integer j, k                                    ! loop indices
   integer endi                                    ! ending longitude index

   if (associated(buf8)) then

      do k=dimind%beg2,dimind%end2
         do j=dimind%beg3,dimind%end3
            endi = nlon(j)
            do i=dimind%beg1,endi
               if (buf8(i,j,k) == -huge(buf8) .or. buf8(i,j,k) == +huge(buf8)) then
                  buf8(i,j,k) = fillvalue
               end if
            end do
         end do
      end do

   else if (associated(buf4)) then

      do k=dimind%beg2,dimind%end2
         do j=dimind%beg3,dimind%end3
            endi = nlon(j)
            do i=dimind%beg1,endi
               if (buf4(i,j,k) == -huge(buf4) .or. buf4(i,j,k) == +huge(buf4)) then
                  buf4(i,j,k) = fillvalue
               end if
            end do
         end do
      end do

   end if

   return
   end subroutine fill_unset

!#######################################################################

   subroutine hbuf_accum_inst (buf4, buf8, field, nacs, dimind, idim, flag_xyfill)
!-----------------------------------------------------------------------
!
! Purpose: Accumulate instantaneous values of field in 2-D hbuf.
!          Set accumulation counter to 1.
!
!-----------------------------------------------------------------------
      real(r4), pointer :: buf4(:,:)    ! 2-D history buffer
      real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
      integer, pointer                 :: nacs(:) ! accumulation counter
      integer, intent(in)              :: idim    ! Longitude dimension of field array
      logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
      real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
      type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
!
! Local indices
!
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer :: ieu, jeu  ! number of elements in each dimension
      integer :: i, k      ! loop indices

      logical :: bad       ! flag indicates input field fillvalues not applied consistently
                           ! with vertical level

      ib = dimind%beg1
      ie = dimind%end1
      jb = dimind%beg2
      je = dimind%end2

      ieu = ie-ib+1
      jeu = je-jb+1

      if (associated(buf8)) then
         do k=1,jeu
            do i=1,ieu
               buf8(i,k) = field(i,k)
            end do
         end do
      else if (associated(buf4)) then
         do k=1,jeu
            do i=1,ieu
               buf4(i,k) = field(i,k)
            end do
         end do
      end if

      if (flag_xyfill) then
         do i=1,ieu
            if (field(i,1) == fillvalue) then
               nacs(i) = 0
            else
               nacs(i) = 1
            end if
            call check_accum (field, idim, ieu, jeu)
         end do
      else
         do i=1,ieu
            nacs(i) = 1
         end do
      end if

      return
   end subroutine hbuf_accum_inst

!#######################################################################

   subroutine check_accum (field, idim, ieu, jeu)
      integer, intent(in)  :: idim
      real(r8), intent(in) :: field(idim,*)   ! real*8 array
      integer, intent(in)  :: ieu,jeu         ! loop ranges

      logical :: bad
      integer :: i,k
!
! For multilevel fields ensure that all levels have fillvalue applied consistently
!
      bad = .false.
      do k=2,jeu
         do i=1,ieu
            if (field(i,1) == fillvalue .and. field(i,k) /= fillvalue .or. &
                field(i,1) /= fillvalue .and. field(i,k) == fillvalue) then
!         write(6,*)'CHECK_ACCUM: inconsistent level values for i,k=',i,k
!		write(6,'(a,2g12.2)')'field(i,1),field(i,k)=',field(i,1),field(i,k)
                bad = .true.
            end if
         end do
      end do

      if (bad) then
         write(6,*)'CHECK_ACCUM: inconsistent level values'
         call endrun ()
      end if

      return
   end subroutine check_accum

!#######################################################################

   subroutine hbuf_accum_add (buf4, buf8, field, nacs, dimind, idim, flag_xyfill)
!-----------------------------------------------------------------------
!
! Purpose: Add the values of field to 2-D hbuf.
!          Increment accumulation counter by 1.
!
!-----------------------------------------------------------------------
      real(r4), pointer :: buf4(:,:)    ! 2-D history buffer
      real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
      integer, pointer                 :: nacs(:) ! accumulation counter
      integer, intent(in) :: idim           ! Longitude dimension of field array
      logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
      real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
      type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
!
! Local indices
!
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer :: ieu, jeu  ! number of elements in each dimension
      integer :: i,k       ! indices

      ib = dimind%beg1
      ie = dimind%end1
      jb = dimind%beg2
      je = dimind%end2

      ieu = ie-ib+1
      jeu = je-jb+1

      if (flag_xyfill) then
         if (associated(buf8)) then
            do k=1,jeu
               do i=1,ieu
                  if (field(i,k) /= fillvalue) then
                     buf8(i,k) = buf8(i,k) + field(i,k)
                  end if
               end do
            end do
         else if (associated(buf4)) then
            do k=1,jeu
               do i=1,ieu
                  if (field(i,k) /= fillvalue) then
                     buf4(i,k) = buf4(i,k) + field(i,k)
                  end if
               end do
            end do
         end if
!
! Ensure input field has fillvalue defined invariant in the z-direction, then increment nacs
!
         call check_accum (field, idim, ieu, jeu)
         do i=1,ieu
            if (field(i,1) /= fillvalue) then
               nacs(i) = nacs(i) + 1
            end if
         end do
      else
         if (associated(buf8)) then
            do k=1,jeu
               do i=1,ieu
                  buf8(i,k) = buf8(i,k) + field(i,k)
               end do
            end do
         else if (associated(buf4)) then
            do k=1,jeu
               do i=1,ieu
                  buf4(i,k) = buf4(i,k) + field(i,k)
               end do
            end do
         end if
         do i=1,ieu
            nacs(i) = nacs(i) + 1
         end do
      end if

      return
   end subroutine hbuf_accum_add

!#######################################################################

   subroutine hbuf_accum_max (buf4, buf8, field, nacs, dimind, idim, flag_xyfill)
!-----------------------------------------------------------------------
!
! Purpose: Accumulate the maximum values of field in 2-D hbuf
!          Set accumulation counter to 1.
!
!-----------------------------------------------------------------------
      real(r4), pointer :: buf4(:,:)    ! 2-D history buffer
      real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
      integer, pointer                 :: nacs(:) ! accumulation counter
      integer, intent(in) :: idim           ! Longitude dimension of field array
      logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
      real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
      type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
!
! Local indices
!
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer :: ieu, jeu  ! number of elements in each dimension
      integer :: i, k

      ib = dimind%beg1
      ie = dimind%end1
      jb = dimind%beg2
      je = dimind%end2

      ieu = ie-ib+1
      jeu = je-jb+1

      if (associated(buf8)) then

         if (flag_xyfill) then
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf8(i,k) = -huge (buf8)
                  end if
                  if (field(i,k) > buf8(i,k) .and. field(i,k) /= fillvalue) then
                     buf8(i,k) = field(i,k)
                  end if
               end do
            end do
         else
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf8(i,k) = -huge (buf8)
                  end if
                  if (field(i,k) > buf8(i,k)) then
                     buf8(i,k) = field(i,k)
                  end if
               end do
            end do
         end if

      else if (associated(buf4)) then

         if (flag_xyfill) then
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf4(i,k) = -huge (buf4)
                  end if
                  if (field(i,k) > buf4(i,k) .and. field(i,k) /= fillvalue) then
                     buf4(i,k) = field(i,k)
                  end if
               end do
            end do
         else
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf4(i,k) = -huge (buf4)
                  end if
                  if (field(i,k) > buf4(i,k)) then
                     buf4(i,k) = field(i,k)
                  end if
               end do
            end do
         end if
      end if

      if (flag_xyfill) then
         call check_accum (field, idim, ieu, jeu)
         do i=1,ieu
            if (field(i,1) /= fillvalue) then
               nacs(i) = 1
            end if
         end do
      else
         do i=1,ieu
            nacs(i) = 1
         end do
      end if

      return
   end subroutine hbuf_accum_max

!#######################################################################

   subroutine hbuf_accum_min (buf4, buf8, field, nacs, dimind, idim, flag_xyfill)
!-----------------------------------------------------------------------
!
! Purpose: Accumulate the minimum values of field in 2-D hbuf
!          Set accumulation counter to 1.
!
!-----------------------------------------------------------------------
      real(r4), pointer :: buf4(:,:)    ! 2-D history buffer
      real(r8), pointer :: buf8(:,:)    ! 2-D history buffer
      integer, pointer                 :: nacs(:) ! accumulation counter
      integer, intent(in) :: idim           ! Longitude dimension of field array
      logical, intent(in)              :: flag_xyfill ! non-applicable xy points flagged with fillvalue
      real(r8),          intent(in ) :: field(idim,*)   ! real*8 array
      type (dim_index_2d), intent(in ) :: dimind  ! 2-D dimension index
!
! Local indices
!
      integer :: ib, ie    ! beginning and ending indices of first dimension
      integer :: jb, je    ! beginning and ending indices of second dimension
      integer :: ieu, jeu  ! number of elements in each dimension
      integer :: i, k

      ib = dimind%beg1
      ie = dimind%end1
      jb = dimind%beg2
      je = dimind%end2

      ieu = ie-ib+1
      jeu = je-jb+1

      if (associated(buf8)) then

         if (flag_xyfill) then
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf8(i,k) = +huge (buf8)
                  end if
                  if (field(i,k) < buf8(i,k) .and. field(i,k) /= fillvalue) then
                     buf8(i,k) = field(i,k)
                  end if
               end do
            end do
         else
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf8(i,k) = +huge (buf8)
                  end if
                  if (field(i,k) < buf8(i,k)) then
                     buf8(i,k) = field(i,k)
                  end if
               end do
            end do
         end if

      else if (associated(buf4)) then

         if (flag_xyfill) then
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf4(i,k) = +huge (buf4)
                  end if
                  if (field(i,k) < buf4(i,k) .and. field(i,k) /= fillvalue) then
                     buf4(i,k) = field(i,k)
                  end if
               end do
            end do
         else
            do k=1,jeu
               do i=1,ieu
                  if (nacs(i) == 0) then
                     buf4(i,k) = +huge (buf4)
                  end if
                  if (field(i,k) < buf4(i,k)) then
                     buf4(i,k) = field(i,k)
                  end if
               end do
            end do
         end if
      end if

      if (flag_xyfill) then
         call check_accum (field, idim, ieu, jeu)
         do i=1,ieu
            if (field(i,1) /= fillvalue) then
               nacs(i) = 1
            end if
         end do
      else
         do i=1,ieu
            nacs(i) = 1
         end do
      end if

      return
   end subroutine hbuf_accum_min

!#######################################################################

   subroutine readin_hbuf (iu, hbuf, numperlat)
!-----------------------------------------------------------------------
!
! Purpose: Read history restart binary file 
!
!-----------------------------------------------------------------------
      use binary_io
      integer, intent(in) :: iu                   ! Logical unit
      integer, intent(in) :: numperlat            ! Length of arr
      type (hbuffer_3d), intent(inout) :: hbuf    ! hbuf to be read

      if (associated(hbuf%buf8)) then
         call readin_r8 (iu, hbuf%buf8, numperlat)
      else if (associated(hbuf%buf4)) then
         call readin_r4 (iu, hbuf%buf4, numperlat)
      end if

   end subroutine readin_hbuf

!#######################################################################

#ifdef SPMD
   subroutine mpigatherv_hbuf (hbuf_send, numsend, mpireal1, hbuf_recv, &
                               numrecv, displs, mpireal2, root, comm)
!-----------------------------------------------------------------------
!
! Purpose: Collects history buffer from each thread on masterproc
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(in   ) :: hbuf_send  ! send buffer
      type (hbuffer_3d), intent(inout) :: hbuf_recv  ! receive buffer
      integer :: numsend                ! number of items to be sent
      integer :: mpireal1               ! MPI real data type for hbuf_send
      integer :: mpireal2               ! MPI real data type for hbuf_recv
      integer :: numrecv(*)             ! number of items to be received
      integer :: displs(*)              ! displacement array
      integer, intent(in) :: root
      integer, intent(in) :: comm
 
      if (associated(hbuf_send%buf8)) then
         call mpigatherv (hbuf_send%buf8, numsend, mpireal1, hbuf_recv%buf8, numrecv, &
                          displs, mpireal2, root, comm)
      else if (associated(hbuf_send%buf4)) then
         call mpigatherv (hbuf_send%buf4, numsend, mpireal1, hbuf_recv%buf4, numrecv, &
                          displs, mpireal2, root, comm)
      end if

   end subroutine mpigatherv_hbuf
#endif

!#######################################################################

   subroutine wrap_put_vara_hbuf (nfid, varid, start, count, hbuf)
!-----------------------------------------------------------------------
!
! Purpose: Output the given portion of the real array.
!
!-----------------------------------------------------------------------
      type (hbuffer_3d), intent(in) :: hbuf    ! buffer to write
      integer, intent(in):: nfid               ! file ID
      integer, intent(in):: varid              ! variable ID
      integer, intent(in):: start(*)           ! array of starting field indices
      integer, intent(in):: count(*)           ! array of count values

      if (associated(hbuf%buf8)) then
         call wrap_put_vara_realx (nfid, varid, start, count, hbuf%buf8)
      else if (associated(hbuf%buf4)) then
         call wrap_put_vara_real  (nfid, varid, start, count, hbuf%buf4)
      end if

   end subroutine wrap_put_vara_hbuf

!#######################################################################

   subroutine gather_chunk_to_field_hbuf (fdim,mdim,ldim, &
                                          nlond,localchunks,globalfield)
!-----------------------------------------------------------------------
!
! Purpose: Reconstruct longitude/latitude field
!          from decomposed chunk data structure
!
!-----------------------------------------------------------------------
      use phys_grid, only: gather_chunk_to_field, gather_chunk_to_field4

      type (hbuffer_3d), intent(in   ) :: localchunks
      type (hbuffer_3d), intent(inout) :: globalfield
      integer, intent(in) :: fdim      ! declared length of first dimension
      integer, intent(in) :: mdim      ! declared length of middle dimension
      integer, intent(in) :: ldim      ! declared length of last dimension
      integer, intent(in) :: nlond     ! declared number of longitudes

      if (associated(localchunks%buf8)) then
         call gather_chunk_to_field (fdim,mdim,ldim, &
                                 nlond,localchunks%buf8,globalfield%buf8)
      else if (associated(localchunks%buf4)) then
         call gather_chunk_to_field4(fdim,mdim,ldim, &
                                 nlond,localchunks%buf4,globalfield%buf4)
      end if

   end subroutine gather_chunk_to_field_hbuf

!#######################################################################

   subroutine scatter_field_to_chunk_hbuf (fdim,mdim,ldim, &
                                      nlond,globalfield,localchunks)
!-----------------------------------------------------------------------
!
! Purpose: Reconstruct longitude/latitude field
!          from decomposed chunk data structure
!
!-----------------------------------------------------------------------
      use phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk4

      type (hbuffer_3d), intent(inout) :: localchunks
      type (hbuffer_3d), intent(in   ) :: globalfield
      integer, intent(in) :: fdim      ! declared length of first dimension
      integer, intent(in) :: mdim      ! declared length of middle dimension
      integer, intent(in) :: ldim      ! declared length of last dimension
      integer, intent(in) :: nlond     ! declared number of longitudes

      if (associated(globalfield%buf8)) then
         call scatter_field_to_chunk (fdim,mdim,ldim, &
                                   nlond,globalfield%buf8,localchunks%buf8)
      else if (associated(globalfield%buf4)) then
         call scatter_field_to_chunk4(fdim,mdim,ldim, &
                                   nlond,globalfield%buf4,localchunks%buf4)
      end if

   end subroutine scatter_field_to_chunk_hbuf

end module history
