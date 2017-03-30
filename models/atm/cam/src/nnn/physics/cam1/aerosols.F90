#include <misc.h>
#include <params.h>

module aerosols
!----------------------------------------------------------------------- 
! 
! Purposes: 
!       read, store, interpolate, and return fields
!         of aerosols to CAM.  The initialization
!         file (mass.nc) is assumed to be a monthly climatology
!         of aerosols from MATCH (on a sigma pressure
!         coordinate system).
!       also provide a "background" aerosol field to correct
!         for any deficiencies in the physical parameterizations
!         This fields is a "tuning" parameter.
!       Public methods:
!       (1) - initialization
!          read aerosol masses from external file
!             also pressure coordinates
!          convert from monthly average values to mid-month values
!       (2) - interpolation (time and vertical)
!          interpolate onto pressure levels of CAM
!          interpolate to time step of CAM
!          return mmr's of aerosols 
!          provide background aerosol based on "tauback"
!       (3) - diagnostics
!          write out various diagnostic fields
! 
! Calling Heirarchy
!   aerosol_initialize (public)
!     add_aer_diag_fields
!   get_aerosol (public)
!     vert_interpolate
!     ramp_sulfate
!     background
!     scale_aerosols
!   aerosol_diagnostics (public)
!   aerosol_indirect (public)
!   get_rf_scales (public)
!   get_int_scales (public)
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid, only: pcols, pver, pverp, begchunk, endchunk
   use pmgrid, only: plon, plat, plev, plevp, masterproc
   use phys_grid, only: get_ncols_p, scatter_field_to_chunk
   use time_manager, only: get_curr_calday
   use infnan, only: inf, bigint

   implicit none

! naer_all is total number of species
! naer is number of species in climatology
! naer_all = naer + 1 (background "species")

  integer, public, parameter :: naer_all = 11
  integer, private, parameter :: naer = 10 

! indices to aerosol array (species portion)

  integer, public, parameter :: &
      idxSUL   =  1, &
      idxSSLT  =  2, &
      idxOCPHO =  7, &
      idxBCPHO =  8, &
      idxOCPHI =  9, &
      idxBCPHI = 10, &
      idxBG    = 11

! indices to sections of array that represent 
! groups of aerosols

  integer, public, parameter :: &
      idxDUSTfirst    = 3, &
      numDUST         = 4, &
      idxCARBONfirst = 7, &
      numCARBON      = 4

! names of aerosols are they are represented in
! the climatology file and for the purposes
! of outputting mmr's.  Appended '_V' indicates field has been vertically summed.

  character(len=8), public, parameter :: aerosol_name(naer_all) =  &
     (/"MSUL_V  "&
      ,"MSSLT_V "&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"&
      ,"Backgrnd"/)
! 
! default values for namelist variables
!

! compute radiative forcing per aerosol
  logical, public :: radforce   = .false. 

! portion of each species group to use in computation
! of relative radiative forcing.

  real(r8), public :: sulscl_rf  = 0._r8 ! 
  real(r8), public :: carscl_rf  = 0._r8
  real(r8), public :: ssltscl_rf = 0._r8
  real(r8), public :: dustscl_rf = 0._r8
  real(r8), public :: bgscl_rf   = 0._r8

! "background" aerosol species mmr.
  real(r8), public :: tauback = 0._r8

! portion of each species group to use in computation
! of aerosol forcing in driving the climate
  real(r8), public :: sulscl  = 1._r8
  real(r8), public :: carscl  = 1._r8
  real(r8), public :: ssltscl = 1._r8
  real(r8), public :: dustscl = 1._r8

  private
  save

  integer :: aernid = -1           ! netcdf id for aerosol file (init to invalid)
  integer :: species_id(naer) = -1 ! netcdf_id of each aerosol species (init to invalid)
  integer :: Mpsid                 ! netcdf id for MATCH PS
  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = bigint       ! index to nxt month in file

  real(r8) :: cdaym = inf          ! calendar day of prv month
  real(r8) :: cdayp = inf          ! calendar day of next month
  real(r8) :: Mid(12)              ! Days into year for mid month date
  data Mid/16.5, 46.0, 75.5, 106.0, 136.5, 167.0, 197.5, 228.5, 259.0, 289.5, 320.0, 350.5 /
  
  public aerosol_initialize  ! read from file, interpolate onto horiz grid
  public get_aerosol         ! interpolate onto pressure levels, and time
  public aerosol_diagnostics ! write out various data about aerosols
  public aerosol_indirect    ! compute change in effective liqw radius from sulfate
  public get_rf_scales       ! compute scale factors for forcing computation
  public get_int_scales      ! compute scale factors for interactive comps
  public aerint              ! read next month info
!
!  values read from file and temporary values used for interpolation
!
!  AEROSOLc is:
!  Cumulative Mass at midpoint of each month
!    on CAM's horizontal grid (col)
!    on MATCH's levels (lev)
!  AEROSOLc(month_ind, col_ind, Match_lev_ind, chunk_ind, species_ind) 
!
  integer, parameter :: paerlev = 28           ! number of levels for aerosol fields (MUST = naerlev)
  integer naerlev                              ! size of level dimension in MATCH data
  real(r8) :: M_hybi(paerlev+1)                ! MATCH hybi
  real(r8) :: M_ps(plon,plat)                  ! surface pressure from MATCH file
  real(r8), allocatable :: AEROSOLc(:,:,:,:,:) ! Aerosol cumulative mass from MATCH
  real(r8), allocatable :: M_ps_cam_col(:,:,:) ! PS from MATCH on Cam Columns
  real(r8), allocatable :: Match_ps_chunk(:,:) ! surface pressure array section at a particular time step

  integer :: mxaerl                            ! Maximum level of background aerosol

contains

subroutine aerosol_initialize
!------------------------------------------------------------------
!  Reads in:
!     file from which to read aerosol Masses on CAM grid. Currently
!        assumed to be MATCH ncep runs, averaged by month.
!     NOTE (Data have been externally interpolated onto CAM grid 
!        and backsolved to provide Mid-month values)
!     
!  Populates:
!     module variables:
!       AEROSOLc(pcols,paerlev+1,begchunk:endchunk,naer,2))
!       AEROSOLc(  column_index
!                , level_index (match levels)
!                , chunk_index 
!                , species_index
!                , month = 1:2 )
!       M_hybi(level_index = Lev_MATCH) = pressure at mid-level.
!       M_ps_cam_col(column,chunk,month) ! PS from MATCH on Cam Columns
!
!  Method:
!    read data from file
!    allocate memory for storage of aerosol data on CAM horizontal grid
!    distribute data to remote nodes
!    populates the module variables
!
!------------------------------------------------------------------
   use ioFileMod, only: getfil
   use filenames, only: bndtvaer
#if ( defined SPMD )
   use mpishorthand
#endif

   include 'netcdf.inc'

! local variables

   integer :: naerlon
   integer :: naerlat
   integer :: naerlev

   integer dateid                       ! netcdf id for date variable
   integer secid                        ! netcdf id for seconds variable
   integer londimid                     ! netcdf id for longitude dimension
   integer latdimid                     ! netcdf id for latitude dimension
   integer levdimid                     ! netcdf id for level dimension
   integer lonid                        ! netcdf id for longitude variable

   integer timesiz                      ! number of time samples (=12) in netcdf file
   integer latid                        ! netcdf id for latitude variable
   integer Mhybiid                      ! netcdf id for MATCH hybi
   integer timeid                       ! netcdf id for time variable
   integer dimids(nf_max_var_dims)      ! variable shape
   integer :: start(4)                  ! start vector for netcdf calls
   integer :: kount(4)                  ! count vector for netcdf calls
   integer mo                           ! month index
   integer m                            ! constituent index
   integer :: n                         ! loop index
   integer :: i,j,k                     ! spatial indices
   integer :: date_aer(12)              ! Date on aerosol dataset (YYYYMMDD)
   integer :: attnum                    ! attribute number
   integer mo_prv                       ! index to previous month

   character(len=256) :: locfn          ! netcdf local filename to open
!
! aerosol_data will be read in from the aerosol boundary dataset, then scattered to chunks
! after filling in the bottom level with zeros
! 
   real(r8) :: aerosol_data(plon,plat,paerlev)    ! aerosol field read in from dataset
   real(r8) :: aerosol_field(plon,paerlev+1,plat) ! aerosol field to be scattered
   real(r8) :: caldayloc                          ! calendar day of current timestep

   call t_startf ('aerosol_init')
   call print_memusage ('Start aerosol_initialize')
!
! Allocate memory for dynamic arrays local to this module
!
   allocate (AEROSOLc(pcols,paerlev+1,begchunk:endchunk,naer,2))
   allocate (M_ps_cam_col(pcols,begchunk:endchunk,2))
   allocate (Match_ps_chunk(pcols,begchunk:endchunk))

! TBH:  HACK to avoid use of uninitialized values when ncols < pcols
  AEROSOLc(:,:,:,:,:) = 0.
  M_ps_cam_col(:,:,:) = 0.
  Match_ps_chunk(:,:) = 0.
! TBH:  END HACK

!
!  Add fields to master list for output
!
   call add_aer_diag_fields ()
!
! compute mxaerl
!
   call setmxaerl ()

   if (masterproc) then
!
! find and open file; abort if fail (getfil(,,0)).
!
      call getfil (bndtvaer, locfn, 0)
      call wrap_open (locfn, 0, aernid)
      write(6,*)'aerosol_initialize: reading aerosol dataset....  If this seems to be taking too'
      write(6,*)'long, perhaps the dataset is being read from an nfs-mounted filesystem.'
!
! First ensure dataset is CAM-ready
!
      if (nf_inq_attid (aernid, nf_global, 'cam-ready', attnum) /= nf_noerr) then
         write(6,*)'aerosol_initialize: global attribute cam-ready not found: aborting'
         write(6,*)'interpaerosols needs to be run to create a cam-ready aerosol dataset'
         call endrun ()
      end if
!
! Get and check dimension info
!
      call wrap_inq_dimid( aernid,  'lon', londimid )
      call wrap_inq_dimid( aernid,  'lev', levdimid )
      call wrap_inq_dimid( aernid, 'time',   timeid )
      call wrap_inq_dimid( aernid,  'lat', latdimid )

      call wrap_inq_dimlen( aernid, londimid, naerlon )
      call wrap_inq_dimlen( aernid, levdimid, naerlev )
      call wrap_inq_dimlen( aernid, latdimid, naerlat )
      call wrap_inq_dimlen( aernid,   timeid, timesiz )

      if (naerlon /= plon .or. naerlat /= plat .or. naerlev /= paerlev .or. timesiz /= 12) then
         write(6,*)'aerosol_initialize: grid mismatch'
         write(6,*)'Model:   plon,    plat,    naerlev, 12     =', plon, plat, naerlev, 12
         write(6,*)'Dataset: naerlon, naerlat, paerlev, timesiz=', naerlon, naerlat, naerlev, timesiz
         call endrun ()
      end if

      call wrap_inq_varid( aernid, 'date',   dateid )
      call wrap_inq_varid( aernid, 'datesec', secid )
      do m = 1, naer
         call wrap_inq_varid( aernid, TRIM(aerosol_name(m)), species_id(m))
      end do

      call wrap_inq_varid( aernid, 'lon', lonid   )
      call wrap_inq_varid( aernid, 'lat', latid   )
!
! quick sanity check on one field
!
      call wrap_inq_vardimid (aernid, species_id(1), dimids)
      if ( (dimids(4) /= timeid) .or. &
           (dimids(3) /= levdimid) .or. &
           (dimids(2) /= latdimid) .or. &
           (dimids(1) /= londimid) ) then
         write(6,*)'AEROSOL READ: Data must be ordered time, lev, lat, lon'
         write(6,*)'data are       ordered as', dimids(4), dimids(3), dimids(2), dimids(1)
         write(6,*)'data should be ordered as', timeid, levdimid, latdimid, londimid
         call endrun ()
      end if
!
! use hybi,PS from MATCH
!
      call wrap_inq_varid( aernid, 'hybi', Mhybiid   )
      call wrap_inq_varid( aernid, 'PS', Mpsid   )
!
! check dimension order for MATCH's surface pressure
!
      call wrap_inq_vardimid (aernid, Mpsid, dimids)
      if ( (dimids(3) /= timeid) .or. &
           (dimids(2) /= latdimid) .or. &
           (dimids(1) /= londimid) ) then
         write(6,*)'AEROSOL READ: Pressure must be ordered time, lat, lon'
         write(6,*)'data are       ordered as', dimids(3), dimids(2), dimids(1)
         write(6,*)'data should be ordered as', timeid, levdimid, latdimid, londimid
         call endrun ()
      end if
! 
! read in hybi from MATCH
!
      call wrap_get_var_realx (aernid, Mhybiid, M_hybi)
!
! Retrieve date and sec variables.
!
      call wrap_get_var_int (aernid, dateid, date_aer)
      if (timesiz < 12) then
         write(6,*)'AEROSOL READ: When cycling aerosols, dataset must have 12 consecutive ', &
                   'months of data starting with Jan'
         write(6,*)'Current dataset has only ',timesiz,' months'
         call endrun ()
      end if
      do mo = 1,12
         if (mod(date_aer(mo),10000)/100 /= mo) then
            write(6,*)'AEROSOL READ: When cycling aerosols, dataset must have 12 consecutive ', &
                      'months of data starting with Jan'
            write(6,*)'Month ',mo,' of dataset says date=',date_aer(mo)
            call endrun ()
         end if
      end do
   end if          ! masterproc
!
! broadcast hybi to nodes
!
#if ( defined SPMD )
   call mpibcast (M_hybi, paerlev+1, mpir8, 0, mpicom)
#endif

   caldayloc = get_curr_calday ()
  
   if (caldayloc < Mid(1)) then
      mo_prv = 12
      mo_nxt =  1
   else if (caldayloc >= Mid(12)) then
      mo_prv = 12
      mo_nxt =  1
   else
      do i = 2 , 12
         if (caldayloc < Mid(i)) then
            mo_prv = i-1
            mo_nxt = i
            exit
         end if
      end do
   end if
!
! Set initial calendar day values
!
   cdaym = Mid(mo_prv)
   cdayp = Mid(mo_nxt)
!
! Retrieve Aerosol Masses (kg/m^2 in each layer), transpose to model order (lon,lev,lat),
! then scatter to slaves.
!
   if (nm /= 1 .or. np /= 2) call endrun ()
   do n=nm,np
      if (n == 1) then
         mo = mo_prv
      else
         mo = mo_nxt
      end if

      do m=1,naer
         if (masterproc) then
            start(:) = (/1,1,1,mo/)
            kount(:) = (/plon,plat,paerlev,1/)
            call wrap_get_vara_realx (aernid, species_id(m), start, kount, aerosol_data(1,1,1))
            do j=1,plat
               do k=1,paerlev
                  aerosol_field(:,k,j) = aerosol_data(:,j,k)
               end do
               aerosol_field(:,paerlev+1,j) = 0.   ! value at bottom
            end do
         end if
         call scatter_field_to_chunk (1, paerlev+1, 1, plon, aerosol_field, AEROSOLc(1,1,begchunk,m,n))
      end do
!
! Retrieve PS from Match
!
      if (masterproc) then
         start(:) = (/1,1,mo,-1/)
         kount(:) = (/plon,plat,1,-1/)
         call wrap_get_vara_realx (aernid, Mpsid, start, kount, M_ps(1,1))
      end if
      call scatter_field_to_chunk (1, 1, 1, plon, M_ps(1,1), M_ps_cam_col(1,begchunk,n))
   end do     ! n=nm,np (=1,2)

   call print_memusage('End aerosol_initialize')
   call t_stopf ('aerosol_init')

   return
end subroutine aerosol_initialize

subroutine get_aerosol(c, ncol1, ncol, pint, AEROSOLt, scale)
!------------------------------------------------------------------
!
!  Input:
!     time at which aerosol mmrs are needed (get_curr_calday())
!     chunk index
!     CAM's vertical grid (pint)
!
!  Output:
!     values for Aerosol Mass Mixing Ratios at specified time
!     on vertical grid specified by CAM (AEROSOLt)
!
!  Method:
!     first determine which indexs of aerosols are the bounding data sets
!     interpolate both onto vertical grid aerm(),aerp().
!     from those two, interpolate in time.
!
!------------------------------------------------------------------

   implicit none
!
! aerosol fields interpolated to current time step
!   on pressure levels of this time step.
! these should be made read-only for other modules
! Is allocation done correctly here?
!
   integer, intent(in) :: c                   ! Chunk Id.
   integer, intent(in) :: ncol1,ncol          ! Column index (first & last)
   real(r8), intent(in) :: pint(pcols,pverp)  ! midpoint pres.
   real(r8), intent(in) :: scale(naer_all)    ! scale each aerosol by this amount

   real(r8), intent(out) :: AEROSOLt(pcols, pver, naer_all) ! aerosols
!
! Local workspace
!
   real(r8) caldayloc                     ! calendar day of current timestep
   real(r8) deltat, fact2, fact1          ! time interpolation factors

   integer i, k, j                        ! spatial indices
   integer m                              ! constituent index
   integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
!   integer ncol                           ! number of columns
   
   real(r8) speciesmin(naer)              ! minimal value for each species
!
! values before current time step "the minus month"
! aerosolm(pcols,pver) is value of preceeding month's aerosol mmr
! aerosolp(pcols,pver) is value of next month's aerosol mmr
!  (think minus and plus or values to left and right of point to be interpolated)
!
   real(r8) AEROSOLm(pcols,pver,naer) ! aerosol mmr from MATCH in column at previous (minus) month
!
! values beyond (or at) current time step "the plus month"
!
   real(r8) AEROSOLp(pcols,pver,naer) ! aerosol mmr from MATCH in column at next (plus) month 
!
! Determine time interpolation factor.  Account for December-January 
! interpolation if cycling dataset.  Code matches that for SST, ozone
!
   caldayloc = get_curr_calday ()

   if (mo_nxt == 1) then                    ! Dec-Jan interpolation
      deltat = cdayp + 365. - cdaym
      if (caldayloc > cdayp) then           ! We're in December
         fact1 = (cdayp + 365. - caldayloc)/deltat
         fact2 = (caldayloc - cdaym)/deltat
      else                                  ! We're in January
         fact1 = (cdayp - caldayloc)/deltat
         fact2 = (caldayloc + 365. - cdaym)/deltat
      end if
   else
      deltat = cdayp - cdaym
      fact1 = (cdayp - caldayloc)/deltat
      fact2 = (caldayloc - cdaym)/deltat
   end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
           fact1 > 1.000001 .or. &
           fact1 < -1.e-6 .or. &
           fact2 > 1.000001 .or. &
           fact2 < -1.e-6) then
      write(6,*)'GET_AEROSOL: Bad fact1 and/or fact2=',fact1,fact2
      call endrun ()
   end if
!
! interpolate (prv and nxt month) bounding datasets onto cam vertical grid.
! compute mass mixing ratios on CAMS's pressure coordinate
!  for both the "minus" and "plus" months
!
!   ncol = get_ncols_p(c)

   call vert_interpolate (M_ps_cam_col(1,c,nm), pint, nm, AEROSOLm, ncol1, ncol, c)
   call vert_interpolate (M_ps_cam_col(1,c,np), pint, np, AEROSOLp, ncol1, ncol, c)
!
! Time interpolate.
!
   do m=1,naer
      do k=1,pver
         do i=ncol1,ncol
            AEROSOLt(i,k,m) = AEROSOLm(i,k,m)*fact1 + AEROSOLp(i,k,m)*fact2
         end do
      end do
   end do

   do i=ncol1,ncol
      Match_ps_chunk(i,c) = M_ps_cam_col(i,c,nm)*fact1 + M_ps_cam_col(i,c,np)*fact2
   end do
!
! get background aerosol (tuning) field
!
   call background (c, ncol1, ncol, pint, AEROSOLt(:, :, idxBG))
!
! exit if mmr is negative (we have previously set
!  cumulative mass to be a decreasing function.)
!
   speciesmin(:) = 0. ! speciesmin(m) = 0 is minimum mmr for each species
 
   do m=1,naer
      do k=1,pver
         do i=ncol1,ncol
            if (AEROSOLt(i, k, m) < speciesmin(m)) then
               write(6,*) 'AEROSOL_INTERPOLATE: negative mass mixing ratio, exiting'
               write(6,*) 'm, column, pver',m, i, k ,AEROSOLt(i, k, m)
               call endrun ()
            end if
         end do
      end do
   end do

   call ramp_sulfate (AEROSOLt)
!
! scale any AEROSOLS as required
!
   call scale_aerosols (AEROSOLt, ncol1, ncol, c, scale)

   return
end subroutine get_aerosol

subroutine vert_interpolate (Match_ps, pint, n, AEROSOL_mmr, ncol1, ncol, c)
!--------------------------------------------------------------------
! Input: match surface pressure, cam interface pressure, 
!        month index, number of columns, chunk index
! 
! Output: Aerosol mass mixing ratio (AEROSOL_mmr)
!
! Method:
!         interpolate column mass (cumulative) from match onto
!           cam's vertical grid (pressure coordinate)
!         convert back to mass mixing ratio
!
!--------------------------------------------------------------------

   use physconst,     only: gravit

   real(r8), intent(out) :: AEROSOL_mmr(pcols,pver,naer)  ! aerosol mmr from MATCH
   real(r8), intent(in) :: Match_ps(pcols)                ! surface pressure at a particular month
   real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from CAM

   integer, intent(in) :: ncol,ncol1,c                          ! chunk index and number of columns
   integer, intent(in) :: n                               ! prv or nxt month index
!
! Local workspace
!
   integer m                           ! index to aerosol species
   integer kupper(pcols)               ! last upper bound for interpolation
   integer i, k, kk, kkstart, kount    ! loop vars for interpolation
   integer isv, ksv, msv               ! loop indices to save

   logical bad                         ! indicates a bad point found

   real(r8) AEROSOL(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper 
                                       ! interface of level in column at particular month
   real(r8) dpl, dpu                   ! lower and upper intepolation factors
   real(r8) v_coord                    ! vertical coordinate
   real(r8) m_to_mmr                   ! mass to mass mixing ratio conversion factor
   real(r8) AER_diff                   ! temp var for difference between aerosol masses

   call t_startf ('vert_interpolate')
!
! Initialize index array 
!
   do i=ncol1,ncol
      kupper(i) = 1
   end do
!
! assign total mass to topmost level
!
   AEROSOL(:,1,:) = AEROSOLc(:,1,c,:,n)
!
! At every pressure level, interpolate onto that pressure level
!
   do k=2,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = paerlev+1
      do i=ncol1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
! for the pressure interpolation should be comparing
! pint(column,lev) with M_hybi(lev)*M_ps_cam_col(month,column,chunk)
!
      do kk=kkstart,paerlev
         do i=ncol1,ncol
            v_coord = pint(i,k)
            if (M_hybi(kk)*Match_ps(i) .lt. v_coord .and. v_coord .le. M_hybi(kk+1)*Match_ps(i)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
! Interpolate in pressure.
!
         if (kount.eq.ncol-ncol1+1) then
            do i=ncol1,ncol
               dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
               dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
               AEROSOL(i,k,:) = &
                    (AEROSOLc(i,kupper(i)  ,c,:,n)*dpl + &
                     AEROSOLc(i,kupper(i)+1,c,:,n)*dpu)/(dpl + dpu)
            enddo !i
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top pressure level for at least some
! of the longitude points.
!
      do i=ncol1,ncol
         if (pint(i,k) .lt. M_hybi(1)*Match_ps(i)) then
            AEROSOL(i,k,:) =  AEROSOLc(i,1,c,:,n)
         else if (pint(i,k) .gt. M_hybi(paerlev+1)*Match_ps(i)) then
            AEROSOL(i,k,:) = 0.0
         else
            dpu = pint(i,k) - M_hybi(kupper(i))*Match_ps(i)
            dpl = M_hybi(kupper(i)+1)*Match_ps(i) - pint(i,k)
            AEROSOL(i,k,:) = &
                 (AEROSOLc(i,kupper(i)  ,c,:,n)*dpl + &
                  AEROSOLc(i,kupper(i)+1,c,:,n)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol-ncol1+1) then
         write(6,*)'vert_interpolate: Bad data: non-monotonicity suspected in dependent variable'
         call endrun ()
      end if
35    continue
   end do

   call t_startf ('vi_checks')
!
! aerosol mass beneath lowest interface (pverp) must be 0
!
   AEROSOL(ncol1:ncol,pverp,:) = 0.
!
! Set mass in layer to zero whenever it is less than 
!   1.e-40 kg/m^2 in the layer
!
   do m = 1, naer
      do k = 1, pver
         do i = ncol1, ncol
            if (AEROSOL(i,k,m) < 1.e-40) AEROSOL(i,k,m) = 0.
         end do
      end do
   end do
!
! Set mass in layer to zero whenever it is less than 
!   10^-15 relative to column total mass
! convert back to mass mixing ratios. 
! exit if mmr is negative
!
   do m = 1, naer
      do k = 1, pver
         do i = ncol1, ncol
            AER_diff = AEROSOL(i,k,m) - AEROSOL(i,k+1,m)
            if( abs(AER_diff) < 1e-15*AEROSOL(i,1,m)) then
               AER_diff = 0.
            end if
            m_to_mmr = gravit / (pint(i,k+1)-pint(i,k))
            AEROSOL_mmr(i,k,m)= AER_diff * m_to_mmr
            if (AEROSOL_mmr(i,k,m) < 0) then
               write(6,*)'vert_interpolate: mmr < 0, m, col, lev, mmr',m, i, k, AEROSOL_mmr(i,k,m)
               write(6,*)'vert_interpolate: aerosol(k),(k+1)',AEROSOL(i,k,m),AEROSOL(i,k+1,m)
               write(6,*)'vert_interpolate: pint(k+1),(k)',pint(i,k+1),pint(i,k)
               write(6,*)'n,c',n,c
               call endrun()
            end if
         end do
      end do
   end do

   call t_stopf ('vi_checks')
   call t_stopf ('vert_interpolate')

   return
end subroutine vert_interpolate

subroutine scale_aerosols(AEROSOLt, ncol1, ncol, lchnk, scale)
!-----------------------------------------------------------------
! scale each species as determined by scale factors
!-----------------------------------------------------------------
  integer, intent(in) :: ncol,ncol1, lchnk ! number of columns and chunk index
  real(r8), intent(in) :: scale(naer_all) ! scale each aerosol by this amount
  real(r8), intent(inout) :: AEROSOLt(pcols, pver, naer_all) ! aerosols
  integer m

  do m = 1, naer_all
     AEROSOLt(ncol1:ncol, :, m) = scale(m)*AEROSOLt(ncol1:ncol, :, m)
  end do

  return
end subroutine scale_aerosols

subroutine background(lchnk, ncol1, ncol, pint, mmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set global mean tropospheric aerosol background (or tuning) field
! 
! Method: 
! Specify aerosol mixing ratio.
! Aerosol mass mixing ratio
! is specified so that the column visible aerosol optical depth is a
! specified global number (tauback). This means that the actual mixing
! ratio depends on pressure thickness of the lowest three atmospheric
! layers near the surface.
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use aer_optics, only: kbg,idxVIS
   use physconst, only: gravit
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <ptrrgrid.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol,ncol1                  ! number of atmospheric columns

   real(r8), intent(in) :: pint(pcols,pverrp)   ! Interface pressure (mks)
!
! Output arguments
!
   real(r8), intent(out) :: mmr(pcols,pverr)    ! "background" aerosol mass mixing ratio
!
!---------------------------Local variables-----------------------------
!
   integer i          ! Longitude index
   integer k          ! Level index
!
   real(r8) mass2mmr  ! Factor to convert mass to mass mixing ratio
   real(r8) mass      ! Mass of "background" aerosol as specified by tauback
!
!-----------------------------------------------------------------------
!
   do i=ncol1,ncol
      mass2mmr =  gravit / (pint(i,pverrp)-pint(i,pverrp-mxaerl))
      do k=1,pverr
!
! Compute aerosol mass mixing ratio for specified levels (1.e3 factor is
! for units conversion of the extinction coefficiant from m2/g to m2/kg)
!
         if ( k >= pverrp-mxaerl ) then
! kaervs is not consistent with the values in aer_optics
! this ?should? be changed.
! rhfac is also implemented differently
            mass = tauback / (1.e3 * kbg(idxVIS))
            mmr(i,k) = mass2mmr*mass
         else
            mmr(i,k) = 0._r8
         endif
!
      enddo
   enddo
!
   return
end subroutine background


subroutine aerosol_diagnostics(state, AEROSOL_mmr)
!-----------------------------------------------------------------
! write out the mass of aerosol in each layer
!-----------------------------------------------------------------

  use physics_types, only: physics_state
  use physconst,     only: rga
  use history,       only: outfld

  type(physics_state), intent(in) :: state     ! state from physics code
  real(r8), intent(in) :: AEROSOL_mmr(pcols, pver, naer_all) ! aerosols

  integer lchnk,ncol  ! chunk index and number of columns
  real(r8) :: MAerosol(pcols, pver, naer_all)  !mass of aerosols in each layer
! column totals of masses
  real(r8) :: TMAllAerosols(pcols), TMCarbon(pcols), TMDust(pcols) 
  real(r8) :: TMAerosol(pcols,naer_all)
  integer k,m         ! level, constituent indices

  TMAllAerosols = 0.0
  TMCarbon = 0.0
  TMDust = 0.0
  TMAerosol = 0.0

  ncol = state%ncol
  lchnk = state%lchnk
!
! output surface pressure from MATCH
!
  call outfld('PS_match',Match_ps_chunk(1,lchnk),pcols,lchnk)
!
! compute column mass of aerosols (and groups)
!
! first, mass in each layer
!
  do m = 1, naer_all
    MAerosol(:ncol,:,m) = AEROSOL_mmr(:ncol,:,m)*state%pdel(:ncol,:)*rga
  enddo
!
! now accumulate mass below upper interface bounding each layer
!
  do k = pver-1,1,-1
    MAerosol(:ncol,k,:) = MAerosol(:ncol,k+1,:) + MAerosol(:ncol,k,:)
  enddo
!
! total mass in column is given by cumulative mass and level = 1
!
  TMAerosol(:ncol,:) = MAerosol(:ncol,1,:)

  do m = idxCARBONfirst, idxCARBONfirst+numCARBON-1
    TMCarbon(:ncol) = TMCarbon(:ncol) + TMAerosol(:ncol,m)
  enddo

  do m = idxDUSTfirst, idxDUSTfirst+numDUST-1
    TMDust(:ncol) = TMDust(:ncol) + TMAerosol(:ncol,m)
  enddo

  TMAllAerosols(:ncol)=TMCarbon(:ncol)            &
                      +TMDust(:ncol)              &
                      +TMAerosol(:ncol,idxSUL )   &
                      +TMAerosol(:ncol,idxSSLT)   
!
! write out column totals of aerosol mmr's and groups
!
  call outfld('MSUL_V  ', MAerosol(:,:,idxSUL) ,pcols, lchnk)
  call outfld('MSSLT_V ', MAerosol(:,:,idxSSLT) ,pcols, lchnk)
  call outfld('MDUST1_V', MAerosol(:,:,idxDUSTfirst) ,pcols, lchnk)
  call outfld('MDUST2_V', MAerosol(:,:,idxDUSTfirst+1) ,pcols, lchnk)
  call outfld('MDUST3_V', MAerosol(:,:,idxDUSTfirst+2) ,pcols, lchnk)
  call outfld('MDUST4_V', MAerosol(:,:,idxDUSTfirst+3) ,pcols, lchnk)
  call outfld('MOCPHO_V', MAerosol(:,:,idxOCPHO) ,pcols, lchnk)
  call outfld('MOCPHI_V', MAerosol(:,:,idxOCPHI) ,pcols, lchnk)
  call outfld('MBCPHO_V', MAerosol(:,:,idxBCPHO) ,pcols, lchnk)
  call outfld('MBCPHI_V', MAerosol(:,:,idxBCPHI) ,pcols, lchnk)
  call outfld('MBG_V   ', MAerosol(:,:,idxBG) ,pcols, lchnk)

  return
end subroutine aerosol_diagnostics

subroutine aerosol_indirect(ncol1, ncol, lchnk,landfrac,pmid,t,qm1,cld,zm,rel)
!--------------------------------------------------------------
! Compute effect of sulfate on effective liquid water radius
!  Method of Martin et. al.
!--------------------------------------------------------------

  use constituents, only: ppcnst, cnst_get_ind
  use history, only: outfld

#include <comctl.h>

  integer, intent(in) :: ncol1,ncol                  ! number of atmospheric columns
  integer, intent(in) :: lchnk                 ! chunk identifier

  real(r8), intent(in) :: landfrac(pcols)      ! land fraction
  real(r8), intent(in) :: pmid(pcols,pver)     ! Model level pressures
  real(r8), intent(in) :: t(pcols,pver)        ! Model level temperatures
  real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! Specific humidity and tracers
  real(r8), intent(in) :: cld(pcols,pver)      ! Fractional cloud cover
  real(r8), intent(in) :: zm(pcols,pver)       ! Height of midpoints (above surface)
  real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
!
! local variables
!
  real(r8) locrhoair(pcols,pver)  ! dry air density            [kg/m^3 ]
  real(r8) lwcwat(pcols,pver)     ! in-cloud liquid water path [kg/m^3 ]
  real(r8) sulfmix(pcols,pver)    ! sulfate mass mixing ratio  [kg/kg  ]
  real(r8) so4mass(pcols,pver)    ! sulfate mass concentration [g/cm^3 ]
  real(r8) Aso4(pcols,pver)       ! sulfate # concentration    [#/cm^3 ]
  real(r8) Ntot(pcols,pver)       ! ccn # concentration        [#/cm^3 ]
  real(r8) relmod(pcols,pver)     ! effective radius           [microns]

  real(r8) wrel(pcols,pver)       ! weighted effective radius    [microns]
  real(r8) wlwc(pcols,pver)       ! weighted liq. water content  [kg/m^3 ]
  real(r8) cldfrq(pcols,pver)     ! frequency of occurance of...
!                                  ! clouds (cld => 0.01)         [fraction]
  real(r8) locPi                  ! my piece of the pi
  real(r8) Rdryair                ! gas constant of dry air   [J/deg/kg]
  real(r8) rhowat                 ! density of water          [kg/m^3  ]
  real(r8) Acoef                  ! m->A conversion factor; assumes
!                                  ! Dbar=0.10, sigma=2.0      [g^-1    ]
  real(r8) rekappa                ! kappa in evaluation of re(lmod)
  real(r8) recoef                 ! temp. coeficient for calc of re(lmod)
  real(r8) reexp                  ! 1.0/3.0
  real(r8) Ntotb                  ! temp var to hold below cloud ccn
! -- Parameters for background CDNC (from `ambient' non-sulfate aerosols)...
  real(r8) Cmarn                  ! Coef for CDNC_marine         [cm^-3]
  real(r8) Cland                  ! Coef for CDNC_land           [cm^-3]
  real(r8) Hmarn                  ! Scale height for CDNC_marine [m]
  real(r8) Hland                  ! Scale height for CDNC_land   [m]
  parameter ( Cmarn = 50.0, Cland = 100.0 )
  parameter ( Hmarn = 1000.0, Hland = 2000.0 )
  real(r8) bgaer                  ! temp var to hold background CDNC
  integer :: ixcldliq               ! index of liquid cloud water

  integer i,k     ! loop indices
!
! Statement functions
!
  logical land    ! is this a column over land?
  land(i) = nint(landfrac(i)).gt.0.5_r8

  if (indirect) then
 
    print *, 'AEROSOLS:  indirect effect is obsolete'
    call endrun

!   ramping is not yet resolved so sulfmix is 0.
    sulfmix(ncol1:ncol,1:pver) = 0._r8

    locPi = 3.141592654
    Rdryair = 287.04
    rhowat = 1000.0
    Acoef = 1.2930E14
    recoef = 3.0/(4.0*locPi*rhowat)
    reexp = 1.0/3.0

    call cnst_get_ind('CLDLIQ', ixcldliq)
    do k=pver,1,-1
      do i = ncol1,ncol
        locrhoair(i,k) = pmid(i,k)/( Rdryair*t(i,k) )
        lwcwat(i,k) = ( qm1(i,k,ixcldliq)/max(0.01_r8,cld(i,k)) )* &
                      locrhoair(i,k)
!          NOTE: 0.001 converts kg/m3 -> g/cm3
        so4mass(i,k) = sulfmix(i,k)*locrhoair(i,k)*0.001
        Aso4(i,k) = so4mass(i,k)*Acoef

        if (Aso4(i,k) <= 280.0) then
           Aso4(i,k) = max(36.0_r8,Aso4(i,k))
           Ntot(i,k) = -1.15E-3*Aso4(i,k)**2 + 0.963*Aso4(i,k)+5.30
           rekappa = 0.80
        else
           Aso4(i,k) = min(1500.0_r8,Aso4(i,k))
           Ntot(i,k) = -2.10E-4*Aso4(i,k)**2 + 0.568*Aso4(i,k)-27.9
           rekappa = 0.67
        end if
        if (land(i)) then ! Account for local background aerosol;
           bgaer = Cland*exp(-(zm(i,k)/Hland))
           Ntot(i,k) = max(bgaer,Ntot(i,k))
        else
           bgaer = Cmarn*exp(-(zm(i,k)/Hmarn))
           Ntot(i,k) = max(bgaer,Ntot(i,k))
        end if

        if (k == pver) then
           Ntotb = Ntot(i,k)
        else
           Ntotb = Ntot(i,k+1)
        end if

        relmod(i,k) = (( (recoef*lwcwat(i,k))/(rekappa*Ntotb))**reexp)*10000.0
        relmod(i,k) = max(4.0_r8,relmod(i,k))
        relmod(i,k) = min(20.0_r8,relmod(i,k))
        if (cld(i,k) >= 0.01) then
           cldfrq(i,k) = 1.0
        else
           cldfrq(i,k) = 0.0
        end if
        wrel(i,k) = relmod(i,k)*cldfrq(i,k)
        wlwc(i,k) = lwcwat(i,k)*cldfrq(i,k)
      end do
    end do
    if(ncol1.ne.ncol) then
     call outfld('MSO4    ',so4mass,pcols,lchnk)
     call outfld('LWC     ',lwcwat ,pcols,lchnk)
     call outfld('CLDFRQ  ',cldfrq ,pcols,lchnk)
     call outfld('WREL    ',wrel   ,pcols,lchnk)
     call outfld('WLWC    ',wlwc   ,pcols,lchnk)
    end if
    write(6,*)'WARNING: indirect calculation has no effects'
  else
    do k = 1, pver
      do i = ncol1, ncol
        relmod(i,k) = rel(i,k)
      end do
    end do
  endif

  call outfld('REL     ',relmod ,pcols,lchnk)

  return
end subroutine aerosol_indirect

subroutine ramp_sulfate(AEROSOL_mmr)
!-------------------------------------------------------------
! increment anthropogenic sulfate mmr 
!
! when this code is resolved, remove line in radctl setting 
!  sulfmix() = 0 in radctl.F90
!-------------------------------------------------------------

  use so4bnd, only: getso4bnd, so4ramp

! doRamp_so4
#include <comctl.h> 

  real(r8), intent(inout) :: AEROSOL_mmr(pcols, pver, naer) ! aerosols
  real(r8) sulfbio(pcols,pver)    ! biogenic sulfate mmr       [kg/kg  ]
  real(r8) sulfant(pcols,pver)    ! anthropogenic sulfate mmr  [kg/kg  ]
  real(r8) sulfscalef             ! sulfate scale factor
  real(r8) sulfmix(pcols,pver)    ! sulfate mass mixing ratio  [kg/kg  ]

  if ( doRamp_so4 ) then
    write(6,*) 'Ramping of sulfate is not yet implemented'
    call endrun()
  endif

!      if ( doRamp_so4 ) then
!         call getso4bnd( lchnk, ncol, sulfbio, sulfant )
!         sulfscalef = so4ramp()
!         do k = 1, pver
!            do i = 1, ncol
!               sulfmix(i,k) = sulfbio(i,k) + sulfscalef*sulfant(i,k)
!            end do
!         end do
!         call outfld('SULFBIO ',sulfbio,pcols,lchnk)
!         call outfld('SULFANT ',sulfant,pcols,lchnk)
!         call outfld('SULFMMR ',sulfmix,pcols,lchnk)
!      else
!         do k = 1, pver
!            do i = 1, ncol
!               sulfmix(i,k) = 0.
!            end do
!         end do
!      endif

  return
end subroutine ramp_sulfate

subroutine setmxaerl()
!------------------------------------
!
! set levels for background aerosol
!
!-----------------------------------

! need hypm from comhyb.h
#include <comhyb.h>

   integer k ! index through vertical levels
!
! mxaerl = max number of levels (from bottom) for background aerosol
! Limit background aerosol height to regions below 900 mb
!
   mxaerl = 0
   do k=pver,1,-1
      if (hypm(k) >= 9.e4) mxaerl = mxaerl + 1
   end do
   mxaerl = max(mxaerl,1)
   if (masterproc) then
      write(6,*)'AEROSOLS:  Background aerosol will be limited to ', &
                'bottom ',mxaerl,' model interfaces. Top interface is ', &
                hypi(pverp-mxaerl),' pascals'
   end if

   return
end subroutine setmxaerl

subroutine add_aer_diag_fields()
!-------------------------------------------------------------
! make sure aerosols are in the master list so they can be written
!-------------------------------------------------------------

  use history, only: phys_decomp, addfld

! description of aerosol species
  character(len=46), parameter :: aerosol_descript(naer_all) =  &
     (/"Sulfate Mass Mixing Ratio                     "&
      ,"Sea Salt Mass Mixing Ratio                    "&
      ,"Dust Mass Mixing Ratio Bin 1                  "&
      ,"Dust Mass Mixing Ratio Bin 2                  "&
      ,"Dust Mass Mixing Ratio Bin 3                  "&
      ,"Dust Mass Mixing Ratio Bin 4                  "&
      ,"Organic Carbon Mass Mixing Ratio (HydroPhobic)"&
      ,"Black Carbon Mass Mixing Ratio   (HydroPhobic)"&
      ,"Organic Carbon Mass Mixing Ratio (HydroPhilic)"&
      ,"Black Carbon Mass Mixing Ratio   (HydroPhilic)"&
      ,"Background Species Mass Mixing Ratio          "/)
  
  call addfld ('PS_match','N/m^2   ',1, 'I','Surface Pressure from aerosol climatology' ,phys_decomp)
  call addfld ('frc_day ','None    ',1, 'I','Portion of time column is lit' ,phys_decomp)

  call addfld ('MSUL_V  ','Kg/m^2  ',pver,'I','Mass of Sulfate in and below layer',phys_decomp)
  call addfld ('MSSLT_V ','Kg/m^2  ',pver,'I','Mass of Sea Salt in and below layer',phys_decomp)
  call addfld ('MDUST1_V','Kg/m^2  ',pver,'I','Mass of Dust bin 1 in and below layer',phys_decomp)
  call addfld ('MDUST2_V','Kg/m^2  ',pver,'I','Mass of Dust bin 2 in and below layer',phys_decomp)
  call addfld ('MDUST3_V','Kg/m^2  ',pver,'I','Mass of Dust bin 3 in and below layer',phys_decomp)
  call addfld ('MDUST4_V','Kg/m^2  ',pver,'I','Mass of Dust bin 4 in and below layer',phys_decomp)
  call addfld ('MOCPHO_V','Kg/m^2  ',pver,'I','Mass of OCPHO in and below layer',phys_decomp)
  call addfld ('MOCPHI_V','Kg/m^2  ',pver,'I','Mass of OCPHI in and below layer',phys_decomp)
  call addfld ('MBCPHO_V','Kg/m^2  ',pver,'I','Mass of BCPHO in and below layer',phys_decomp)
  call addfld ('MBCPHI_V','Kg/m^2  ',pver,'I','Mass of BCPHI in and below layer',phys_decomp)
  call addfld ('MBG_V   ','Kg/m^2  ',pver,'I','Mass of Background Aerosol in and below layer',phys_decomp)

  call addfld ('SULOD_v ','None    ',1, 'I','Sulfate Optical Depth in visible', phys_decomp)
  call addfld ('SSLTOD_v','None    ',1, 'I','Sea Salt Optical Depth in visible', phys_decomp)
  call addfld ('CAROD_v ','None    ',1, 'I','Carbon Optical Depth in visible', phys_decomp)
  call addfld ('DUSTOD_v','None    ',1, 'I','Dust Optical Depth in visible', phys_decomp)
  call addfld ('BGOD_v  ','None    ',1, 'I','Background Aerosol Optical Depth in visible', phys_decomp)
  call addfld ('AEROD_v ','None    ',1, 'I','Total Aerosol Optical Depth in visible', phys_decomp)
  call addfld ('AERSSA_v','None    ',1, 'I','Total Aerosol Single Scattering Albedo in visible', phys_decomp)
  call addfld ('AERASM_v','None    ',1, 'I','Total Aerosol Asymmetry Parameter in visible', phys_decomp)
  call addfld ('AERFWD_v','None    ',1, 'I','Total Aerosol Forward Scattering in visible', phys_decomp)
!
! addfld calls for aerosol forcing-only calculations
!
  if(radforce) then

    call addfld ('FSNT_RF ','W/m^2   ',1, 'I','Total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNTC_RF','W/m^2   ',1, 'I','Clear sky total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNS_RF ','W/m^2   ',1, 'I','Surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNSC_RF','W/m^2   ',1, 'I','Clear sky surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('QRS_RF  ','K/s     ',pver, 'I','Solar heating rate (radforce)' ,phys_decomp)

  endif
  
  return

end subroutine add_aer_diag_fields

subroutine get_rf_scales(scales)

  real(r8), intent(out)::scales(naer_all)  ! scale aerosols by this amount
    
  integer i                                  ! loop index

  scales(idxBG) = bgscl_rf
  scales(idxSUL) = sulscl_rf
  scales(idxSSLT) = ssltscl_rf

  do i = idxCARBONfirst, idxCARBONfirst+numCARBON-1
    scales(i) = carscl_rf
  enddo

  do i = idxDUSTfirst, idxDUSTfirst+numDUST-1
    scales(i) = dustscl_rf
  enddo

end subroutine get_rf_scales

subroutine get_int_scales(scales)

  real(r8), intent(out)::scales(naer_all)  ! scale each aerosol by this amount
   
  integer i                                  ! index through species

  scales(idxBG) = 1._r8
  scales(idxSUL) = sulscl 
  scales(idxSSLT) = ssltscl 

  do i = idxCARBONfirst, idxCARBONfirst+numCARBON-1
    scales(i) = carscl 
  enddo

  do i = idxDUSTfirst, idxDUSTfirst+numDUST-1
    scales(i) = dustscl 
  enddo

  return
end subroutine get_int_scales

subroutine aerint ()
   implicit none

   integer :: ntmp                                ! used in index swapping
   integer :: start(4)                            ! start vector for netcdf calls
   integer :: kount(4)                            ! count vector for netcdf calls
   integer :: i,j,k                               ! spatial indices
   integer :: m                                   ! constituent index

   real(r8) :: caldayloc                          ! calendar day of current timestep
   real(r8) :: aerosol_data(plon,plat,paerlev)    ! aerosol field read in from dataset
   real(r8) :: aerosol_field(plon,paerlev+1,plat) ! aerosol field to be scattered
!
! determine if need to read in next month data
! also determine time interpolation factors
!
   caldayloc = get_curr_calday ()  
!
! If model time is past current forward ozone timeslice, then
! masterproc reads in the next timeslice for time interpolation.  Messy logic is 
! for interpolation between December and January (mo_nxt == 1).  Just like
! oznint, sstint.
!
   if (caldayloc > cdayp .and. .not. (mo_nxt == 1 .and. caldayloc > cdaym)) then
      mo_nxt = mod(mo_nxt,12) + 1
      cdaym = cdayp
      cdayp = Mid(mo_nxt)
!
! Check for valid date info
!
      if (.not. (mo_nxt == 1 .or. caldayloc <= cdayp)) then
         write(6,*)'AERINT: Non-monotonicity suspected in input aerosol data: stopping.'
         call endrun ()
      end if

      ntmp = nm
      nm = np
      np = ntmp

      do m=1,naer
         if (masterproc) then
            start(:) = (/1,1,1,mo_nxt/)
            kount(:) = (/plon,plat,paerlev,1/)
            call wrap_get_vara_realx (aernid, species_id(m), start, kount, aerosol_data(1,1,1))
            do j=1,plat
               do k=1,paerlev
                  aerosol_field(:,k,j) = aerosol_data(:,j,k)
               end do
               aerosol_field(:,paerlev+1,j) = 0.   ! value at bottom
            end do
         end if
         call scatter_field_to_chunk (1, paerlev+1, 1, plon, aerosol_field, AEROSOLc(1,1,begchunk,m,np))
      end do
!
! Retrieve PS from Match
!
      if (masterproc) then
         start(:) = (/1,1,mo_nxt,-1/)
         kount(:) = (/plon,plat,1,-1/)
         call wrap_get_vara_realx (aernid, Mpsid, start, kount, M_ps(1,1))
         write(6,*)'AERINT: Read aerosols data for julian day', Mid(mo_nxt)
      end if
      call scatter_field_to_chunk (1, 1, 1, plon, M_ps(1,1), M_ps_cam_col(1,begchunk,np))
   end if

   return
end subroutine aerint

end module aerosols
