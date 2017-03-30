#include <misc.h>
#include <params.h>
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: analyses
!
! !DESCRIPTION:	Module to handle analyses that will be used to push model
!
! Public interfaces:
!
!	analyses_ini -- Initialization and reading of dataset.
!	analyses_int -- Interpolate dataset analyses to current time.
!       analyses_nudge -- Nudge model prognostics with analyses
!
! $Id: analyses.F90,v 1.1.4.2 2002/06/06 22:38:25 olson Exp $
!
!----------------------------------------------------------------------- 

module analyses
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,    only: plon, plev, plat
  use ppgrid,    only: pcols, pver, begchunk, endchunk
  use phys_grid, only: scatter_field_to_chunk, get_ncols_p
  use ioFileMod, only: getfil
  use pmgrid,   only: masterproc
  use runtime_opts, only:  ncid_analyses, less_surface_nudging, nudge_dse_not_T
  use mpishorthand

  implicit none
!----------------------------------------------------------------------- 
! PUBLIC: Make default data and interfaces private
!----------------------------------------------------------------------- 
!
! ! PUBLIC MEMBER FUNCTIONS:
!
  public analyses_ini   ! Initialization
  public analyses_int   ! Time interpolation of analyses
  public analyses_nudge ! Nudge model prognostics with analyses

  public t_a, u_a, v_a, q_a, ps_a, s_a

!===============================================================================
!EOP
!===============================================================================
!----------------------------------------------------------------------- 
! PRIVATE: Everthing else is private to this module
!----------------------------------------------------------------------- 
  private   ! By default all data is private to this module
  integer, parameter :: totsz=2000
  real(r8), parameter :: daysperyear = 365.0  ! Number of days in a year

  real(r8), allocatable, dimension(:,:,:,:) :: t_bdy  ! Temp   analyses (6 time slices)
  real(r8), allocatable, dimension(:,:,:,:) :: s_bdy  ! Temp   analyses (6 time slices)
  real(r8), allocatable, dimension(:,:,:,:) :: u_bdy  ! u-wind analyses (6 time slices)
  real(r8), allocatable, dimension(:,:,:,:) :: v_bdy  ! v-wind analyses (6 time slices)
  real(r8), allocatable, dimension(:,:,:,:) :: q_bdy  ! moisture analyses (6 time slices)
  real(r8), allocatable, dimension(:,:,:)   :: ps_bdy ! surface pressure analyses (6 time slices)

  real(r8), allocatable, dimension(:,:,:)   :: t_a,s_a    ! Temp:   interpolated between time slices
  real(r8), allocatable, dimension(:,:,:)   :: u_a    ! u-wind; interpolated between time slices
  real(r8), allocatable, dimension(:,:,:)   :: v_a    ! v-wind; interpolated between time slices
  real(r8), allocatable, dimension(:,:,:)   :: q_a    ! moisture; interpolated between time slices
  real(r8), allocatable, dimension(:,:)     :: ps_a   ! surface press; interpolated between time slices

  real(r8) :: cdays(6)            ! Calendar days for 6 time slices read in
  real(r8) :: deltat              ! time (days) between interpolating data

  integer :: nm2,nm1,n0,np1,np2,np3 ! Array indices for 6 time slice data
  integer :: analyses_id          ! netcdf id for analysis variable
  integer :: lonsiz               ! size of longitude dimension on dataset
  integer :: levsiz               ! size of level dimension on dataset
  integer :: latsiz               ! size of latitude dimension on dataset
  integer :: timesiz              ! size of time dimension on dataset
  integer :: n_current            ! current time index of dataset
  integer :: date_analyses(totsz) ! Date on dataset (YYYYMMDD)
  integer :: sec_analyses(totsz)  ! seconds of date on dataset (0-86399)

!===============================================================================
CONTAINS
!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: analyses_ini
!
! !DESCRIPTION:
!
! Do initial read of time-varying analysis boundary dataset(s), reading three
! consecutive time slices on either side of the current model date.
!
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
subroutine analyses_ini
!
! !USES:
!
  use error_messages, only: alloc_err
  use time_manager, only: get_curr_date, get_curr_calday
  use filenames, only: bndtva, m_analyses, max_analyses
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpiint, mpir8
#endif

  implicit none
!
! EOP
!
!---------------------------Common blocks-------------------------------
!#include <comctl.h>
!#include <comlun.h>
!---------------------------Local variables-----------------------------
  integer cnt3(4)               ! array of counts for each dimension
  integer strt3(4)              ! array of starting indices
  integer  :: yr, mon, day      ! components of a date
  integer  :: ncdate_l          ! current date in integer format [yyyymmdd]
  integer  :: ncsec_l           ! current time of day [seconds]
  integer i,mm,m,n,nn           ! indices
  integer istat                 ! error return
  integer bail_out       ! flag to stop time slice search
  integer ntmp           ! temporary
  real(r8) calday_l      ! current julian day including fraction
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
  real(r8) psvar(plon     ,plat,6) ! work space 
  real(r8) tvar (plon,plev,plat,6) ! work space 
  real(r8) svar (plon,plev,plat,6) ! work space 
  real(r8) uvar (plon,plev,plat,6) ! work space 
  real(r8) vvar (plon,plev,plat,6) ! work space 
  real(r8) qvar (plon,plev,plat,6) ! work space 
  real(r8) zvar (plon,plat,plev  ) ! work space 
  character(len=256) :: locfn    ! netcdf local filename to open 
!
!-----------------------------------------------------------------------
!
! Initialize time indices
!


  nm2      = 1
  nm1      = 2
  n0       = 3
  np1      = 4
  np2      = 5
  np3      = 6

  cdays(:) = -1.
  bail_out = 0
!
! Allocate space for data.
!
  allocate( t_a (pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'analyses_ini', 't_a' ,pcols*pver*(endchunk-begchunk+1) )
  allocate( s_a (pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'analyses_ini', 's_a' ,pcols*pver*(endchunk-begchunk+1) )
  allocate( u_a (pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'u_a' ,pcols*pver*(endchunk-begchunk+1) )
  allocate( v_a (pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'v_a' ,pcols*pver*(endchunk-begchunk+1) )
  allocate( q_a (pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'q_a' ,pcols*pver*(endchunk-begchunk+1) )
  allocate( ps_a(pcols,begchunk:endchunk)     , stat=istat )
  call alloc_err( istat, 'analyses_ini', 'ps_a',pcols*     (endchunk-begchunk+1) )

  allocate( t_bdy (pcols,pver,begchunk:endchunk,6), stat=istat )
  call alloc_err( istat, 'analyses_ini', 't_bdy' ,pcols*pver*(endchunk-begchunk+1)*6 )
  allocate( s_bdy (pcols,pver,begchunk:endchunk,6), stat=istat )
  call alloc_err( istat, 'analyses_ini', 's_bdy' ,pcols*pver*(endchunk-begchunk+1)*6 )
  allocate( u_bdy (pcols,pver,begchunk:endchunk,6), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'u_bdy' ,pcols*pver*(endchunk-begchunk+1)*6 )
  allocate( v_bdy (pcols,pver,begchunk:endchunk,6), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'v_bdy' ,pcols*pver*(endchunk-begchunk+1)*6 )
  allocate( q_bdy (pcols,pver,begchunk:endchunk,6), stat=istat )
  call alloc_err( istat, 'analyses_ini', 'q_bdy' ,pcols*pver*(endchunk-begchunk+1)*6 )
  allocate( ps_bdy(pcols,begchunk:endchunk,6)     , stat=istat )
  call alloc_err( istat, 'analyses_ini', 'ps_bdy',pcols*     (endchunk-begchunk+1)*6 )
!
! SPMD: Master does all the work.
!
  if (masterproc) then
!
! Extract current model date info
!
     call get_curr_date(yr, mon, day, ncsec_l)
     calday_l = get_curr_calday()
     ncdate_l = yr*10000 + mon*100 + day
     caldayloc = calday_l + yr*daysperyear
!
! Kluge for restarts:  In case the previous analysis file contains the lower bracketing
!                      time sample, always read it in first on restart.
!
     if(m_analyses > 1) m_analyses = m_analyses - 1

     do m = m_analyses,max_analyses
        mm = m
        call getfil (bndtva(mm), locfn)
        call wrap_open (locfn, 0, ncid_analyses)
        write(6,*)'ANALYSES_INI: NCOPN returns id ',ncid_analyses,' for file ',trim(locfn)
!
! Get and check dimension info and retrieve entire date and sec variables.
!
        call check_file
!
! Get time slice (insert into "np3" time index; will get cycled backwards)
!
        do n=1,timesiz
           nn = n
           call bnddyi(date_analyses(n), sec_analyses(n), cdays(np3))
           write(6,*)'ANALYSES_INI: Read analysis data for date ',date_analyses(n),sec_analyses(n)

           yr         = date_analyses(n)/10000
           cdays(np3) = cdays(np3) + yr*daysperyear
           if(cdays(np2) >= cdays(np3)) then
              write(6,*)'ANALYSES_INI:  List of analysis time samples not in ascending chronological order'
              call endrun
           end if

           strt3(1)   = 1
           strt3(2)   = 1
           strt3(3)   = 1
           strt3(4)   = n
           cnt3(1)    = lonsiz
           cnt3(2)    = latsiz
           cnt3(3)    = levsiz
           cnt3(4)    = 1
           if (nudge_dse_not_T) then 
             call get_analyses('S'     ,plon    ,plev    ,plat    ,strt3   , &
                             cnt3    ,zvar    ,svar (1,1,1,np3)  )
           endif
           call get_analyses('T'     ,plon    ,plev    ,plat    ,strt3   , &
                             cnt3    ,zvar    ,tvar (1,1,1,np3)  )
           call get_analyses('U'     ,plon    ,plev    ,plat    ,strt3   , &
                             cnt3    ,zvar    ,uvar (1,1,1,np3)  )
           call get_analyses('V'     ,plon    ,plev    ,plat    ,strt3   , &
                             cnt3    ,zvar    ,vvar (1,1,1,np3)  )
           call get_analyses('Q'     ,plon    ,plev    ,plat    ,strt3   , &
                             cnt3    ,zvar    ,qvar (1,1,1,np3)  )
           strt3(1)   = 1
           strt3(2)   = 1
           strt3(3)   = n
           strt3(4)   = 1
           cnt3(1)    = lonsiz
           cnt3(2)    = latsiz
           cnt3(3)    = 1
           cnt3(4)    = 1
           call get_analyses('PS'    ,plon    ,1       ,plat    ,strt3   , &
                             cnt3    ,zvar    ,psvar(1,1,np3)    )
!
! If found three time slices beyond current model time, bail out
!
           if (caldayloc < cdays(np3)) then
              bail_out = bail_out + 1
              if(bail_out == 3) goto 10
           end if
!
! Else, cycle time indices and read next time slice
!
           ntmp = nm2
           nm2  = nm1
           nm1  = n0
           n0   = np1
           np1  = np2
           np2  = np3
           np3  = ntmp
        end do
        write(6,*)'ANALYSES_INI: nf_close(',ncid_analyses,') = ',trim(bndtva(mm))
        call wrap_close (ncid_analyses)
     end do
     write(6,*)'ANALYSES_INI: Failed to find analyses bracketing the date ', &
          ncdate_l, ncsec_l
     call endrun
10   continue

     n_current  = nn
     m_analyses = mm

     do i = 1,6
        if(cdays(i) < 0.) then
           write(6,*)'ANALYSES_INI: Failed to find 3 dates before and 3 dates after ', &
                      ncdate_l, ncsec_l
           call endrun
        endif
     end do

     deltat = cdays(np3) - cdays(np2)
     if ( abs(cdays(np2) - cdays(np1) - deltat) > 1.e-6  .or. &
          abs(cdays(np1) - cdays(n0 ) - deltat) > 1.e-6  .or. &
          abs(cdays(n0 ) - cdays(nm1) - deltat) > 1.e-6  .or. &
          abs(cdays(nm1) - cdays(nm2) - deltat) > 1.e-6) then
        write(6,*)'ANALYSES_INI: Time samples are not evenly spaced in ascending order'
        call endrun
     end if
!
! Compute ln(Ps) for time interpolation
!
     psvar(:,:,:) = log( psvar(:,:,:) )
  end if

  call mpibcast( timesiz      , 1    , mpiint, 0, mpicom )
  call mpibcast( date_analyses, totsz, mpiint, 0, mpicom )
  call mpibcast( sec_analyses , totsz, mpiint, 0, mpicom )
  call mpibcast( m_analyses   , 1    , mpiint, 0, mpicom )
  call mpibcast( n_current    , 1    , mpiint, 0, mpicom )
  call mpibcast( nm2          , 1    , mpiint, 0, mpicom )
  call mpibcast( nm1          , 1    , mpiint, 0, mpicom )
  call mpibcast( n0           , 1    , mpiint, 0, mpicom )
  call mpibcast( np1          , 1    , mpiint, 0, mpicom )
  call mpibcast( np2          , 1    , mpiint, 0, mpicom )
  call mpibcast( np3          , 1    , mpiint, 0, mpicom )
  call mpibcast( deltat       , 1    , mpir8 , 0, mpicom )
  call mpibcast( cdays        , 6    , mpir8 , 0, mpicom )

  call scatter_field_to_chunk(1,plev,6,plon,svar ,s_bdy )
  call scatter_field_to_chunk(1,plev,6,plon,tvar ,t_bdy )
  call scatter_field_to_chunk(1,plev,6,plon,uvar ,u_bdy )
  call scatter_field_to_chunk(1,plev,6,plon,vvar ,v_bdy )
  call scatter_field_to_chunk(1,plev,6,plon,qvar ,q_bdy )
  call scatter_field_to_chunk(1,   1,6,plon,psvar,ps_bdy)

  return
end subroutine analyses_ini

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: analyses_int
!
! !DESCRIPTION:
!
! Time interpolate analyses to current time, reading in new data if
! necessary.
!
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine analyses_int(ztodt)
!
! !USES:
!
  use time_manager, only: get_curr_date, get_curr_calday, get_nstep
  use filenames, only: bndtva, m_analyses, max_analyses
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpiint
#endif
  use runtime_opts, only: analyses_time_interp, nudge_dse_not_t

  implicit none
!
! EOP
!
!---------------------------Common blocks-------------------------------
!#include <comctl.h>
!#include <comlun.h>
! ------------------------ Arguments------------------------------------
  real(r8), intent(in) :: ztodt  ! twice time step unless nstep=0
!---------------------------Local variables-----------------------------
  integer cnt3(4)        ! array of counts for each dimension
  integer strt3(4)       ! array of starting indices
  integer  :: yr, mon, day ! components of a date
  integer  :: ncsec_l    ! current time of day [seconds]
  integer i,j,k,lchnk    ! indices
  integer ncol           ! number of columns in current chunk
  integer ntmp           ! temporary
  integer kstep          ! local time step
  real(r8) denom1        ! |
  real(r8) denom2        ! |
  real(r8) denom3        ! |
  real(r8) denom4        ! | - weights for Lag. time interpolation
  real(r8) denom5        ! |
  real(r8) denom6        ! |
  real(r8) tmp1          ! |
  real(r8) tmp2          ! |
  real(r8) tmp3          ! |
  real(r8) tmp4          ! | -- tmp variables for Lag. time interpolation
  real(r8) tmp5          ! |
  real(r8) tmp6          ! |
  real(r8) coef23        ! |
  real(r8) coef45        ! |
  real(r8) coef123       ! |
  real(r8) coef456       ! |
  real(r8) wgt1          ! weight for time-interpolant
  real(r8) wgt2          ! weight for time-interpolant
  real(r8) wgt3          ! weight for time-interpolant
  real(r8) wgt4          ! weight for time-interpolant
  real(r8) wgt5          ! weight for time-interpolant
  real(r8) wgt6          ! weight for time-interpolant
  real(r8) calday_l      ! current julian day including fraction
  real(r8) fact1, fact2  ! time interpolation factors
  real(r8) rdt           ! Hermite cubic weight
  real(r8) rdt6          ! Hermite cubic weight
  real(r8) hb,ha,dhb,dha ! Hermite cubic weights
  real(r8) fxb, fxa      ! Hermite cubic weights
  real(r8) fac           ! factor applied in Hermite limiter
  real(r8) deli          ! linear derivative
  real(r8) caldayloc     ! calendar day (includes yr if no cycling)
  real(r8) psvar(plon     ,plat) ! work space 
  real(r8) tvar (plon,plev,plat) ! work space 
  real(r8) svar (plon,plev,plat) ! work space 
  real(r8) uvar (plon,plev,plat) ! work space 
  real(r8) vvar (plon,plev,plat) ! work space 
  real(r8) qvar (plon,plev,plat) ! work space 
  real(r8) zvar (plon,plat,plev) ! work space 
  character(len=256) :: locfn    ! netcdf local filename to open 
!
! Use year information only if a multiyear dataset
!
  kstep = get_nstep()
  call get_curr_date(yr, mon, day, ncsec_l)
  calday_l = get_curr_calday()
  caldayloc = calday_l + yr*daysperyear
!
! If model time is past current forward timeslice, read in the next
! timeslice for time interpolation.
!
  if (caldayloc >= cdays(np1)) then
     n_current = n_current + 1
     if (n_current > timesiz) then
        if (masterproc) then

           write(6,*)'ANALYSES_INT: nf_close(',ncid_analyses,') = ',trim(bndtva(m_analyses))
           call wrap_close (ncid_analyses)
           m_analyses = m_analyses + 1

           if(m_analyses .gt. max_analyses) then
              write(6,*)'ANALYSES_INT: Attempt to read past end of analysis dataset(s)'
              call endrun
           end if

           call getfil (bndtva(m_analyses), locfn)
           call wrap_open (locfn, 0, ncid_analyses)
           write(6,*)'ANALYSES_INT: NCOPN returns id ',ncid_analyses,' for file ',trim(locfn)

           call check_file
           n_current = 1
        end if
#if (defined SPMD )
        call mpibcast( timesiz      , 1    , mpiint, 0, mpicom )
        call mpibcast( date_analyses, totsz, mpiint, 0, mpicom )
        call mpibcast( sec_analyses , totsz, mpiint, 0, mpicom )
        call mpibcast( m_analyses   , 1    , mpiint, 0, mpicom )
        call mpibcast( n_current    , 1    , mpiint, 0, mpicom )
#endif
     end if
!
! Cycle time indices and read next time slice
!
     ntmp = nm2
     nm2  = nm1
     nm1  = n0
     n0   = np1
     np1  = np2
     np2  = np3
     np3  = ntmp
     call bnddyi(date_analyses(n_current), sec_analyses(n_current), cdays(np3))
     yr         = date_analyses(n_current)/10000
     cdays(np3) = cdays(np3) + yr*daysperyear

     if ( abs(cdays(np3) - cdays(np2) - deltat) > 1.e-6) then
        if (masterproc) then
           write(6,*)'ANALYSES_INT: Input analysis for date',date_analyses(n_current), &
                ' sec ',sec_analyses(n_current), 'not evenly spaced from previous analysis dates'
        endif
        call endrun
     endif
     if (masterproc) then
        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        strt3(4) = n_current
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = levsiz
        cnt3(4)  = 1
        call get_analyses('T'     ,plon    ,plev    ,plat    ,strt3   , &
                          cnt3    ,zvar    ,tvar    )
        if (nudge_dse_not_T) then 
        call get_analyses('S'     ,plon    ,plev    ,plat    ,strt3   , &
                          cnt3    ,zvar    ,svar    )
        end if
       call get_analyses('U'     ,plon    ,plev    ,plat    ,strt3   , &
                          cnt3    ,zvar    ,uvar    )
        call get_analyses('V'     ,plon    ,plev    ,plat    ,strt3   , &
                          cnt3    ,zvar    ,vvar    )
        call get_analyses('Q'     ,plon    ,plev    ,plat    ,strt3   , &
                          cnt3    ,zvar    ,qvar    )
        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = n_current
        strt3(4) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1
        cnt3(4)  = 1
        call get_analyses('PS'    ,plon    ,1       ,plat    ,strt3   , &
                          cnt3    ,zvar    ,psvar   )
!
! Compute ln(Ps) for time interpolation
!
        psvar(:,:) = log( psvar(:,:) )

        write(6,*)'ANALYSES_INT: Read analyses for date (yyyymmdd) ',   &
                   date_analyses(n_current), &
                  ' sec ',sec_analyses(n_current)
     endif
    if (nudge_dse_not_T) then
    call scatter_field_to_chunk(1,plev,1,plon,svar ,s_bdy (1,1,begchunk,np3))
    end if
    call scatter_field_to_chunk(1,plev,1,plon,tvar ,t_bdy (1,1,begchunk,np3))
     call scatter_field_to_chunk(1,plev,1,plon,uvar ,u_bdy (1,1,begchunk,np3))
     call scatter_field_to_chunk(1,plev,1,plon,vvar ,v_bdy (1,1,begchunk,np3))
     call scatter_field_to_chunk(1,plev,1,plon,qvar ,q_bdy (1,1,begchunk,np3))
     call scatter_field_to_chunk(1,   1,1,plon,psvar,ps_bdy(1,  begchunk,np3))
  end if

  fact1  = (cdays(np1) - caldayloc)/deltat
  fact2  = (caldayloc  - cdays(n0))/deltat

  if(analyses_time_interp .eq. 'LINEAR') then
!
! Time interpolation:  Linear
!
     wgt1   = 0.
     wgt2   = 0.
     wgt3   = fact1
     wgt4   = fact2
     wgt5   = 0.
     wgt6   = 0.

  elseif(analyses_time_interp .eq. 'CUBIC') then
!
! Time interpolation:  Lag. Cubic
!
!
! Initialize Lag cubic weights
!
     denom2 = -1._r8/6._r8
     denom3 =  0.5_r8
     denom4 = -0.5_r8
     denom5 =  1._r8/6._r8

     tmp2   =  fact2 + 1.
     tmp5   =  fact2 - 2.
     coef23 = -fact1*tmp5
     coef45 =  fact2*tmp2
     wgt1   =  0.
     wgt2   =  denom2*coef23*fact2
     wgt3   =  denom3*coef23*tmp2
     wgt4   =  denom4*coef45*tmp5
     wgt5   = -denom5*coef45*fact1
     wgt6   =  0.

  elseif(analyses_time_interp .eq. 'QUINTIC') then
!
! Time interpolation:  Lag. quintic
!
! Initialize Lag. quintic weights
!
     denom1 = -1._r8/120._r8
     denom2 =  1._r8/ 24._r8
     denom3 = -1._r8/ 12._r8
     denom4 =  1._r8/ 12._r8
     denom5 = -1._r8/ 24._r8
     denom6 =  1._r8/120._r8

     tmp1   =  fact2 + 2.
     tmp2   =  fact2 + 1.
     tmp5   =  fact2 - 2.
     tmp6   =  fact2 - 3.
     coef123= -fact1*tmp5*tmp6
     coef456=  fact2*tmp1*tmp2
     wgt1   =  denom1*coef123*tmp2*fact2
     wgt2   =  denom2*coef123*tmp1*fact2
     wgt3   =  denom3*coef123*tmp1*tmp2
     wgt4   =  denom4*coef456*tmp5*tmp6
     wgt5   = -denom5*coef456*tmp6*fact1
     wgt6   = -denom6*coef456*tmp5*fact1

  else
     if (masterproc) then
        write(6,*) 'ANALYSES_INT: Namelist variable ANALYSES_TIME_INTERP must be', &
                   ' set to "LINEAR", "CUBIC", or "QUINTIC" only.'
        write(6,*) 'Currently set to:  ',analyses_time_interp
     endif
     call endrun
  end if

  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     do k = 1,pver
        do i = 1,ncol
           if (nudge_dse_not_T) then 
              s_a(i,k,lchnk) = s_bdy(i,k,lchnk,nm2)*wgt1 + s_bdy(i,k,lchnk,nm1)*wgt2 + &
                            s_bdy(i,k,lchnk,n0 )*wgt3 + s_bdy(i,k,lchnk,np1)*wgt4 + &
                            s_bdy(i,k,lchnk,np2)*wgt5 + s_bdy(i,k,lchnk,np3)*wgt6
           end if
           t_a(i,k,lchnk) = t_bdy(i,k,lchnk,nm2)*wgt1 + t_bdy(i,k,lchnk,nm1)*wgt2 + &
                            t_bdy(i,k,lchnk,n0 )*wgt3 + t_bdy(i,k,lchnk,np1)*wgt4 + &
                            t_bdy(i,k,lchnk,np2)*wgt5 + t_bdy(i,k,lchnk,np3)*wgt6
           u_a(i,k,lchnk) = u_bdy(i,k,lchnk,nm2)*wgt1 + u_bdy(i,k,lchnk,nm1)*wgt2 + &
                            u_bdy(i,k,lchnk,n0 )*wgt3 + u_bdy(i,k,lchnk,np1)*wgt4 + &
                            u_bdy(i,k,lchnk,np2)*wgt5 + u_bdy(i,k,lchnk,np3)*wgt6
           v_a(i,k,lchnk) = v_bdy(i,k,lchnk,nm2)*wgt1 + v_bdy(i,k,lchnk,nm1)*wgt2 + &
                            v_bdy(i,k,lchnk,n0 )*wgt3 + v_bdy(i,k,lchnk,np1)*wgt4 + &
                            v_bdy(i,k,lchnk,np2)*wgt5 + v_bdy(i,k,lchnk,np3)*wgt6
!!         q_a(i,k,lchnk) = q_bdy(i,k,lchnk,nm2)*wgt1 + q_bdy(i,k,lchnk,nm1)*wgt2 + &
!!                          q_bdy(i,k,lchnk,n0 )*wgt3 + q_bdy(i,k,lchnk,np1)*wgt4 + &
!!                          q_bdy(i,k,lchnk,np2)*wgt5 + q_bdy(i,k,lchnk,np3)*wgt6
        end do
     end do
     do i = 1,ncol
        ps_a(i,lchnk) = ps_bdy(i,lchnk,nm2)*wgt1 + ps_bdy(i,lchnk,nm1)*wgt2 + &
                        ps_bdy(i,lchnk,n0 )*wgt3 + ps_bdy(i,lchnk,np1)*wgt4 + &
                        ps_bdy(i,lchnk,np2)*wgt5 + ps_bdy(i,lchnk,np3)*wgt6
     end do
  end do
!
! Alternative cubic interpolation for "q":  Hermite cubic
!
  if(analyses_time_interp .eq. 'CUBIC') then
!
! Compute weights for q and its time derivatives
!
     fac  = 3.*(1. - 10.*epsilon(fac))
     rdt  = 1./deltat
     rdt6 = 1./(6.*deltat)
  
     hb  = ( 3.0 - 2.0*fact1 )*fact1**2
     ha  = ( 3.0 - 2.0*fact2 )*fact2**2
     dhb = -deltat*( fact1 - 1. )*fact1**2
     dha =  deltat*( fact2 - 1. )*fact2**2
  
     do lchnk=begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        do k = 1,pver
           do i = 1,ncol
!
! "Before", "after" time derivatives of q
!
              fxb  = (   - 2.*q_bdy(i,k,lchnk,nm1) &
                         - 3.*q_bdy(i,k,lchnk,n0 ) &
                         + 6.*q_bdy(i,k,lchnk,np1) &
                         -    q_bdy(i,k,lchnk,np2) )*rdt6
              fxa  = (        q_bdy(i,k,lchnk,nm1) &
                         - 6.*q_bdy(i,k,lchnk,n0 ) &
                         + 3.*q_bdy(i,k,lchnk,np1) &
                         + 2.*q_bdy(i,k,lchnk,np2) )*rdt6
!
! Apply limiter to derivatives
!
              deli = (        q_bdy(i,k,lchnk,np1) - &
                              q_bdy(i,k,lchnk,n0 ) )*rdt
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxb   .le. 0.0  ) fxb = 0.
              if( deli*fxa   .le. 0.0  ) fxa = 0.
              if( abs( fxb ) .gt. tmp2 ) fxb = tmp1
              if( abs( fxa ) .gt. tmp2 ) fxa = tmp1
!
! Compute result
!
              q_a(i,k,lchnk) = q_bdy(i,k,lchnk,n0 )*hb + q_bdy(i,k,lchnk,np1)*ha + &
                               fxb*dhb + fxa*dha
           end do
        end do
     end do
  endif

  return
end subroutine analyses_int

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: analyses_nudge
!
! !DESCRIPTION: 
! Nudge model prognostics with analyeses
!
! Author: J. Olson
!
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine analyses_nudge(ztodt   ,tau     ,ncol    ,nver    ,analysis, &
                          field   ,tend    ,l_tend   )
  implicit none
!
! !INPUT PARAMETERS:
!
  real(r8), intent(in)    :: ztodt                ! 2 delta t (model time increment)
  real(r8), intent(in)    :: tau                  ! time scale for nudging model  with analyses (seconds)
  integer , intent(in)    :: ncol                 ! number of atmospheric columns
  integer , intent(in)    :: nver                 ! levels in column
  real(r8), intent(in)    :: analysis(pcols,nver) ! analysis data
  real(r8), intent(in)    :: field(pcols,nver)    ! model prognostic
!
! OUTPUT PARAMETERS:
!
  real(r8), intent(out)   :: tend (pcols,nver)    ! model prognostic
  logical , intent(out)   :: l_tend               ! logical flag indicating tendency is being
!                                                 ! computed and passed to calling routine.
!                                                 ! (for use in "physics_update")

!
!-----------------------------------------------------------------------
!
! EOP
!
!---------------------------Local variables-----------------------------
  integer  i,k               ! indices
  real(r8) denom, x, y
!-----------------------------------------------------------------------
!
  l_tend = .true.

  do k = 1,nver
     denom = tau + ztodt
     if (less_surface_nudging) then
       x = k-(nver-3) ! Add a ramping function to nudging timescale as approach surface.
       y = max(0.,x) ! y = 0 above z(end-3:end), i.e. above bottom four levels. y increases linearly beneath.
       denom = denom + (exp(5.*y) - 1.) ! The bracketed term is zero above the 4th lowest layer but  increases exponentially to 800 days (6.5e7 s) as we approach the bottom layer. . 
     endif 
     do i = 1,ncol
        tend (i,k) = (analysis(i,k) - field(i,k))/denom
     end do
  end do

  return
end subroutine analyses_nudge


!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: analyses_nudge
!
! !DESCRIPTION: 
! Get analysis data from netCDF file
!
! Author: J. Olson
!
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine get_analyses(name    ,nlon    ,nlev    ,nlat    ,strt3   , &
                        cnt3    ,zvar    ,xvar    )

  implicit none
!
! !INPUT PARAMETERS:
!
  character*(*), intent(in)  :: name     ! variable name
  integer      , intent(in)  :: nlon     ! number of lons
  integer      , intent(in)  :: nlev     ! number of levs
  integer      , intent(in)  :: nlat     ! number of lats
  integer      , intent(in)  :: strt3(4) ! array of starting indices
  integer      , intent(in)  :: cnt3 (4) ! array of counts for each dimension
!
! OUTPUT PARAMETERS:
!
  real(r8)     , intent(out) :: zvar (nlon,nlat,nlev) ! netCDF variable
  real(r8)     , intent(out) :: xvar (nlon,nlev,nlat) ! netCDF variable with
!                                                     ! dimensions transposed
!
!-----------------------------------------------------------------------
!
! EOP
!
!
!---------------------------Common blocks-------------------------------
!#include <comlun.h>
!---------------------------Local variables-----------------------------
!
  integer i,j,k             ! indices
!
!-----------------------------------------------------------------------
!
  if(.not. masterproc) then
     write(6,*)'GET_ANALYSES: This routine must only be called from "masterproc"'
     call endrun
  endif
  call wrap_inq_varid( ncid_analyses, name, analyses_id )
  call wrap_get_vara_realx (ncid_analyses,analyses_id,strt3,cnt3,zvar)
  do j = 1,nlat
     do k = 1,nlev
        do i = 1,nlon
           xvar(i,k,j) = zvar(i,j,k)
        end do
     end do
  end do

  return
end subroutine get_analyses

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: check_file
!
! !DESCRIPTION:
!
! Do initial read of time-varying analysis boundary dataset
!
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
subroutine check_file
!
! !USES:
!

  use rgrid, only: nlon

#include <comctl.h>

!
! EOP
!
!---------------------------Common blocks-------------------------------
!#include <comctl.h>
!#include <comlun.h>
!---------------------------Local variables-----------------------------
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
  integer levdimid              ! netcdf id for level variable
  integer lonid                 ! netcdf id for longitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer nlonid                ! netcdf id for nlon variable (rgrid)
  integer nlon_analyses(plat)   ! number of lons per lat on bdy dataset
  integer j                     ! index
!-----------------------------------------------------------------------
!
! SPMD: Master does all the work.
!
  if(.not. masterproc) then
     write(6,*)'CHECK_FILE: This routine must only be called from "masterproc"'
     call endrun
  endif
!
! Get and check dimension info
!
  call wrap_inq_dimid( ncid_analyses, 'time', timeid  )
  call wrap_inq_dimid( ncid_analyses, 'lon' , londimid)
  call wrap_inq_dimid( ncid_analyses, 'lat' , latdimid)
  call wrap_inq_dimid( ncid_analyses, 'lev' , levdimid)

  call wrap_inq_dimlen( ncid_analyses, londimid, lonsiz   )
  if (lonsiz /= plon) then
     write(6,*)'CHECK_FILE: lonsiz=',lonsiz,' must = plon=',plon
     call endrun
  end if
  call wrap_inq_dimlen( ncid_analyses, latdimid, latsiz   )
  if (latsiz /= plat) then
     write(6,*)'CHECK_FILE: latsiz=',latsiz,' must = plat=',plat
     call endrun
  end if
  call wrap_inq_dimlen( ncid_analyses, levdimid, levsiz   )
  if (levsiz /= plev) then
     write(6,*)'CHECK_FILE: levsiz=',levsiz,' must = plev=',plev
     call endrun
  end if
  call wrap_inq_dimlen( ncid_analyses, timeid, timesiz   )
!
! Check to make sure space allocated for time variables is sufficient
!
  if (timesiz > totsz) then
     write(6,*)'CHECK_FILE:  Allocated space for analysis data is insufficient.'
     write(6,*)'Please increase parameter totsz to',timesiz,' and recompile.'
     call endrun
  end if
!
! Check to ensure reduced or not grid of dataset matches that of model
!
  print*,'fullgrid=',fullgrid
  if (fullgrid) then
     call wrap_inq_varid( ncid_analyses, 'lon', lonid   )
  else
     call wrap_inq_varid (ncid_analyses, 'nlon', nlonid)
     call wrap_get_var_int (ncid_analyses, nlonid, nlon_analyses)
     do j=1,plat
        if (nlon_analyses(j) /= nlon(j)) then
           write(6,*)'CHECK_FILE: model grid does not match dataset grid'
           call endrun
        end if
     end do
  end if
  
  call wrap_inq_varid( ncid_analyses, 'date'   , dateid      )
  call wrap_inq_varid( ncid_analyses, 'datesec', secid       )
  call wrap_inq_varid( ncid_analyses, 'lat'    , latid       )
!
! Retrieve entire date and sec variables.
!
  call wrap_get_var_int (ncid_analyses,dateid,date_analyses)
  call wrap_get_var_int (ncid_analyses,secid,sec_analyses)

  return
end subroutine check_file

end module analyses
