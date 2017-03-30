#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! !MODULE: comsrf
!
! !DESCRIPTION:	Module to handle surface fluxes for the subcomponents of cam/csm
!
! Public interfaces:
!
!	update       
!	init           
!	zero           
!
!-----------------------------------------------------------------------
module comsrf
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents, only: pcnst, pnats
  use ppgrid, only: pcols, begchunk, endchunk, pvermx
  use phys_grid, only: get_ncols_p
  use infnan, only: inf, uninit_r8
  use pmgrid, only: plon, plat

  implicit none

!----------------------------------------------------------------------- 
! PUBLIC: Make default data and interfaces private
!----------------------------------------------------------------------- 
!
! ! PUBLIC MEMBER FUNCTIONS:
!
  public initialize_comsrf          ! Set the surface temperature and sea-ice fraction

  public   ! By default all data is public to this module

  public surface_state
  public srfflx_state
  public srfflx_parm

! Public interfaces

  public srfflx_parm_reset
  public srfflx_state_reset
  public update_ocnice
  public update_srf_fluxes

#ifdef FLUXDAMP
  public :: surface_state
  public :: srfflx_state
  public :: srfflx_parm
#endif

  integer, parameter :: plevmx = 4       ! number of subsurface levels

  character*8 tsnam(plevmx)              ! names of sub-surface temperature fields

  real(r8), allocatable:: landm(:,:)     ! land/ocean/sea ice flag
  real(r8), allocatable:: sgh(:,:)       ! land/ocean/sea ice flag
  real(r8), allocatable:: sicthk(:,:)    ! cam sea-ice thickness (m)
  real(r8), allocatable:: snowhice(:,:)  ! snow depth (liquid water) ovr ice
  real(r8), allocatable:: snowhland(:,:) !snow depth (liquid water) ovr lnd
  real(r8), allocatable:: fv(:,:)        ! needed for dry dep velocities (over land)
  real(r8), allocatable:: ram1(:,:)      ! needed for dry dep velocities (over land)
  real(r8), allocatable:: fsns(:,:)      ! surface absorbed solar flux
  real(r8), allocatable:: fsnt(:,:)      ! Net column abs solar flux at model top
  real(r8), allocatable:: flns(:,:)      ! Srf longwave cooling (up-down) flux
  real(r8), allocatable:: flnt(:,:)      ! Net outgoing lw flux at model top
  real(r8), allocatable:: srfrpdel(:,:)  ! 1./(pint(k+1)-pint(k))
  real(r8), allocatable:: psm1(:,:)      ! surface pressure
  real(r8), allocatable:: absorb(:,:)    ! cam surf absorbed solar flux (W/m2)
  real(r8), allocatable:: prcsnw(:,:)    ! cam tot snow precip
  real(r8), allocatable:: landfrac(:,:)  ! fraction of sfc area covered by land
  real(r8), allocatable:: landfrac_field(:,:)  ! fraction of sfc area covered by land (needed globally)
  real(r8), allocatable:: ocnfrac(:,:)   ! frac of sfc area covered by open ocean
  real(r8), allocatable:: icefrac(:,:)   ! fraction of sfc area covered by seaice
  real(r8), allocatable:: previcefrac(:,:) ! fraction of sfc area covered by seaice
  real(r8), allocatable:: trefmxav(:,:)  ! diagnostic: tref max over the day
  real(r8), allocatable:: trefmnav(:,:)  ! diagnostic: tref min over the day

  real(r8), allocatable:: asdirice(:,:)
  real(r8), allocatable:: aldirice(:,:)
  real(r8), allocatable:: asdifice(:,:)
  real(r8), allocatable:: aldifice(:,:)

  real(r8), allocatable:: asdirlnd(:,:)
  real(r8), allocatable:: aldirlnd(:,:)
  real(r8), allocatable:: asdiflnd(:,:)
  real(r8), allocatable:: aldiflnd(:,:)

  real(r8), allocatable:: asdirocn(:,:)
  real(r8), allocatable:: aldirocn(:,:)
  real(r8), allocatable:: asdifocn(:,:)
  real(r8), allocatable:: aldifocn(:,:)

  real(r8), allocatable:: lwuplnd(:,:)
  real(r8), allocatable:: lwupocn(:,:)
  real(r8), allocatable:: lwupice(:,:)

  real(r8), allocatable:: tsice(:,:)
  real(r8), allocatable:: tsice_rad(:,:) ! Equivalent LW_up temperature

!JR Renamed sst to tsocn to avoid potential conflict with other CAM variable of the same name

  real(r8), allocatable:: tsocn(:,:)
  real(r8), allocatable:: Focn(:,:)
  real(r8), allocatable:: frzmlt(:,:)
  real(r8), allocatable:: aice(:,:)      ! CSIM ice fraction


!---------------------------------------------------------------------------
  type surface_state
     real(r8) :: tbot(pcols)         ! bot level temperature
     real(r8) :: zbot(pcols)         ! bot level height above surface
     real(r8) :: ubot(pcols)         ! bot level u wind
     real(r8) :: vbot(pcols)         ! bot level v wind
     real(r8) :: qbot(pcols)         ! bot level specific humidity
     real(r8) :: pbot(pcols)         ! bot level pressure
     real(r8) :: flwds(pcols)        ! 
     real(r8) :: precsc(pcols)       !
     real(r8) :: precsl(pcols)       !
     real(r8) :: precc(pcols)        ! 
     real(r8) :: precl(pcols)        ! 
     real(r8) :: soll(pcols)         ! 
     real(r8) :: sols(pcols)         ! 
     real(r8) :: solld(pcols)        !
     real(r8) :: solsd(pcols)        !
     real(r8) :: srfrad(pcols)       !
     real(r8) :: thbot(pcols)        ! 
     real(r8) :: tssub(pcols,plevmx) ! cam surface/subsurface temperatures
  end type surface_state

  type srfflx_state
     integer :: lchnk     ! chunk index
     integer :: ncol      ! number of active columns

     real(r8) :: asdir(pcols)            ! albedo: shortwave, direct
     real(r8) :: asdif(pcols)            ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols)            ! albedo: longwave, direct
     real(r8) :: aldif(pcols)            ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols)             ! longwave up radiative flux
     real(r8) :: lhf(pcols)              ! latent heat flux
     real(r8) :: shf(pcols)              ! sensible heat flux
     real(r8) :: wsx(pcols)              ! surface u-stress (N)
     real(r8) :: wsy(pcols)              ! surface v-stress (N)
     real(r8) :: tref(pcols)             ! ref height surface air temp
     real(r8) :: qref(pcols)             ! ref height specific humidity
                                         ! (diagnostic field, used only if coupled)
     real(r8) :: ts(pcols)               ! sfc temp (merged w/ocean if coupled)
     real(r8) :: sst(pcols)              ! sea sfc temp
     real(r8) :: cflx(pcols,pcnst+pnats) ! constituent flux (evap)
  end type srfflx_state

!---------------------------------------------------------------------------
! This is for surface fluxes returned from individual parameterizations
!---------------------------------------------------------------------------

  type srfflx_parm
     character(len=24) :: name    ! name of parameterization which produced tendencies

     real(r8) :: asdir(pcols)     ! albedo: shortwave, direct
     real(r8) :: asdif(pcols)     ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols)     ! albedo: longwave, direct
     real(r8) :: aldif(pcols)     ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols)      ! longwave up radiative flux
     real(r8) :: lhf(pcols)       ! latent heat flux
     real(r8) :: shf(pcols)       ! sensible heat flux
     real(r8) :: wsx(pcols)       ! surface u-stress (N)
     real(r8) :: wsy(pcols)       ! surface v-stress (N)
     real(r8) :: tref(pcols)      ! ref height surface air temp
     real(r8) :: ts(pcols)        ! sfc temp (merged w/ocean if coupled)
     real(r8) :: cflx(pcols,pcnst+pnats)     ! pcnst+pnats constituent flux (evap)
  end type srfflx_parm

  type (surface_state), allocatable :: surface_state2d(:)
  type (srfflx_state), allocatable :: srfflx_state2d(:)
  type (srfflx_parm), allocatable :: srfflx_parm2d(:)
  type (srfflx_parm), allocatable :: srfflx_parm2d_ocn(:)


! Private module data

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
! !IROUTINE: initialize_comsrf
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying ice boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine initialize_comsrf
!-----------------------------------------------------------------------
!
! Purpose:
! Initialize surface data
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
    implicit none

    integer k,c      ! level, constituent indices

    allocate (landm   (pcols,begchunk:endchunk))
    allocate (sgh     (pcols,begchunk:endchunk))
    allocate (sicthk  (pcols,begchunk:endchunk))
    allocate (snowhice(pcols,begchunk:endchunk))
    allocate (snowhland(pcols,begchunk:endchunk))
    allocate (fv(pcols,begchunk:endchunk))
    allocate (ram1(pcols,begchunk:endchunk))
    allocate (fsns    (pcols,begchunk:endchunk))         
    allocate (fsnt    (pcols,begchunk:endchunk))         
    allocate (flns    (pcols,begchunk:endchunk))         
    allocate (flnt    (pcols,begchunk:endchunk))         
    allocate (lwupocn    (pcols,begchunk:endchunk))         
    allocate (lwuplnd    (pcols,begchunk:endchunk))         
    allocate (lwupice    (pcols,begchunk:endchunk))         
    allocate (srfrpdel(pcols,begchunk:endchunk))
    allocate (psm1    (pcols,begchunk:endchunk))
    allocate (absorb  (pcols,begchunk:endchunk))
    allocate (prcsnw  (pcols,begchunk:endchunk))
    allocate (landfrac(pcols,begchunk:endchunk))
    allocate (ocnfrac (pcols,begchunk:endchunk))
    allocate (icefrac (pcols,begchunk:endchunk))
    allocate (previcefrac (pcols,begchunk:endchunk))
    allocate (trefmxav(pcols,begchunk:endchunk))
    allocate (trefmnav(pcols,begchunk:endchunk))

    allocate (tsice_rad(pcols,begchunk:endchunk))
    allocate (asdirice(pcols,begchunk:endchunk))
    allocate (aldirice(pcols,begchunk:endchunk))
    allocate (asdifice(pcols,begchunk:endchunk))
    allocate (aldifice(pcols,begchunk:endchunk))

    allocate (asdirlnd(pcols,begchunk:endchunk))
    allocate (aldirlnd(pcols,begchunk:endchunk))
    allocate (asdiflnd(pcols,begchunk:endchunk))
    allocate (aldiflnd(pcols,begchunk:endchunk))

    allocate (asdirocn(pcols,begchunk:endchunk))
    allocate (aldirocn(pcols,begchunk:endchunk))
    allocate (asdifocn(pcols,begchunk:endchunk))
    allocate (aldifocn(pcols,begchunk:endchunk))

    allocate (tsice(pcols,begchunk:endchunk))
    allocate (tsocn(pcols,begchunk:endchunk))
    allocate (Focn(pcols,begchunk:endchunk))
    allocate (frzmlt(pcols,begchunk:endchunk))
    allocate (aice(pcols,begchunk:endchunk))

    allocate (srfflx_state2d(begchunk:endchunk))
    allocate (srfflx_parm2d(begchunk:endchunk))
    allocate (surface_state2d(begchunk:endchunk))
    allocate (srfflx_parm2d_ocn(begchunk:endchunk))

!
! Initialize to NaN or Inf
!
    landm    (:,:) = inf
    sgh      (:,:) = inf
    sicthk   (:,:) = inf
!
! snowhice and snowhland are only defined for a single surface type
! elements of the array outside valid surface points must be set to
! zero if these fields are to be written to netcdf history files.
!
    snowhice (:,:) = 0.
    snowhland(:,:) = 0.
    fsns     (:,:) = inf
    fsnt     (:,:) = inf
    flns     (:,:) = inf
    flnt     (:,:) = inf
    lwupocn  (:,:) = 0.
    lwuplnd  (:,:) = 0.
    lwupice  (:,:) = 0.
    srfrpdel (:,:) = inf
    psm1     (:,:) = inf
    absorb   (:,:) = inf
    prcsnw   (:,:) = inf
    landfrac (:,:) = inf
    ocnfrac  (:,:) = inf
    icefrac  (:,:) = inf
    previcefrac  (:,:) = uninit_r8
    trefmxav (:,:) = inf
    trefmnav (:,:) = inf

    asdirice (:,:) = 0.
    aldirice (:,:) = 0.
    asdifice (:,:) = 0.
    aldifice (:,:) = 0.

    asdirlnd (:,:) = 0.
    aldirlnd (:,:) = 0.
    asdiflnd (:,:) = 0.
    aldiflnd (:,:) = 0.

    asdirocn (:,:) = 0.
    aldirocn (:,:) = 0.
    asdifocn (:,:) = 0.
    aldifocn (:,:) = 0.

    tsice     (:,:) = inf
    tsice_rad (:,:) = inf

    tsocn (:,:) = inf
    Focn (:,:) = inf
    frzmlt (:,:) = inf
    aice (:,:) = inf
!
! Sub-surface temperatures
!
    if (plevmx > 9) then
       write(6,*)'INITIALIZE_COMSRF: Cannot handle more than 9 subsurface levels'
       call endrun ()
    endif
    do k=1,plevmx
       write(unit=tsnam(k),fmt='(''TS'',i1,''     '')') k
    end do

    do c = begchunk,endchunk
       srfflx_state2d(c)%asdir  (:) = 0.
       srfflx_state2d(c)%asdif  (:) = 0.
       srfflx_state2d(c)%aldir  (:) = 0.
       srfflx_state2d(c)%aldif  (:) = 0.
       srfflx_state2d(c)%lwup   (:) = 0.
       srfflx_state2d(c)%lhf    (:) = 0.
       srfflx_state2d(c)%shf    (:) = 0.
       srfflx_state2d(c)%cflx   (:,:) = 0.
       srfflx_state2d(c)%wsx    (:) = 0.
       srfflx_state2d(c)%wsy    (:) = 0.
       srfflx_state2d(c)%tref   (:) = 0.
       srfflx_state2d(c)%ts     (:) = 0.
       srfflx_state2d(c)%sst    (:) = 0.

       surface_state2d(c)%tbot  (:) = 0.
       surface_state2d(c)%zbot  (:) = 0.
       surface_state2d(c)%ubot  (:) = 0.
       surface_state2d(c)%vbot  (:) = 0.
       surface_state2d(c)%qbot  (:) = 0.
       surface_state2d(c)%pbot  (:) = 0.
       surface_state2d(c)%thbot (:) = 0.
       surface_state2d(c)%tssub (:,:) = 0.

       call srfflx_parm_reset (srfflx_parm2d(c))
       call srfflx_parm_reset (srfflx_parm2d_ocn(c))
    end do

    allocate (landfrac_field(plon,plat))
    landfrac_field(:,:) = inf

    return
  end subroutine initialize_comsrf

!===========================================================================
  subroutine srfflx_parm_reset(parm)
!---------------------------------------------------------------------------
!
! Purpose:
! Zero fluxes that are update by land,ocn,ice sub processes
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------
    implicit none
    type(srfflx_parm), intent(inout) :: parm ! Parameterization tendencies

    parm%asdir(:) = 0.
    parm%asdif(:) = 0.
    parm%aldir(:) = 0.
    parm%aldif(:) = 0.
    parm%lwup (:) = 0.
    parm%lhf  (:) = 0.
    parm%shf  (:) = 0.
    parm%cflx (:,:) = 0.
    parm%wsx  (:) = 0.
    parm%wsy  (:) = 0.
    parm%tref (:) = 0.
    parm%ts   (:) = 0.

    return
  end subroutine srfflx_parm_reset

  subroutine srfflx_state_reset (state)
!-----------------------------------------------------------------------
!
! Purpose:
! Zero fluxes that are update by land,ocn,ice sub processes
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------
    implicit none
    type(srfflx_state), intent(inout)  :: state   ! srfflx state variables

    state%asdir(:) = 0.
    state%asdif(:) = 0.
    state%aldir(:) = 0.
    state%aldif(:) = 0.
    state%lwup(:) = 0.
    state%lhf(:) = 0.
    state%shf(:) = 0.
    state%cflx(:,:) = 0.
    state%wsx(:) = 0.
    state%wsy(:) = 0.
    state%tref(:) = 0.
    state%ts(:) = 0.
    state%sst(:) = 0.

    return
  end subroutine srfflx_state_reset

  subroutine update_srf_fluxes (state, parm, frac, ncol)
!-----------------------------------------------------------------------
!
! Purpose:
! update surface fluxes
!
! Method:
!
! Author: John Truesdale
!
!------------------------------Arguments--------------------------------

    use physconst, only: stebol
    use phys_grid, only: get_ncols_p

    implicit none

    integer, intent(in)               :: ncol        ! number of columns
    type(srfflx_parm), intent(inout)  :: parm
    type(srfflx_state), intent(inout) :: state
    real(r8) :: frac(pcols)
!
! Local workspace
!
    integer :: i,c,m     ! longitude, chunk, constituent indices

    do i=1,ncol
       if (frac(i) > 0.) then
          state%asdir(i) = state%asdir(i) + parm%asdir(i) * frac(i)
          state%asdif(i) = state%asdif(i) + parm%asdif(i) * frac(i)
          state%aldir(i) = state%aldir(i) + parm%aldir(i) * frac(i)
          state%aldif(i) = state%aldif(i) + parm%aldif(i) * frac(i)
          state%lwup(i)  = state%lwup(i)  + parm%lwup(i)  * frac(i)
          state%lhf(i)   = state%lhf(i)   + parm%lhf(i)   * frac(i)
          state%shf(i)   = state%shf(i)   + parm%shf(i)   * frac(i)
          state%wsx(i)   = state%wsx(i)   + parm%wsx(i)   * frac(i)
          state%wsy(i)   = state%wsy(i)   + parm%wsy(i)   * frac(i)
          state%tref(i)  = state%tref(i)  + parm%tref(i)  * frac(i) 
!
! if we are calculating ts for a non-fractional grid box (ie all land or
! all ocean or all ice then use the ts given by the parameterization) 
! otherwise calculate ts based on the grid averaged lwup 
!
!jt pull this code out after bit for bit testing
!jt
          if (frac(i) == 1.) then  
             state%ts(i) = state%ts(i) + parm%ts(i) * frac(i) 
          else
             state%ts(i) = sqrt(sqrt(state%lwup(i)/stebol))
          end if

          do m=1,pcnst+pnats
             state%cflx(i,m) = state%cflx(i,m) + parm%cflx(i,m) * frac(i)
          end do
       end if
    end do
!
! zero srfflx parameterization
!
    call srfflx_parm_reset(parm)

    return
  end subroutine update_srf_fluxes

  subroutine update_ocnice (c)
!-----------------------------------------------------------------------
!
! Purpose: update surface fractions.  This routine MUST NOT be called when
!          running through the coupler.  
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------
     use phys_grid, only: get_ncols_p
     use time_manager, only: is_first_step

     implicit none
!
! Arguments
!
     integer, intent(in) :: c  ! chunk index
!
! Local workspace
!
     integer :: i, ncol

#ifdef COUP_CSM
     write(6,*) 'UPDATE_OCNICE: This routine should NEVER be called in coupled mode.'
     write(6,*) 'Fractions need to come from coupler.'
     call endrun ()
#endif

     call t_startf ('update_ocnice')
!
! Convert from non-land ice fraction (aice) to gridbox ice fraction (icefrac)
! Set previcefrac so new sea ice can be flagged when it forms.
! The array previcefrac is irrelevant in SOM mode.
!
     ncol = get_ncols_p (c)
     do i=1,ncol
        previcefrac(i,c) = icefrac(i,c)
        icefrac(i,c) = aice(i,c)*(1. - landfrac(i,c))
     end do

     do i=1,ncol
        ocnfrac(i,c) = 1.0 - landfrac(i,c) - icefrac(i,c)
     end do
!
! Ensure that surface fractions add up to 1
!
     call verify_fractions (c, ncol)
     call t_stopf ('update_ocnice')
     
     return
  end subroutine update_ocnice

  subroutine verify_fractions (c, ncol)
!-----------------------------------------------------------------------
!
! Purpose: Ensure that surface fractions are valid
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------
     use phys_grid, only: get_ncols_p

     implicit none

! Arguments

     integer, intent(in) :: c     ! chunk index
     integer, intent(in) :: ncol  ! number of columns
!
! Local Workspace
!
     integer :: i                      ! loop index
     logical :: bad                    ! flag
     real(r8) :: icesv, ocnsv, lndsv   ! saved values
     real(r8) :: sfcfrac               ! icefrac + ocnfrac + landfrac
     real(r8) :: delta                 ! 1. - sum of sfc fractions

     call t_startf ('verify_fractions')
!
! Buffer previous fractions for output of flns fsns surface fluxes for individual
! components
!
     call bounding (icefrac(1,c), c, ncol, 'ICE')
     call bounding (ocnfrac(1,c), c, ncol, 'OCN')
     call bounding (landfrac(1,c), c, ncol, 'LND')

     bad = .false.
     do i=1,ncol
        sfcfrac = icefrac(i,c) + ocnfrac(i,c) + landfrac(i,c)
        delta   = 1. - sfcfrac
        if (abs (delta) > 10.*epsilon (1._r8)) then
           bad = .true.
           icesv = icefrac(i,c)
           ocnsv = ocnfrac(i,c)
           lndsv = landfrac(i,c)
        end if
     end do

     if (bad) then
        write(6,*)'VERIFY_FRACTIONS: total surface fraction: ', icesv, ocnsv, lndsv
        call endrun ()
     end if

     call t_stopf ('verify_fractions')

     return
  end subroutine verify_fractions

  subroutine bounding (frac, c, ncol, string)
!-----------------------------------------------------------------------
!
! Purpose: Ensure that input array frac is bounded between 0 and 1 to
!          within roundoff.  Only SOM currently requires the epsilon ()
!          application.
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------
     implicit none
!
! Arguments
!
     real(r8), intent(in) :: frac(pcols)      ! surface fraction (land, ocn, or ice)
     integer, intent(in) :: c                 ! chunk index
     integer, intent(in) :: ncol              ! number of columns
     character(len=*), intent(in) :: string   ! string for error print
!
! Local Workspace
!
     integer :: i, isave   ! loop index
     integer :: csave      ! chunk index of bad point
     logical :: bad        ! flag indicates some bad points found
     real(r8) :: fracsave  ! bad fraction value

     bad = .false.
     do i=1,ncol
        if (frac(i) < -10.*epsilon (1._r8) .or. frac(i) > 1. + 10.*epsilon (1._r8)) then
           bad = .true.
           fracsave = frac(i)
           isave = i
           csave = c
        end if
     end do

     if (bad) then
        write(6,*)'BOUNDING: bad frac=',fracsave,' at i,c=', isave, csave, ' ', string
        call endrun ()
     end if

     return
  end subroutine bounding

end module comsrf
