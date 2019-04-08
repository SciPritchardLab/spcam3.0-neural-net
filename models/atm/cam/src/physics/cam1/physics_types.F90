#include <misc.h>
! CLOUDBRAIN doesn't actually do anyting here but import some things
#define CLOUDBRAIN
!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  use constituents, only: ppcnst, qmin, cnst_name
  use geopotential, only: geopotential_dse
  use physconst,    only: zvir, gravit, cpair, rair
  use dycore,       only: dycore_is

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend

! Public interfaces

  public physics_update
  public physics_ptend_reset
  public physics_ptend_init
  public physics_dme_adjust  ! adjust dry mass and energy for change in water
                             ! cannot be applied to eul or sld dycores
 
!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,   &! chunk index
          ncol      ! number of active columns
     real(r8), dimension(pcols)                  :: &
          ps,      &! surface pressure
          phis      ! surface geopotential
     real(r8), dimension(pcols,pver)             :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s) 
          pmid,    &! midpoint pressure (Pa) 
          pdel,    &! layer thickness (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          lnpmid,  &! ln(pmid)
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
#ifdef PVBUDGET
          pv,      &  ! Potential vorticity (added by Mike Pritchard)
#endif
#ifdef CLOUDBRAIN
          tap, qap, & ! T,Q after physics.
          qcap, qiap, & ! QC and QI after physics
          dtv, vd01, & ! vertical diffusion of T, Q (for energy checker)
          vap, &      ! SR, also add v-component of wind
#endif
          zm      ! geopotential height above surface at midpoints (m)
     real(r8), dimension(pcols,pver,ppcnst) :: &
          q         ! constituent mixing ratio (kg/kg moist air)

     real(r8), dimension(pcols,pver+1)           :: &
          pint,    &! interface pressure (Pa)
          lnpint,  &! ln(pint)
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(pcols)                  :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     
  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend
     real(r8), dimension(pcols,pver)             :: dtdt, dudt, dvdt
     real(r8), dimension(pcols     )             :: flx_net
     real(r8), dimension(pcols)                  :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend
     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls,               &! true if dsdt is returned
          lu,               &! true if dudt is returned
          lv,               &! true if dvdt is returned
          lq(ppcnst)         ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(pcols,pver)             :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,ppcnst) :: &
          q                  ! consituent tendencies (kg/kg/s)
  end type physics_ptend


!===============================================================================
contains
!===============================================================================

!===============================================================================
  subroutine physics_update(state, tend, ptend, dt)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use geopotential, only: geopotential_dse
    use physconst,    only: cpair, gravit, rair, zvir
    use constituents, only: cnst_get_ind
!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables
    type(physics_tend ), intent(inout)  :: tend    ! Physics tendencies

    real(r8), intent(in) :: dt                     ! time step
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3
!-----------------------------------------------------------------------
    ncol = state%ncol

! Update u,v fields
    if(ptend%lu) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%u  (i,k) = state%u  (i,k) + ptend%u(i,k) * dt
             tend%dudt(i,k) = tend%dudt(i,k) + ptend%u(i,k)
          end do
       end do
    end if

    if(ptend%lv) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%v  (i,k) = state%v  (i,k) + ptend%v(i,k) * dt
             tend%dvdt(i,k) = tend%dvdt(i,k) + ptend%v(i,k)
          end do
       end do
    end if

! Update dry static energy
    if(ptend%ls) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%s(i,k)   = state%s(i,k)   + ptend%s(i,k) * dt
             tend%dtdt(i,k) = tend%dtdt(i,k) + ptend%s(i,k)/cpair
          end do
       end do
    end if

! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    do m = 1, ppcnst
       if(ptend%lq(m)) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
             end do
          end do
! now test for mixing ratios which are too small
          name = trim(ptend%name) // '/' // trim(cnst_name(m))
          call qneg3(trim(name), state%lchnk, ncol, pcols, pver, 1, qmin(m), state%q(1,1,m))
       end if
    end do

! special test for cloud water
    if(ptend%lq(ixcldliq)) then
       if (ptend%name == 'cldwat') then
#ifdef PERGRO
          where (state%q(:ncol,:pver,ixcldliq) < 1.e-12)
             state%q(:ncol,:pver,ixcldliq) = 0.
          end where
#endif
       else if (ptend%name == 'convtran1') then
          where (state%q(:ncol,:pver,ixcldliq) < 1.e-36)
             state%q(:ncol,:pver,ixcldliq) = 0.
          end where
       end if
    end if
    if(ptend%lq(ixcldice)) then
       if (ptend%name == 'cldwat') then
#ifdef PERGRO
          where (state%q(:ncol,:pver,ixcldice) < 1.e-12)
             state%q(:ncol,:pver,ixcldice) = 0.
          end where
#endif
       else if (ptend%name == 'convtran1') then
          where (state%q(:ncol,:pver,ixcldice) < 1.e-36)
             state%q(:ncol,:pver,ixcldice) = 0.
          end where
       end if
    end if

! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(                                                                    &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(1,1,1),state%phis , rair        , gravit      , cpair        , &
            zvir        , state%t     , state%zi    , state%zm    , ncol         )
    end if

! Reset all parameterization tendency flags to false
    call physics_ptend_reset(ptend)

  end subroutine physics_update

!===============================================================================
  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) ptend%s = 0.
    if(ptend%lu) ptend%u = 0.
    if(ptend%lv) ptend%v = 0.
    do m = 1, ppcnst
       if(ptend%lq(m)) ptend%q(:,:,m) = 0.
    end do

    ptend%name  = "none"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .FALSE.
    ptend%lu    = .FALSE.
    ptend%lv    = .FALSE.

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend)
!-----------------------------------------------------------------------
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    ptend%name  = "none"
    ptend%lq(:) = .true.
    ptend%ls    = .true.
    ptend%lu    = .true.
    ptend%lv    = .true.

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================
  subroutine physics_dme_adjust(state, tend, qini, dt)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    ! 
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer 
    !         interfaces and midpoints to the surface pressure. The result is no longer in 
    !         the original hybrid coordinate. 
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    ! 
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    ! 
    !-----------------------------------------------------------------------

    implicit none
    !
    ! Arguments
    !
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    real(r8),            intent(in   ) :: qini(pcols,pver)    ! initial specific humidity
    real(r8),            intent(in   ) :: dt                  ! model physics timestep
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: lchnk         ! chunk identifier
    integer  :: ncol          ! number of atmospheric columns
    integer  :: i,k,m         ! Longitude, level indices
    real(r8) :: fdq(pcols)    ! mass adjustment factor
    real(r8) :: te(pcols)     ! total energy in a layer
    real(r8) :: utmp(pcols)   ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)   ! temp variable for recalculating the initial v values
    !
    !-----------------------------------------------------------------------
    ! verify that the dycore is not sld or eul
    if (dycore_is('SLD') .or. dycore_is('EUL')) return

    lchnk = state%lchnk
    ncol  = state%ncol

    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    do k = 1, pver

       ! adjusment factor is just change in water vapor
       fdq(:ncol) = 1. + state%q(:ncol,k,1) - qini(:ncol,k)

       ! adjust constituents to conserve mass in each layer
       do m = 1, ppcnst
          state%q(:ncol,k,m) = state%q(:ncol,k,m) / fdq(:ncol)
       end do

       if (adjust_te) then
          ! compute specific total energy of unadjusted state (J/kg)
          te(:ncol) = state%s(:ncol,k) + 0.5*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2) 

          ! recompute initial u,v from the new values and the tendencies
          utmp(:ncol) = state%u(:ncol,k) - dt * tend%dudt(:ncol,k)
          vtmp(:ncol) = state%v(:ncol,k) - dt * tend%dvdt(:ncol,k)
          ! adjust specific total energy and specific momentum (velocity) to conserve each
          te     (:ncol)   = te     (:ncol)     / fdq(:ncol)
          state%u(:ncol,k) = state%u(:ncol,k  ) / fdq(:ncol)
          state%v(:ncol,k) = state%v(:ncol,k  ) / fdq(:ncol)
          ! compute adjusted u,v tendencies
          tend%dudt(:ncol,k) = (state%u(:ncol,k) - utmp(:ncol)) / dt
          tend%dvdt(:ncol,k) = (state%v(:ncol,k) - vtmp(:ncol)) / dt

          ! compute adjusted static energy
          state%s(:ncol,k) = te(:ncol) - 0.5*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
       end if

! compute new total pressure variables
       state%pdel  (:ncol,k  ) = state%pdel(:ncol,k  ) * fdq(:ncol)
       state%pint  (:ncol,k+1) = state%pint(:ncol,k  ) + state%pdel(:ncol,k)
       state%lnpint(:ncol,k+1) = log(state%pint(:ncol,k+1))
       state%rpdel (:ncol,k  ) = 1./ state%pdel(:ncol,k  )
    end do

! compute new T,z from new s,q,dp
    if (adjust_te) then
       call geopotential_dse(state%lnpint, state%lnpmid  , state%pint ,  &
            state%pmid  , state%pdel    , state%rpdel,  &
            state%s     , state%q(1,1,1), state%phis , rair, gravit, cpair, zvir, &
            state%t     , state%zi      , state%zm   , ncol      )
    end if

  end subroutine physics_dme_adjust
!-----------------------------------------------------------------------

end module physics_types
