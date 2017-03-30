!======================================================================
! This is an interface for test tracers: 
!      test1 = radon, test2 = unit, test3 = ozone
!
! It uses the rnozunit module to initialize 
!   mixing ratios, fluxes and calculate tendencies
!
! Author D.Bundy Oct 2002
!
!======================================================================

#include <misc.h>
#include <params.h>

module test_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

! Public interfaces
  public test_tracers_register                  ! register constituent
  public test_tracers_implements_cnst           ! true if named constituent is implemented by this package
  public test_tracers_init_cnst                 ! initialize constituent
  public test_tracers_timestep_tend             ! calculate tendencies

! Data from namelist variables
  logical, public :: trace_test1  = .false.     ! true => turn on test tracer code with 1 tracer
  logical, public :: trace_test2  = .false.     ! true => turn on test tracer code with 2 tracers
  logical, public :: trace_test3  = .false.     ! true => turn on test tracer code with 3 tracers

! Private module data
  integer, parameter :: ncnst=3                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'TT_RN', 'TT_UN', 'TT_OZ'/)

  integer :: ixtrct=-999                             ! index of 1st constituent

contains
!======================================================================
subroutine test_tracers_register()
!----------------------------------------------------------------------- 
!
! Purpose: register advected tracers
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

   use physconst,    only: mwdry, cpair
   use constituents, only: cnst_add, advected
   
   implicit none
!---------------------------Local workspace-----------------------------
   integer :: m                                   ! dummy
!-----------------------------------------------------------------------

   if (trace_test1 .or. trace_test2 .or. trace_test3) &
      call cnst_add(cnst_names(1), advected, mwdry, cpair, 0._r8, ixtrct, &
                    readiv=.false.)
   if (trace_test2 .or. trace_test3) &
      call cnst_add(cnst_names(2), advected, mwdry, cpair, 0._r8, m, &
                    readiv=.false.)
   if (trace_test3) &
      call cnst_add(cnst_names(3), advected, mwdry, cpair, 0._r8, m, &
                    readiv=.false.)

end subroutine test_tracers_register
!======================================================================

function test_tracers_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------Arguments---------------------------------

   character(len=*), intent(in) :: name   ! constituent name
   logical :: test_tracers_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
   integer :: m
!-----------------------------------------------------------------------

   test_tracers_implements_cnst = .false.
   do m = 1, ncnst
      if (name == cnst_names(m)) then
         test_tracers_implements_cnst = .true.
         return
      end if
   end do
end function test_tracers_implements_cnst

!===============================================================================
subroutine test_tracers_init_cnst(name, q)

!----------------------------------------------------------------------- 
!
! Purpose: initialize test tracers
!-----------------------------------------------------------------------

  use pmgrid,     only: plon, plev, plat
  use rnozunit,   only: init_rn, init_unit, init_oz

  implicit none

  character(len=*), intent(in) :: name
  real(r8), intent(out), dimension(plon,plev,plat) :: q    ! kg tracer/kg dry air

  if (name == cnst_names(1)) then
     call init_rn(q)
  else if (name == cnst_names(2)) then
     call init_unit(q)
  else if (name == cnst_names(3)) then
     call init_oz(q)
  end if
end subroutine test_tracers_init_cnst

!======================================================================

subroutine test_tracers_timestep_tend(state, ptend, cflx, landfrac, deltat)

!----------------------------------------------------------------------- 
!
! Purpose: compute test tracer mixing ratio tendencies and surface fluxes.
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
  use ppgrid,        only: pcols, pver
  use constituents,  only: ppcnst
  use rnozunit,      only: flux_rn, tend_rn

  implicit none

! Arguments
   type(physics_state), intent(in)  :: state          ! state variables
   type(physics_ptend), intent(out) :: ptend          ! package tendencies
   real(r8),            intent(in)  :: deltat         ! timestep
   real(r8),            intent(in)  :: landfrac(pcols) ! Land fraction
   real(r8), intent(inout) :: cflx(pcols,ppcnst) ! Surface constituent flux (kg/m^2/s)

! Local variables
   integer  :: lchnk           ! chunk identifier
   integer  :: ncol            ! number of atmospheric columns in chunk
   real(r8) :: rlat(pcols)     ! current latitudes(radians)
   real(r8) :: rlon(pcols)     ! current longitudes(radians)
!-----------------------------------------------------------------------

! Initialize output tendency structure
   call physics_ptend_init(ptend)
   ptend%name  = 'test_tracers'

! Initialize chunk id and size
   lchnk = state%lchnk
   ncol  = state%ncol

! Get lat/lon values in radians
   call get_rlat_all_p(lchnk, ncol, rlat)
   call get_rlon_all_p(lchnk, ncol, rlon)

! calculate flux and tendency for tracer 1 (radon)
   if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
      call flux_rn(ncol, rlat, rlon, landfrac, cflx(:,ixtrct))
      ! return RN mixing ratio tendency
      ptend%lq(ixtrct) = .true.
      call tend_rn(ncol, state%q(:,:,ixtrct), deltat, ptend%q(:,:,ixtrct))
   endif

! calculate tendency for tracer 3 (ozone)
   if ( trace_test3) then
      ! return pseudo-O3 mixing ratio tendency required to zero the bottom layer.
      ptend%lq(ixtrct+2)    = .true.
      ptend%q(1:ncol,pver,ixtrct+2) = -state%q(1:ncol,pver,ixtrct+2)/deltat
   endif

end subroutine test_tracers_timestep_tend
!======================================================================

end module test_tracers
