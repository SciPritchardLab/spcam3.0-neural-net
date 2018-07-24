#include <misc.h>
#include <params.h>

module rnozunit

! Implements radon, unit, and psuedo-ozone test tracers.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plev, plat
  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

! Public interfaces
  public init_rn    ! initialize RN distribution
  public init_unit  ! initialize unit tracer
  public init_oz    ! initialize psuedo-ozone distribution
  public tend_rn    ! tendency due to radon decay
  public flux_rn    ! radon surface flux

contains

!======================================================================
subroutine init_rn(q)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 1: radon
!    1) Radon, init to zero, surface fluxes from WCRP95, 5.5
!       day e-folding decay.
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! copied from test_tracers.F90:initesttr()  D.Bundy Oct 2002
! 
!-----------------------------------------------------------------------
   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
!-----------------------------------------------------------------------
!
! Initialize radon tracer to zero.
!
   q = 0.0
end subroutine init_rn

!======================================================================
subroutine init_unit(q)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 2.
!    2) conserved unit tracer
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! copied from test_tracers.F90:initesttr()  D.Bundy Oct 15 2002
!-----------------------------------------------------------------------
   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
!-----------------------------------------------------------------------
! Initialize conserved unit tracer.

   q = 1.0
end subroutine init_unit

!======================================================================
subroutine init_oz(q)

!----------------------------------------------------------------------- 
! Purpose:
! Initialize test tracer 3.
!    3) ozone-like tracer, init to 1.e-9 above ~100mb, zero
!       elsewhere, re-zero the bottom level at each timestep.
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! copied from test_tracers.F90:initesttr()  D.Bundy Oct 15 2002
!-----------------------------------------------------------------------
   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air

! Local variables
   integer :: k
!-----------------------------------------------------------------------
! Initialize strat tracer to 1.e-9 above 100mb
!
   q = 0.0
   do k = 1, 5
      q(:,k,:) = 1.e-9
   end do
end subroutine init_oz

!======================================================================

subroutine tend_rn(ncol, rn, deltat, rnsnk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 1 (radon)
! 
! Method:
!-------------------------Code History----------------------------------
!
! test_tracers.F90:rndecay()
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! D. Bundy Oct 2002
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: rn(pcols,pver)    ! radon mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: rnsnk(pcols,pver) ! conversion rate
                                              ! (kg rn /(s kg moist air))
!--------------------------Local Variables------------------------------
!
   real(r8), parameter :: a = 2.10e-6  ! lifetime
!-----------------------------------------------------------------------
!
!   calculate tendencies using Euler Backward
!
   rnsnk(1:ncol,:) = -rn(1:ncol,:)*a / (1. + a*deltat)

end subroutine tend_rn

!======================================================================

subroutine flux_rn(ncol, rlat, rlon, landfrac, flux )
!----------------------------------------------------------------------- 
!
! Purpose: Set surface fluxes for radon for WCRP95 RN-PB simulation.
! 
! Method:
!  The flux is specified non-zero over land between 60S - 70N, except
!  exclude Greenland.
!
!  Flux strength:
!  60S - 60N:  3.69e-21 kg/m^2/s
!  60N - 70N:  (3.69e-21)/2 kg/m^2/s
!
!  This land source is has been adjusted so that the total radon flux is
!  15 kg/yr for a T42 grid.
!
!-------------------------Code History----------------------------------
! test_tracers.F90:rnsfwcrp()
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! Written as pkg_testtrac.F90:flux(): D. Bundy Oct 2002
!
!-----------------------------------------------------------------------

!  use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid, only: pcols, pver
!   use ppgrid
!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments-------------------------------------

   integer,  intent(in)  :: ncol            ! number of atmospheric columns
   real(r8), intent(in)  :: rlat(pcols)     ! current latitudes(radians)
   real(r8), intent(in)  :: rlon(pcols)     ! current longitudes(radians)
   real(r8), intent(in)  :: landfrac(pcols) ! landfraction
   real(r8), intent(out) :: flux(pcols)     ! specified radon flux in kg/m^2/s

!--------------------------Local Variables------------------------------

   real(r8), parameter :: rad2deg=360./6.283185308  ! convert radians to degrees

   integer  :: i
   real(r8) :: landflx   ! land flux
   real(r8) :: landflxn  ! (land flux)/2
   real(r8) :: latdeg    ! latitude in degrees

!--------------------------Statement functions--------------------------

   logical land
   land(i) = nint(landfrac(i)) == 1
!-----------------------------------------------------------------------
!
   landflx = 3.7796e-21   ! scaled so total flux is 15 kg/yr on T42 grid
   landflxn = landflx/2.

   do i = 1, ncol
      flux(i) = 0.
      latdeg = rlat(i) * rad2deg
      if ( latdeg .ge. -60.  .and.  latdeg .le. 60. ) then    ! 60S - 60N
         if ( land(i) ) flux(i) = landflx
      else if ( latdeg .gt. 60. .and. latdeg .le. 70 ) then   ! 60N - 70N
         if (rlon(i)*rad2deg .le. 300.0) then            ! 0 - 300E excludes Greenland
            if ( land(i) ) flux(i) = landflxn
         end if
      end if
   end do
end subroutine flux_rn

end module rnozunit
