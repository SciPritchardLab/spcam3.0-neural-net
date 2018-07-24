#include <misc.h>
#include <params.h>
subroutine srfsice(lchnk   ,ncol    , icefrac ,snowh   ,ubot    ,&
                   vbot    ,tbot    , qbot    ,thbot   ,zbot    ,&
                   pmidm1  ,srfrad  , tssub   ,qflx    ,taux    ,&
                   tauy    ,ts      , shflx   ,lhflx   ,lwup    ,&
                   tref    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute sea ice to atmosphere surface fluxes; then compute
! sea ice temperature change.
!
! Method: 
! Temperatures over sea-ice surfaces are specified in 'plevmx' layers of
! fixed thickness and thermal properties.  The forecast temperatures are
! determined from a backward/implicit diffusion calculation using
! linearized sensible/latent heat fluxes. The bottom ocean temperature
! is fixed at -2C, allowing heat flux exchange with underlying ocean.
! Temperature over sea ice is not allowed to exceed melting temperature.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: srfsice.F90,v 1.2.2.2 2002/11/22 02:11:15 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use comsrf, only: plevmx
  use constituents, only: pcnst, pnats
  implicit none
#include <comtsc.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: lchnk                ! chunk identifier
  integer , intent(in) :: ncol                 ! number of atmospheric columns

  real(r8), intent(in) :: icefrac(pcols)       ! seaice fraction
  real(r8), intent(in) :: snowh(pcols)         ! Snow depth (liquid water equivalent)
  real(r8), intent(in) :: ubot(pcols)          ! Bottom level u wind
  real(r8), intent(in) :: vbot(pcols)          ! Bottom level v wind
  real(r8), intent(in) :: tbot(pcols)          ! Bottom level temperature
  real(r8), intent(in) :: qbot(pcols)          ! Bottom level specific humidity
  real(r8), intent(in) :: thbot(pcols)         ! Bottom level potential temperature
  real(r8), intent(in) :: zbot(pcols)          ! Bottom level height above surface
  real(r8), intent(in) :: pmidm1(pcols)        ! Bottom level pressure
  real(r8), intent(in) :: srfrad(pcols)        ! Srf solar abs flux plus down longwave
  real(r8), intent(inout):: tssub(pcols,plevmx) ! Surface/sub-surface temperat
  real(r8), intent(inout):: qflx(pcols,pcnst+pnats)    ! Constituent flux (kg/m2/s)
  real(r8), intent(inout):: taux(pcols)          ! X surface stress (N/m2)
  real(r8), intent(inout):: tauy(pcols)          ! Y surface stress (N/m2)
  real(r8), intent(inout):: ts(pcols)            ! surface temperature (K)
  real(r8), intent(inout):: shflx(pcols)         ! Surface sensible heat flux (J/m2/s)
  real(r8), intent(inout):: lhflx(pcols)         ! Surface latent   heat flux (J/m2/s)
  real(r8), intent(inout):: lwup(pcols)          ! surface longwave up flux (W/m2)
  real(r8), intent(inout):: tref(pcols)          ! 2m reference temperature
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer indx(pcols)       ! column index array (land)
  integer m                 ! constituent index
  integer npts              ! Number of land points
  integer i,ii              ! column indices
  integer k                 ! Sub-surface level index
  integer isrftyp(pcols)    ! Integer type for surface

  real(r8) tsbsf(pcols,plevmx)  ! Non-adjusted srfc/sub-srfc temperatures
  real(r8) fnt(pcols)           ! Net surface flux for input conditions
  real(r8) dfntdt(pcols)        ! ts partial derivative of net srf flx
  real(r8) ltheat(pcols)        ! Latent heat for given sfc conditions
!-----------------------------------------------------------------------
!
! Find sea ice points
!
  npts = 0
  do i=1,ncol
     if (icefrac(i) > 0.) then
        npts = npts + 1
        indx(npts) = i
     end if
  end do
  if (npts==0) return
!
! Set latent heat for evaporation sea ice surface
!
  do ii=1,npts
     i = indx(ii)
     isrftyp(i) = 2
     ltheat(i) = latvap + latice
  end do
!
! Compute surface fluxes, derviatives, and exchange coefficiants
!
  call flxsice(indx    ,npts    ,pmidm1  ,ubot    ,vbot    , &
               tbot    ,qbot    ,thbot   ,zbot    ,srfrad  , &
               tssub(1,1),ltheat,fnt     ,dfntdt  ,shflx   , &
               lhflx   ,taux    ,tauy    ,tref    )
!
! Initialize surface/subsurface temperatures for srftsb
!
  do k=1,plevmx
     do ii=1,npts
        i = indx(ii)
        tsbsf(i,k) = tssub(i,k)
     end do
  end do
!
! Diffusion calculation for temperature
!
  call srftsb(isrftyp ,indx    ,npts    ,fnt     ,dfntdt  , &
              snowh   ,tsbsf   )
! 
! Modification to updated surface temperatures
! Reset temperature to melting point.
!
  do ii=1,npts
     i = indx(ii)
     do k=1,plevmx
        tsbsf(i,k) = min(tsbsf(i,k),tmelt)
     end do
  end do
!
! Update surface and sub-surface temperatures
!
  do k=1,plevmx
     do ii=1,npts
        i = indx(ii)
        tssub(i,k) = tsbsf(i,k)
     end do
  end do
  do ii=1,npts
     i = indx(ii)
     ts(i) = tssub(i,1)
     lwup(i) = stebol * ts(i)**4
  end do
!
! Evaluate constituent fluxes
!
  do ii=1,npts
     i = indx(ii)
     qflx(i,1) = lhflx(i)/ltheat(i)
  end do
!
! Set non-water constituent fluxes to zero
!
  do m=2,pcnst+pnats
     do ii=1,npts
        i = indx(ii)
        qflx(i,m) = 0.
     end do
  end do
!    
  return
end subroutine srfsice

