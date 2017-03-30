#include <misc.h>
#include <params.h>

subroutine flxsice(indx    ,npts    ,pmidm1  ,ubot    ,vbot    , &
                   tbot    ,qbot    ,thbot   ,zbot    ,srfrad  , &
                   ts      ,ltheat  ,fnt     ,dfntdt  ,shf     , &
                   lhf     ,taux    ,tauy    ,tref    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute sea ice to atmosphere surface fluxes of sensible, latent heat
! and stress components
!
! Method: 
! Follows the same basic parameterizations as for ocean surfaces
! 
! Author: Bill Large/M.Vertenstein, Sep. 1995
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use physconst, only: rair, cpair, cpvir, zvir, gravit, stebol

  implicit none

#include <parpbl.h>

!------------------------------Arguments--------------------------------

  integer , intent(in) :: indx(pcols)    ! Column index array (land)
  integer , intent(in) :: npts           ! Number of land points
  real(r8), intent(in) :: pmidm1(pcols)  ! Bottom level pressure
  real(r8), intent(in) :: ubot(pcols)    ! Bottom level u wind
  real(r8), intent(in) :: vbot(pcols)    ! Bottom level v wind
  real(r8), intent(in) :: tbot(pcols)    ! Bottom level temperature
  real(r8), intent(in) :: qbot(pcols)    ! Bottom level specific humidity
  real(r8), intent(in) :: thbot(pcols)   ! Bottom level potential temperature
  real(r8), intent(in) :: zbot(pcols)    ! Bottom level height above surface
  real(r8), intent(in) :: srfrad(pcols)  ! Solar absorbed plus down longwave flux
  real(r8), intent(in) :: ts(pcols)      ! Surface temperature
  real(r8), intent(in) :: ltheat(pcols)  ! Latent heat for given srf conditions
  real(r8), intent(out)::  fnt(pcols)     ! Net surface flux for input conditions (W/m2)
  real(r8), intent(out):: dfntdt(pcols)  ! Net surface flux ts partial derivative (W/m2)

  real(r8), intent(inout) :: shf(pcols)  ! Initial sensible heat flux (W/m2)
  real(r8), intent(inout) :: lhf(pcols)  ! Initial latent heat flux (W/m2)
  real(r8), intent(inout) :: taux(pcols) ! X surface stress (N/m2)
  real(r8), intent(inout) :: tauy(pcols) ! Y surface stress (N/m2)
  real(r8), intent(inout) :: tref(pcols) ! 2m reference temperature
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,ii            ! Column indices
  real(r8) ssq            ! Surface saturation specific humidity
  real(r8) ustar          ! ustar
  real(r8) tstar          ! tstar
  real(r8) qstar          ! qstar
  real(r8) vmag           ! Surface wind magnitude
  real(r8) thvbot         ! Bottom level virtual potential temperature
  real(r8) delt           ! potential T difference (k)
  real(r8) delq           ! specific humidity difference (kg/kg)
  real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
  real(r8) rhn            ! sqrt of neutral exchg coeff (heat & tracers)
  real(r8) ren            ! sqrt of neutral exchange coeff (water)
  real(r8) rd             ! sqrt of exchange coefficient (momentum)
  real(r8) rh             ! sqrt of exchange coefficient(heat & tracers)
  real(r8) re             ! sqrt of exchange coefficient (water) 
  real(r8) hol            ! Ref height (10m) / monin-obukhov length
  real(r8) xsq            ! Temporary variable
  real(r8) xqq            ! Temporary variable
  real(r8) alz            ! ln(zbot/10)
  real(r8) cp             ! Specific heat of moist air
  real(r8) tau            ! Reference height stress
  real(r8) psimh          ! Stability function at ref lev (momentum)
  real(r8) psixh          ! Stability function at ref lev (heat & tracers) 
  real(r8) stable         ! Stability factor
  real(r8) rbot(pcols)    ! Density at bottom model level
  real(r8) dssqdt         ! Derivate of qs wrt surface temperature
  real(r8) dshf(pcols)    ! Ts partial derivative for sensible heat flux
  real(r8) dlhf(pcols)    ! Ts partial derivative for latent heat flux
  real(r8) lwup(pcols)    ! Long wave up surface flux
  real(r8) bn             ! exchange coef funct for interpolation
  real(r8) bh             ! exchange coef funct for interpolation
  real(r8) fac            ! interpolation factor
  real(r8) ln0            ! log factor for interpolation
  real(r8) ln3            ! log factor for interpolation
  real(r8), parameter :: ztref = 2.0 ! reference height for air temperature
!-----------------------------------------------------------------------

!------------------------------Functions--------------------------------
  real(r8) psimhu         ! Unstable part of psimh
  real(r8) psixhu         ! Unstable part of psixh
  real(r8) qsat           ! Saturation humidty of air
  real(r8) dqsatdt        ! Derivative of qsat wrt surface temperature
  real(r8) xd             ! Dummy argument
  real(r8) Tk             ! Temperature (K)
!
  qsat(Tk)    = 640380. / exp(5107.4/Tk)
  dqsatdt(Tk) = (5107.4 / Tk**2) * 640380. / exp(5107.4/Tk)
  psimhu(xd)  = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) + 1.571
  psixhu(xd)  = 2. * log((1. + xd*xd)/2.)
!-----------------------------------------------------------------------
!
! Loop over ice points
!
  do ii=1,npts
     i = indx(ii)
!
!-----------------------------------------------------------------------
! Determine some necessary variables
!-----------------------------------------------------------------------
!
     rbot(i)= pmidm1(i)/(rair*tbot(i))
     vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
     thvbot = thbot(i) * (1.0 + zvir*qbot(i))
     ssq    = 0.98 * qsat(ts(i)) / rbot(i)
     dssqdt = 0.98 * dqsatdt(ts(i)) / rbot(i) 
     delt   = thbot(i) - ts(i)
     delq   = qbot(i) - ssq 
     alz    = log(zbot(i)/zref) 
     cp     = cpair*(1. + cpvir*ssq) 
!
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Determine roots of neutral exchange coefficients
!
     rdn = xkar/log(zref/zzsice)
     rhn = rdn
     ren = rdn
!
! Determine initial guess of ustar,tstar and qstar
!
     ustar = rdn*vmag
     tstar = rhn*delt
     qstar = ren*delq
!
! Compute stability and evaluate all stability functions
! Stable if (thbot > ts or hol > 0 )
!
     hol = xkar * gravit * zbot(i) * (tstar/thvbot + qstar/(1./zvir+qbot(i))) / ustar**2
     hol = sign( min(abs(hol),10._r8), hol )
     stable = 0.5 + sign(0.5_r8 , hol)
     xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
     xqq   = sqrt(xsq)
     psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
     psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift all coeffs to measurement height and stability
!
     rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
     rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
     re = ren / (1.+ren/xkar*(alz-psixh))
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
     ustar = rd * vmag 
     tstar = rh * delt 
     qstar = re * delq 
!
!---------------------------------------------------------------
! Second iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Recompute stability & evaluate all stability functions  
! Stable if (thbot > ts or hol > 0 )
! 
     hol = xkar * gravit * zbot(i) * (tstar/thvbot + qstar/(1./zvir+qbot(i))) / ustar**2
     hol = sign( min(abs(hol),10._r8), hol )
     stable = 0.5 + sign(0.5_r8 , hol)
     xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
     xqq   = sqrt(xsq)
     psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
     psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift all coeffs to measurement height and stability
!
     rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
     rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
     re = ren / (1.+ren/xkar*(alz-psixh)) 
!
!---------------------------------------------------------------
! Compute the fluxes
!---------------------------------------------------------------
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
     ustar = rd * vmag 
     tstar = rh * delt 
     qstar = re * delq 
!
! Compute surface stress components
!
     tau     =  rbot(i) * ustar * ustar 
     taux(i) = -tau * ubot(i) / vmag 
     tauy(i) = -tau * vbot(i) / vmag 
!
! Compute heat flux components at current surface temperature
! (Define positive latent and sensible heat as upwards into the atm)
!
     shf(i) = -cp * tau * tstar / ustar 
     lhf(i) = -ltheat(i) * tau * qstar / ustar
     lwup(i) = stebol * ts(i)**4 
!
! Compute net surface flux surface temperature derivative at the current
! surface temperature (ignore the variation of the exchange coefficients
! with temperature).
!
     dshf(i) = cp * rbot(i) * rd*rh * vmag
     dlhf(i) = ltheat(i) * rbot(i) * rd*re * vmag * dssqdt
!
! Compute net surface flux at current surface temperature
! (Define positive net flux as downwards into surface)
!
     fnt(i) = srfrad(i) - lwup(i) - shf(i) - lhf(i)
!
! Compute derivate of net surface flux (ignore changes due to radiation)
!
     dfntdt(i) = -(dshf(i) + dlhf(i)) - stebol * 4.*ts(i)**3 
!
!---------------------------------------------------------------
! Following Geleyn(1988), interpolate ts to fixed height zref
!---------------------------------------------------------------
!
! Compute function of exchange coefficients. Assume that 
! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
! 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh 
!
     bn = xkar/rdn
     bh = xkar/rh
!
! Interpolation factor for stable and unstable cases
!
     ln0 = log(1.0 + (ztref/zbot(i))*(exp(bn) - 1.0))
     ln3 = log(1.0 + (ztref/zbot(i))*(exp(bn - bh) - 1.0))
     fac = (ln0 - ztref/zbot(i)*(bn - bh))/bh * stable + (ln0 - ln3)/bh * (1.-stable)
     fac = min(max(fac,0._r8),1._r8)
!
! Actual interpolation
!
     tref(i) = ts(i) + (tbot(i) - ts(i))*fac

  end do
!
  return
end subroutine flxsice

