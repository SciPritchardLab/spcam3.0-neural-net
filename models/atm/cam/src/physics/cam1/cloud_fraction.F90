#include <misc.h>
#include <params.h>
module cloud_fraction
  !
  ! Module for cloud fraction routines.
  !
  ! $Id: cloud_fraction.F90,v 1.1.4.4 2003/12/22 21:10:47 jmccaa Exp $
  !
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  private
  save
  !
  ! Public interfaces
  !
  public cldfrc_init ! Inititialization of cloud_fraction run-time parameters
  public cldfrc   !  Computation of cloud fraction
  !
  ! Private data
  !
  real(r8) rhminl                ! minimum rh for low stable clouds
  real(r8) rhminh                ! minimum rh for high stable clouds
  real(r8) sh1,sh2               ! parameters for shallow convection cloud fraction
  real(r8) dp1,dp2               ! parameters for deep convection cloud fraction
  real(r8) premit                ! top pressure bound for mid level cloud
  !
contains  

  subroutine cldfrc_init()
    !
    ! Purpose:
    ! Initialize cloud fraction run-time parameters
    !
    ! Author: J. McCaa
    !    
    use dycore, only: dycore_is, get_resolution

    if ( dycore_is ('LR') ) then
       rhminl = .91
       rhminh = .80
       sh1 = 0.04
       sh2 = 500.0
       dp1 = 0.10
       dp2 = 500.0
       premit = 750.e2  ! top of area defined to be mid-level cloud
    else
       if ( get_resolution() == 'T85' ) then
          rhminl = .91
          rhminh = .70
          sh1 = 0.07
          sh2 = 500.0
          dp1 = 0.14
          dp2 = 500.0
          premit = 250.e2  ! top of area defined to be mid-level cloud
       else
          rhminl = .90
          rhminh = .80
          sh1 = 0.07
          sh2 = 500.0
          dp1 = 0.14
          dp2 = 500.0
          premit = 750.e2  ! top of area defined to be mid-level cloud
       endif
    endif

  end subroutine cldfrc_init

  subroutine cldfrc(lchnk   ,ncol    , &
       pmid    ,temp    ,q       ,omga    , phis, &
       cldtop  ,cldbot  ,cloud   ,clc     ,pdel    , &
       cmfmc   ,cmfmc2  ,landfrac,snowh   ,concld  ,cldst   , &
       ts      ,sst     ,ps      ,zdu     ,ocnfrac ,&
       rhu00   ,relhum  ,dindex )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud fraction 
    ! 
    ! 
    ! Method: 
    ! This calculate cloud fraction using a relative humidity threshold
    ! The threshold depends upon pressure, and upon the presence or absence 
    ! of convection as defined by a reasonably large vertical mass flux 
    ! entering that layer from below.
    ! 
    ! Author: Many. Last modified by Jim McCaa
    ! 
    !-----------------------------------------------------------------------
    use ppgrid
    use physconst, only: cappa, gravit, rair, tmelt
    use cldconst
    use wv_saturation, only: aqsat
    use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
    use dycore,        only: dycore_is, get_resolution

    implicit none

    real(r8), parameter :: pnot = 1.e5       ! reference pressure
    real(r8), parameter :: lapse = 6.5e-3    ! U.S. Standard Atmsophere lapse rate
    real(r8), parameter :: premib = 750.e2   ! bottom pressure bound of middle cloud
    real(r8), parameter :: pretop = 1.0e2    ! pressure bounding high cloud
    !
    ! Arguments
    !
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: dindex                 ! 0 or 1 to perturb rh

    real(r8), intent(in) :: pmid(pcols,pver)      ! midpoint pressures
    real(r8), intent(in) :: temp(pcols,pver)      ! temperature
    real(r8), intent(in) :: q(pcols,pver)         ! specific humidity
    real(r8), intent(in) :: omga(pcols,pver)      ! vertical pressure velocity
    real(r8), intent(in) :: cldtop(pcols)         ! top level of convection
    real(r8), intent(in) :: cldbot(pcols)         ! bottom level of convection
    real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
    real(r8), intent(in) :: cmfmc2(pcols,pverp)   ! shallow convective mass flux--m sub c
    real(r8), intent(in) :: snowh(pcols)          ! snow depth (liquid water equivalent)
    real(r8), intent(in) :: pdel(pcols,pver)      ! pressure depth of layer
    real(r8), intent(in) :: landfrac(pcols)       ! Land fraction
    real(r8), intent(in) :: ocnfrac(pcols)        ! Ocean fraction
    real(r8), intent(in) :: ts(pcols)             ! surface temperature
    real(r8), intent(in) :: sst(pcols)            ! sea surface temperature
    real(r8), intent(in) :: ps(pcols)             ! surface pressure
    real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection
    real(r8), intent(in) :: phis(pcols)           ! surface geopotential

    !
    ! Output arguments
    !
    real(r8), intent(out) :: cloud(pcols,pver)     ! cloud fraction
    real(r8), intent(out) :: clc(pcols)            ! column convective cloud amount
    real(r8), intent(out) :: cldst(pcols,pver)     ! cloud fraction
    real(r8), intent(out) :: rhu00(pcols,pver)     ! RH threshold for cloud
    real(r8), intent(out) :: relhum(pcols,pver)    ! RH 
    !      real(r8) dmudp                 ! measure of mass detraining in a layer
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) concld(pcols,pver)    ! convective cloud cover
    real(r8) cld                   ! intermediate scratch variable (low cld)
    real(r8) dthdpmn(pcols)         ! most stable lapse rate below 750 mb
    real(r8) dthdp                 ! lapse rate (intermediate variable)
    real(r8) es(pcols,pver)        ! saturation vapor pressure
    real(r8) qs(pcols,pver)        ! saturation specific humidity
    real(r8) rhwght                ! weighting function for rhlim transition
    real(r8) rh(pcols,pver)        ! relative humidity
    real(r8) rhdif                 ! intermediate scratch variable
    real(r8) strat                 ! intermediate scratch variable
    real(r8) theta(pcols,pver)     ! potential temperature
    real(r8) rhlim                 ! local rel. humidity threshold estimate
    real(r8) coef1                 ! coefficient to convert mass flux to mb/d
    real(r8) clrsky(pcols)         ! temporary used in random overlap calc
    real(r8) rpdeli(pcols,pver-1) ! 1./(pmid(k+1)-pmid(k))
    real(r8) rhpert                !the specified perturbation to rh
    real(r8) deepcu                ! deep convection cloud fraction
    real(r8) shallowcu             ! shallow convection cloud fraction

    logical cldbnd(pcols)          ! region below high cloud boundary

    integer i,k                    ! longitude, level indices
    integer kp1
    integer kdthdp(pcols)
    integer numkcld                ! number of levels in which to allow clouds

    real(r8) thetas(pcols)                    ! ocean surface potential temperature
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    !
    ! Statement functions
    !
    logical land

    land(i) = nint(landfrac(i)) == 1
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)

    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and 
    ! stratocumulus clouds do not compete with convective cloud for which one creates 
    ! the most cloud.  They contribute collectively to the total grid-box average cloud 
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated 
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:,:) = 2.0
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01 

    !
    ! Evaluate potential temperature and relative humidity
    !
    call aqsat(temp    ,pmid    ,es      ,qs      ,pcols   , &
         1    ,ncol    ,pver    ,1       ,pver    )
    do k=1,pver
       do i=1,ncol
          theta(i,k)  = temp(i,k)*(pnot/pmid(i,k))**cappa
          rh(i,k)     = q(i,k)/qs(i,k)*(1.0+float(dindex)*rhpert)
          !
          !  record relhum, rh itself will later be modified related with concld
          !
          relhum(i,k) = rh(i,k)
          cloud(i,k)  = 0.
          cldst(i,k)  = 0.
          concld(i,k) = 0.
       end do
    end do
    !
    ! Initialize other temporary variables
    !
    do i=1,ncol
       ! Adjust thetas(i) in the presence of non-zero ocean heights.
       ! This reduces the temperature for positive heights according to a standard lapse rate.
       if(ocnfrac(i).gt.0.01) thetas(i)  = &
            ( sst(i) - lapse * phis(i) / gravit) * (pnot/ps(i))**cappa
       if(ocnfrac(i).gt.0.01.and.sst(i).lt.260.) &
            write(6,*) 'COLDSST: encountered in cldfrc:', lchnk,i,ocnfrac(i),sst(i)
       clc(i) = 0.0
    end do
    coef1 = gravit*864.0    ! conversion to millibars/day

    do k=1,pver-1
       do i=1,ncol
          rpdeli(i,k) = 1./(pmid(i,k+1) - pmid(i,k))
       end do
    end do

    !
    ! Estimate of local convective cloud cover based on convective mass flux
    ! Modify local large-scale relative humidity to account for presence of 
    ! convective cloud when evaluating relative humidity based layered cloud amount
    !
    do k=1,pver
       do i=1,ncol
          concld(i,k) = 0.0
       end do
    end do
    !
    ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !   
#ifndef PERGRO
    do k=1,pver-1
       do i=1,ncol
          shallowcu = max(0.0,min(sh1*log(1.0+sh2*cmfmc2(i,k+1)),0.30))
          deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60))
          concld(i,k) = min(shallowcu + deepcu,0.80)
          rh(i,k) = (rh(i,k) - concld(i,k))/(1.0 - concld(i,k))
       end do
    end do
#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH 
    !====================================================================
    !
    numkcld = pver
    do k=2,numkcld
       kp1 = min(k + 1,pver)
       do i=1,ncol
          !
          cldbnd(i) = pmid(i,k).ge.pretop
          !
          if ( pmid(i,k).ge.premib ) then
             !==============================================================
             ! This is the low cloud (below premib) block
             !==============================================================
             ! enhance low cloud activation over land with no snow cover
             if (land(i) .and. (snowh(i) <= 0.000001)) then
                rhlim = rhminl - 0.10
             else
                rhlim = rhminl
             endif
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else if ( pmid(i,k).lt.premit ) then
             !==============================================================
             ! This is the high cloud (above premit) block
             !==============================================================
             !
             rhlim = rhminh
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else
             !==============================================================
             ! This is the middle cloud block
             !==============================================================
             !
             !       linear rh threshold transition between thresholds for low & high cloud
             !
             rhwght = (premib-(max(pmid(i,k),premit)))/(premib-premit)
             
             if (land(i) .and. (snowh(i) <= 0.000001)) then
                rhlim = rhminh*rhwght + (rhminl - 0.10)*(1.0-rhwght)
             else
                rhlim = rhminh*rhwght + rhminl*(1.0-rhwght)
             endif
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          end if
          !==================================================================================
          ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
          !==================================================================================
          !      !
          !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
          !      !
          rhu00(i,k)=rhlim
          !==================================================================================
       end do
       !
       ! Final evaluation of layered cloud fraction
       !
    end do
    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON 
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND 
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.  
    !
    !==================================================================================
    !
    ! Find most stable level below 750 mb for evaluating stratus regimes
    !
    do i=1,ncol
       ! Nothing triggers unless a stability greater than this minimum threshold is found
       dthdpmn(i) = -0.125
       kdthdp(i) = 0
    end do
    !
    do k=2,pver
       do i=1,ncol
          if (pmid(i,k) >= premib .and. ocnfrac(i).gt. 0.01) then
             ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
             dthdp = 100.0*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
             if (dthdp < dthdpmn(i)) then
                dthdpmn(i) = dthdp
                kdthdp(i) = k     ! index of interface of max inversion
             end if
          end if
       end do
    end do

    ! Also check between the bottom layer and the surface
    ! Only perform this check if the criteria were not met above

    do i = 1,ncol
       if ( kdthdp(i) .eq. 0 .and. ocnfrac(i).gt.0.01) then
          dthdp = 100.0 * (thetas(i) - theta(i,pver)) / (ps(i)-pmid(i,pver))
          if (dthdp < dthdpmn(i)) then
             dthdpmn(i) = dthdp
             kdthdp(i) = pver     ! index of interface of max inversion
          endif
       endif
    enddo

    do i=1,ncol
       if (kdthdp(i) /= 0) then
          k = kdthdp(i)
          kp1 = min(k+1,pver)
          ! Note: strat will be zero unless ocnfrac > 0.01
          strat = min(1._r8,max(0._r8, ocnfrac(i) * ((theta(i,k700)-thetas(i))*.057-.5573) ) )
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
       end if
    end do
    !
    ! AGGREGATE CLOUD CONTRIBUTIONS (cldst should be zero everywhere except at level kdthdp(i))
    !
    do k=1,pver
       do i=1,ncol
          !
          !       which is greater; standard layered cloud amount or stratocumulus diagnosis
          !
          cloud(i,k) = max(cloud(i,k),cldst(i,k))
          !
          !       add in the contributions of convective cloud (determined separately and accounted
          !       for by modifications to the large-scale relative humidity.
          !
          cloud(i,k) = min(cloud(i,k)+concld(i,k), 1.0_r8)
       end do
    end do

    !
    return
  end subroutine cldfrc
end module cloud_fraction
