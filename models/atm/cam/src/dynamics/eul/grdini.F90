#include <misc.h>
#include <params.h>

subroutine grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
                  lam     ,phi     ,dphi    ,gw      ,sinlam  , &
                  coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
                  detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize model and extended grid parameters
! Initialize weights for Lagrange cubic derivative estimates
! Initialize weights for Lagrange cubic interpolant
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id: grdini.F90,v 1.1.44.1 2002/06/15 13:47:43 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use rgrid
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
!
! Input arguments
!
   integer, intent(in) :: pmap              ! dimension of artificial vert. grid
!
   real(r8), intent(in) :: etamid(plev)         ! full-level model vertical grid
   real(r8), intent(in) :: etaint(plevp)        ! half-level model vertical grid
   real(r8), intent(in) :: gravit               ! gravitational constant
!
! Output arguments
!
   real(r8), intent(out) :: dlam(platd)          ! longitudinal grid interval (radians)
   real(r8), intent(out) :: lam   (plond,platd)  ! longitudinal coords of extended grid
   real(r8), intent(out) :: phi   (platd)        ! latitudinal  coords of extended grid
   real(r8), intent(out) :: dphi  (platd)        ! latitude intervals (radians)
   real(r8), intent(out) :: gw    (plat)         ! Gaussian weights
   real(r8), intent(out) :: sinlam(plond,platd)  ! sin(lam) model domain only
   real(r8), intent(out) :: coslam(plond,platd)  ! cos(lam) model domain only
   real(r8), intent(out) :: lbasdy(4,2,platd)    ! latitude derivative weights
   real(r8), intent(out) :: lbasdz(4,2,plev)     ! vertical (full levels) deriv weights
   real(r8), intent(out) :: lbassd(4,2,plevp)    ! vertical (half levels) deriv weights
   real(r8), intent(out) :: lbasiy(4,2,platd)    ! Lagrange cubic interp weights (lat.)
   real(r8), intent(out) :: detam (plev)         ! intervals between vertical full levs.
   real(r8), intent(out) :: detai (plevp)        ! intervals between vertical half levs.
!
   integer, intent(out) :: kdpmpf(pmap)      ! artificial full vertical grid indices
   integer, intent(out) :: kdpmph(pmap)      ! artificial half vertical grid indices
!
   real(r8), intent(out) :: cwava(plat)          ! weight applied to global integrals
!
!-----------------------------------------------------------------------
!
!  pmap    Dimension of artificial evenly spaced vertical grid arrays
!  etamid  Full-index hybrid-levels in vertical grid.
!  etaint  Half-index hybrid-levels from sig(1/2) = etaint(1) = 0. to
!          sig(plev+1/2) = etaint(plevp) = 1.
!  gravit  Gravitational constant.
!  dlam    Length of increment in longitude grid.
!  lam     Longitude values in the extended grid.
!  phi     Latitude values in the extended grid.
!  dphi    Interval between latitudes in the extended grid
!  gw      Gauss weights for latitudes in the global grid.  (These sum
!          to 2.0.)
!  sinlam  Sine of longitudes in global grid (no extension points).
!  coslam  Cosine of longitudes in global grid (no extension points).
!  lbasdy  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced latitude grid
!  lbasdz  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (corresponding to model
!          full levels).
!  lbassd  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (corresponding to model
!          half levels).
!  lbasiy  Weights for Lagrange cubic interpolation on the
!          unequally spaced latitude grid
!  detam   Increment between model mid-levels ("full" levels)
!  detai   Increment between model interfaces ("half" levels).
!  kdpmpf  Array of indicies of the model full levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  kdpmph  Array of indicies of the model half levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  cwava   1./(plon*gravit)
!
!---------------------------Local variables-----------------------------
!
   integer j                 ! index
   integer k                 ! index
!
   real(r8) etamln(plev)         ! log(etamid)
   real(r8) etailn(plevp)        ! log(etaint)
   real(r8) detamln(plev)        ! dlog(etamid)
   real(r8) detailn(plevp)       ! dlog(etaint)
!
!-----------------------------------------------------------------------
!
! Initialize extended horizontal grid coordinates.
!
   call grdxy(dlam    ,lam     ,phi     ,gw      ,sinlam  , &
      coslam  )
!
! Basis functions for computing Lagrangian cubic derivatives
! on unequally spaced latitude and vertical grids.
!
   call basdy(phi     ,lbasdy  )
   call basdz(plev    ,etamid  ,lbasdz  )
   call basdz(plevp   ,etaint  ,lbassd  )
!
! Basis functions for computing weights for Lagrangian cubic
! interpolation on unequally spaced latitude grids.
!
   call basiy(phi     ,lbasiy  )
!
! Compute interval lengths in latitudinal grid
!
   do j = 1,platd-1
      dphi(j) = phi(j+1) - phi(j)
   end do
!
! Compute interval lengths in vertical grids.
!
   do k = 1,plev
      etamln(k) = log(etamid(k))
   end do
   do k = 1,plevp
      etailn(k) = log(etaint(k))
   end do
   do k = 1,plev-1
      detam  (k) = etamid(k+1) - etamid(k)
      detamln(k) = etamln(k+1) - etamln(k)
   end do
   do k = 1,plev
      detai  (k) = etaint(k+1) - etaint(k)
      detailn(k) = etailn(k+1) - etailn(k)
   end do
!
! Build artificial evenly spaced vertical grid for use in determining
! vertical position of departure point.
! Build one grid for full model levels and one for half levels.
!
   call vrtmap(plev    ,pmap    ,etamln  ,detamln ,kdpmpf  )
   call vrtmap(plevp   ,pmap    ,etailn  ,detailn ,kdpmph  )
!
! Compute moisture integration constant
!
   do j=1,plat
      cwava(j) = 1./(nlon(j)*gravit)
   end do
!
   return
end subroutine grdini
