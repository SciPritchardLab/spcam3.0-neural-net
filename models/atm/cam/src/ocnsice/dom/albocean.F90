#include <misc.h>
#include <params.h>

subroutine albocean(lchnk   ,ncol    , &
                    coszrs  ,asdir,aldir, &
                    asdif,aldif)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute surface albedos
!
! Method: 
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
! Ocean           Uses solar zenith angle to compute albedo for direct
!                 radiation; diffuse radiation values constant; albedo
!                 independent of spectral interval and other physical
!                 factors such as ocean surface wind speed.
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: albocean.F90,v 1.3.8.4 2003/02/27 00:58:19 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  use comsrf,       only: ocnfrac

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: lchnk            ! chunk identifier
  integer , intent(in) :: ncol             ! number of atmospheric columns

  real(r8), intent(in) :: coszrs(pcols)    ! Cosine solar zenith angle
  real(r8), intent(inout) :: asdir(pcols)  ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r8), intent(inout) :: aldir(pcols)  ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r8), intent(inout) :: asdif(pcols)  ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r8), intent(inout) :: aldif(pcols)  ! Srf alb for diffuse rad  0.7-5.0 micro-ms
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! Longitude index
!-----------------------------------------------------------------------
!
! ocean albedos function of solar zenith angle only, and
! independent of spectral interval:
!
  do i=1,ncol
     if (ocnfrac(i,lchnk) > 0. .and. coszrs(i)>0.) then
        aldir(i)  = (.026/(coszrs(i)**1.7 + .065)) + &
             (.15*(coszrs(i) - 0.10)*(coszrs(i) - 0.50)*(coszrs(i) - 1.00))
        asdir(i)  = aldir(i)
        aldif(i) = 0.06
        asdif(i) = 0.06
     else
!JR These should be inf not zero but accumulating inf in history barfs.
        asdir(i) = 0.
        aldir(i) = 0.
        asdif(i) = 0.
        aldif(i) = 0.
     end if
  end do
!
  return
end subroutine albocean

