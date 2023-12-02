#include <misc.h>
#include <params.h>

subroutine albice(lchnk   ,ncol    ,snowh   ,coszrs  ,asdir  ,&
                  aldir   ,asdif   ,aldif   )

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
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: albice.F90,v 1.2.2.1 2002/06/15 13:48:48 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use comsrf, only :icefrac
  implicit none
#include <albedo.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: lchnk            ! chunk identifier
  integer , intent(in) :: ncol             ! number of atmospheric columns

  real(r8), intent(in) :: snowh(pcols)     ! Snow depth (liquid water equivalent)
  real(r8), intent(in) :: coszrs(pcols)    ! Cosine solar zenith angle
  real(r8), intent(out):: asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r8), intent(out):: aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r8), intent(out):: asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r8), intent(out):: aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! Longitude index
  real(r8) frsnow           ! Horizontal fraction of snow cover
  real(r8) snwhgt           ! Physical snow height
  real(r8) rghsnw           ! Roughness for horizontal snow cover fractn
  real(r8) sasdir(pcols)    ! Snow alb for direct rad  0.2-0.7 micro-ms
  real(r8) saldir(pcols)    ! Snow alb for direct rad  0.7-5.0 micro-ms
  real(r8) sasdif(pcols)    ! Snow alb for diffuse rad  0.2-0.7 micro-ms
  real(r8) saldif(pcols)    ! Snow alb for diffuse rad  0.7-5.0 micro-ms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Initialize all sea ice surface albedos to zero
!
  do i=1,ncol
     if (icefrac(i,lchnk) > 0.) then
        asdir(i) = 0.
        aldir(i) = 0.
        asdif(i) = 0.
        aldif(i) = 0.
     end if
  end do

  do i=1,ncol
     if (icefrac(i,lchnk) > 0. .and. coszrs(i)>0.0) then
        asdir(i)  = sices
        aldir(i)  = sicel
        asdif(i) = asdir(i)
        aldif(i) = aldir(i)
        sasdif(i) = snws
        saldif(i) = snwl
     end if
  end do

  do i=1,ncol
     if (icefrac(i,lchnk) > 0.) then
        if (snowh(i)>0. .and. coszrs(i)>0.) then
           if (coszrs(i)<0.5) then
!
! Zenith angle regime 1 ( coszrs < 0.5 ).
! Set direct snow albedos (limit to 0.98 max)
!
              sasdir(i) = min(0.98_r8,sasdif(i) + (1. - sasdif(i))*0.5* &
                   ((3./(1. + 4.*coszrs(i))) - 1.))
              saldir(i) = min(0.98_r8,saldif(i) + (1. - saldif(i))*0.5* &
                   ((3./(1. + 4.*coszrs(i))) - 1.))
           else
!
! Zenith angle regime 2 ( coszrs >= 0.5 )
!
              sasdir(i) = snws
              saldir(i) = snwl
           end if
!
! Compute both diffuse and direct total albedos
!
           snwhgt = 20.*snowh(i)
           rghsnw = 0.25
           frsnow = snwhgt/(rghsnw + snwhgt)
           asdir(i)  = asdir(i) *(1. - frsnow) + sasdir(i) *frsnow
           aldir(i)  = aldir(i) *(1. - frsnow) + saldir(i) *frsnow
           asdif(i) = asdif(i)*(1. - frsnow) + sasdif(i)*frsnow
           aldif(i) = aldif(i)*(1. - frsnow) + saldif(i)*frsnow
        end if
     end if
  end do
!
  return
end subroutine albice

