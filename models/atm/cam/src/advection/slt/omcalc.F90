#include <misc.h>
#include <params.h>

subroutine omcalc(rcoslat ,d       ,u       ,v       ,dpsl    , &
                  dpsm    ,pmid    ,pdel    ,rpmid   ,pbot    , &
                  omga    ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate vertical pressure velocity (omga = dp/dt)
! 
! Method: 
! First evaluate the expressions for omega/p, then rescale to omega at
! the end.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: omcalc.F90,v 1.1.2.1 2002/06/15 13:46:58 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  implicit none
#include <comhyb.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                 ! lonitude dimension
  real(r8), intent(in) :: rcoslat(nlon)        ! 1 / cos(lat)
  real(r8), intent(in) :: d(plond,plev)        ! divergence
  real(r8), intent(in) :: u(plond,plev)        ! zonal wind * cos(lat)
  real(r8), intent(in) :: v(plond,plev)        ! meridional wind * cos(lat)
  real(r8), intent(in) :: dpsl(plond)          ! longitudinal component of grad ln(ps)
  real(r8), intent(in) :: dpsm(plond)          ! latitudinal component of grad ln(ps)
  real(r8), intent(in) :: pmid(plond,plev)     ! mid-level pressures
  real(r8), intent(in) :: pdel(plond,plev)     ! layer thicknesses (pressure)
  real(r8), intent(in) :: rpmid(plond,plev)    ! 1./pmid
  real(r8), intent(in) :: pbot(plond)          ! bottom interface pressure
  real(r8), intent(out):: omga(plond,plev)     ! vertical pressure velocity
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k               ! longitude, level indices
  real(r8) hkk(plond)       ! diagonal element of hydrostatic matrix
  real(r8) hlk(plond)       ! super diagonal element
  real(r8) suml(plond)      ! partial sum over l = (1, k-1)
  real(r8) vgpk             ! v dot grad ps
  real(r8) tmp              ! vector temporary
!-----------------------------------------------------------------------
!
! Zero partial sum
!
  do i=1,nlon
     suml(i) = 0.
  end do
!
! Pure pressure part: top level
!
  do i=1,nlon
     hkk(i) = 0.5*rpmid(i,1)
     omga(i,1) = -hkk(i)*d(i,1)*pdel(i,1)
     suml(i) = suml(i) + d(i,1)*pdel(i,1)
  end do
!
! sum(k)(v(j)*ps*grad(lnps)*db(j)) part. Not normally invoked since 
! the top layer is normally a pure pressure layer.
!
  if (1>=nprlev) then
     do i=1,nlon
        vgpk = rcoslat(i)*(u(i,1)*dpsl(i) + v(i,1)*dpsm(i))*pbot(i)
        tmp = vgpk*hybd(1)
        omga(i,1) = omga(i,1) + hybm(1)*rpmid(i,1)*vgpk - hkk(i)*tmp
        suml(i) = suml(i) + tmp
     end do
  end if
!
! Integrals to level above bottom
!
  do k=2,plev-1
!
! Pure pressure part
!
     do i=1,nlon
        hkk(i) = 0.5*rpmid(i,k)
        hlk(i) =     rpmid(i,k)
        omga(i,k) = -hkk(i)*d(i,k)*pdel(i,k) - hlk(i)*suml(i)
        suml(i) = suml(i) + d(i,k)*pdel(i,k)
     end do
!
! v(j)*grad(lnps) part
!
     if (k>=nprlev) then
        do i=1,nlon
           vgpk = rcoslat(i)*(u(i,k)*dpsl(i) + v(i,k)*dpsm(i))*pbot(i)
           tmp = vgpk*hybd(k)
           omga(i,k) = omga(i,k) + hybm(k)*rpmid(i,k)*vgpk - hkk(i)*tmp
           suml(i) = suml(i) + tmp
        end do
     end if
  end do
!
! Pure pressure part: bottom level
!
  do i=1,nlon
     hkk(i) = 0.5*rpmid(i,plev)
     hlk(i) =     rpmid(i,plev)
     omga(i,plev) = -hkk(i)*d(i,plev)*pdel(i,plev) - hlk(i)*suml(i)
  end do
!
! v(j)*grad(lnps) part. Normally invoked, but omitted if the model is
! running in pure pressure coordinates throughout (e.g. stratospheric 
! mechanistic model).
!
  if (plev>=nprlev) then
     do i=1,nlon
        vgpk = rcoslat(i)*(u(i,plev)*dpsl(i) + v(i,plev)*dpsm(i))* pbot(i)
        omga(i,plev) = omga(i,plev) + hybm(plev)*rpmid(i,plev)*vgpk - &
             hkk(i)*vgpk*hybd(plev)
     end do
  end if
!
! The above expressions give omega/p. Rescale to omega.
!
  do k=1,plev
     do i=1,nlon
        omga(i,k) = omga(i,k)*pmid(i,k)
     end do
  end do
!
  return
end subroutine omcalc

