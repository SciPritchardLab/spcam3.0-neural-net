#include <misc.h>
#include <params.h>

subroutine difcor(klev    ,ztodt   ,delps   ,u       ,v       , &
                  qsave   ,pdel    ,pint    ,t       ,tdif    , &
                  udif    ,vdif    ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Add correction term to t and q horizontal diffusions and
! determine the implied heating rate due to momentum diffusion.
! 
! Method: 
! 1. Add correction term to t and q horizontal diffusions. This term
! provides a partial correction of horizontal diffusion on hybrid (sigma)
! surfaces to horizontal diffusion on pressure surfaces. The appropriate
! function of surface pressure (delps, which already contains the diffusion
! coefficient and the time step) is computed during the transform
! from  spherical harmonic coefficients to grid point values. This term
! can only be applied in the portion of the vertical domain in which
! biharmonic horizontal diffusion is employed. In addition, the term is
! unnecessary on pure pressure levels.
!
! 2. Determine the implied heating rate due to momentum diffusion in order
! to conserve total energy and add it to the temperature.
! Reduce complex matrix (ac) to upper Hessenburg matrix (ac)
!
! Author: D. Williamson
! 
!-----------------------------------------------------------------------
!
! $Id: difcor.F90,v 1.1.2.1 2002/06/15 13:46:56 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use physconst, only: cpair, cpvir

  implicit none

#include <comctl.h>
#include <comhyb.h>

!------------------------------Arguments--------------------------------

  integer , intent(in) :: klev                  ! k-index of top hybrid level
  integer , intent(in) :: nlon                  ! longitude dimension
  real(r8), intent(in) :: ztodt                 ! twice time step unless nstep = 0
  real(r8), intent(in) :: delps(plond)          ! srf press function for correction
  real(r8), intent(in) :: u(plond,plev)         ! u-wind
  real(r8), intent(in) :: v(plond,plev)         ! v-wind
  real(r8), intent(in) :: qsave(plond,plev)     ! moisture fm prv fcst
  real(r8), intent(in) :: pdel(plond,plev)      ! pdel(k) = pint(k+1) - pint(k)
  real(r8), intent(in) :: pint(plond,plevp)     ! pressure at model interfaces
  real(r8), intent(inout) :: t(plond,plev)      ! temperature
  real(r8), intent(inout) :: tdif(plond,plev)   ! initial/final temperature diffusion
  real(r8), intent(inout) :: udif(plond,plev)   ! initial/final u-momentum diffusion
  real(r8), intent(inout) :: vdif(plond,plev)   ! initial/final v-momentum diffusion

!---------------------------Local workspace-----------------------------

  integer i,k                ! longitude, level indices
  real(r8) tcor(plond,plev)  ! temperature correction term
!-----------------------------------------------------------------------
!
! Compute the pressure surface correction term for horizontal diffusion of
! temperature. 
!
  do k=klev,plev
     if (k==1) then
        do i=1,nlon
           tcor(i,k) = delps(i)*0.5/pdel(i,k)*(hybi(k+1)*(t(i,k+1)-t(i,k)))*pint(i,plevp)
        end do
     else if (k==plev) then
        do i=1,nlon
           tcor(i,k) = delps(i)*0.5/pdel(i,k)*(hybi(k)*(t(i,k)-t(i,k-1)))*pint(i,plevp)
        end do
     else
        do i=1,nlon
           tcor(i,k) = delps(i)*0.5/pdel(i,k)*(hybi(k+1)*(t(i,k+1)-t(i,k)) + &
                       hybi(k  )*(t(i,k)-t(i,k-1)))*pint(i,plevp)
        end do
     end if
  end do
!
! Add the temperture diffusion correction to the diffusive heating term 
! and to the temperature.
!
  if (.not.adiabatic .and. .not.ideal_phys) then
     do k=klev,plev
        do i=1,nlon
           tdif(i,k) = tdif(i,k) + tcor(i,k)/ztodt
           t(i,k) = t(i,k) + tcor(i,k)
        end do
     end do
!
! Convert momentum diffusion tendencies to heating rates in order to 
! conserve internal energy. Add the heating to the temperature and to 
! diffusive heating term.
!
     do k=1,plev
        do i=1,nlon
           t(i,k) = t(i,k) - ztodt * (u(i,k)*udif(i,k) + v(i,k)*vdif(i,k)) / &
                                     (cpair*(1. + cpvir*qsave(i,k)))
           tdif(i,k) = tdif(i,k) - (u(i,k)*udif(i,k) + v(i,k)*vdif(i,k)) / &
                                   (cpair*(1. + cpvir*qsave(i,k)))
        end do
     end do
  end if
 
  return
end subroutine difcor

