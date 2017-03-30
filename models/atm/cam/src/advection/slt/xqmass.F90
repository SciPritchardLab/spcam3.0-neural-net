#include <misc.h>
#include <params.h>

subroutine xqmass(cwava   ,etamid  ,w       ,qo      ,qn      , &
                  xo      ,xn      ,pdela   ,pdelb   ,hwxal   , &
                  hwxbl   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute comtribution of current latitude to global integrals necessary
! to compute the fixer for the non-water constituents.
! 
! Method: 
! 
! Author: J. Olson, March 1994
! 
!-----------------------------------------------------------------------
!
! $Id: xqmass.F90,v 1.1.2.3 2002/11/22 02:11:08 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use constituents, only: pcnst

  implicit none

!---------------------------Arguments-----------------------------------
  real(r8), intent(in)  :: cwava                ! normalization factor
  real(r8), intent(in)  :: etamid(plev)         ! vertical coords at midpoints 
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: qo(plond,plev      ) ! q old            (pre -SLT)
  real(r8), intent(in)  :: qn(plond,plev      ) ! q new            (post-SLT)
  real(r8), intent(in)  :: xo(plond,plev,pcnst) ! old constituents (pre -SLT)
  real(r8), intent(in)  :: xn(plond,plev,pcnst) ! new constituents (post-SLT)
  real(r8), intent(in)  :: pdela(plond,plev)    ! pressure diff between interfaces
  integer , intent(in) :: nlon                  ! number of longitudes
                                                ! based pure pressure part of hybrid grid
  real(r8), intent(in)  :: pdelb(plond,plev)    ! pressure diff between interfaces
                                                ! based sigma part of hybrid grid
  real(r8), intent(inout) :: hwxal(pcnst,4)     ! partial integrals (weighted by pure
                                                ! pressure part of hybrid pressures)
  real(r8), intent(inout) :: hwxbl(pcnst,4)     ! partial integrals (weighted by sigma
                                                ! part of hybrid pressures)
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                     ! longitude index
  integer k                     ! level index
  integer m                     ! constituent index
  integer n                     ! index for partial integral
  real(r8) a                    ! integral constant
  real(r8) xdx,xq1,xqdq,xdxq1   ! work elements
  real(r8) xdxqdq               ! work elements
  real(r8) hwak(4),hwbk(4)      ! work arrays
  real(r8) q1 (plond,plev)      ! work array
  real(r8) qdq(plond,plev)      ! work array
  real(r8) hwalat(4)            ! partial integrals (weighted by pure
!                               ! pressure part of hybrid pressures)
  real(r8) hwblat(4)            ! partial integrals (weighted by sigma
!                               ! part of hybrid pressures)
  real(r8) etamsq(plev)         ! etamid*etamid
!-----------------------------------------------------------------------
!
  a = cwava*w*0.5
  do k = 1,plev
     etamsq(k) = etamid(k)*etamid(k)
  end do
!
! Compute terms involving water vapor mixing ratio
!
  do k = 1,plev
     do i = 1,nlon
        q1 (i,k) = 1. - qn(i,k)
        qdq(i,k) = qn(i,k)*abs(qn(i,k) - qo(i,k))
     end do
  end do
!
! Compute partial integrals for non-water constituents
!
  do m = 2,pcnst
     do n = 1,4
        hwalat(n) = 0.
        hwblat(n) = 0.
     end do
     do k = 1,plev
        do n = 1,4
           hwak(n) = 0.
           hwbk(n) = 0.
        end do

        do i = 1,nlon
           xdx    = xn(i,k,m)*abs(xn(i,k,m) - xo(i,k,m))
           xq1    = xn(i,k,m)*q1 (i,k)
           xqdq   = xn(i,k,m)*qdq(i,k)
           xdxq1  = xdx      *q1 (i,k)
           xdxqdq = xdx      *qdq(i,k)

           hwak(1) = hwak(1) + xq1   *pdela(i,k)
           hwbk(1) = hwbk(1) + xq1   *pdelb(i,k)
           hwak(2) = hwak(2) + xqdq  *pdela(i,k)
           hwbk(2) = hwbk(2) + xqdq  *pdelb(i,k)
           hwak(3) = hwak(3) + xdxq1 *pdela(i,k)
           hwbk(3) = hwbk(3) + xdxq1 *pdelb(i,k)
           hwak(4) = hwak(4) + xdxqdq*pdela(i,k)
           hwbk(4) = hwbk(4) + xdxqdq*pdelb(i,k)
        end do

        hwalat(1) = hwalat(1) + hwak(1)
        hwblat(1) = hwblat(1) + hwbk(1)
        hwalat(2) = hwalat(2) + hwak(2)*etamid(k)
        hwblat(2) = hwblat(2) + hwbk(2)*etamid(k)
        hwalat(3) = hwalat(3) + hwak(3)*etamid(k)
        hwblat(3) = hwblat(3) + hwbk(3)*etamid(k)
        hwalat(4) = hwalat(4) + hwak(4)*etamsq(k)
        hwblat(4) = hwblat(4) + hwbk(4)*etamsq(k)
     end do
!
! The 0.5 factor arises because Gaussian weights sum to 2
!
     do n = 1,4
        hwxal(m,n) = hwxal(m,n) + hwalat(n)*a
        hwxbl(m,n) = hwxbl(m,n) + hwblat(n)*a
     end do
  end do

  return
end subroutine xqmass






