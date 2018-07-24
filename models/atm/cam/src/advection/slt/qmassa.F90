#include <misc.h>
#include <params.h>

subroutine qmassa(cwava   ,w       ,q3      ,pdel    ,hw1lat  , &
                  nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to mass of constituents
! being advected by slt.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: qmassa.F90,v 1.1.2.2 2002/11/22 02:11:07 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use constituents, only: pcnst, pnats

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: q3(plond,plev,pcnst) ! constituents
  real(r8), intent(in)  :: pdel(plond,plev)     ! pressure diff between interfaces
  real(r8), intent(out) :: hw1lat(pcnst)        ! accumulator
!-----------------------------------------------------------------------
!
!---------------------------Local variables-----------------------------
  integer i,k,m             ! longitude, level, constituent indices
  real(r8) const            ! temporary constant
!-----------------------------------------------------------------------
!
! Integration factor (the 0.5 factor arises because gaussian weights sum to 2)
!
  const = cwava*w*0.5
  do m=1,pcnst
     hw1lat(m) = 0.
  end do
!
! Compute mass integral for water ONLY
!
  do k=1,plev
     do i=1,nlon
        hw1lat(1) = hw1lat(1) + q3(i,k,1)*pdel(i,k)
     end do
  end do
!
! Compute mass integral for non-water constituents (on a DRY basis)
!
  do m=2,pcnst
     do k=1,plev
        do i=1,nlon
           hw1lat(m) = hw1lat(m) + q3(i,k,m)*(1. - q3(i,k,1))*pdel(i,k)
        end do
     end do
  end do

  do m = 1,pcnst
     hw1lat(m) = hw1lat(m)*const
  end do

  return
end subroutine qmassa

