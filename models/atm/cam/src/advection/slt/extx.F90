#include <misc.h>
#include <params.h>

subroutine extx (pkcnst, pkdim, fb, kloop)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Copy data to the longitude extensions of the extended array
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: extx.F90,v 1.1.2.3 2002/06/15 13:46:57 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use rgrid
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pkcnst    ! dimension construct for 3-D arrays
  integer , intent(in) :: pkdim     ! vertical dimension
  real(r8), intent(inout) :: fb(plond,pkdim*pkcnst,beglatex:endlatex) ! constituents
  integer,  intent(in) :: kloop ! Limit extent of loop of pkcnst
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! longitude index
  integer j                 ! latitude  index
  integer k                 ! vertical  index
  integer nlond             ! extended longitude dim
  integer i2pi              ! start of eastern long. extension
  integer pk                ! k extent to loop over
!-----------------------------------------------------------------------
!
! Fill west edge points.
!
  pk = pkdim*kloop
  if(nxpt >= 1) then
     do j=beglatex,endlatex
        do k=1,pk
           do i=1,nxpt
              fb(i,k,j) = fb(i+nlonex(j),k,j)
           end do
        end do
     end do
  end if
!
! Fill east edge points
!
  do j=beglatex,endlatex
     i2pi = nxpt + nlonex(j) + 1
     nlond = nlonex(j) + 1 + 2*nxpt
     do k=1,pk
        do i=i2pi,nlond
           fb(i,k,j) = fb(i-nlonex(j),k,j)
        end do
     end do
  end do

  return
end subroutine extx
