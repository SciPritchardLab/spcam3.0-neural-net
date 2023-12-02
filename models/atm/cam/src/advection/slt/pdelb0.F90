#include <misc.h>
#include <params.h>

subroutine pdelb0(ps      ,pdelb   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the pressure intervals between the interfaces for the "B"
! (surface pressure dependent) portion of the hybrid grid only. 
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: pdelb0.F90,v 1.1.2.1 2002/06/15 13:46:58 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none
#include <comhyb.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon               ! longitude dimension
  real(r8), intent(in) :: ps(plond)          ! surface Pressure
  real(r8), intent(out):: pdelb(plond,plev)  ! pressure difference between interfaces
                                             ! (pressure defined using the "B" part  
                                             ! of the hybrid grid only)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k               ! longitude, level indices
!-----------------------------------------------------------------------
!
! Compute del P(B)
!
  do k = 1,plev
     do i = 1,nlon
        pdelb(i,k) = hybd(k)*ps(i)
     end do
  end do

  return
end subroutine pdelb0

