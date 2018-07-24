#include <misc.h>
#include <params.h>

subroutine vertinterp_onecol(nlev, pmid, pout, arrin, arrout)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Vertically interpolate input array to output pressure level
! Copy values at boundaries.
! 
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: nlev              ! vertical dimension
  real(r8), intent(in)  :: pmid(nlev)  ! input level pressure levels 
  real(r8), intent(in)  :: pout              ! output pressure level 
  real(r8), intent(in)  :: arrin(nlev) ! input  array
  real(r8), intent(out) :: arrout     ! output value (interpolated)
!--------------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer k               ! indices
  integer kupper     ! Level indices for interpolation
  real(r8) dpu              ! upper level pressure difference
  real(r8) dpl              ! lower level pressure difference
  logical found      ! true if input levels found
  logical error             ! error flag 
!-----------------------------------------------------------------
!
! Initialize index array and logical flags
!
     found  = .false.
     kupper = 1
     error = .false.
!     
! Store level indices for interpolation. 
! If all indices for this level have been found, 
! do the interpolation 
!     
  do k=1,nlev-1
        if ((.not. found) .and. pmid(k)<pout .and. pout<=pmid(k+1)) then
           found = .true.
           kupper = k
        end if
  end do
!
! If we've fallen through the k=1,nlev-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top data level for at least some
! of the longitude points.
!
     if (pout <= pmid(1)) then
        arrout = arrin(1)
     else if (pout >= pmid(nlev)) then
        arrout = arrin(nlev)
     else if (found) then
        dpu = pout - pmid(kupper)
        dpl = pmid(kupper+1) - pout
        arrout = (arrin(kupper)*dpl + arrin(kupper+1)*dpu)/(dpl + dpu)
     else
        error = .true.
     end if
!     
! Error check
!
  if (error) then
     write(6,*)'VERTINTERP: ERROR FLAG'
     call endrun
  end if

  return
end subroutine vertinterp_onecol
