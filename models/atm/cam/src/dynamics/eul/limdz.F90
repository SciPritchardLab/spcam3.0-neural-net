#include <misc.h>
#include <params.h>

subroutine limdz(f       ,dsig    ,fst     ,fsb     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Apply SCMO limiter to vertical derivative estimates on a vertical
! slice.
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
! $Id: limdz.F90,v 1.3.4.1 2002/06/15 13:47:49 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use constituents, only: pcnst
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   integer plevm1
   parameter( plevm1 = plev - 1 )
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: f(plon,plev,pcnst)   ! input field
   real(r8), intent(in) :: dsig(plev)           ! size of vertical interval

   integer, intent(in) :: nlon
!
! Input/output arguments
!
   real(r8), intent(inout) :: fst(plon,plev,pcnst) ! z-derivative at top of interval
   real(r8), intent(inout) :: fsb(plon,plev,pcnst) ! z-derivative at bot of interval
!
!-----------------------------------------------------------------------
!
!  f       Field values used to compute the discrete differences for
!          each interval in the vertical grid.
!  dsig    Increment in the sigma-coordinate value for each interval.
!  fst     Limited derivative at the top of each interval.
!  fsb     Limited derivative at the bottom of each interval.
!
!---------------------------Local variables-----------------------------
!
   integer i                 ! longitude   index
   integer k                 ! vertical    index
   integer m                 ! constituent index
!
   real(r8) rdsig                ! 1./dsig
   real(r8) deli(plon)           ! simple linear derivative
!
!------------------------------Externals--------------------------------
!
   external scm0
!
!-----------------------------------------------------------------------
!
! Loop over fields.
!
   do m = 1,pcnst
      do k = 1,plev-1
         rdsig = 1.0/dsig(k)
         do i = 1,nlon
            deli(i) = ( f(i,k+1,m) - f(i,k,m) )*rdsig
         end do
         call scm0(nlon,deli,fst(1,k,m),fsb(1,k,m) )
      end do
   end do
!
   return
end subroutine limdz
