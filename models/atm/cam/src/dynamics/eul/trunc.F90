#include <misc.h>
#include <params.h>
subroutine trunc
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: trunc.F90,v 1.1.44.3 2003/08/13 19:16:13 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
!---------------------------Local variables-----------------------------
!
#if ( defined PVP )
   integer n              ! Loop index over diagonals
   integer ik2            ! K+2
#else
   integer m              ! loop index
#endif
!
!-----------------------------------------------------------------------
!
! trunc first evaluates truncation parameters for a general pentagonal 
! truncation for which the following parameter relationships are true
!
! 0 .le. |m| .le. ptrm
!
! |m| .le. n .le. |m|+ptrn for |m| .le. ptrk-ptrn
!
! |m| .le. n .le. ptrk     for (ptrk-ptrn) .le. |m| .le. ptrm
!
! Most commonly utilized truncations include:
!  1: triangular  truncation for which ptrk=ptrm=ptrn
!  2: rhomboidal  truncation for which ptrk=ptrm+ptrn
!  3: trapezoidal truncation for which ptrn=ptrk .gt. ptrm
!
! Simple sanity check
! It is necessary that ptrm .ge. ptrk-ptrn .ge. 0
!
   if (ptrm.lt.(ptrk-ptrn)) then
      write(6,*)'TRUNC: Error in truncation parameters'
      write(6,*)'       ntrm.lt.(ptrk-ptrn)'
      call endrun
   end if
   if (ptrk.lt.ptrn) then
      write(6,*)'TRUNC: Error in truncation parameters'
      write(6,*)'       ptrk.lt.ptrn'
      call endrun
   end if
!
! Evaluate pointers and displacement info based on truncation params
!
! The following ifdef logic seems to have something do with SPMD 
! implementation, although it's not clear how this info is used
! Dave, can you check this with JR?
!
#if ( defined PVP )
   ncoefi(1) = 1
   ik2 = ptrk + 2
   do n=1,pmax
      ncoefi(n+1) = ncoefi(n) + min0(pmmax,ik2-n)
      nalp(n) = ncoefi(n) - 1
      nco2(n) = ncoefi(n)*2
      nm(n) = ncoefi(n+1) - ncoefi(n)
   end do
#else
   nstart(1) = 0
   nlen(1) = ptrn + 1
   do m=2,pmmax
      nstart(m) = nstart(m-1) + nlen(m-1)
      nlen(m) = min0(ptrn+1,ptrk+2-m)
   end do
!      write(6,*)'Starting index  length'
!      do m=1,ptrm+1
!         write(6,'(1x,i14,i8)')nstart(m),nlen(m)
!      end do
#endif
!
! Define break-even point for vector lengths in GRCALC.  Don't implement
! for non-PVP machines
!
#if ( defined PVP )
   ncutoff = pmax
   if (2*nm(1).lt.plev) ncutoff = 0
   do n=2,pmax,2
      if (2*nm(n).lt.plev) then
         ncutoff = n
         goto 100
      end if
   end do
100 continue
   write(6,*)'TRUNC: n cutoff for GRCALC vectorization = ',ncutoff
#endif
!
! Assign wavenumbers  and spectral offsets if not SPMD
!
#if ( ! defined SPMD )
   do m=1,pmmax
      locm(m,0) = m
#if ( ! defined PVP )
      lnstart(m) = nstart(m)
#endif
   enddo
#endif
!
   return
end subroutine trunc
