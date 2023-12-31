#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs. non-PVP.
! This is due to the fact that spectral coefficients are stored consecutively
! along diagonals of M-N wavenumber space when the target architecture is 
! PVP (optimal for vectorization), and along total wavenumber N otherwise
! (optimal for message-passing).
#if ( defined PVP )
subroutine tstep(n       ,zdt     ,ztdtsq  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! the n(th) diagonal. (Note, n is distinct from the two dimensional
! wavenumber which is also often denoted n.) The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: tstep.F90,v 1.2.28.3 2003/08/13 19:16:13 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------
   use pmgrid
   use pspect
   use comspe
!-----------------------------------------------------------------------
   implicit none
!------------------------------Commons----------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
use commap
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: n            ! index of spectral diagonal being calculated
!                            this call (not two dimensional wavenumber)
   real(r8), intent(in) :: zdt             ! timestep, dt (seconds)
   real(r8), intent(in) :: ztdtsq(2*pnmax) ! dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
   real(r8) z(2*pmmax,plev) ! workspace for computation of spectral array d
   real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
   real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
   real(r8) ztemp           ! temporary workspace

   integer ne           ! index into ztdtsq
   integer m            ! diagonal element (index) of complex array
   integer k,kk         ! level indices
   integer irh          ! index into levels of spectral arrays
   integer irhr,irhi    ! index into real, imaginary coefficients
   integer isp          ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
! Set offsets for beginning of diagonal being calculated this call
!
   isp = nco2(n) - 2
   ne = 2*(n-1)
   do k=1,plev
!
! Coefficients for diagonal terms
!
      hhref = 0.5*href(k,k)
      hbps = 0.5*bps(k)
!
! Loop along current diagonal (in spectral space)
! Add lnps and diagonal (vertical space) T terms to d(t-1)
!
      do m=1,2*nm(n)
         d(isp+m,k) = d(isp+m,k) + &
            ztdtsq(ne+m)*(hhref*t(isp+m,k) + hbps*alps(isp+m))
      end do
      if (k.lt.plev) then
         do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to d(t-1)
!
            hhref = 0.5*href(kk,k)
            do m=1,2*nm(n)
               d(isp+m,k) = d(isp+m,k) + ztdtsq(ne+m)*hhref*t(isp+m,kk)
            end do
         end do
      end if
   end do                    ! k=1,plev (calculation level)
!
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
   do k=1,plev
      do m=1,2*pmmax
         z(m,k) = 0.
      end do
   end do
   do k=1,plev
!
! Initialize offset for diagonals (inner loop over levels)
! Start inner loop over levels (for matrix multiply)
!
      irhr = nco2(n) - 3
      irhi = irhr + 1
      do kk=1,plev
!
! Multiply right hand side by inverse matrix
!
!DIR$ IVDEP
         do m=1,nm(n)
            z(2*m-1,k) = z(2*m-1,k) + bm1(kk,k,m+n-1)*d(irhr+2*m,kk)
            z(2*m  ,k) = z(2*m  ,k) + bm1(kk,k,m+n-1)*d(irhi+2*m,kk)
         end do
      end do                  ! inner loop over levels
   end do                    ! outer loop over levels
!
! Move solution for divergence to d
!
   irh = nco2(n) - 2
   do k=1,plev
      do m=1,2*nm(n)
         d(irh+m,k) = z(m,k)
      end do
   end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
   do k=1,plev
      ztemp = zdt*hypd(k)/hypi(plevp)
      do m=1,2*nm(n)
         alps(isp+m) = alps(isp+m) - ztemp*d(isp+m,k)
      end do
   end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
   do k=1,plev
      do kk=1,plev
         ztemp = zdt*tau(kk,k)
         do m=1,2*nm(n)
            t(isp+m,k) = t(isp+m,k) - ztemp*d(isp+m,kk)
         end do
      end do
   end do
!
   return
#else
   subroutine tstep(lm      ,zdt     ,ztdtsq  )
!-----------------------------------------------------------------------
!
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: tstep.F90,v 1.2.28.3 2003/08/13 19:16:13 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use pspect
      use comspe
      use commap

      implicit none

!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: lm            ! local Fourier wavenumber index

      real(r8), intent(in) :: zdt             ! timestep, dt (seconds)
      real(r8), intent(in) :: ztdtsq(pnmax)   ! dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
      real(r8) z(2*pnmax,plev) ! workspace for computation of spectral array d
      real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
      real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
      real(r8) ztemp           ! temporary workspace

      integer m            ! global wavenumber index
      integer n,j          ! 2-d wavenumber index
      integer k,kk         ! level indices
      integer lmr,lmc      ! real and imaginary spectral indices
      integer ir,ii        ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
!
      m  = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do k=1,plev
!
! Coefficients for diagonal terms
!
         hhref = 0.5*href(k,k)
         hbps = 0.5*bps(k)
!
! Loop along total wavenumber index (in spectral space)
! Add lnps and diagonal (vertical space) T terms to d(t-1)
!
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            d(ir,k) = d(ir,k) + ztdtsq(n+m-1)*(hhref*t(ir,k) + hbps*alps(ir))
            d(ii,k) = d(ii,k) + ztdtsq(n+m-1)*(hhref*t(ii,k) + hbps*alps(ii))
         end do
         if (k.lt.plev) then
            do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to d(t-1)
!
               hhref = 0.5*href(kk,k)
               do n=1,nlen(m)
                  ir = lmc + 2*n - 1
                  ii = ir + 1
                  d(ir,k) = d(ir,k) + ztdtsq(n+m-1)*hhref*t(ir,kk)
                  d(ii,k) = d(ii,k) + ztdtsq(n+m-1)*hhref*t(ii,kk)
               end do
            end do
         end if
      end do                    ! k=1,plev (calculation level)
!
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
      do k=1,plev
         do j=1,2*pnmax
            z(j,k) = 0.
         end do
      end do
      do k=1,plev
         do kk=1,plev
!
! Multiply right hand side by inverse matrix
!
!DIR$ IVDEP
            do n=1,nlen(m)
               ir = lmc + 2*n - 1
               ii = ir + 1
               z(2*n-1,k) = z(2*n-1,k) + bm1(kk,k,m+n-1)*d(ir,kk)
               z(2*n  ,k) = z(2*n  ,k) + bm1(kk,k,m+n-1)*d(ii,kk)
            end do
         end do                  ! inner loop over levels
      end do                    ! outer loop over levels
!
! Move solution for divergence to d
!
      do k=1,plev
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            d(ir,k) = z(2*n-1,k)
            d(ii,k) = z(2*n  ,k)
         end do
      end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
      do k=1,plev
         ztemp = zdt*hypd(k)/hypi(plevp)
         do n=1,nlen(m)
            ir = lmc + 2*n - 1
            ii = ir + 1
            alps(ir) = alps(ir) - ztemp*d(ir,k)
            alps(ii) = alps(ii) - ztemp*d(ii,k)
         end do
      end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
      do k=1,plev
         do kk=1,plev
            ztemp = zdt*tau(kk,k)
            do n=1,nlen(m)
               ir = lmc + 2*n - 1
               ii = ir + 1
               t(ir,k) = t(ir,k) - ztemp*d(ir,kk)
               t(ii,k) = t(ii,k) - ztemp*d(ii,kk)
            end do
         end do
      end do
!
      return
#endif
   end subroutine tstep

