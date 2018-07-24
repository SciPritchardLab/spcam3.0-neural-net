#include <misc.h>
#include <params.h>

! Note that this routine has 2 complete blocks of code for PVP vs. non-PVP.
! Make sure to make appropriate coding changes where necessary.
#if ( defined PVP )

subroutine grcalcs (irow    ,ztodt   ,grts    ,grths   ,grds    ,&
                    grzs    ,grus    ,gruhs   ,grvs    ,grvhs   ,&
                    grpss   ,grdpss  ,grpms   ,grpls   )
!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.6.6 2003/12/15 18:52:50 hender Exp $
! $Author: hender $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra, ez
   use comhd

   implicit none

!
! Input arguments
!
   integer, intent(in) :: irow         ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
!
! Output arguments: symmetric fourier coefficients
!
   real(r8), intent(out) :: grts(plond,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grths(plond,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grds(plond,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grzs(plond,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grus(plond,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruhs(plond,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1)) 
   real(r8), intent(out) :: grvs(plond,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvhs(plond,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpss(plond)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpss(plond)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpms(plond)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpls(plond)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1s(plond,plev)   ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1s(plond,plev)  ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1s(plond,plev)   ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1s(plond,plev)  ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) zdfac(2*pnmax,plev) ! horiz. diffusion factor (vort,div) (complex)
   real(r8) tqfac(2*pnmax,plev) ! horiz. diffusion factor (t,q) (complex)
   real(r8) alp2(2*pspt)        ! Legendre functions (complex)
   real(r8) dalp2(2*pspt)       ! derivative of Legendre functions (complex)
   real(r8) alpn2(2*pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn2(2*pspt)      ! (a/(n(n+1)))*derivative of Legendre functions (complex)
   real(r8) dlpnz               ! (a/(n(n+1)))*H(0,1) for conversion bet abs and rel vort
   real(r8) zurcor              ! conversion term relating abs. & rel. vort.

   integer k                ! level index
   integer m                ! diagonal element(index) of spectral array
   integer n                ! meridional wavenumber index
   integer ne               ! index into spectral arrays
   integer mn               ! index into spectral arrays
   integer mnc              ! index into spectral arrays
   integer mnev             ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
! Expand polynomials and derivatives to complex form to allow largest 
! possible vector length and multiply by appropriate factors
!
   do n=1,pmax
      ne = n - 1
!dir$ ivdep
      do m=1,nmreduced(n,irow)
         mnc = 2*(m+nalp(n))
         mn = m + nalp(n)
         alp2(mnc-1) = alp(mn,irow)
         alp2(mnc  ) = alp(mn,irow)
         dalp2(mnc-1) = dalp(mn,irow)*ra
         dalp2(mnc  ) = dalp(mn,irow)*ra
         alpn2(mnc-1) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         alpn2(mnc  ) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         dalpn2(mnc-1) = dalp(mn,irow)*(rsq(m+ne)*ra)
         dalpn2(mnc  ) = dalp(mn,irow)*(rsq(m+ne)*ra)
      end do
   end do
   dlpnz = dalpn2(2*nalp(2)+1)
   zurcor = ez*dlpnz
!
! Initialize sums
!
   grzs(:,:)  = 0.
   grds(:,:)  = 0.
   gruhs(:,:) = 0.
   grvhs(:,:) = 0.
   grths(:,:) = 0.
   grpss(:)   = 0.
   grus(:,:)  = 0.
   grvs(:,:)  = 0.
   grts(:,:)  = 0.
   grpls(:)   = 0.
   grpms(:)   = 0.
   grdpss(:)   = 0.

   do k=1,plev
!
! Diffusion factors: expand for longest possible vectors
!
!dir$ ivdep
      do n=1,pnmax
         zdfac(n*2-1,k) = -hdifzd(n,k)
         zdfac(n*2  ,k) = -hdifzd(n,k)
         tqfac(n*2-1,k) = -hdiftq(n,k)
         tqfac(n*2  ,k) = -hdiftq(n,k)
      end do

      gru1s(:,k) = 0.
      gruh1s(:,k) = 0.
      grv1s(:,k) = 0.
      grvh1s(:,k) = 0.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n=1,ncutoff,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            grts (m,k) = grts (m,k) + t(mnev,k)*alp2(mn)
            grths(m,k) = grths(m,k) + t(mnev,k)*alp2(mn)*tqfac(m+ne,k)
            grds(m,k) = grds(m,k) + d(mnev,k)*alp2(mn)
            grzs(m,k) = grzs(m,k) + vz(mnev,k)*alp2(mn)
            gru1s (m,k) = gru1s (m,k) + d(mnev,k)*alpn2(mn)
            gruh1s(m,k) = gruh1s(m,k) + d(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
            grv1s (m,k) = grv1s (m,k) + vz(mnev,k)*alpn2(mn)
            grvh1s(m,k) = grvh1s(m,k) + vz(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
         end do
      end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n=2,ncutoff,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            grus (m,k) = grus (m,k) + vz(mnev,k)*dalpn2(mn)
            gruhs(m,k) = gruhs(m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            grvs (m,k) = grvs (m,k) - d(mnev,k)*dalpn2(mn)
            grvhs(m,k) = grvhs(m,k) - d(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
         end do
      end do
   end do
!
! For short diagonals, repeat above loops with vectorization in vertical,
! instead of along diagonals, to keep vector lengths from getting too short.
!
   if (ncutoff.lt.pmax) then
      do n=ncutoff+1,pmax,2   ! ncutoff guaranteed even
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            do k=1,plev
               grts (m,k) = grts (m,k) + t(mnev,k)*alp2(mn)
               grths(m,k) = grths(m,k) + t(mnev,k)*alp2(mn)*tqfac(m+ne,k)
               grds(m,k) = grds(m,k) + d(mnev,k)*alp2(mn)
               grzs(m,k) = grzs(m,k) + vz(mnev,k)*alp2(mn)
               gru1s (m,k) = gru1s (m,k) + d(mnev,k)*alpn2(mn)
               gruh1s(m,k) = gruh1s(m,k) + d(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
               grv1s (m,k) = grv1s (m,k) + vz(mnev,k)*alpn2(mn)
               grvh1s(m,k) = grvh1s(m,k) + vz(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do

      do n=ncutoff+2,pmax,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            do k=1,plev
               grus (m,k) = grus (m,k) + vz(mnev,k)*dalpn2(mn)
               gruhs(m,k) = gruhs(m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
               grvs (m,k) = grvs (m,k) - d(mnev,k)*dalpn2(mn)
               grvhs(m,k) = grvhs(m,k) - d(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do
   end if                      ! ncutoff.lt.pmax

   do k=1,plev
!
! Combine the two parts of u(m) and v(m)
!
!dir$ ivdep
      do m=1,nmmax(irow)
         grus (2*m-1,k) = grus (2*m-1,k) + gru1s (2*m  ,k)
         gruhs(2*m-1,k) = gruhs(2*m-1,k) + gruh1s(2*m  ,k)
         grus (2*m  ,k) = grus (2*m  ,k) - gru1s (2*m-1,k)
         gruhs(2*m  ,k) = gruhs(2*m  ,k) - gruh1s(2*m-1,k)
         grvs (2*m-1,k) = grvs (2*m-1,k) + grv1s (2*m  ,k)
         grvhs(2*m-1,k) = grvhs(2*m-1,k) + grvh1s(2*m  ,k)
         grvs (2*m  ,k) = grvs (2*m  ,k) - grv1s (2*m-1,k)
         grvhs(2*m  ,k) = grvhs(2*m  ,k) - grvh1s(2*m-1,k)
      end do
!
! Remove Coriolis contribution to absolute vorticity from u(m)
! Correction for u:zeta=vz-ez=(zeta+f)-f
!
      grus(1,k) = grus(1,k) - zurcor
   end do
!
!-----------------------------------------------------------------------
! Computation for single level variables.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=1,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpss (m) = grpss (m) + alps(mnev)*alp2(mn)
         grdpss(m) = grdpss(m) + alps(mnev)*alp2(mn)*hdfst4(ne+(m+1)/2)*ztodt
      end do
   end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=2,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpms(m) = grpms(m) + alps(mnev)*dalp2(mn)
      end do
   end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
   do m=1,nmmax(irow)
      grpls(2*m-1) = -grpss(2*m  )*ra*xm(m)
      grpls(2*m  ) =  grpss(2*m-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalcs


subroutine grcalca (irow    ,ztodt   ,grta    ,grtha   ,grda    ,&
                    grza    ,grua    ,gruha   ,grva    ,grvha   ,&
                    grpsa   ,grdpsa  ,grpma   ,grpla   )

!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.6.6 2003/12/15 18:52:50 hender Exp $
! $Author: hender $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra, ez
   use comhd

   implicit none

!
! Input arguments
!
   integer, intent(in) :: irow         ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
!
! Output arguments: antisymmetric fourier coefficients
!
   real(r8), intent(out) :: grta(plond,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grtha(plond,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grda(plond,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grza(plond,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grua(plond,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruha(plond,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grva(plond,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvha(plond,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpsa(plond)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpsa(plond)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpma(plond)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpla(plond)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1a(plond,plev)   ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1a(plond,plev)  ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1a(plond,plev)   ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1a(plond,plev)  ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) zdfac(2*pnmax,plev) ! horiz. diffusion factor (vort,div) (complex)
   real(r8) tqfac(2*pnmax,plev) ! horiz. diffusion factor (t,q) (complex)
   real(r8) alp2(2*pspt)        ! Legendre functions (complex)
   real(r8) dalp2(2*pspt)       ! derivative of Legendre functions (complex)
   real(r8) alpn2(2*pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn2(2*pspt)      ! (a/(n(n+1)))*derivative of Legendre functions (complex)
   real(r8) dlpnz               ! (a/(n(n+1)))*H(0,1) for conversion bet abs and rel vort
   real(r8) zurcor              ! conversion term relating abs. & rel. vort.

   integer k                ! level index
   integer m                ! diagonal element(index) of spectral array
   integer n                ! meridional wavenumber index
   integer ne               ! index into spectral arrays
   integer mn               ! index into spectral arrays
   integer mnc              ! index into spectral arrays
   integer mnev             ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
! Expand polynomials and derivatives to complex form to allow largest 
! possible vector length and multiply by appropriate factors
!
   do n=1,pmax
      ne = n - 1
!dir$ ivdep
      do m=1,nmreduced(n,irow)
         mnc = 2*(m+nalp(n))
         mn = m + nalp(n)
         alp2(mnc-1) = alp(mn,irow)
         alp2(mnc  ) = alp(mn,irow)
         dalp2(mnc-1) = dalp(mn,irow)*ra
         dalp2(mnc  ) = dalp(mn,irow)*ra
         alpn2(mnc-1) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         alpn2(mnc  ) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         dalpn2(mnc-1) = dalp(mn,irow)*(rsq(m+ne)*ra)
         dalpn2(mnc  ) = dalp(mn,irow)*(rsq(m+ne)*ra)
      end do
   end do

   dlpnz = dalpn2(2*nalp(2)+1)
   zurcor = ez*dlpnz
!
! Initialize sums
!
   grza(:,:)  = 0.
   grda(:,:)  = 0.
   gruha(:,:) = 0.
   grvha(:,:) = 0.
   grtha(:,:) = 0.
   grpsa(:)   = 0.
   grua(:,:)  = 0.
   grva(:,:)  = 0.
   grta(:,:)  = 0.
   grpla(:)   = 0.
   grpma(:)   = 0.
   grdpsa(:)   = 0.

   do k=1,plev
!
! Diffusion factors: expand for longest possible vectors
!
!dir$ ivdep
      do n=1,pnmax
         zdfac(n*2-1,k) = -hdifzd(n,k)
         zdfac(n*2  ,k) = -hdifzd(n,k)
         tqfac(n*2-1,k) = -hdiftq(n,k)
         tqfac(n*2  ,k) = -hdiftq(n,k)
      end do
      
      gru1a(:,k) = 0.
      gruh1a(:,k) = 0.
      grv1a(:,k) = 0.
      grvh1a(:,k) = 0.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n=1,ncutoff,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            grua (m,k) = grua (m,k) + vz(mnev,k)*dalpn2(mn)
            gruha(m,k) = gruha(m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            grva (m,k) = grva (m,k) - d(mnev,k)*dalpn2(mn)
            grvha(m,k) = grvha(m,k) - d(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
         end do
      end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n=2,ncutoff,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            grta (m,k) = grta (m,k) + t(mnev,k)*alp2(mn)
            grtha(m,k) = grtha(m,k) + t(mnev,k)*alp2(mn)*tqfac(m+ne,k)
            grda(m,k) = grda(m,k) + d(mnev,k)*alp2(mn)
            grza(m,k) = grza(m,k) + vz(mnev,k)*alp2(mn)
            gru1a (m,k) = gru1a (m,k) + d(mnev,k)*alpn2(mn)
            gruh1a(m,k) = gruh1a(m,k) + d(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
            grv1a (m,k) = grv1a (m,k) + vz(mnev,k)*alpn2(mn)
            grvh1a(m,k) = grvh1a(m,k) + vz(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
         end do
      end do
   end do
!
! For short diagonals, repeat above loops with vectorization in vertical,
! instead of along diagonals, to keep vector lengths from getting too short.
!
   if (ncutoff.lt.pmax) then
      do n=ncutoff+1,pmax,2   ! ncutoff guaranteed even
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            do k=1,plev
               grua (m,k) = grua (m,k) + vz(mnev,k)*dalpn2(mn)
               gruha(m,k) = gruha(m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
               grva (m,k) = grva (m,k) - d(mnev,k)*dalpn2(mn)
               grvha(m,k) = grvha(m,k) - d(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do

      do n=ncutoff+2,pmax,2
         ne = 2*(n-1)
         do m=1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn = m + 2*nalp(n)
            do k=1,plev
               grta (m,k) = grta (m,k) + t(mnev,k)*alp2(mn)
               grtha(m,k) = grtha(m,k) + t(mnev,k)*alp2(mn)*tqfac(m+ne,k)
               grda(m,k) = grda(m,k) + d(mnev,k)*alp2(mn)
               grza(m,k) = grza(m,k) + vz(mnev,k)*alp2(mn)
               gru1a (m,k) = gru1a (m,k) + d(mnev,k)*alpn2(mn)
               gruh1a(m,k) = gruh1a(m,k) + d(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
               grv1a (m,k) = grv1a (m,k) + vz(mnev,k)*alpn2(mn)
               grvh1a(m,k) = grvh1a(m,k) + vz(mnev,k)*alpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do
   end if                      ! ncutoff.lt.pmax

   do k=1,plev
!
! Combine the two parts of u(m) and v(m)
!
!dir$ ivdep
      do m=1,nmmax(irow)
         grua (2*m-1,k) = grua (2*m-1,k) + gru1a (2*m  ,k)
         gruha(2*m-1,k) = gruha(2*m-1,k) + gruh1a(2*m  ,k)
         grua (2*m  ,k) = grua (2*m  ,k) - gru1a (2*m-1,k)
         gruha(2*m  ,k) = gruha(2*m  ,k) - gruh1a(2*m-1,k)
         grva (2*m-1,k) = grva (2*m-1,k) + grv1a (2*m  ,k)
         grvha(2*m-1,k) = grvha(2*m-1,k) + grvh1a(2*m  ,k)
         grva (2*m  ,k) = grva (2*m  ,k) - grv1a (2*m-1,k)
         grvha(2*m  ,k) = grvha(2*m  ,k) - grvh1a(2*m-1,k)
      end do
   end do
!
!-----------------------------------------------------------------------
! Computation for single level variables.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=1,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpma(m) = grpma(m) + alps(mnev)*dalp2(mn)
      end do
   end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=2,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpsa (m) = grpsa (m) + alps(mnev)*alp2(mn)
         grdpsa(m) = grdpsa(m) + alps(mnev)*alp2(mn)*hdfst4(ne+(m+1)/2)*ztodt
      end do
   end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
   do m=1,nmmax(irow)
      grpla(2*m-1) = -grpsa(2*m  )*ra*xm(m)
      grpla(2*m  ) =  grpsa(2*m-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalca

#else

subroutine grcalcs (irow    ,ztodt   ,grts    ,grths   ,grds    ,&
                    grzs    ,grus    ,gruhs   ,grvs    ,grvhs   ,&
                    grpss   ,grdpss  ,grpms   ,grpls   )
!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.6.6 2003/12/15 18:52:50 hender Exp $
! $Author: hender $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra, ez
   use comhd

   implicit none

!
! Input arguments
!
   integer, intent(in) :: irow         ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
!
! Output arguments: symmetric fourier coefficients
!
   real(r8), intent(out) :: grts(2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grths(2*maxm,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grds(2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grzs(2*maxm,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grus(2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruhs(2*maxm,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1)) 
   real(r8), intent(out) :: grvs(2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvhs(2*maxm,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpss(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpss(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpms(2*maxm)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpls(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1s(2*maxm)        ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1s(2*maxm)       ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1s(2*maxm)        ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1s(2*maxm)       ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) alpn(lpspt)          ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn(lpspt)         ! (a/(n(n+1)))*derivative of Legendre functions (complex)
   real(r8) zurcor               ! conversion term relating abs. & rel. vort.

   integer k                ! level index
   integer lm, m            ! local and global Fourier wavenumber indices of spectral array
   integer mlength          ! number of local wavenumbers
   integer n                ! meridional wavenumber index
   integer ir,ii            ! spectral indices
   integer lmr,lmc          ! spectral indices
   integer lmwave0          ! local index for wavenumber 0
   integer lmrwave0         ! local offset for wavenumber 0
   real(r8) tmp,raxm  ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
   lmwave0 = -1
   lmrwave0 = 0
   dalpn(2) = 0.0
   mlength = numm(iam)
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      if (m .eq. 1) then
         lmwave0 = lm
         lmrwave0 = lmr
      endif
      raxm = ra*xm(m)
      do n=1,nlen(m)
         alpn(lmr+n) = lalp(lmr+n,irow)*rsq(m+n-1)*raxm
         dalpn(lmr+n) = ldalp(lmr+n,irow)*rsq(m+n-1)*ra
      end do
   end do
   zurcor = ez*dalpn(lmrwave0 + 2)
!
! Initialize sums
!
   grzs(:,:)  = 0.
   grds(:,:)  = 0.
   gruhs(:,:) = 0.
   grvhs(:,:) = 0.
   grths(:,:) = 0.
   grpss(:)   = 0.
   grus(:,:)  = 0.
   grvs(:,:)  = 0.
   grts(:,:)  = 0.
   grpls(:)   = 0.
   grpms(:)   = 0.
   grdpss(:)   = 0.

   do k=1,plev
      gru1s(:) = 0.
      gruh1s(:) = 0.
      grv1s(:) = 0.
      grvh1s(:) = 0.
!
! Loop over n for t,q,d,and end of u and v
!
      do lm=1,mlength
         m = locm(lm,iam)
         lmr = lnstart(lm)
         lmc = 2*lmr
         do n=1,nlen(m),2
            ir = lmc + 2*n - 1
            ii = ir + 1
            grts (2*lm-1,k) = grts (2*lm-1,k) + t(ir,k)*lalp(lmr+n,irow)
            grts (2*lm  ,k) = grts (2*lm  ,k) + t(ii,k)*lalp(lmr+n,irow)
!
            tmp = lalp(lmr+n,irow)*hdiftq(n+m-1,k)
            grths(2*lm-1,k) = grths(2*lm-1,k) - t(ir,k)*tmp
            grths(2*lm  ,k) = grths(2*lm  ,k) - t(ii,k)*tmp
!
            grds(2*lm-1,k) = grds(2*lm-1,k) + d(ir,k)*lalp(lmr+n,irow)
            grds(2*lm  ,k) = grds(2*lm  ,k) + d(ii,k)*lalp(lmr+n,irow)
!
            grzs(2*lm-1,k) = grzs(2*lm-1,k) + vz(ir,k)*lalp(lmr+n,irow)
            grzs(2*lm  ,k) = grzs(2*lm  ,k) + vz(ii,k)*lalp(lmr+n,irow)
!
            gru1s (2*lm-1) = gru1s (2*lm-1) + d(ir,k)*alpn(lmr+n)
            gru1s (2*lm  ) = gru1s (2*lm  ) + d(ii,k)*alpn(lmr+n)
!
            tmp = alpn(lmr+n)*hdifzd(n+m-1,k)
            gruh1s(2*lm-1) = gruh1s(2*lm-1) - d(ir,k)*tmp
            gruh1s(2*lm  ) = gruh1s(2*lm  ) - d(ii,k)*tmp
!
            grv1s (2*lm-1) = grv1s (2*lm-1) + vz(ir,k)*alpn(lmr+n)
            grv1s (2*lm  ) = grv1s (2*lm  ) + vz(ii,k)*alpn(lmr+n)
!
            grvh1s(2*lm-1) = grvh1s(2*lm-1) - vz(ir,k)*tmp
            grvh1s(2*lm  ) = grvh1s(2*lm  ) - vz(ii,k)*tmp
         end do
      end do

      do lm=1,mlength
         m = locm(lm,iam)
         lmr = lnstart(lm)
         lmc = 2*lmr
         do n=2,nlen(m),2
            ir = lmc + 2*n - 1
            ii = ir + 1
!
            grus (2*lm-1,k) = grus (2*lm-1,k) + vz(ir,k)*dalpn(lmr+n)
            grus (2*lm  ,k) = grus (2*lm  ,k) + vz(ii,k)*dalpn(lmr+n)
!
            tmp = dalpn(lmr+n)*hdifzd(n+m-1,k)
            gruhs(2*lm-1,k) = gruhs(2*lm-1,k) - vz(ir,k)*tmp
            gruhs(2*lm  ,k) = gruhs(2*lm  ,k) - vz(ii,k)*tmp
!
            grvs (2*lm-1,k) = grvs (2*lm-1,k) - d(ir,k)*dalpn(lmr+n)
            grvs (2*lm  ,k) = grvs (2*lm  ,k) - d(ii,k)*dalpn(lmr+n)
!
            grvhs(2*lm-1,k) = grvhs(2*lm-1,k) + d(ir,k)*tmp
            grvhs(2*lm  ,k) = grvhs(2*lm  ,k) + d(ii,k)*tmp
         end do
      end do
!
! Combine the two parts of u(m) and v(m)
!
      do lm=1,mlength
         grus (2*lm-1,k) = grus (2*lm-1,k) + gru1s (2*lm  )
         gruhs(2*lm-1,k) = gruhs(2*lm-1,k) + gruh1s(2*lm  )
         grus (2*lm  ,k) = grus (2*lm  ,k) - gru1s (2*lm-1)
         gruhs(2*lm  ,k) = gruhs(2*lm  ,k) - gruh1s(2*lm-1)
         grvs (2*lm-1,k) = grvs (2*lm-1,k) + grv1s (2*lm  )
         grvhs(2*lm-1,k) = grvhs(2*lm-1,k) + grvh1s(2*lm  )
         grvs (2*lm  ,k) = grvs (2*lm  ,k) - grv1s (2*lm-1)
         grvhs(2*lm  ,k) = grvhs(2*lm  ,k) - grvh1s(2*lm-1)
      end do
!
! Remove Coriolis contribution to absolute vorticity from u(m)
! Correction for u:zeta=vz-ez=(zeta+f)-f
!
      if (lmwave0 .ne. -1) then
!        grus(1,k) = grus(1,k) - zurcor
         grus(2*lmwave0-1,k) = grus(2*lmwave0-1,k) - zurcor
      endif
   end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do n=1,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!         
         grpss (2*lm-1) = grpss (2*lm-1) + alps(ir)*lalp(lmr+n,irow)
         grpss (2*lm  ) = grpss (2*lm  ) + alps(ii)*lalp(lmr+n,irow)
!
         grdpss(2*lm-1) = grdpss(2*lm-1) + alps(ir)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
         grdpss(2*lm  ) = grdpss(2*lm  ) + alps(ii)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
      end do
   end do

   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do n=2,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!
         grpms(2*lm-1) = grpms(2*lm-1) + alps(ir)*ldalp(lmr+n,irow)*ra
         grpms(2*lm  ) = grpms(2*lm  ) + alps(ii)*ldalp(lmr+n,irow)*ra
      end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
      grpls(2*lm-1) = -grpss(2*lm  )*ra*xm(m)
      grpls(2*lm  ) =  grpss(2*lm-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalcs


subroutine grcalca (irow    ,ztodt   ,grta    ,grtha   ,grda    ,&
                    grza    ,grua    ,gruha   ,grva    ,grvha   ,&
                    grpsa   ,grdpsa  ,grpma   ,grpla   )

!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.6.6 2003/12/15 18:52:50 hender Exp $
! $Author: hender $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra
   use comhd

   implicit none

!
! Input arguments
!
   integer, intent(in) :: irow         ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
!
! Output arguments: antisymmetric fourier coefficients
!
   real(r8), intent(out) :: grta(2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grtha(2*maxm,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grda(2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grza(2*maxm,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grua(2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruha(2*maxm,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grva(2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvha(2*maxm,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpsa(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpsa(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpma(2*maxm)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpla(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1a(2*maxm)        ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1a(2*maxm)       ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1a(2*maxm)        ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1a(2*maxm)       ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) alpn(lpspt)          ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn(lpspt)         ! (a/(n(n+1)))*derivative of Legendre functions (complex)

   integer k                ! level index
   integer lm, m            ! local and global Fourier wavenumber indices of spectral array
   integer mlength          ! number of local wavenumbers
   integer n                ! meridional wavenumber index
   integer ir,ii            ! spectral indices
   integer lmr,lmc          ! spectral indices

   real(r8) tmp,raxm  ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
   mlength = numm(iam)
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      raxm = ra*xm(m)
      do n=1,nlen(m)
         alpn(lmr+n) = lalp(lmr+n,irow)*rsq(m+n-1)*raxm
         dalpn(lmr+n) = ldalp(lmr+n,irow)*rsq(m+n-1)*ra
      end do
   end do
!
! Initialize sums
!
   grza(:,:)  = 0.
   grda(:,:)  = 0.
   gruha(:,:) = 0.
   grvha(:,:) = 0.
   grtha(:,:) = 0.
   grpsa(:)   = 0.
   grua(:,:)  = 0.
   grva(:,:)  = 0.
   grta(:,:)  = 0.
   grpla(:)   = 0.
   grpma(:)   = 0.
   grdpsa(:)   = 0.

   do k=1,plev
      gru1a(:) = 0.
      gruh1a(:) = 0.
      grv1a(:) = 0.
      grvh1a(:) = 0.
!
! Loop over n for t,q,d,and end of u and v
!
      do lm=1,mlength
         m = locm(lm,iam)
         lmr = lnstart(lm)
         lmc = 2*lmr
         do n=1,nlen(m),2
            ir = lmc + 2*n - 1
            ii = ir + 1
!
            grua (2*lm-1,k) = grua (2*lm-1,k) + vz(ir,k)*dalpn(lmr+n)
            grua (2*lm  ,k) = grua (2*lm  ,k) + vz(ii,k)*dalpn(lmr+n)
!
            tmp = dalpn(lmr+n)*hdifzd(n+m-1,k)
            gruha(2*lm-1,k) = gruha(2*lm-1,k) - vz(ir,k)*tmp
            gruha(2*lm  ,k) = gruha(2*lm  ,k) - vz(ii,k)*tmp
!
            grva (2*lm-1,k) = grva (2*lm-1,k) - d(ir,k)*dalpn(lmr+n)
            grva (2*lm  ,k) = grva (2*lm  ,k) - d(ii,k)*dalpn(lmr+n)
!
            grvha(2*lm-1,k) = grvha(2*lm-1,k) + d(ir,k)*tmp
            grvha(2*lm  ,k) = grvha(2*lm  ,k) + d(ii,k)*tmp
         end do
      end do

      do lm=1,mlength
         m = locm(lm,iam)
         lmr = lnstart(lm)
         lmc = 2*lmr
         do n=2,nlen(m),2
            ir = lmc + 2*n - 1
            ii = ir + 1
            grta (2*lm-1,k) = grta (2*lm-1,k) + t(ir,k)*lalp(lmr+n,irow)
            grta (2*lm  ,k) = grta (2*lm  ,k) + t(ii,k)*lalp(lmr+n,irow)
!
            tmp = lalp(lmr+n,irow)*hdiftq(n+m-1,k)
            grtha(2*lm-1,k) = grtha(2*lm-1,k) - t(ir,k)*tmp
            grtha(2*lm  ,k) = grtha(2*lm  ,k) - t(ii,k)*tmp
!
            grda(2*lm-1,k) = grda(2*lm-1,k) + d(ir,k)*lalp(lmr+n,irow)
            grda(2*lm  ,k) = grda(2*lm  ,k) + d(ii,k)*lalp(lmr+n,irow)
!
            grza(2*lm-1,k) = grza(2*lm-1,k) + vz(ir,k)*lalp(lmr+n,irow)
            grza(2*lm  ,k) = grza(2*lm  ,k) + vz(ii,k)*lalp(lmr+n,irow)
!
            gru1a (2*lm-1) = gru1a (2*lm-1) + d(ir,k)*alpn(lmr+n)
            gru1a (2*lm  ) = gru1a (2*lm  ) + d(ii,k)*alpn(lmr+n)
!
            tmp = alpn(lmr+n)*hdifzd(n+m-1,k)
            gruh1a(2*lm-1) = gruh1a(2*lm-1) - d(ir,k)*tmp
            gruh1a(2*lm  ) = gruh1a(2*lm  ) - d(ii,k)*tmp
!
            grv1a (2*lm-1) = grv1a (2*lm-1) + vz(ir,k)*alpn(lmr+n)
            grv1a (2*lm  ) = grv1a (2*lm  ) + vz(ii,k)*alpn(lmr+n)
!
            grvh1a(2*lm-1) = grvh1a(2*lm-1) - vz(ir,k)*tmp
            grvh1a(2*lm  ) = grvh1a(2*lm  ) - vz(ii,k)*tmp
         end do
      end do
!
! Combine the two parts of u(m) and v(m)
!
      do lm=1,mlength
         grua (2*lm-1,k) = grua (2*lm-1,k) + gru1a (2*lm  )
         gruha(2*lm-1,k) = gruha(2*lm-1,k) + gruh1a(2*lm  )
         grua (2*lm  ,k) = grua (2*lm  ,k) - gru1a (2*lm-1)
         gruha(2*lm  ,k) = gruha(2*lm  ,k) - gruh1a(2*lm-1)
         grva (2*lm-1,k) = grva (2*lm-1,k) + grv1a (2*lm  )
         grvha(2*lm-1,k) = grvha(2*lm-1,k) + grvh1a(2*lm  )
         grva (2*lm  ,k) = grva (2*lm  ,k) - grv1a (2*lm-1)
         grvha(2*lm  ,k) = grvha(2*lm  ,k) - grvh1a(2*lm-1)
      end do
   end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do n=1,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1

         grpma(2*lm-1) = grpma(2*lm-1) + alps(ir)*ldalp(lmr+n,irow)*ra
         grpma(2*lm  ) = grpma(2*lm  ) + alps(ii)*ldalp(lmr+n,irow)*ra
      end do
   end do

   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do n=2,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!
         grpsa (2*lm-1) = grpsa (2*lm-1) + alps(ir)*lalp(lmr+n,irow)
         grpsa (2*lm  ) = grpsa (2*lm  ) + alps(ii)*lalp(lmr+n,irow)
!
         grdpsa(2*lm-1) = grdpsa(2*lm-1) + alps(ir)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
         grdpsa(2*lm  ) = grdpsa(2*lm  ) + alps(ii)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
      end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
      grpla(2*lm-1) = -grpsa(2*lm  )*ra*xm(m)
      grpla(2*lm  ) =  grpsa(2*lm-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalca

#endif
