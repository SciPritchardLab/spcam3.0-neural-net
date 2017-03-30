#include <misc.h>
#include <params.h>

module spetru

!----------------------------------------------------------------------- 
! 
! Purpose: Spectrally truncate initial data fields.
!
! Method: Truncate one or a few fields at a time, to minimize the 
!         memory  requirements
! 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified to implement processing of subsets of fields: P. Worley, May 2003
! 
!-----------------------------------------------------------------------

contains

subroutine spetru_phis(phis    ,phis_hires)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate PHIS input field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,   only: plon, plond, plev, plat
   use pspect
   use comspe
   use rgrid,    only: nlon, nmmax
   use commap,   only: w, xm
   use dynconst, only: rearth

   implicit none

#include <comctl.h>
#include <comfft.h>
!
! Input/Output arguments
!
   real(r8), intent(inout) :: phis(plond,plat)       ! Fourier -> spec. coeffs. for sfc geo.
   logical, intent(in) :: phis_hires          ! true => PHIS came from hi res topo file
!
!---------------------------Local workspace-----------------------------
!
   real(r8) phi(2,psp/2)       ! used in spectral truncation of phis
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) phialpr,phialpi    ! phi*alp (real and imaginary)
   real(r8) zwalp              ! zw*alp
   real(r8) zw                 ! w**2
   real(r8) filtlim            ! filter function
   real(r8) ft                 ! filter multiplier for spectral coefficients
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer i                   ! longitude index
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index (non-PVP)
!                                   along-diagonal index (PVP)
   integer n                   ! latitudinal wavenumber index (non-PVP)
!                                   diagonal index (PVP)
   integer nspec
#if ( defined PVP )              
   integer ic                  ! complex coeff. index
   integer ialp                ! index into legendre polynomials
#else
   integer mr                  ! spectral indices
#endif
!
!-----------------------------------------------------------------------
!
! Zero spectral array
!
   phi(:,:) = 0.
!
! Transform grid -> fourier
!
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      call fft991(phis(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,-1)
   end do                    ! lat=1,plat
!
! Loop over latitude pairs
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
      zw = w(irow)*2.
      do i=1,2*nmmax(irow)
!
! Compute symmetric and antisymmetric components
!
         tmp1 = 0.5*(phis(i,latm) - phis(i,latp))
         tmp2 = 0.5*(phis(i,latm) + phis(i,latp))
         phis(i,latm) = tmp1
         phis(i,latp) = tmp2
      end do
!     
! Compute phi*mn
!
#if ( defined PVP )
      do n=1,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            zwalp = zw*alp(ialp+m,irow)
            phi(1,ic+m) = phi(1,ic+m) + zwalp*phis(2*m-1,latp)
            phi(2,ic+m) = phi(2,ic+m) + zwalp*phis(2*m  ,latp)
         end do
      end do
!
      do n=2,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            zwalp = zw*alp(ialp+m,irow)
            phi(1,ic+m) = phi(1,ic+m) + zwalp*phis(2*m-1,latm)
            phi(2,ic+m) = phi(2,ic+m) + zwalp*phis(2*m  ,latm)
         end do
      end do
#else
      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=1,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            phi(1,mr+n) = phi(1,mr+n) + zwalp*phis(2*m-1,latp)
            phi(2,mr+n) = phi(2,mr+n) + zwalp*phis(2*m  ,latp)
         end do

         do n=2,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            phi(1,mr+n) = phi(1,mr+n) + zwalp*phis(2*m-1,latm)
            phi(2,mr+n) = phi(2,mr+n) + zwalp*phis(2*m  ,latm)
         end do
      end do
#endif
   enddo                  ! irow=1,plat/2
!
   if (phis_hires) then
!
! Apply spectral filter to phis
!     filter is a function of n 
!        if n < filter limit then
!           spectral_coeff = spectral_coeff * (1. - (float(n)/filtlim)**2)
!        else         
!           spectral_coeff = 0.
!        endif
!     where filter limit = 1.4*PTRN
!     
      filtlim = float(int(1.4*float(ptrn)))
#if ( defined PVP )
      do n=1,pmax
         ic = ncoefi(n) - 1
         do m=1,nm(n)
            nspec=m-1+n
            ft = 1. - (float(nspec)/filtlim)**2
            if (float(nspec) .ge. filtlim) ft = 0.
            phi(1,ic+m) = phi(1,ic+m)*ft
            phi(2,ic+m) = phi(2,ic+m)*ft
         end do
      end do
#else   
      do m=1,pmmax
         mr = nstart(m)
         do n=1,nlen(m)
            nspec=m-1+n
            ft = 1. - (float(nspec)/filtlim)**2
            if (float(nspec) .ge. filtlim) ft = 0.
            phi(1,mr+n) = phi(1,mr+n)*ft 
            phi(2,mr+n) = phi(2,mr+n)*ft 
         end do
      end do
#endif   
      call hordif1(rearth,phi)
   end if
!
! Compute grid point values of phi*.
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
!
! Zero fourier fields
!
      phis(:,latm) = 0.
      phis(:,latp) = 0.
!
! Compute(phi*)m
!
#if ( defined PVP )
      do n=1,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            phialpr = phi(1,ic+m)*alp(ialp+m,irow)
            phialpi = phi(2,ic+m)*alp(ialp+m,irow)
!
            phis(2*m-1,latm) = phis(2*m-1,latm) + phialpr
            phis(2*m  ,latm) = phis(2*m  ,latm) + phialpi
         end do
      end do
!
      do n=2,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            phialpr = phi(1,ic+m)*alp(ialp+m,irow)
            phialpi = phi(2,ic+m)*alp(ialp+m,irow)
!
            phis(2*m-1,latp) = phis(2*m-1,latp) + phialpr
            phis(2*m  ,latp) = phis(2*m  ,latp) + phialpi
!
         end do
      end do
#else
      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=1,nlen(m),2
            phialpr = phi(1,mr+n)*alp(mr+n,irow)
            phialpi = phi(2,mr+n)*alp(mr+n,irow)
!     
            phis(2*m-1,latm) = phis(2*m-1,latm) + phialpr
            phis(2*m  ,latm) = phis(2*m  ,latm) + phialpi
         end do
      end do

      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=2,nlen(m),2
            phialpr = phi(1,mr+n)*alp(mr+n,irow)
            phialpi = phi(2,mr+n)*alp(mr+n,irow)
!     
            phis(2*m-1,latp) = phis(2*m-1,latp) + phialpr
            phis(2*m  ,latp) = phis(2*m  ,latp) + phialpi
!
         end do
      end do
#endif
!
! Recompute real fields from symmetric and antisymmetric parts
!
      do i=1,nlon(latm)+2
         tmp1 = phis(i,latm) + phis(i,latp)
         tmp2 = phis(i,latm) - phis(i,latp)
         phis(i,latm) = tmp1
         phis(i,latp) = tmp2
      end do

   enddo                 ! irow=1,plat/2
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(phis(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,+1)
   enddo

   return
end subroutine spetru_phis

!************************************************************************
subroutine spetru_ps(ps      ,dpsl    ,dpsm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate PS input field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,   only: plon, plond, plev, plat
   use pspect
   use comspe
   use rgrid,    only: nlon, nmmax
   use commap,   only: w, xm
   use dynconst, only: ra

   implicit none

#include <comctl.h>
#include <comfft.h>
!
! Input/Output arguments
!
   real(r8), intent(inout) :: ps(plond,plat)         ! Fourier -> spec. coeffs. for ln(ps)
!
! Output arguments
!
   real(r8), intent(out) :: dpsl(plond,plat)         ! Spectrally trunc d(ln(ps))/d(longitude)
   real(r8), intent(out) :: dpsm(plond,plat)         ! Spectrally trunc d(ln(ps))/d(latitude)

!
!---------------------------Local workspace-----------------------------
!
   real(r8) alps_tmp(psp)      ! used in spectral truncation of phis
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) zwalp              ! zw*alp
   real(r8) psdalpr,psdalpi    ! alps (real and imaginary)*dalp
   real(r8) psalpr,psalpi      ! alps (real and imaginary)*alp
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index (non-PVP)
!                                   along-diagonal index (PVP)
   integer n                   ! latitudinal wavenumber index (non-PVP)
!                                   diagonal index (PVP)
   integer nspec
#if ( defined PVP )              
   integer ic                  ! complex coeff. index
   integer ialp                ! index into legendre polynomials
#else
   integer mr,mc               ! spectral indices
#endif
!
!-----------------------------------------------------------------------
!
! Zero spectral array
!
   alps_tmp(:) = 0.
!
! Compute the 2D quantities which are transformed to spectral space:
!  ps= ln(ps). 
!
   do lat=1,plat
     irow = lat
     if (lat.gt.plat/2) irow = plat - lat + 1
     do i=1,nlon(lat)
         ps(i,lat) = log(ps(i,lat))
      end do
!
! Transform grid -> fourier
!
      call fft991(ps(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,-1)

   end do                    ! lat=1,plat
!
! Loop over latitude pairs
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
      zw = w(irow)*2.
      do i=1,2*nmmax(irow)
!
! Compute symmetric and antisymmetric components
!
         tmp1 = 0.5*(ps(i,latm) - ps(i,latp))
         tmp2 = 0.5*(ps(i,latm) + ps(i,latp))
         ps(i,latm) = tmp1
         ps(i,latp) = tmp2

      end do
!     
! Compute ln(p*)mn
!
#if ( defined PVP )
      do n=1,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            zwalp = zw*alp(ialp+m,irow)
            ir = 2*(ic+m) - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*ps(2*m-1,latp)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*ps(2*m  ,latp)
         end do
      end do
!
      do n=2,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            zwalp = zw*alp(ialp+m,irow)
            ir = 2*(ic+m) - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*ps(2*m-1,latm)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*ps(2*m  ,latm)
         end do
      end do
#else
      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=1,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            ir = mc + 2*n - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*ps(2*m-1,latp)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*ps(2*m  ,latp)
         end do

         do n=2,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            ir = mc + 2*n - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*ps(2*m-1,latm)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*ps(2*m  ,latm)
         end do
      end do
#endif
   enddo                  ! irow=1,plat/2
!
! Compute grid point values of:ln(p*) and grad(ln(p*)).
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
!
! Zero fourier fields
!
      ps(:,latm) = 0.
      ps(:,latp) = 0.

      dpsl(:,latm) = 0.
      dpsl(:,latp) = 0.

      dpsm(:,latm) = 0.
      dpsm(:,latp) = 0.

!
! Compute(ln(p*),grad(ln(p*)))m
!
#if ( defined PVP )
      do n=1,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            ir = 2*(ic+m) - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(ialp+m,irow)
            psalpi = alps_tmp(ii)*alp(ialp+m,irow)
!
            ps(2*m-1,latm) = ps(2*m-1,latm) + psalpr
            ps(2*m  ,latm) = ps(2*m  ,latm) + psalpi
            dpsl(2*m-1,latm) = dpsl(2*m-1,latm) - psalpi*ra
            dpsl(2*m  ,latm) = dpsl(2*m  ,latm) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(ialp+m,irow)
            psdalpi = alps_tmp(ii)*dalp(ialp+m,irow)
!
            dpsm(2*m-1,latp) = dpsm(2*m-1,latp) + psdalpr*ra
            dpsm(2*m  ,latp) = dpsm(2*m  ,latp) + psdalpi*ra
         end do
      end do
!
      do n=2,pmax,2
         ic = ncoefi(n) - 1
         ialp = nalp(n)
         do m=1,nmreduced(n,irow)
            ir = 2*(ic+m) - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(ialp+m,irow)
            psalpi = alps_tmp(ii)*alp(ialp+m,irow)
!
            ps(2*m-1,latp) = ps(2*m-1,latp) + psalpr
            ps(2*m  ,latp) = ps(2*m  ,latp) + psalpi
            dpsl(2*m-1,latp) = dpsl(2*m-1,latp) - psalpi*ra
            dpsl(2*m  ,latp) = dpsl(2*m  ,latp) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(ialp+m,irow)
            psdalpi = alps_tmp(ii)*dalp(ialp+m,irow)
!     
            dpsm(2*m-1,latm) = dpsm(2*m-1,latm) + psdalpr*ra
            dpsm(2*m  ,latm) = dpsm(2*m  ,latm) + psdalpi*ra
         end do
      end do
!
#else
      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=1,nlen(m),2
            ir = mc + 2*n - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(mr+n,irow)
            psalpi = alps_tmp(ii)*alp(mr+n,irow)
!     
            ps(2*m-1,latm) = ps(2*m-1,latm) + psalpr
            ps(2*m  ,latm) = ps(2*m  ,latm) + psalpi
            dpsl(2*m-1,latm) = dpsl(2*m-1,latm) - psalpi*ra
            dpsl(2*m  ,latm) = dpsl(2*m  ,latm) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(mr+n,irow)
            psdalpi = alps_tmp(ii)*dalp(mr+n,irow)
!
            dpsm(2*m-1,latp) = dpsm(2*m-1,latp) + psdalpr*ra
            dpsm(2*m  ,latp) = dpsm(2*m  ,latp) + psdalpi*ra
         end do
      end do

      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=2,nlen(m),2
            ir = mc + 2*n - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(mr+n,irow)
            psalpi = alps_tmp(ii)*alp(mr+n,irow)
!     
            ps(2*m-1,latp) = ps(2*m-1,latp) + psalpr
            ps(2*m  ,latp) = ps(2*m  ,latp) + psalpi
            dpsl(2*m-1,latp) = dpsl(2*m-1,latp) - psalpi*ra
            dpsl(2*m  ,latp) = dpsl(2*m  ,latp) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(mr+n,irow)
            psdalpi = alps_tmp(ii)*dalp(mr+n,irow)
!
            dpsm(2*m-1,latm) = dpsm(2*m-1,latm) + psdalpr*ra
            dpsm(2*m  ,latm) = dpsm(2*m  ,latm) + psdalpi*ra
         end do
      end do

#endif
      do m=1,nmmax(irow)
         dpsl(2*m-1,latm) = xm(m)*dpsl(2*m-1,latm)
         dpsl(2*m  ,latm) = xm(m)*dpsl(2*m  ,latm)
         dpsl(2*m-1,latp) = xm(m)*dpsl(2*m-1,latp)
         dpsl(2*m  ,latp) = xm(m)*dpsl(2*m  ,latp)
      end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
      do i=1,nlon(latm)+2
!
         tmp1 = ps(i,latm) + ps(i,latp)
         tmp2 = ps(i,latm) - ps(i,latp)
         ps(i,latm) = tmp1
         ps(i,latp) = tmp2
!
         tmp1 = dpsl(i,latm) + dpsl(i,latp)
         tmp2 = dpsl(i,latm) - dpsl(i,latp)
         dpsl(i,latm) = tmp1
         dpsl(i,latp) = tmp2
!
         tmp1 = dpsm(i,latm) + dpsm(i,latp)
         tmp2 = dpsm(i,latm) - dpsm(i,latp)
         dpsm(i,latm) = tmp1
         dpsm(i,latp) = tmp2
      end do
!
   enddo                 ! irow=1,plat/2
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(ps(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,+1)
      call fft991(dpsl(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,+1)
      call fft991(dpsm(1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),1,+1)
!
! Convert from ln(ps) to ps
!
      do i=1,nlon(lat)
         ps(i,lat) = exp(ps(i,lat))
      end do
!
   enddo

   return
end subroutine spetru_ps

!************************************************************************

subroutine spetru_t(t3)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate T input field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,   only: plon, plond, plev, plat
   use pspect
   use comspe
   use rgrid,    only: nlon, nmmax
   use commap,   only: w, xm

   implicit none

#include <comctl.h>
#include <comfft.h>
!
! Input/Output arguments
!
   real(r8), intent(inout) :: t3(plond,plev,plat)    ! Fourier -> spec. coeffs. for temperature
!
!---------------------------Local workspace-----------------------------
!
   real(r8) t_tmp(psp)         ! used in spectral truncation of t
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) tmpr               ! vector temporary (real)
   real(r8) tmpi               ! vector temporary (imaginary)
   real(r8) zwalp              ! zw*alp
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index (non-PVP)
!                                   along-diagonal index (PVP)
   integer n                   ! latitudinal wavenumber index (non-PVP)
!                                   diagonal index (PVP)
   integer nspec
#if ( defined PVP )              
   integer ic                  ! complex coeff. index
   integer ialp                ! index into legendre polynomials
#else
   integer mr,mc               ! spectral indices
#endif
!
!-----------------------------------------------------------------------
!
! Transform grid -> fourier
!
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      call fft991(t3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,-1)
   end do                    ! lat=1,plat
!
! Loop over vertical levels
!
   do k=1,plev
!
! Zero spectral array
!
      t_tmp(:) = 0.
!
! Loop over latitude pairs
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
         zw = w(irow)*2.
!
! Multi-level field: T
!
         do i=1,2*nmmax(irow)
            tmp1 = 0.5*(t3(i,k,latm) - t3(i,k,latp))
            tmp2 = 0.5*(t3(i,k,latm) + t3(i,k,latp))
            t3(i,k,latm) = tmp1
            t3(i,k,latp) = tmp2
         end do
!     
! Compute tmn
!
#if ( defined PVP )
         do n=1,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
            do m=1,nmreduced(n,irow)
               zwalp  = zw*alp (ialp+m,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*t3(2*m-1,k,latp) 
               t_tmp(ii) = t_tmp(ii) + zwalp*t3(2*m  ,k,latp)
            end do
         end do
!
         do n=2,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
            do m=1,nmreduced(n,irow)
               zwalp  = zw*alp (ialp+m,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*t3(2*m-1,k,latm)
               t_tmp(ii) = t_tmp(ii) + zwalp*t3(2*m  ,k,latm)
            end do
         end do
#else
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*t3(2*m-1,k,latp)
               t_tmp(ii) = t_tmp(ii) + zwalp*t3(2*m  ,k,latp)
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*t3(2*m-1,k,latm)
               t_tmp(ii) = t_tmp(ii) + zwalp*t3(2*m  ,k,latm)
            end do
         end do
#endif
      enddo                ! irow=1,plat/2
!
! Compute grid point values of:t.
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
!
! Zero fourier fields
!
         t3(:,k,latm) = 0.
         t3(:,k,latp) = 0.
#if ( defined PVP )
         do n=1,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
!DIR$ IVDEP
            do m=1,nmreduced(n,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
!
               tmpr = t_tmp(ir)*alp(ialp+m,irow)
               tmpi = t_tmp(ii)*alp(ialp+m,irow)
               t3(2*m-1,k,latm) = t3(2*m-1,k,latm) + tmpr
               t3(2*m  ,k,latm) = t3(2*m  ,k,latm) + tmpi
            end do
         end do
!
         do n=2,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
!DIR$ IVDEP
            do m=1,nmreduced(n,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
!
               tmpr = t_tmp(ir)*alp(ialp+m,irow)
               tmpi = t_tmp(ii)*alp(ialp+m,irow)
               t3(2*m-1,k,latp) = t3(2*m-1,k,latp) + tmpr
               t3(2*m  ,k,latp) = t3(2*m  ,k,latp) + tmpi
            end do
         end do
#else
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!
               tmpr = t_tmp(ir)*alp(mr+n,irow)
               tmpi = t_tmp(ii)*alp(mr+n,irow)
               t3(2*m-1,k,latm) = t3(2*m-1,k,latm) + tmpr
               t3(2*m  ,k,latm) = t3(2*m  ,k,latm) + tmpi
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!
               tmpr = t_tmp(ir)*alp(mr+n,irow)
               tmpi = t_tmp(ii)*alp(mr+n,irow)
               t3(2*m-1,k,latp) = t3(2*m-1,k,latp) + tmpr
               t3(2*m  ,k,latp) = t3(2*m  ,k,latp) + tmpi
            end do
         end do
#endif
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = t3(i,k,latm) + t3(i,k,latp)
            tmp2 = t3(i,k,latm) - t3(i,k,latp)
            t3(i,k,latm) = tmp1
            t3(i,k,latp) = tmp2
         end do
      enddo                ! irow=1,plat/2
   enddo                   ! k=1,plev
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(t3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,+1)
   enddo

   return
end subroutine spetru_t

!***********************************************************************

subroutine spetru_uv(u3      ,v3      ,vort    ,div)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate U, V input fields.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,   only: plon, plond, plev, plat
   use pspect
   use comspe
   use rgrid,    only: nlon, nmmax
   use commap,   only: w, xm, rsq, cs
   use dynconst, only: ez, ra

   implicit none

#include <comctl.h>
#include <comfft.h>
!
! Input/Output arguments
!
   real(r8), intent(inout) :: u3(plond,plev,plat)    ! Fourier -> spec. coeffs. for u-wind
   real(r8), intent(inout) :: v3(plond,plev,plat)    ! Fourier -> spec. coeffs. for v-wind
!
! Output arguments
!
   real(r8), intent(out) :: vort(plond,plev,plat)    ! Spectrally truncated vorticity
   real(r8), intent(out) :: div(plond,plev,plat)     ! Spectrally truncated divergence

!
!---------------------------Local workspace-----------------------------
!
   real(r8) d_tmp(psp)         ! used in spectral truncation of div
   real(r8) vz_tmp(psp)        ! used in spectral truncation of vort
   real(r8) alpn(pspt)         ! alp*rsq*xm*ra
   real(r8) dalpn(pspt)        ! dalp*rsq*ra
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) tmpr               ! vector temporary (real)
   real(r8) tmpi               ! vector temporary (imaginary)
   real(r8) zcor               ! correction for absolute vorticity
   real(r8) zwalp              ! zw*alp
   real(r8) zwdalp             ! zw*dalp
   real(r8) zrcsj              ! ra/(cos**2 latitude)
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif
   real(r8) zsqcs

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index (non-PVP)
!                                   along-diagonal index (PVP)
   integer n                   ! latitudinal wavenumber index (non-PVP)
!                                   diagonal index (PVP)
   integer nspec
#if ( defined PVP )              
   integer ne                  ! index into rsq
   integer ic                  ! complex coeff. index
   integer ialp                ! index into legendre polynomials
#else
   integer mr,mc               ! spectral indices
#endif

   call print_memusage ('spetru_uv')
!
!-----------------------------------------------------------------------
!
! Compute the quantities which are transformed to spectral space:
!   1. u = u*sqrt(1-mu*mu),   u * cos(phi)
!   2. v = v*sqrt(1-mu*mu),   v * cos(phi)
!
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      zsqcs = sqrt(cs(irow))
      do k=1,plev
         do i=1,nlon(lat)
            u3(i,k,lat) = u3(i,k,lat)*zsqcs
            v3(i,k,lat) = v3(i,k,lat)*zsqcs
         end do
      end do
!
! Transform grid -> fourier
! 1st transform: U,V,T: note contiguity assumptions
! 2nd transform: LN(PS).  3rd transform: surface geopotential
!
      call fft991(u3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,-1)
      call fft991(v3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,-1)

   end do                    ! lat=1,plat
!
! Multi-level fields: U, V
!
   do k=1,plev
!
! Zero spectral arrays
!
      vz_tmp(:) = 0.
      d_tmp(:) = 0.
!
! Loop over latitude pairs
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
         zrcsj = ra/cs(irow)
         zw = w(irow)*2.
         do i=1,2*nmmax(irow)

            tmp1 = 0.5*(u3(i,k,latm) - u3(i,k,latp))
            tmp2 = 0.5*(u3(i,k,latm) + u3(i,k,latp))
            u3(i,k,latm) = tmp1
            u3(i,k,latp) = tmp2

            tmp1 = 0.5*(v3(i,k,latm) - v3(i,k,latp))
            tmp2 = 0.5*(v3(i,k,latm) + v3(i,k,latp))
            v3(i,k,latm) = tmp1
            v3(i,k,latp) = tmp2

         end do
!     
! Compute vzmn and dmn
!
#if ( defined PVP )
         do n=1,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
            do m=1,nmreduced(n,irow)
               zwdalp = zw*dalp(ialp+m,irow)
               zwalp  = zw*alp (ialp+m,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v3(2*m-1,k,latm) + &
                  xm(m)*zwalp*u3(2*m  ,k,latp))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v3(2*m  ,k,latm) - &
                  xm(m)*zwalp*u3(2*m-1,k,latp))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u3(2*m-1,k,latm) - &
                  xm(m)*zwalp*v3(2*m  ,k,latp))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u3(2*m  ,k,latm) + &
                  xm(m)*zwalp*v3(2*m-1,k,latp))*zrcsj
            end do
         end do
!
         do n=2,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
            do m=1,nmreduced(n,irow)
               zwdalp = zw*dalp(ialp+m,irow)
               zwalp  = zw*alp (ialp+m,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v3(2*m-1,k,latp) + &
                  xm(m)*zwalp*u3(2*m  ,k,latm))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v3(2*m  ,k,latp) - &
                  xm(m)*zwalp*u3(2*m-1,k,latm))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u3(2*m-1,k,latp) - &
                  xm(m)*zwalp*v3(2*m  ,k,latm))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u3(2*m  ,k,latp) + &
                  xm(m)*zwalp*v3(2*m-1,k,latm))*zrcsj
            end do
         end do
#else
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               zwdalp = zw*dalp(mr+n,irow)
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v3(2*m-1,k,latm) + &
                  xm(m)*zwalp*u3(2*m  ,k,latp))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v3(2*m  ,k,latm) - &
                  xm(m)*zwalp*u3(2*m-1,k,latp))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u3(2*m-1,k,latm) - &
                  xm(m)*zwalp*v3(2*m  ,k,latp))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u3(2*m  ,k,latm) + &
                  xm(m)*zwalp*v3(2*m-1,k,latp))*zrcsj
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               zwdalp = zw*dalp(mr+n,irow)
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v3(2*m-1,k,latp) + &
                  xm(m)*zwalp*u3(2*m  ,k,latm))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v3(2*m  ,k,latp) - &
                  xm(m)*zwalp*u3(2*m-1,k,latm))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u3(2*m-1,k,latp) - &
                  xm(m)*zwalp*v3(2*m  ,k,latm))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u3(2*m  ,k,latp) + &
                  xm(m)*zwalp*v3(2*m-1,k,latm))*zrcsj
            end do
         end do
#endif
      enddo               ! irow=1,plat/2
!
! Compute grid point values of:u,v,vz, and d.
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
#if ( defined PVP )
         zcor = ez*alp(nalp(2)+1,irow)
#else
         zcor = ez*alp(2,irow)
#endif
!
! Compute(u,v,vz,d)m
!
#if ( defined PVP )
         do n=1,pmax
            ne = n - 1
            ialp = nalp(n)
            do m=1,nmreduced(n,irow)
               alpn (ialp+m) =  alp(ialp+m,irow)*rsq(ne+m)*xm(m)*ra
               dalpn(ialp+m) = dalp(ialp+m,irow)*rsq(ne+m)      *ra
            end do
         end do
#else
         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=1,nlen(m)
!
! These statements will likely not be bfb since xm*ra is now a scalar
!
               alpn (mr+n) =  alp(mr+n,irow)*rsq(n+m-1)*xm(m)*ra
               dalpn(mr+n) = dalp(mr+n,irow)*rsq(n+m-1)      *ra
            end do
         end do
#endif
!
! Zero fourier fields
!
         u3(:,k,latm) = 0.
         u3(:,k,latp) = 0.

         v3(:,k,latm) = 0.
         v3(:,k,latp) = 0.

         vort(:,k,latm) = 0.
         vort(:,k,latp) = 0.

         div(:,k,latm) = 0.
         div(:,k,latp) = 0.

#if ( defined PVP )
         do n=1,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
!DIR$ IVDEP
            do m=1,nmreduced(n,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
!
               tmpr = d_tmp(ir)*alpn(ialp+m)
               tmpi = d_tmp(ii)*alpn(ialp+m)
               u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpi
               u3(2*m  ,k,latm) = u3(2*m  ,k,latm) - tmpr
!
               tmpr = d_tmp(ir)*dalpn(ialp+m)
               tmpi = d_tmp(ii)*dalpn(ialp+m)
               v3(2*m-1,k,latp) = v3(2*m-1,k,latp) - tmpr
               v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(ialp+m)
               tmpi = vz_tmp(ii)*dalpn(ialp+m)
               u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpr
               u3(2*m  ,k,latp) = u3(2*m  ,k,latp) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(ialp+m)
               tmpi = vz_tmp(ii)*alpn(ialp+m)
               v3(2*m-1,k,latm) = v3(2*m-1,k,latm) + tmpi
               v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpr
!
               tmpr = d_tmp(ir)*alp(ialp+m,irow)
               tmpi = d_tmp(ii)*alp(ialp+m,irow)
               div(2*m-1,k,latm) = div(2*m-1,k,latm) + tmpr
               div(2*m  ,k,latm) = div(2*m  ,k,latm) + tmpi
!
               tmpr = vz_tmp(ir)*alp(ialp+m,irow)
               tmpi = vz_tmp(ii)*alp(ialp+m,irow)
               vort(2*m-1,k,latm) = vort(2*m-1,k,latm) + tmpr
               vort(2*m  ,k,latm) = vort(2*m  ,k,latm) + tmpi
            end do
         end do
!
         do n=2,pmax,2
            ic = ncoefi(n) - 1
            ialp = nalp(n)
!DIR$ IVDEP
            do m=1,nmreduced(n,irow)
               ir = 2*(ic+m) - 1
               ii = ir + 1
!
               tmpr = d_tmp(ir)*alpn(ialp+m)
               tmpi = d_tmp(ii)*alpn(ialp+m)
               u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpi
               u3(2*m  ,k,latp) = u3(2*m  ,k,latp) - tmpr
!     
               tmpr = d_tmp(ir)*dalpn(ialp+m)
               tmpi = d_tmp(ii)*dalpn(ialp+m)
               v3(2*m-1,k,latm) = v3(2*m-1,k,latm) - tmpr
               v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(ialp+m)
               tmpi = vz_tmp(ii)*dalpn(ialp+m)
               u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpr
               u3(2*m  ,k,latm) = u3(2*m  ,k,latm) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(ialp+m)
               tmpi = vz_tmp(ii)*alpn(ialp+m)
               v3(2*m-1,k,latp) = v3(2*m-1,k,latp) + tmpi
               v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpr
!
               tmpr = d_tmp(ir)*alp(ialp+m,irow)
               tmpi = d_tmp(ii)*alp(ialp+m,irow)
               div(2*m-1,k,latp) = div(2*m-1,k,latp) + tmpr
               div(2*m  ,k,latp) = div(2*m  ,k,latp) + tmpi
!
               tmpr = vz_tmp(ir)*alp(ialp+m,irow)
               tmpi = vz_tmp(ii)*alp(ialp+m,irow)
               vort(2*m-1,k,latp) = vort(2*m-1,k,latp) + tmpr
               vort(2*m  ,k,latp) = vort(2*m  ,k,latp) + tmpi
            end do
         end do
!
! Correction to get the absolute vorticity.
!
         vort(1,k,latp) = vort(1,k,latp) + zcor
#else
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!
               tmpr = d_tmp(ir)*alpn(mr+n)
               tmpi = d_tmp(ii)*alpn(mr+n)
               u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpi
               u3(2*m  ,k,latm) = u3(2*m  ,k,latm) - tmpr
!
               tmpr = d_tmp(ir)*dalpn(mr+n)
               tmpi = d_tmp(ii)*dalpn(mr+n)
               v3(2*m-1,k,latp) = v3(2*m-1,k,latp) - tmpr
               v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(mr+n)
               tmpi = vz_tmp(ii)*dalpn(mr+n)
               u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpr
               u3(2*m  ,k,latp) = u3(2*m  ,k,latp) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(mr+n)
               tmpi = vz_tmp(ii)*alpn(mr+n)
               v3(2*m-1,k,latm) = v3(2*m-1,k,latm) + tmpi
               v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpr
!
               tmpr = d_tmp(ir)*alp(mr+n,irow)
               tmpi = d_tmp(ii)*alp(mr+n,irow)
               div(2*m-1,k,latm) = div(2*m-1,k,latm) + tmpr
               div(2*m  ,k,latm) = div(2*m  ,k,latm) + tmpi
!
               tmpr = vz_tmp(ir)*alp(mr+n,irow)
               tmpi = vz_tmp(ii)*alp(mr+n,irow)
               vort(2*m-1,k,latm) = vort(2*m-1,k,latm) + tmpr
               vort(2*m  ,k,latm) = vort(2*m  ,k,latm) + tmpi
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!     
               tmpr = d_tmp(ir)*alpn(mr+n)
               tmpi = d_tmp(ii)*alpn(mr+n)
               u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpi
               u3(2*m  ,k,latp) = u3(2*m  ,k,latp) - tmpr
!
               tmpr = d_tmp(ir)*dalpn(mr+n)
               tmpi = d_tmp(ii)*dalpn(mr+n)
               v3(2*m-1,k,latm) = v3(2*m-1,k,latm) - tmpr
               v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(mr+n)
               tmpi = vz_tmp(ii)*dalpn(mr+n)
               u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpr
               u3(2*m  ,k,latm) = u3(2*m  ,k,latm) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(mr+n)
               tmpi = vz_tmp(ii)*alpn(mr+n)
               v3(2*m-1,k,latp) = v3(2*m-1,k,latp) + tmpi
               v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpr
!
               tmpr = d_tmp(ir)*alp(mr+n,irow)
               tmpi = d_tmp(ii)*alp(mr+n,irow)
               div(2*m-1,k,latp) = div(2*m-1,k,latp) + tmpr
               div(2*m  ,k,latp) = div(2*m  ,k,latp) + tmpi
!
               tmpr = vz_tmp(ir)*alp(mr+n,irow)
               tmpi = vz_tmp(ii)*alp(mr+n,irow)
               vort(2*m-1,k,latp) = vort(2*m-1,k,latp) + tmpr
               vort(2*m  ,k,latp) = vort(2*m  ,k,latp) + tmpi
            end do
         end do
!
! Correction to get the absolute vorticity.
!     
         vort(1,k,latp) = vort(1,k,latp) + zcor
#endif
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = u3(i,k,latm) + u3(i,k,latp)
            tmp2 = u3(i,k,latm) - u3(i,k,latp)
            u3(i,k,latm) = tmp1
            u3(i,k,latp) = tmp2
!
            tmp1 = v3(i,k,latm) + v3(i,k,latp)
            tmp2 = v3(i,k,latm) - v3(i,k,latp)
            v3(i,k,latm) = tmp1
            v3(i,k,latp) = tmp2
!
            tmp1 = vort(i,k,latm) + vort(i,k,latp)
            tmp2 = vort(i,k,latm) - vort(i,k,latp)
            vort(i,k,latm) = tmp1
            vort(i,k,latp) = tmp2
!
            tmp1 = div(i,k,latm) + div(i,k,latp)
            tmp2 = div(i,k,latm) - div(i,k,latp)
            div(i,k,latm) = tmp1
            div(i,k,latp) = tmp2
         end do
      enddo               ! irow=1,plat/2
   enddo                  ! k=1,plev
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(u3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,+1)
      call fft991(v3(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,+1)
      call fft991(vort(1,1,lat),work,trig(1,irow),ifax(1,irow),1, &
                  plond,nlon(lat),plev,+1)
      call fft991(div(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plond, &
                  nlon(lat),plev,+1)
!
! Convert U,V to u,v
!
      zsqcs = sqrt(cs(irow))
      do k=1,plev
         do i=1,nlon(lat)
            u3(i,k,lat) = u3(i,k,lat)/zsqcs
            v3(i,k,lat) = v3(i,k,lat)/zsqcs
         end do
      end do
   enddo

   return
end subroutine spetru_uv

end module spetru
