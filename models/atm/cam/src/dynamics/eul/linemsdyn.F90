#include <misc.h>
#include <params.h>

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Control non-linear dynamical terms, FFT and combine terms
! in preparation for Fourier -> spectral quadrature.
! 
! Method: 
! The naming convention is as follows:
!  - prefix gr contains grid point values before FFT and Fourier
!     coefficients after
!  - t, q, d, z and ps refer to temperature, specific humidity,
!     divergence, vorticity and surface pressure
!  - "1" suffix to an array => symmetric component current latitude pair
!  - "2" suffix to an array => antisymmetric component.
!
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------

subroutine linemsdyn_bft(                                         &
                     lat     ,nlon    ,psm1    ,psm2    ,u3m1    , &
                     u3m2    ,v3m1    ,v3m2    ,t3m1    ,t3m2    , &
                     q3m1    ,etadot  ,etamid  ,                   &
                     ztodt   , vcour   ,vmax   ,vmaxt   ,          &
                     detam   ,t2      ,fu      ,fv      ,          &
                     divm1   ,vortm2  ,divm2   ,vortm1  ,phis    , &
                     dpsl    ,dpsm    ,omga    ,cwava   ,flx_net , &
                     fftbuf             )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Control non-linear dynamical terms and fill FFT buffer 
! in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
!
!-----------------------------------------------------------------------
!
! $Id: linemsdyn.F90,v 1.12.2.4 2003/06/13 15:50:41 hender Exp $
! $Author: hender $

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use constituents, only: pcnst, pnats
   use pspect
   use comslt
   use commap
   use history, only: outfld
   use time_manager, only: get_step_size

   implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>
!
! Input arguments
!     
   integer lat               ! latitude index for S->N storage
   integer nlon

   real(r8), intent(in) :: psm1(plond)        ! surface pressure (time n)
   real(r8), intent(in) :: psm2(plond)      ! surface pressure (time n-1)
   real(r8), intent(in) :: u3m1(plond,plev)   ! u-wind (time n)
   real(r8), intent(in) :: u3m2(plond,plev) ! u-wind (time n-1)
   real(r8), intent(in) :: v3m1(plond,plev)   ! v-wind (time n)
   real(r8), intent(in) :: v3m2(plond,plev) ! v-wind (time n-1)
   real(r8), intent(in) :: t3m1(plond,plev)   ! temperature (time n)
   real(r8), intent(in) :: q3m1(plond,plev,pcnst+pnats)   ! constituent conc(time n: h2o first)
   real(r8), intent(inout) :: etadot(plon,plevp) ! vertical motion (3-d used by slt)
   real(r8), intent(in) :: etamid(plev)     ! midpoint values of eta (a+b)
   real(r8), intent(in) :: ztodt            ! 2*timestep unless nstep = 0
   real(r8), intent(in) :: detam(plev)      ! maximum Courant number in vert.
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: t2(plond,plev)   ! t tend
   real(r8), intent(inout) :: fu(plond,plev)   ! nonlinear term - u momentum eqn.
   real(r8), intent(inout) :: fv(plond,plev)   ! nonlinear term - v momentum eqn.
   real(r8), intent(inout) :: divm1(plond,plev)
   real(r8), intent(inout) :: vortm2(plond,plev)
   real(r8), intent(inout) :: divm2(plond,plev)
   real(r8), intent(inout) :: vortm1(plond,plev)
   real(r8), intent(inout) :: phis(plond)
   real(r8), intent(inout) :: dpsl(plond)
   real(r8), intent(inout) :: dpsm(plond)
   real(r8), intent(inout) :: omga(plond,plev)
   real(r8), intent(inout) :: t3m2(plond,plev) ! temperature (time n-1)
   real(r8), intent(in)    :: cwava                  ! weight for global water vapor int.
   real(r8), intent(in)    :: flx_net(plond)         ! net flux from physics
!     
! Output arguments
!     
   real(r8), intent(out) :: fftbuf(plond,plev,9) ! buffer used for in-place FFTs
   real(r8), intent(out) :: vcour(plev)      ! maximum Courant number in vert.
   real(r8), intent(out) :: vmax(plev)       ! maximum wind speed squared (m^2/s^2)
   real(r8), intent(out) :: vmaxt(plev)      ! maximum truncated wind speed (m^2/s^2)
!     
!---------------------------Local workspace-----------------------------
!     
   real(r8) :: dtime          ! timestep size
   real(r8) :: bpstr(plond)   ! 
   real(r8) pmid(plond,plev)  ! pressure at model levels (time n)
   real(r8) rpmid(plond,plev) ! 1./pmid
   real(r8) pint(plond,plevp) ! pressure at model interfaces (n  )
   real(r8) pdel(plond,plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) rpdel(plond,plev) ! 1./pdel
   real(r8) tdyn(plond,plev)   ! temperature for dynamics
   real(r8) logpsm1(plond)     ! log(psm1)
   real(r8) logpsm2(plond)   ! log(psm2)
   real(r8) engy(plond,plev) ! kinetic energy
   real(r8) ut(plond,plev)   ! (u*T) - heat flux - zonal
   real(r8) vt(plond,plev)   ! (v*T) - heat flux - meridional
   real(r8) drhs(plond,plev) ! RHS of divergence eqn. (del^2 term)
   real(r8) lvcour           ! local vertical courant number
   real(r8) dtdz             ! dt/detam(k)
   real(r8) ddivdt(plond,plev) ! temporary workspace
   real(r8) ddpn(plond)      ! complete sum of d*delta p
   real(r8) vpdsn(plond)     ! complete sum V dot grad(ln(ps)) delta b
   real(r8) dpslat(plond,plev) ! Pressure gradient term 
   real(r8) dpslon(plond,plev) ! Pressure gradient term 
   real(r8) coslat           ! cosine(latitude)
   real(r8) rcoslat          ! 1./cosine(latitude)
   real(r8) rhypi            ! 1./hypi(plevp)

   real(r8) wind             ! u**2 + v**2 (m/s)
   real(r8) utfac            ! asymmetric truncation factor for courant calculation
   real(r8) vtfac            ! asymmetric truncation factor for courant calculation

   real(r8) tmp                  ! accumulator
   integer i,k,kk            ! longitude,level,constituent indices
   integer, parameter :: tdyndex = 1     ! indices into fftbuf 
   integer, parameter :: fudex = 2
   integer, parameter :: fvdex = 3
   integer, parameter :: utdex = 4
   integer, parameter :: vtdex = 5
   integer, parameter :: drhsdex = 6
   integer, parameter :: vortdyndex = 7
   integer, parameter :: divdyndex = 8
   integer, parameter :: bpstrdex = 9
!
! This group of arrays are glued together via equivalence to exbuf for
! communication from LINEMSBC.
!
!
!-----------------------------------------------------------------------
!
!
! Compute maximum wind speed this latitude (used in Courant number estimate)
!
   if (ptrm .lt. ptrn) then
      utfac = float(ptrm)/float(ptrn)
      vtfac = 1.
   else if (ptrn .lt. ptrm) then
      utfac = 1.
      vtfac = float(ptrn)/float(ptrm) 
   else if (ptrn .eq. ptrm) then
      utfac = 1.
      vtfac = 1.
   end if
   do k=1,plev
      vmax(k) = 0.
      vmaxt(k) = 0.
      do i=1,nlon
         wind = u3m2(i,k)**2 + v3m2(i,k)**2
         vmax(k) = max(wind,vmax(k))
!
! Change to Courant limiter for non-triangular truncations.
!
         wind = utfac*u3m2(i,k)**2 + vtfac*v3m2(i,k)**2
         vmaxt(k) = max(wind,vmaxt(k))
      end do
   end do
!
! Variables needed in tphysac
!
   coslat = cos(clat(lat))
   rcoslat = 1./coslat
!
! Set current time pressure arrays for model levels etc.
!
   call plevs0(nlon,plond,plev,psm1,pint,pmid,pdel)
   do k=1,plev
      do i=1,nlon
         rpmid(i,k) = 1./pmid(i,k)
         rpdel(i,k) = 1./pdel(i,k)
      end do
   end do
!
! Accumulate statistics for diagnostic print
!
   call stats(lat,     pint,    pdel,      psm1,   &
              vortm1,  divm1,   t3m1,      q3m1(:,:,1), nlon  )
!
! Compute log(surface pressure) for use by grmult and when adding tendency.
!
   do i=1,nlon
      logpsm1(i) = log(psm1(i))
      logpsm2(i) = log(psm2(i))
   end do
!     
! Compute integrals
!     
   call plevs0(nlon,plond,plev,psm2,pint,pmid,pdel)
   call engy_te (cwava,w(lat),t3m2,u3m2,v3m2,phis    ,pdel, tmp  ,nlon)
   engy1lat(lat) = tmp 
   call plevs0(nlon,plond,plev,psm1,pint,pmid,pdel)
!
! Include top/bottom flux integral to energy integral
!
   call flxint  (w(lat) ,flx_net ,tmp  ,nlon )
   engy1lat(lat) = engy1lat(lat) + tmp *ztodt
!
! Calculate non-linear terms in tendencies
!
   if (adiabatic) t2(:,:) = 0.
   call grmult(rcoslat ,divm1     ,q3m1(1,1,1),t3m1   ,u3m1    , &
               v3m1    ,vortm1    ,t3m2    ,phis    ,dpsl    , &
               dpsm    ,omga    ,pdel    ,pint(1,plevp),logpsm2, &
               logpsm1 ,rpmid   ,rpdel   ,fu      ,fv      , &
               t2      ,ut      ,vt      ,drhs    ,pmid    , &
               etadot  ,etamid  ,engy    ,ddpn    ,vpdsn   , &
               dpslon  ,dpslat  ,nlon    )
!
! Add tendencies to previous timestep values of surface pressure,
! temperature, and (if spectral transport) moisture.  Store *log* surface
! pressure in bpstr array for transform to spectral space.
!
   do i=1,nlon
      bpstr(i) = logpsm2(i) - ztodt*(vpdsn(i)+ddpn(i))/psm1(i)
   end do

   rhypi = 1./hypi(plevp)
   do k=1,plev
      do i=1,nlon
         ddivdt(i,k) = ztodt*(0.5*divm2(i,k) - divm1(i,k))
         bpstr(i) = bpstr(i) - ddivdt(i,k)*hypd(k)*rhypi
         tdyn(i,k) = t3m2(i,k) + ztodt*t2(i,k)
      end do
   end do

   do k=1,plev
      do kk=1,plev
         do i=1,nlon
#ifdef HADVTEST
!
!jr Remove semi-implicit contribution for advection test
!
!jr              tdyn(i,k) = tdyn(i,k) - ddivdt(i,kk)*tau(kk,k)
#else
            tdyn(i,k) = tdyn(i,k) - ddivdt(i,kk)*tau(kk,k)
#endif
         end do
      end do
   end do

!
! Compute maximum vertical Courant number this latitude.
!
   dtime = get_step_size()
   vcour(:) = 0.
   do k=2,plev
      dtdz = dtime/detam(k-1)
      do i=1,nlon
         lvcour = abs(etadot(i,k))*dtdz
         vcour(k) = max(lvcour,vcour(k))
      end do
   end do

   call outfld('ETADOT  ',etadot,plon,lat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Apply cos(lat) to momentum terms before fft
!
   do k=1,plev
      do i=1,nlon
         fu(i,k) = coslat*fu(i,k)
         fv(i,k) = coslat*fv(i,k)
         ut(i,k) = coslat*ut(i,k)
         vt(i,k) = coslat*vt(i,k)
      end do
   end do

!
! Copy fields into FFT buffer
!
   do k=1,plev
      do i=1,nlon
!
! undifferentiated terms
         fftbuf(i,k,tdyndex) = tdyn(i,k)
! longitudinally and latitudinally differentiated terms
         fftbuf(i,k,fudex)   = fu(i,k)
         fftbuf(i,k,fvdex)   = fv(i,k)
         fftbuf(i,k,utdex)   = ut(i,k)
         fftbuf(i,k,vtdex)   = vt(i,k)
         fftbuf(i,k,drhsdex) = drhs(i,k)
! vort,div
         fftbuf(i,k,vortdyndex) = vortm2(i,k)
         fftbuf(i,k,divdyndex)  = divm2(i,k)
!
      enddo
   enddo
! ps
   do i=1,nlon
      fftbuf(i,1,bpstrdex) = bpstr(i)
   enddo

   return
end subroutine linemsdyn_bft


subroutine linemsdyn_fft(fftbuf,fftbuf2)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute FFT of non-linear dynamical terms
! in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, September 2002
!
!-----------------------------------------------------------------------
!
! $Id: linemsdyn.F90,v 1.12.2.4 2003/06/13 15:50:41 hender Exp $
! $Author: hender $

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use rgrid
#if (defined SPMD)
# if ( defined TIMING_BARRIERS )
   use mpishorthand, only: mpicom
# endif
   use comspe
#endif

   implicit none

#include <comfft.h>
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: fftbuf(plond,plev,9,beglat:endlat) 
                            ! buffer used for in-place FFTs
!     
! Output arguments
!     
#if (defined SPMD)
   real(r8), intent(out) :: fftbuf2(2*maxm,plev,9,plat) 
                            ! buffer for returning reorderd Fourier coefficients
#else
   real(r8), intent(in) :: fftbuf2(1) 
                            ! buffer unused
#endif
!     
!---------------------------Local workspace-----------------------------
!     
! The "work" array has a different size requirement depending upon whether
! the proprietary Cray assembly language version of the FFT library
! routines, or the all-Fortran version, is being used.
!     
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev*9)
#else 
   real(r8) work((plon+1)*pcray) ! workspace array for fft991
#endif
   integer lat               ! latitude index
   integer inc               ! increment for fft991
   integer isign             ! flag indicates transform direction
   integer ntr               ! number of transforms to perform
   integer ifld, k, i
!
   inc = 1
   isign = -1
   ntr = 8*plev + 1
!$OMP PARALLEL DO PRIVATE (LAT,WORK)
   do lat=beglat,endlat
      call fft991(fftbuf(1,1,1,lat)     ,work    ,trig(1,lat),ifax(1,lat),inc     ,&
                  plond   ,nlon(lat)    ,ntr     ,isign   )
   enddo
!
#if ( defined SPMD )
!
!  reorder Fourier coefficients
!
#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_realloc4a')
   call mpibarrier (mpicom)
   call t_stopf ('sync_realloc4a')
#endif
   call t_startf('realloc4a')
   call realloc4a(fftbuf, fftbuf2)
   call t_stopf('realloc4a')
#endif

   return
end subroutine linemsdyn_fft

subroutine linemsdyn_aft(                                          &
                     irow    ,fftbufs ,fftbufn ,                   &
                     grlps1  ,grt1    ,grz1    ,grd1    ,          &
                     grfu1   ,grfv1   ,grut1   ,grvt1   ,grrh1   , &
                     grlps2  ,grt2    ,grz2    ,grd2    ,grfu2   , &
                     grfv2   ,grut2   ,grvt2   ,grrh2              )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Combine terms in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, September 2002
!
!-----------------------------------------------------------------------
!
! $Id: linemsdyn.F90,v 1.12.2.4 2003/06/13 15:50:41 hender Exp $
! $Author: hender $

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
#if (defined SPMD)
   use comspe, only: numm, maxm
#else
   use comspe, only: maxm
   use rgrid, only: nmmax
#endif

   implicit none
!
! Input arguments
!     
   integer, intent(in) :: irow                ! latitude pair index

#if (defined SPMD)
   real(r8), intent(in) :: fftbufs(2*maxm,plev,9) ! southern latitude Fourier coefficients
   real(r8), intent(in) :: fftbufn(2*maxm,plev,9) ! northern latitude Fourier coefficients
#else
   real(r8), intent(in) :: fftbufs(plond,plev,9) ! southern latitude Fourier coefficients
   real(r8), intent(in) :: fftbufn(plond,plev,9) ! northern latitude Fourier coefficients
#endif
!     
! Output arguments
!     
   real(r8), intent(out) :: grlps1(2*maxm)  ! sym. undiff. term in lnps eqn.
   real(r8), intent(out) :: grlps2(2*maxm)  ! antisym undiff. term in lnps eqn.
   real(r8), intent(out) :: grt1(2*maxm,plev) ! sym. undiff. term in t eqn.
   real(r8), intent(out) :: grt2(2*maxm,plev) ! antisym. undiff. term in t eqn.
   real(r8), intent(out) :: grz1(2*maxm,plev) ! sym. undiff. term in z eqn.
   real(r8), intent(out) :: grz2(2*maxm,plev) ! antisym. undiff. term in z eqn.
   real(r8), intent(out) :: grd1(2*maxm,plev) ! sym. undiff. term in d eqn.
   real(r8), intent(out) :: grd2(2*maxm,plev) ! antisym. undiff. term in d eqn.
   real(r8), intent(out) :: grfu1(2*maxm,plev) ! sym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfu2(2*maxm,plev) ! antisym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfv1(2*maxm,plev) ! sym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grfv2(2*maxm,plev) ! antisym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grut1(2*maxm,plev) ! sym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grut2(2*maxm,plev) ! antisym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grvt1(2*maxm,plev) ! sym. mu derivative term in t eqn.
   real(r8), intent(out) :: grvt2(2*maxm,plev) ! antisym. mu deriv. term in t eqn.
   real(r8), intent(out) :: grrh1(2*maxm,plev) ! sym. del**2 term in d eqn.
   real(r8), intent(out) :: grrh2(2*maxm,plev) ! antisym. del**2 term in d eqn.
!     
!---------------------------Local workspace-----------------------------
!     
   integer i,k            ! longitude,level indices
   integer mlength        ! number of wavenumbers
   integer, parameter :: tdyndex = 1     ! indices into fftbuf 
   integer, parameter :: fudex = 2
   integer, parameter :: fvdex = 3
   integer, parameter :: utdex = 4
   integer, parameter :: vtdex = 5
   integer, parameter :: drhsdex = 6
   integer, parameter :: vortdyndex = 7
   integer, parameter :: divdyndex = 8
   integer, parameter :: bpstrdex = 9
!
#if (defined SPMD)
   mlength = numm(iam)
#else
   mlength = nmmax(irow)
#endif
   do k=1,plev
      do i=1,2*mlength

         grt1(i,k) = 0.5*(fftbufn(i,k,tdyndex)+fftbufs(i,k,tdyndex))
         grt2(i,k) = 0.5*(fftbufn(i,k,tdyndex)-fftbufs(i,k,tdyndex))

         grz1(i,k) = 0.5*(fftbufn(i,k,vortdyndex)+fftbufs(i,k,vortdyndex))
         grz2(i,k) = 0.5*(fftbufn(i,k,vortdyndex)-fftbufs(i,k,vortdyndex))

         grd1(i,k) = 0.5*(fftbufn(i,k,divdyndex)+fftbufs(i,k,divdyndex))
         grd2(i,k) = 0.5*(fftbufn(i,k,divdyndex)-fftbufs(i,k,divdyndex))

         grfu1(i,k) = 0.5*(fftbufn(i,k,fudex)+fftbufs(i,k,fudex))
         grfu2(i,k) = 0.5*(fftbufn(i,k,fudex)-fftbufs(i,k,fudex))

         grfv1(i,k) = 0.5*(fftbufn(i,k,fvdex)+fftbufs(i,k,fvdex))
         grfv2(i,k) = 0.5*(fftbufn(i,k,fvdex)-fftbufs(i,k,fvdex))

         grut1(i,k) = 0.5*(fftbufn(i,k,utdex)+fftbufs(i,k,utdex))
         grut2(i,k) = 0.5*(fftbufn(i,k,utdex)-fftbufs(i,k,utdex))

         grvt1(i,k) = 0.5*(fftbufn(i,k,vtdex)+fftbufs(i,k,vtdex))
         grvt2(i,k) = 0.5*(fftbufn(i,k,vtdex)-fftbufs(i,k,vtdex))

         grrh1(i,k) = 0.5*(fftbufn(i,k,drhsdex)+fftbufs(i,k,drhsdex))
         grrh2(i,k) = 0.5*(fftbufn(i,k,drhsdex)-fftbufs(i,k,drhsdex))

      end do
   end do

   do i=1,2*mlength
      grlps1(i) = 0.5*(fftbufn(i,1,bpstrdex)+fftbufs(i,1,bpstrdex))
      grlps2(i) = 0.5*(fftbufn(i,1,bpstrdex)-fftbufs(i,1,bpstrdex))
   end do

   return
end subroutine linemsdyn_aft

