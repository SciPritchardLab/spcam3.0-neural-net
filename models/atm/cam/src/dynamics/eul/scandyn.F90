#include <misc.h>
#include <params.h>

subroutine scandyn (ztodt,   etadot,  etamid,  grlps1,  grt1,   &
                    grz1,    grd1,    grfu1,   grfv1,   grut1,  &
                    grvt1,   grrh1,   grlps2,  grt2,    grz2,   &
                    grd2,    grfu2,   grfv2,   grut2,   grvt2,  &
                    grrh2,   vcour,   vmax2d,  vmax2dt, detam,  &
                    cwava,   flx_net, t2,      fu,      fv)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! "After coupling" gaussian latitude scan for which some of the physics
! and nonlinear dynamics calculations are completed.  The main loop over
! latitude in this routine is multitasked.
!
! Note: the "ifdef" constructs in this routine are associated with the
! message-passing version of CAM.  Messages are sent  which
! have no relevance to the shared-memory case.  
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------
!
! $Id: scandyn.F90,v 1.8.6.3 2003/06/13 15:50:49 hender Exp $
! $Author: hender $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use prognostics
   use rgrid
   use comspe, only: maxm
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                       ! two delta t unless nstep =0
   real(r8), intent(inout) :: etadot(plon,plevp,beglat:endlat)     ! vertical motion (slt)
   real(r8), intent(in) :: etamid(plev)                ! hybrd coord value at levels
   real(r8), intent(in) :: detam(plev)      
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSDYN and and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
   real(r8), intent(in)    :: cwava(plat)           ! weight applied to global integrals
   real(r8), intent(in)    :: flx_net(plond,beglat:endlat)         ! net flx from physics
   real(r8), intent(inout) :: t2(plond,plev,beglat:endlat)         ! tot dT/dt to to physics
   real(r8), intent(inout) :: fu(plond,plev,beglat:endlat)         ! u wind tend
   real(r8), intent(inout) :: fv(plond,plev,beglat:endlat)         ! v wind tend
!
! Output arguments
!
   real(r8), intent(out) :: grlps1(2*maxm,plat/2)      ! sym. undiff. term in lnps eqn.
   real(r8), intent(out) :: grlps2(2*maxm,plat/2)      ! antisym undiff. term in lnps eqn.
   real(r8), intent(out) :: grt1(2*maxm,plev,plat/2)   ! sym. undiff. term in t eqn.
   real(r8), intent(out) :: grt2(2*maxm,plev,plat/2)   ! antisym. undiff. term in t eqn.
   real(r8), intent(out) :: grz1(2*maxm,plev,plat/2)   ! sym. undiff. term in z eqn.
   real(r8), intent(out) :: grz2(2*maxm,plev,plat/2)   ! antisym. undiff. term in z eqn.
   real(r8), intent(out) :: grd1(2*maxm,plev,plat/2)   ! sym. undiff. term in d eqn.
   real(r8), intent(out) :: grd2(2*maxm,plev,plat/2)   ! antisym. undiff. term in d eqn.
   real(r8), intent(out) :: grfu1(2*maxm,plev,plat/2)  ! sym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfu2(2*maxm,plev,plat/2)  ! antisym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfv1(2*maxm,plev,plat/2)  ! sym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grfv2(2*maxm,plev,plat/2)  ! antisym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grut1(2*maxm,plev,plat/2)  ! sym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grut2(2*maxm,plev,plat/2)  ! antisym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grvt1(2*maxm,plev,plat/2)  ! sym. mu derivative term in t eqn.
   real(r8), intent(out) :: grvt2(2*maxm,plev,plat/2)  ! antisym. mu deriv. term in t eqn.
   real(r8), intent(out) :: grrh1(2*maxm,plev,plat/2)  ! sym. del**2 term in d eqn.
   real(r8), intent(out) :: grrh2(2*maxm,plev,plat/2)  ! antisym. del**2 term in d eqn.
   real(r8), intent(out) :: vcour(plev,plat)           ! maximum Courant number in vert.
   real(r8), intent(out) :: vmax2d(plev,plat)          ! max. wind at each level, latitude
   real(r8), intent(out) :: vmax2dt(plev,plat)         ! max. truncated wind at each lvl,lat

! Local variables

   integer irow              ! latitude pair index
   integer lat,j,latn,lats   ! latitude indices
!
! FFT buffers
!
   real(r8), allocatable:: fftbuf_in(:,:,:,:)          ! fftbuf_in(plond,plev,9,beglat:endlat) 
   real(r8), allocatable:: fftbuf_out(:,:,:,:)         ! fftbuf_out(2*maxm,plev,9,plat)
!
   allocate(fftbuf_in(plond,plev,9,beglat:endlat))
#if ( defined SPMD )
   allocate(fftbuf_out(2*maxm,plev,9,plat))
#else
   allocate(fftbuf_out(1,1,1,1))
#endif
!
!$OMP PARALLEL DO PRIVATE (J, LAT)
   do lat=beglat,endlat

      j = j1 - 1 + lat
      call linemsdyn_bft (lat, nlon(lat), ps(1,lat,n3m1), ps(1,lat,n3m2), u3(i1,1,j,n3m1), &
                      u3(i1,1,j,n3m2), v3(i1,1,j,n3m1), v3(i1,1,j,n3m2), t3(i1,1,j,n3m1), t3(i1,1,j,n3m2), &
                      q3(i1,1,1,j,n3m1), etadot(1,1,lat), etamid, &
                      ztodt, vcour(1,lat), vmax2d(1,lat), vmax2dt(1,lat),       &
                      detam, t2(1,1,lat), fu(1,1,lat), fv(1,1,lat),                     &
                      div(1,1,lat,n3m1), vort(1,1,lat,n3m2), div(1,1,lat,n3m2), vort(1,1,lat,n3m1), &
                      phis(1,lat), dpsl(1,lat), dpsm(1,lat), omga(1,1,lat), &
                      cwava(lat), flx_net(1,lat), fftbuf_in(1,1,1,lat) )
   end do

   call linemsdyn_fft (fftbuf_in,fftbuf_out)

!$OMP PARALLEL DO PRIVATE (IROW, LATN, LATS)
   do irow=1,plat/2

      lats = irow
      latn = plat - irow + 1
#if ( defined SPMD )
      call linemsdyn_aft (irow, fftbuf_out(1,1,1,lats), fftbuf_out(1,1,1,latn), &
                      grlps1(1,irow), grt1(1,1,irow), grz1(1,1,irow), grd1(1,1,irow), &
                      grfu1(1,1,irow),  grfv1(1,1,irow),   &
                      grut1(1,1,irow), grvt1(1,1,irow), grrh1(1,1,irow), grlps2(1,irow),grt2(1,1,irow),    &
                      grz2(1,1,irow), grd2(1,1,irow), grfu2(1,1,irow), grfv2(1,1,irow),  grut2(1,1,irow),  &
                      grvt2(1,1,irow), grrh2(1,1,irow) )
#else
      call linemsdyn_aft (irow, fftbuf_in(1,1,1,lats), fftbuf_in(1,1,1,latn), &
                      grlps1(1,irow), grt1(1,1,irow), grz1(1,1,irow), grd1(1,1,irow), &
                      grfu1(1,1,irow),  grfv1(1,1,irow),   &
                      grut1(1,1,irow), grvt1(1,1,irow), grrh1(1,1,irow), grlps2(1,irow),grt2(1,1,irow),    &
                      grz2(1,1,irow), grd2(1,1,irow), grfu2(1,1,irow), grfv2(1,1,irow),  grut2(1,1,irow),  &
                      grvt2(1,1,irow), grrh2(1,1,irow) )
#endif
   end do
!
   deallocate(fftbuf_in)
   deallocate(fftbuf_out)

   return
end subroutine scandyn

