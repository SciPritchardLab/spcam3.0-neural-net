#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  Make sure to make appropriate coding changes where
! necessary.
#if ( defined PVP )
subroutine dyndrv(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
                  grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
                  grz2    ,grd2    ,grfu2   ,grfv2   ,vmax2d  , &
                  vmax2dt ,vcour   )
!-----------------------------------------------------------------------
!
! Purpose:
! Driving routine for Gaussian quadrature, semi-implicit equation
! solution and linear part of horizontal diffusion.
! The need for this interface routine is to have a multitasking
! driver for the spectral space routines it invokes.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: dyndrv.F90,v 1.5.2.3 2003/02/05 21:53:24 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe, only: maxm
  use commap
  use time_manager, only: get_step_size

  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: grlps1(2*maxm,plat/2)      ! ----------------------------
  real(r8), intent(in)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grz1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grd1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grfu1 (2*maxm,plev,plat/2) ! | see quad for definitions:
  real(r8), intent(in)   :: grfv1 (2*maxm,plev,plat/2) ! | these variables are
  real(r8), intent(in)   :: grlps2(2*maxm,plat/2)      ! | declared here for data
  real(r8), intent(in)   :: grt2  (2*maxm,plev,plat/2) ! | scoping
  real(r8), intent(in)   :: grq2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grz2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grd2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(in)   :: grfv2 (2*maxm,plev,plat/2) ! ----------------------------
  real(r8), intent(inout)   :: vmax2d (plev,plat)      ! max. wind at each level, lat
  real(r8), intent(inout)   :: vmax2dt(plev,plat)      ! max. truncated wind at each lvl, lat
  real(r8), intent(inout)   :: vcour  (plev,plat)      ! maximum Courant number in slice
!
!---------------------------Local workspace-----------------------------
!
  real(r8) ztdtsq(2*pnmax)  ! dt*(n(n+1)/a^2)
  real(r8) zdt              ! dt
  real(r8) ztdt             ! dt
  integer irow              ! latitude pair index
  integer n                 ! meridional wavenumber index
  integer k                 ! level index
!
!-----------------------------------------------------------------------
!
  call t_startf('dyn')

!$OMP PARALLEL DO PRIVATE (IROW)

  do irow = 1,plat/2
     call dyn(irow    ,grlps1(1,irow),grt1(1,1,irow),grq1(1,1,irow),grz1(1,1,irow)   , &
              grd1(1,1,irow),grfu1(1,1,irow),grfv1(1,1,irow),grlps2(1,irow), &
                 grt2(1,1,irow) , &
              grq2(1,1,irow),grz2 (1,1,irow),grd2 (1,1,irow),grfu2(1,1,irow), &
                 grfv2(1,1,irow)  )
  end do
  call t_stopf('dyn')
!
!-----------------------------------------------------------------------
!
! Build expanded vector with del^2 response function
!
  zdt  = get_step_size()*0.5
  ztdt = 2.*zdt
!DIR$ IVDEP
  do n=1,pnmax
     ztdtsq(2*n-1) = ztdt*sq(n)
     ztdtsq(2*n  ) = ztdt*sq(n)
  end do
!
! Perform Gaussian quadrature (multitasked loop)
!
  call t_startf('quad_tstep')

!$OMP PARALLEL DO PRIVATE (N)

  do n=1,pmax
     call quad(n       ,grlps1  ,grlps2  ,grt1    ,grq1    , &
               grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!
! Transform spectral terms into vertical normal modes
!
#ifdef HADVTEST
!
!jr Turn off semi-implicit so T equation is horizontal advection only
!        call tstep(n       ,ztdtsq  )
#else
     call tstep(n       ,ztdtsq)
  end do
  call t_stopf('quad_tstep')
!
! Solve vertically de-coupled semi-implicit divergence/vorticity system
! in normal mode space
!
  call t_startf('vertnm')
!$OMP PARALLEL DO PRIVATE(K)
  do k=1,plev
     call vertnm(k)
  end do
  call t_stopf('vertnm')
!
! Transform divergence and vorticity back from normal mode space and
! complete time advance.
!
  call t_startf('tstep1')
!$OMP PARALLEL DO PRIVATE(N)
  do n=1,pmax
     call tstep1(n       ,zdt     )
  end do
#endif
  call t_stopf('tstep1')
!
! Find out if courant limit has been exceeded.  If so, the limiter will
! be applied in HORDIF
!
  call t_startf('courlim')
  call courlim(vmax2d  ,vmax2dt ,vcour   )
  call t_stopf('courlim')
!
! Linear part of horizontal diffusion (multitasked loop)
!
  call t_startf('hordif')
!$OMP PARALLEL DO PRIVATE(K)
  do k=1,plev
     call hordif(k       ,ztdt    )
  end do
  call t_stopf('hordif')
!
  return
end subroutine dyndrv

#else

subroutine dyndrv(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
                  grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
                  grz2    ,grd2    ,grfu2   ,grfv2   ,vmax2d  , &
                  vmax2dt ,vcour   )
!-----------------------------------------------------------------------
!
! Purpose:
! Driving routine for Gaussian quadrature, semi-implicit equation
! solution and linear part of horizontal diffusion.
! The need for this interface routine is to have a multitasking
! driver for the spectral space routines it invokes.
!
! Author:  J. Rosinski
! Modified: P. Worley, October 2002
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use commap
  use time_manager, only: get_step_size
#if ( defined SPMD ) && ( defined TIMING_BARRIERS )
  use mpishorthand
#endif

  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(inout)   :: grlps1(2*maxm,plat/2)      ! ----------------------------
  real(r8), intent(inout)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grz1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grd1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfu1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfv1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grlps2(2*maxm,plat/2)      ! | definitions: these
  real(r8), intent(inout)   :: grt2  (2*maxm,plev,plat/2) ! | variables are declared here
  real(r8), intent(inout)   :: grq2  (2*maxm,plev,plat/2) ! | for data scoping
  real(r8), intent(inout)   :: grz2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grd2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfv2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: vmax2d (plev,plat)         ! max. wind at each level, latitude
  real(r8), intent(inout)   :: vmax2dt(plev,plat)         ! max. truncated wind at each lvl,lat
  real(r8), intent(inout)   :: vcour  (plev,plat)         ! maximum Courant number in slice
!
!---------------------------Local workspace-----------------------------
!
  real(r8) ztdtsq(pnmax)    ! dt*(n(n+1)/a^2)
  real(r8) zdt              ! dt
  real(r8) ztdt             ! dt
  integer irow              ! latitude pair index
  integer lm                ! local longitudinal wavenumber index
  integer n                 ! spectral index
  integer k                 ! level index
!
!-----------------------------------------------------------------------
!
  call t_startf('dyn')

!$OMP PARALLEL DO PRIVATE (IROW)

  do irow=1,plat/2
     call dyn(irow    ,grlps1(1,irow),grt1(1,1,irow),grq1(1,1,irow),grz1(1,1,irow), &
              grd1(1,1,irow),grfu1(1,1,irow),grfv1(1,1,irow),grlps2(1,irow)  , &
                 grt2(1,1,irow), &
              grq2(1,1,irow),grz2 (1,1,irow),grd2 (1,1,irow),grfu2 (1,1,irow), &
                 grfv2(1,1,irow)  )
  end do

  call t_stopf('dyn')
!
!-----------------------------------------------------------------------
!
! Build vector with del^2 response function
!
  zdt  = get_step_size()*0.5
  ztdt = 2.*zdt
  do n=1,pnmax
     ztdtsq(n) = ztdt*sq(n)
  end do

  call t_startf('quad')

!$OMP PARALLEL DO PRIVATE(LM)

  do lm=1,numm(iam)
!
! Perform Gaussian quadrature
!
     call quad(lm      ,grlps1  ,grlps2  ,grt1    ,grq1    , &
               grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
   end do

   call t_stopf('quad')
   call t_startf('tstep')

!$OMP PARALLEL DO PRIVATE(LM)

   do lm=1,numm(iam)
!
! Complete time advance, solve vertically coupled semi-implicit system
!
#ifdef HADVTEST
!
!jr Turn off semi-implicit so T equation is horizontal advection only
!        call tstep(lm      ,ztdtsq  )
#else
     call tstep(lm      ,ztdtsq  )
!
! Solve vertically de-coupled semi-implicit divergence/vorticity system
! in normal mode space
!
     call vertnm(lm)
!
! Transform divergence and vorticity back from normal mode space and
! complete time advance.
!
     call tstep1(lm      ,zdt     )
#endif
  end do
  call t_stopf('tstep')
!
! Find out if courant limit has been exceeded.  If so, the limiter will
! be applied in HORDIF
!
  call t_startf('courlim')
  call courlim(vmax2d  ,vmax2dt ,vcour   )
  call t_stopf('courlim')
!
! Linear part of horizontal diffusion
!
  call t_startf('hordif')

!$OMP PARALLEL DO PRIVATE(K)

  do k=1,plev
     call hordif(k       ,ztdt    )
  end do
  call t_stopf('hordif')

  return
end subroutine dyndrv
#endif
