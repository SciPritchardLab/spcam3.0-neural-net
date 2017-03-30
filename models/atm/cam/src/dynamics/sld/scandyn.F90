#include <misc.h>
#include <params.h>

subroutine scandyn (ztodt   ,pmap    ,kdpmpf  ,kdpmph  ,lam     , &
                    phi     ,dphi    ,sinlam  ,coslam  ,lbasdy  , &
                    lbasdz  ,lbasiy  ,lbasiz  ,lbassi  ,detam   , &
                    detai   ,dlam    ,cwava   ,etamid  ,etaint  , &
                    grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                    grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                    grfu    ,grfv    ,t2      ,flx_net , &
                    vcour   ,vmax2d  ,vmax2dt )
!-----------------------------------------------------------------------
!
! Purpose:
! Driving routine for semi-lagrangian transport and SLD dynamics.
! Set up  FFT and combine terms in preparation for Fourier -> spectral
! quadrature.
! 
! The latitude loop in this routine is multitasked.
! The naming convention is as follows:
!  - prefix gr contains grid point values before FFT and Fourier
!    coefficients after
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: scandyn.F90,v 1.14.2.4 2003/02/05 21:53:34 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use prognostics
  use comslt
  use rgrid
  use commap
  use physconst, only:
  use comspe, only: maxm
#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none

#include <comfft.h>
#include <comhyb.h>

!-----------------------------------------------------------------------
!
  integer, parameter :: iter=1  ! number of iterations to be used in dep pt calc.
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                     ! twice the timestep unless nstep=0
  integer , intent(in)   :: pmap                      ! dimension of artificial array
                                                      ! used to locate vertical interval
                                                      ! in which departure point falls
  integer , intent(in)   :: kdpmpf (pmap)             ! mapping from artificial array to
                                                      ! model levels
  integer , intent(in)   :: kdpmph (pmap)             ! mapping from artificial array to
                                                      ! model interfaces
  real(r8), intent(in)   :: lam    (plond,platd)      ! long. coordinates of model grid
  real(r8), intent(in)   :: phi    (platd)            ! lat.  coordinates of model grid
  real(r8), intent(in)   :: dphi   (platd)            ! latitudinal grid increments
  real(r8), intent(in)   :: sinlam (plond,platd)      ! sine of longitude
  real(r8), intent(in)   :: coslam (plond,platd)      ! cosine of longitude
  real(r8), intent(in)   :: lbasdy (4,2,platd)        ! basis functions for lat deriv est.
  real(r8), intent(in)   :: lbasdz (4,2,plev)         ! basis functions for vert deriv est
  real(r8), intent(in)   :: lbasiy (4,2,platd)        ! basis functions for Lagrange intrp
  real(r8), intent(in)   :: lbasiz (4,2,plev)         ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: lbassi (4,2,plevp)        ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: detam  (plev)             ! delta eta at levels
  real(r8), intent(in)   :: detai  (plevp)            ! delta eta at interfaces
  real(r8), intent(in)   :: dlam   (platd)            ! longitudinal grid increment
  real(r8), intent(in)   :: cwava  (plat)             ! weight for global water vapor int.
  real(r8), intent(in)   :: etamid (plev)             ! eta at levels
  real(r8), intent(in)   :: etaint (plevp)            ! eta at interfaces

  real(r8), intent(out)   :: grlps1(2*maxm,plat/2)      ! ------------------------------
  real(r8), intent(out)   :: grlps2(2*maxm,plat/2)      ! |
  real(r8), intent(out)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grt2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq2  (2*maxm,plev,plat/2) ! |- see quad for definitions
  real(r8), intent(out)   :: grfu1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv2 (2*maxm,plev,plat/2) ! ------------------------------
  real(r8), intent(inout) :: grfu  (plond,plev,beglat:endlat)    ! nonlinear term - u momentum eqn
  real(r8), intent(inout) :: grfv  (plond,plev,beglat:endlat)    ! nonlinear term - v momentum eqn
  real(r8), intent(inout) :: t2    (plond,plev,beglat:endlat)    ! tot dT/dt to to physics
  real(r8), intent(in)   :: flx_net(plond,beglat:endlat)         ! net flx from physics
  real(r8), intent(out)  :: vcour  (plev,plat)          ! maximum Courant number in vert.
  real(r8), intent(out)  :: vmax2d (plev,plat)          ! max. wind at each level, latitude
  real(r8), intent(out)  :: vmax2dt(plev,plat)          ! max. truncated wind at each lvl,lat
!
!---------------------------Local workspace-----------------------------
!
  integer m                 ! constituent index
  integer irow              ! latitude pair index
  integer lat,j,latn,lats   ! latitude indices
  real(r8) onepeps          ! 1 + epssld
  real(r8) dtr              ! 1/dt
  real(r8) fftbuf_in(plond,plev,5,beglat:endlat) ! buffer used for FFTs
#if ( defined SPMD )
  real(r8) fftbuf_out(2*maxm,plev,5,plat) ! buffer used for FFTs
#else
  real(r8) fftbuf_out(1) ! buffer unused
#endif
!
!-----------------------------------------------------------------------
!

!$OMP PARALLEL DO PRIVATE (J, LAT)

  do lat=beglat,endlat
     j = j1 - 1 + lat

     call linemsdyn(lat                     ,ps     (1,lat,n3)   ,u3   (i1,1,j,n3)  , &
                       u3     (i1,1,j,n3m1) ,v3     (i1,1,j,n3)  ,                    &
                    v3        (i1,1,j,n3m1) ,t3     (i1,1  ,j,n3),t3   (i1,1,j,n3m1), &
                       q3     (i1,1,1,j,n3) ,etadot (i1,1,j,n3)  ,                    &
                    etadot    (i1,1,j,n3m1) ,etamid              ,ztodt             , &
                       vcour  (1,lat)       ,vmax2d(1,lat)       ,vmax2dt(1,lat)    , &
                       detam               ,                                          &
                    ed1       (1,1,lat)     ,grfu   (1,1,lat)    ,grfv (1,1,lat)    , &
                       lnpssld(i1,1,j)      ,prhssld(i1,1,j)     ,                    &
                    tarrsld   (i1,1,j)      ,parrsld(i1,1,j)     ,t2   (1,1,lat)    , &
                       div    (1,1,lat,n3)  ,tl     (1,1,lat)    ,                    &
                    tm        (1,1,lat)     ,ql     (1,1,lat)    ,qm   (1,1,lat)    , &
                       dpsl   (1,lat)       ,dpsm   (1,lat)      ,                    &
                    phis      (1,lat)       ,phisl  (1,lat)      ,phism(1,lat)      , &
                       omga   (1,1,lat)     ,                                         &
                    u3sld  (i1,1,j)     ,v3sld(i1,1,j)     ,                          &
                       urhs   (1,1,lat)     ,vrhs   (1,1,lat)    ,                    &
#ifdef QVORTDAMP
                    u3aux   (i1,1,j,n3), u3aux(i1,1,j,n3m1), v3aux(i1,1,j,n3), v3aux (i1,1,j,n3m1),&
                    u3sldaux (i1,1,j), v3sldaux(i1,1,j), &
#endif

                    trhs      (1,1,lat)     ,prhs   (1,1,lat)    ,nlon (lat),         &
                       cwava(lat), flx_net(1,lat))

  end do

#if ( defined SPMD )
#ifdef TIMING_BARRIERS
  call t_startf ('sync_bndexch')
  call mpibarrier (mpicom)
  call t_stopf ('sync_bndexch')
#endif
!
! Communicate boundary information 
!
  call t_startf ('bndexch')
  call bndexch
  call t_stopf ('bndexch')
#endif

  onepeps = 1. + epssld
  dtr     = 2./(ztodt*onepeps)
!
! Fill latitude/longitude extensions of constituents and dynamics terms
!
  call extys(pcnst+pnats,plev    ,q3(1,1,1,beglatex,n3), pcnst)
  call extyv(1       ,plev    ,coslam  ,sinlam  ,u3sld  (1,1,beglatex     ), &
                                                 v3sld  (1,1,beglatex     ))
  call extys(1       ,plevp                     ,etadot (1,1,beglatex,n3m1), 1)
!
  call extyv(1       ,plev    ,coslam  ,sinlam  ,u3     (1,1,beglatex,n3m1), &
                                                 v3     (1,1,beglatex,n3m1))
  call extys(1       ,plev                      ,t3     (1,1,beglatex,n3m1), 1)
  call extys(1       ,plev                      ,lnpssld(1,1,beglatex     ), 1)
  call extys(1       ,plev                      ,prhssld(1,1,beglatex     ), 1)
!      
  call extx (pcnst+pnats,plev                 ,q3     (1,1,1,beglatex,n3), pcnst)
  call extx (1       ,plev                      ,u3sld  (1,1,beglatex     ), 1)
  call extx (1       ,plev                      ,v3sld  (1,1,beglatex     ), 1)
  call extx (1       ,plevp                     ,etadot (1,1,beglatex,n3m1), 1)
!
  call extx (1       ,plev                      ,u3     (1,1,beglatex,n3m1), 1)
  call extx (1       ,plev                      ,v3     (1,1,beglatex,n3m1), 1)
  call extx (1       ,plev                      ,t3     (1,1,beglatex,n3m1), 1)
  call extx (1       ,plev                      ,lnpssld(1,1,beglatex     ), 1)
  call extx (1       ,plev                      ,prhssld(1,1,beglatex     ), 1)

#ifdef QVORTDAMP
  call extyv(1       ,plev    ,coslam  ,sinlam  ,u3sldaux  (1,1,beglatex     ), &
                                                 v3sldaux  (1,1,beglatex     ))
  call extyv(1       ,plev    ,coslam  ,sinlam  ,u3aux  (1,1,beglatex,n3m1), &
                                                 v3aux  (1,1,beglatex,n3m1))
  call extx (1       ,plev                      ,u3sldaux  (1,1,beglatex     ), 1)
  call extx (1       ,plev                      ,v3sldaux  (1,1,beglatex     ), 1)
  call extx (1       ,plev                      ,u3aux  (1,1,beglatex,n3m1), 1)
  call extx (1       ,plev                      ,v3aux  (1,1,beglatex,n3m1), 1)
#endif
!
! Begin SLT interpolation
!

!$OMP PARALLEL DO PRIVATE(LAT)

  do lat = beglat,endlat
     call scanslt_bft(ztodt   ,lat     ,dtr     ,iter    ,pmap    , &
                      kdpmpf  ,kdpmph  ,lam     ,phi     ,dphi    , &
                      lbasdy  ,lbasdz  ,lbasiy  ,lbasiz  ,lbassi  , &
                      detam   ,detai   ,dlam    ,cwava   ,etamid  , &
                      etaint  ,grfu    ,grfv    ,ps      ,u3      , &
                      v3      ,t3      ,q3      ,lnpssld ,prhssld , &
                      tarrsld ,parrsld ,n3      ,n3m1    ,u3sld   , &
                      v3sld   , &
#if ( defined QVORTDAMP )
                        u3aux,v3aux,u3sldaux,v3sldaux, &
#endif
                      etadot  ,nlon(lat), fftbuf_in(1,1,1,lat) )
  end do                    ! end lat-loop
!
  call scanslt_fft(fftbuf_in,fftbuf_out)
!
!$OMP PARALLEL DO PRIVATE (IROW, LATN, LATS)

  do irow=1,plat/2

      lats = irow
      latn = plat - irow + 1
#if ( defined SPMD )
      call scanslt_aft(irow, fftbuf_out(1,1,1,lats), fftbuf_out(1,1,1,latn), &
                       grlps1(1,irow)  ,grlps2(1,irow) , &
                       grt1(1,1,irow)  ,grt2(1,1,irow) , &
                       grq1(1,1,irow)  ,grq2(1,1,irow) , &
                       grfu1(1,1,irow) ,grfu2(1,1,irow), &
                       grfv1(1,1,irow) ,grfv2(1,1,irow)  )
#else
      call scanslt_aft(irow, fftbuf_in(1,1,1,lats), fftbuf_in(1,1,1,latn), &
                       grlps1(1,irow)  ,grlps2(1,irow) , &
                       grt1(1,1,irow)  ,grt2(1,1,irow) , &
                       grq1(1,1,irow)  ,grq2(1,1,irow) , &
                       grfu1(1,1,irow) ,grfu2(1,1,irow), &
                       grfv1(1,1,irow) ,grfv2(1,1,irow)  )
#endif
  end do                    ! end irow-loop
!
!
  return
end subroutine scandyn
