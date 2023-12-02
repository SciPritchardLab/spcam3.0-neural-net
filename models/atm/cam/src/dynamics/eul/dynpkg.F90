#include <misc.h>
#include <params.h>

subroutine dynpkg (t2      ,fu      ,fv      ,etamid  ,etaint  , &
                   cwava   ,detam   ,dlam    ,lam     ,phi     , &
                   dphi    ,sinlam  ,coslam  ,lbasdy  ,lbasdz  , &
                   lbassd  ,lbasiy  ,detai   ,kdpmpf  ,kdpmph  , &
                   flx_net ,ztodt   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routines for dynamics and transport.
! 
! Method: 
! 
! Author: 
! Original version:  CCM3
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   integer pmap   ! max dimension of evenly spaced vert. grid used 
!                    ! by SLT code to map the departure pts into true 
!                    ! model levels.
   parameter ( pmap = 20000 )
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(inout) :: t2(plond,plev,beglat:endlat)         ! temp tendency
   real(r8), intent(inout) :: fu(plond,plev,beglat:endlat)         ! u wind tendency
   real(r8), intent(inout) :: fv(plond,plev,beglat:endlat)         ! v wind tendency

   real(r8), intent(in) :: etamid(plev)                ! vertical coords at midpoints 
   real(r8), intent(in) :: etaint(plevp)               ! vertical coords at interfaces
   real(r8), intent(in) :: cwava(plat)                 ! weight applied to global integrals
   real(r8), intent(in) :: detam(plev)                 ! intervals between vert full levs.
   real(r8), intent(in) :: dlam(platd)                 ! longitudinal grid interval (radians)
   real(r8), intent(in) :: lam(plond,platd)            ! longitude coords of extended grid
   real(r8), intent(in) :: phi(platd)                  ! latitude  coords of extended grid
   real(r8), intent(in) :: dphi(platd)                 ! latitude intervals (radians)
   real(r8), intent(in) :: sinlam(plond,platd)         ! sin(lam) model domain only           
   real(r8), intent(in) :: coslam(plond,platd)         ! cos(lam) model domain only           
   real(r8), intent(in) :: lbasdy(4,2,platd)           ! latitude derivative weights          
   real(r8), intent(in) :: lbasdz(4,2,plev)            ! vert (full levels) deriv wghts 
   real(r8), intent(in) :: lbassd(4,2,plevp)           ! vert (half levels) deriv wghts 
   real(r8), intent(in) :: lbasiy(4,2,platd)           ! Lagrange cubic interp wghts (lat.) 
   real(r8), intent(in) :: detai(plevp)                ! intervals between vert half levs.
   integer, intent(in) :: kdpmpf(pmap)                 ! artificial full vert grid indices
   integer, intent(in) :: kdpmph(pmap)                 ! artificial half vert grid indices
   real(r8), intent(in) :: flx_net(plond,beglat:endlat)  ! net flux from physics
   real(r8), intent(in) :: ztodt                       ! twice time step unless nstep=0
!
!---------------------------Local workspace-----------------------------
!
   real(r8) etadot(plon,plevp,beglat:endlat)     ! Vertical motion (slt)
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSAC and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
   real(r8) grlps1(2*maxm,plat/2)      ! ------------------------------
   real(r8) grlps2(2*maxm,plat/2)      ! |
   real(r8) grt1(2*maxm,plev,plat/2)   ! |
   real(r8) grt2(2*maxm,plev,plat/2)   ! |
   real(r8) grz1(2*maxm,plev,plat/2)   ! |
   real(r8) grz2(2*maxm,plev,plat/2)   ! |
   real(r8) grd1(2*maxm,plev,plat/2)   ! |
   real(r8) grd2(2*maxm,plev,plat/2)   ! |
   real(r8) grfu1(2*maxm,plev,plat/2)  ! |- see quad for definitions
   real(r8) grfu2(2*maxm,plev,plat/2)  ! | 
   real(r8) grfv1(2*maxm,plev,plat/2)  ! |
   real(r8) grfv2(2*maxm,plev,plat/2)  ! |
   real(r8) grut1(2*maxm,plev,plat/2)  ! |
   real(r8) grut2(2*maxm,plev,plat/2)  ! |
   real(r8) grvt1(2*maxm,plev,plat/2)  ! |
   real(r8) grvt2(2*maxm,plev,plat/2)  ! |
   real(r8) grrh1(2*maxm,plev,plat/2)  ! |
   real(r8) grrh2(2*maxm,plev,plat/2)  ! ------------------------------
   real(r8) :: vcour(plev,plat)        ! maximum Courant number in slice
   real(r8) :: vmax2d(plev,plat)       ! max. wind at each level, latitude
   real(r8) :: vmax2dt(plev,plat)      ! max. truncated wind at each lvl,lat
!
!----------------------------------------------------------
! SCANDYN Dynamics scan
!----------------------------------------------------------
!
   call t_startf('scandyn')
   call scandyn(ztodt   ,etadot  ,etamid  ,grlps1  ,grt1    ,  &
                grz1    ,grd1    ,grfu1   ,grfv1   ,grut1   ,  &
                grvt1   ,grrh1   ,grlps2  ,grt2    ,grz2    ,  &
                grd2    ,grfu2   ,grfv2   ,grut2   ,grvt2   ,  &
                grrh2   ,vcour   ,vmax2d,  vmax2dt ,detam   ,  &
                cwava   ,flx_net ,t2      ,fu      ,fv      )
   call t_stopf('scandyn')
!
!----------------------------------------------------------
! SLT scan from south to north
!----------------------------------------------------------
!
   call t_startf('scanslt')
   call scanslt(ztodt   ,pmap    ,etadot  ,kdpmpf  ,kdpmph  ,  &
                lam     ,phi     ,dphi    ,sinlam  ,coslam  ,  &
                lbasdy  ,lbasdz  ,lbassd  ,lbasiy  ,detam   ,  &
                detai   ,dlam    ,cwava   ,etamid  ,etaint  )
   call t_stopf('scanslt')
!
!----------------------------------------------------------
! Accumulate spectral coefficients
!----------------------------------------------------------
!
   call t_startf('dyndrv')
   allocate( vz  (2*lpspt,plev) )
   allocate( d   (2*lpspt,plev) )
   allocate( t   (2*lpspt,plev) )
   allocate( alps(2*lpspt) )

   call dyndrv(grlps1  ,grt1    ,grz1    ,grd1    ,grfu1   ,  &
               grfv1   ,grut1   ,grvt1   ,grrh1   ,grlps2  ,  &
               grt2    ,grz2    ,grd2    ,grfu2   ,grfv2   ,  &
               grut2   ,grvt2   ,grrh2   ,vmax2d  ,vmax2dt ,  &
               vcour   )
   call t_stopf('dyndrv')
!
!----------------------------------------------------------
! Second gaussian scan (spectral -> grid)
!----------------------------------------------------------
!
   call t_startf('scan2')
   call scan2 (ztodt,   cwava,   etamid)

   deallocate( vz )
   deallocate( d )
   deallocate( t )
   deallocate( alps )
   call t_stopf('scan2')

   return
end subroutine dynpkg

