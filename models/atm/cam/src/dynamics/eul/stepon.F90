#include <misc.h>
#include <params.h>

subroutine stepon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics, dynamics, 
! transport
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Restructured:      J. Truesdale, May 1999
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup, write_inithist, inithist
   use pmgrid
   use pspect
   use comslt
   use rgrid
   use prognostics
   use restart, only: write_restart
#if (defined COUP_CSM)
   use ccsm_msg, only: csmstop, ccsmfin
#endif

   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use phys_buffer,    only: pbuf
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   use commap
   use physconst, only: gravit
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day, get_curr_date
   implicit none

   integer pmap   ! max dimension of evenly spaced vert. grid used 
!                    ! by SLT code to map the departure pts into true 
!                    ! model levels.
   parameter ( pmap = 20000 )
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comqfl.h>
!-----------------------------------------------------------------------
!
   integer kdpmpf(pmap)             ! artificial full vert grid indices
   integer kdpmph(pmap)             ! artificial half vert grid indices

   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) hyad (plev)             ! del (A)
   real(r8) lam(plond,platd)        ! longitude coords of extended grid
   real(r8) phi(platd)              ! latitude  coords of extended grid
   real(r8) dphi(platd)             ! latitude intervals (radians)
   real(r8) gw(plat)                ! Gaussian weights
   real(r8) sinlam(plond,platd)     ! sin(lam) model domain only           
   real(r8) coslam(plond,platd)     ! cos(lam) model domain only           
   real(r8) lbasdy(4,2,platd)       ! latitude derivative weights          
   real(r8) lbasdz(4,2,plev)        ! vert (full levels) deriv wghts 
   real(r8) lbassd(4,2,plevp)       ! vert (half levels) deriv wghts 
   real(r8) lbasiy(4,2,platd)       ! Lagrange cubic interp wghts (lat.) 
   real(r8) detam(plev)             ! intervals between vert full levs.
   real(r8) detai(plevp)            ! intervals between vert half levs.
   real(r8) dlam(platd)             ! longitudinal grid interval (radians)
   real(r8) cwava(plat)             ! weight applied to global integrals
   real(r8) etamid(plev)            ! vertical coords at midpoints 
   real(r8) etaint(plevp)           ! vertical coords at interfaces
   real(r8), allocatable :: t2(:,:,:) ! temp tendency
   real(r8), allocatable :: fu(:,:,:) ! u wind tendency
   real(r8), allocatable :: fv(:,:,:) ! v wind tendency
   real(r8) flx_net(plond,beglat:endlat)       ! net flux from physics
   real(r8) coslat(plond)
   real(r8) rcoslat(plond)
   real(r8) rpmid(plond,plev)
   real(r8) pdel(plond,plev)
   real(r8) pint(plond,plevp)
   real(r8) pmid(plond,plev)
   real(r8) dtime               ! timestep size
   real(r8) ztodt               ! twice time step unless nstep=0
   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep
!
   integer i,k,lat,j,begj       ! longitude,level,latitude indices
   integer iter
   integer :: yr, mon, day    ! year, month, and day components of a date
   integer :: ncsec           ! current time of day [seconds]
   logical l_write_inithist
!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files
!
!-----------------------------------------------------------------------
   call t_startf ('stepon_startup')
   dtime = get_step_size()
!
! Define eta coordinates: Used for calculation etadot vertical velocity 
! for slt.
!
   do k=1,plev
      etamid(k) = hyam(k) + hybm(k)
   end do
   do k=1,plevp
      etaint(k) = hyai(k) + hybi(k)
   end do
!
! Set slt common block variables
!
   call grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
               lam     ,phi     ,dphi    ,gw      ,sinlam  , &
               coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
               detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )
!
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time 
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
   if (is_first_step()) then
      do lat=beglat,endlat
         j = j1 - 1 + lat
         do k=1,plev
            do i=1,nlon(lat)
               lammp(i,k,lat) = float(i-1)*dlam(j1-1+lat)
               phimp(i,k,lat) = clat(lat)
               sigmp(i,k,lat) = etamid(k)
            end do
         end do
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1./coslat(i)
         end do
!     
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plond, plev, ps(1,lat,n3), pint, pmid, pdel)
!
         do k=1,plev
            do i=1,nlon(lat)
               rpmid(i,k) = 1./pmid(i,k)
            end do
         end do
!
! Calculate vertical motion field
!
         call omcalc (rcoslat, div(1,1,lat,n3), u3(i1,1,j,n3), v3(i1,1,j,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid   ,pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
      end do
   end if
!
! Compute pdel from "A" portion of hybrid vertical grid
!
   do k=1,plev
      hyad(k) = hyai(k+1) - hyai(k)
   end do
   do k=1,plev
      do i=1,plon
         pdela(i,k) = hyad(k)*ps0
      end do
   end do

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))
   allocate(t2(plond,plev,beglat:endlat))
   allocate(fu(plond,plev,beglat:endlat))
   allocate(fv(plond,plev,beglat:endlat))
!
! Beginning of basic time step loop
!
   call t_stopf ('stepon_startup')

! Begin time loop.

   do

      call t_startf('stepon_st')
      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

      ztodt = 2.0*dtime
!
! If initial time step adjust dt
!
      if (is_first_step()) ztodt = dtime
!
! adjust hydrostatic matrices if the time step has changed.  This only
! happens on transition from time 0 to time 1. 
! The CMIC$ DO ALL ... construct is a "phony loop" to fool the low level
! Cray matrix library utilities into *not* multitasking, since these 
! utilities give DIFFERENT answers for different values of $NCPUS.  Useful 
! work is done only for iter=1.
      if (get_nstep() == 1) then
!CMIC$ DO ALL SHARED (dtime) PRIVATE (ITER)
         do iter=1,2
            call settau(dtime, iter)
         end do
      end if
!
!----------------------------------------------------------
! PHYSPKG  Call the Physics package
!----------------------------------------------------------
!
      begj = beglatex + numbnd

      call t_stopf('stepon_st')
      call t_startf('d_p_coupling')
      call d_p_coupling (ps(1,beglat,n3m2), t3(i1,1,begj,n3m2), u3(i1,1,begj,n3m2), &
                         v3(i1,1,begj,n3m2), q3(i1,1,1,begj,n3m2), &
                         omga, phis, phys_state, phys_tend, pbuf)
      call t_stopf('d_p_coupling')

      call t_startf('phys_driver')
      if (ideal_phys) then
         call phys_idealized(phys_state, phys_tend, ztodt, etamid)
      else if (adiabatic) then
         call phys_adiabatic(phys_state, phys_tend)
      else
         call physpkg(phys_state, gw, ztodt, phys_tend, pbuf)
      end if
      call t_stopf('phys_driver')
   
      call t_startf('p_d_coupling')
      call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net, &
                         qminus(i1,1,1,begj),  q3(i1,1,1,begj,n3))
      call t_stopf('p_d_coupling')

!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------

      call t_startf('dynpkg')
      call dynpkg (t2      ,fu      ,fv      ,etamid  ,etaint  , &
                   cwava   ,detam   ,dlam    ,lam     ,phi     , &
                   dphi    ,sinlam  ,coslam  ,lbasdy  ,lbasdz  , &
                   lbassd  ,lbasiy  ,detai   ,kdpmpf  ,kdpmph  , &
                   flx_net ,ztodt   )
      call t_stopf('dynpkg')

      call t_startf('stepon_st')
      if (is_first_step() .or. is_first_restart_step()) then
         call print_memusage ('stepon after dynpkg')
      end if

! Set end of run flag.

#if ( ! defined COUP_CSM )
      if (is_last_step()) nlend = .true.
#else
      if (csmstop) then
         if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
         if (is_end_curr_day()) nlend = .true.
      end if
#endif
!
!----------------------------------------------------------
! History and restart logic: Write and/or dispose history tapes if required
!----------------------------------------------------------
!
      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Write restart file
!
      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if
!
! Dispose necessary files
!
      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if
!
! Advance timestep before returning to top of loop
!
      call advance_timestep()
      call get_curr_date(yr, mon, day, ncsec)
      
! Write initial file if requested.

      l_write_inithist = .false.
      if (inithist == '6-HOURLY'  ) then
         if (mod(ncsec, 21600) == 0) l_write_inithist = .true.
      end if
      if (ncsec == 0                               .and. inithist == 'DAILY'  ) l_write_inithist = .true.
      if (ncsec == 0 .and. day == 1                .and. inithist == 'MONTHLY') l_write_inithist = .true.
      if (ncsec == 0 .and. day == 1 .and. mon == 1 .and. inithist == 'YEARLY' ) l_write_inithist = .true.
      if (l_write_inithist) call write_inithist

      call t_stopf('stepon_st')
!
! Check for end of run
!
      if (nlend) then
         call print_memusage ('End stepon')
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(t2)
         deallocate(fu)
         deallocate(fv)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

   end do  ! End of timestep loop

end subroutine stepon
