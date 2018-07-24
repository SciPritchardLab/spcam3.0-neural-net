#include <misc.h>
#include <params.h>

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon --- Loop over time, perform physics and dynamics
!
! !INTERFACE:
subroutine stepon

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup, write_inithist, inithist
   use pmgrid
   use constituents, only: pcnst, pnats
   use comsrf, only: surface_state2d
   use rgrid
   use prognostics
   use commap, only: w
   use dynamics_vars, only: cosp, sinp, cose, sine, coslon, sinlon,    &
                      ak, bk, ks, ptop
   use fv_prints, only: fv_out
   use restart, only: write_restart
   use dynconst, only: omega, rearth
   use physconst, only: gravit, rair

#if (defined COUP_CSM)
   use ccsm_msg, only: csmstop, ccsmfin
#endif

   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use phys_buffer,    only: pbuf
   use dp_coupling
   use physconst, only: zvir, cappa, rga, cpair
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           get_curr_date, is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day

#if defined( SPMD )
   use spmd_dyn, only: comm_z, ijk_yz_to_xy
   use parutilitiesmodule, only: pargatherreal, parcollective3d, sumop
   use mod_irreg, only: mp_sendirr, mp_recvirr
#endif

   implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>
#include <comsta.h>

!
! !DESCRIPTION:
!
!   Loop over time, calling driving routines for physics, dynamics, 
!   transport
!
!   Warning: The Lin-Rood dynamical core ghost points are hidden from
!   the user. If plon is not equal to plond there may be unexpected 
!   results (especially since many physics arrays are defined with 
!   plond). (Glenn Grant)

! !REVISION HISTORY: 
!   00.01.15    Grant      Creation
!   00.08.23    Lin        Modifications
!   01.03.26    Sawyer     Added ProTeX documentation
!   01.06.12    Sawyer     Use dynamics_vars for LR-specific static variables
!   01.06.28    Mirin      Many changes for multi-2D decomposition;
!                          see comments near beginning of main time loop.
!   01.07.13    Mirin      Hid multi-2D decomposition coding
!   01.12.20    Mirin      Changes to fv_out interface
!   02.02.13    Eaton      Pass surface_state2d through fv_out interface.
!   03.08.05    Sawyer     Removed Rayleigh Friction, related variables
!   03.08.13    Sawyer     Ghost zone in u3sxy removed.
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i,k,j,m             ! longitude, level, latitude and tracer indices
   integer t                   ! history tape indices
   integer :: ncdate           ! current date in integer format [yyyymmdd]
   integer :: ncsec            ! time of day relative to current date [seconds]
   integer :: yr, mon, day     ! year, month, day components of a date
   logical l_write_inithist

   real(r8) te0          ! Total energy before dynamics

! Velocity tendencies on a-grid
   real(r8), allocatable :: dudt(:,:,:)
   real(r8), allocatable :: dvdt(:,:,:)
! Work array
   real(r8), allocatable :: dummy3(:,:,:)

!-----------------------------------------------------------------------
! The following arrays are for secondary 2D x-y decomposition
!-----------------------------------------------------------------------

   real(r8), allocatable :: phisxy(:,:)    ! Surface geopotential
   real(r8), allocatable ::   psxy(:,:)    ! Surface pressure
   real(r8), allocatable :: omgaxy(:,:,:)  ! vertical pressure velocity
   real(r8), allocatable ::  u3sxy(:,:,:)  ! Staggered grid winds, latitude
   real(r8), allocatable ::  v3sxy(:,:,:)  ! Satggered grid winds, longitude
   real(r8), allocatable :: delpxy(:,:,:)  ! delta pressure
   real(r8), allocatable ::   ptxy(:,:,:)  ! virtual potential temperature
   real(r8), allocatable ::   q3xy(:,:,:,:)! Moisture and constituents
   real(r8), allocatable ::  pkzxy(:,:,:)  ! finite-volume mean pk
   real(r8), allocatable ::   pkxy(:,:,:)  ! pe**cappa
   real(r8), allocatable ::   pexy(:,:,:)  ! edge pressure
   real(r8), allocatable :: pilnxy(:,:,:)  ! ln(pe)
   real(r8), allocatable :: dudtxy(:,:,:)
   real(r8), allocatable :: dvdtxy(:,:,:)
   real(r8), allocatable :: dummy3xy(:,:,:)

!-----------------------------------------------------------------------

! Other local variables

   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep

   logical full_phys    ! flag to convert pt output from fvcore: on output
                        ! pt will be virtual temp (deg K) if full_phys is .T.
                        ! virtual potential temperature if false.
   integer pdt          ! Physics time step
   real(r8) :: dtime    ! Physics time step

! for fv_out
   integer freq_diag          ! Output frequency in seconds
   logical fv_monitor         ! Monitor Mean/Max/Min fields every time step
   data freq_diag  / 21600 /  ! time interval (sec) for calling fv_out
   data fv_monitor / .true. / ! This is CPU-time comsuming; set it to false for
                              ! production runs

   real (r8), allocatable :: delpz(:,:,:), pez(:,:,:)

!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files

!-----------------------------------------------------------------------

   if ( use_eta ) then
!-----------------------------------------------
! Use ak and bk from the internal set_eta routine
!-----------------------------------------------
      do k = 1, plev+1
         if (iam == 0) write(6,*) k, (hyai(k)-ak(k)*1.e-5),(hybi(k)-bk(k))
         hyai(k) = ak(k) * 1.e-5
         hybi(k) = bk(k)
      end do
   else
!-----------------------------------------
! Use ak and bk as specified by CAM IC
!-----------------------------------------
      do k = 1, plev+1
         ak(k) = hyai(k) * 1.e5
         bk(k) = hybi(k)
         if( bk(k) == 0. ) ks = k-1
      end do
      ptop = ak(1)
      if ( iam == 0 ) then
         write(*,*) 'Using hyai & hybi from IC:', 'KS=',ks,' PTOP=',ptop
      endif
   endif

!----------------------------------------------------------
! Lin-Rood dynamical core initialization
!----------------------------------------------------------

   pdt = get_step_size()    ! Physics time step
   dtime = pdt

   if (plon .ne. plond) then
      print *, "STEPON: PLOND must be set equal to PLON when using"
      print *, "the Lin-Rood dynamical core. Stopping."
      print *, "plond = ",plond
      print *, "plon  = ",plon
      call endrun
   endif

#if (!defined STAGGERED)
   print *,"STEPON: pre-processor variable STAGGERED must be set"
   print *,"in you misc.h file. Enter: #define STAGGERED"
   print *,"Then recompile CAM. Quitting."
   call endrun
#endif

! full_phys is set to true, the output from fvcore is virtual temperature
! (deg K), as needed by the physics.  If full_phys is false, the output
! is virtual potential temperature, as needed in the idealized case.
! Note: zvir=0 in idealized case. Therefore, pt is simply the scaled
! potential temp.

   full_phys = .true.

   if ( ideal_phys .or. adiabatic )  then
       full_phys = .false.
       zvir = 0.
   endif

   if (myid_z .eq. 0) then

!$omp parallel do private(i,j)

      do j = beglat, endlat
         do i=1, plon
            pe(i,1,j) = ptop
         enddo
      enddo

   endif

   if ( nlres ) then
!
! Do not recalculate delta pressure (delp) if this is a restart run.
! Re. SJ Lin: The variable "delp" (pressure thikness for a Lagrangian
! layer) must be in the restart file. This is because delp will be
! modified "after" the physics update (to account for changes in water
! vapor), and it can not be reproduced by surface pressure and the
! ETA coordinate's a's and b's.

      if (npr_z .eq. 1) then

!$omp parallel do private(i,j,k)
         do j = beglat, endlat
            do k=1, plev
               do i=1, plon
                  pe(i,k+1,j) = pe(i,k,j) + delp(i,j,k)
               enddo
            enddo
         enddo

      else

#if defined( SPMD )
         allocate (delpz(plon, beglat:endlat, plev))
         allocate (pez(plon, plevp, beglat:endlat))
         delpz(:,:,:) = 0.
         pez(:,:,:) = 0.
!$omp parallel do private(i,j,k)
         do j = beglat, endlat
            do k=1, plevp
               do i=1, plon
                  pez(i,k,j) = 0.
               enddo
            enddo
         enddo
!$omp parallel do private(i,j,k)
         do k=1, plev
            do j = beglat, endlat
               do i=1, plon
                  delpz(i,j,k) = 0.
               enddo
            enddo
         enddo
         call pargatherreal(comm_z, 0, delp, strip3zaty, delpz)
         call parcollective3d(comm_z, sumop, plon, endlat-beglat+1, plev, delpz)

!$omp parallel do private(i,j)
   do j = beglat, endlat
      do i=1, plon
         pez(i,1,j) = ptop
      enddo
   enddo

!$omp parallel do private(i,j,k)
         do j = beglat, endlat
            do k=1, plev
               do i=1, plon
                  pez(i,k+1,j) = pez(i,k,j) + delpz(i,j,k)
               enddo
            enddo
         enddo

!$omp parallel do private(i,j,k)
         do j = beglat, endlat
            do k=beglev, endlevp1
               do i=1, plon
                  pe(i,k,j) = pez(i,k,j)
               enddo
            enddo
         enddo

#endif
      endif

   else
 
! Initial run --> generate pe and delp from the surface pressure
 
!$omp parallel do private(i,j,k)
         do j = beglat, endlat
            do k=max(2,beglev),endlevp1
               do i=1, plon
                  pe(i,k,j) = ak(k) + bk(k) * ps(i,j)
               enddo
            enddo
         enddo

!$omp parallel do private(i,j,k)
         do k = beglev,endlev
            do j = beglat, endlat
               do i=1, plon
                  delp(i,j,k) = pe(i,k+1,j) - pe(i,k,j)
               enddo
            enddo
         enddo
   endif

!----------------------------------------------------------
! Check total dry air mass; set to 982.22 mb if initial run
! Print out diagnostic message if restart run
!----------------------------------------------------------

  if ( full_phys ) then
     call dryairm( plon,   plat,  plev, beglat, endlat,       &
                   ng_d,   beglev,endlev,                     &
                   .true., ptop,  ps,   q3,     pcnst+pnats,  &
                   pcnst,  delp,  pe,   nlres )
  endif

! Initialize pk, edge pressure to the cappa power.

!$omp parallel do private(i,j,k)
      do k = beglev, endlevp1
         do j = beglat, endlat
            do i = 1, plon
               pk(i,j,k) = pe(i,k,j)**cappa
            enddo
         enddo
      enddo

! Generate pkz, the conversion factor betw pt and t3

   call pkez(1,      plon,   plev,   beglat,  endlat,              &
             beglev, endlev,  1,      plon,    pe,    pk,  cappa,  &
             ks, piln,   pkz,  .false. )

   if ( .not. nlres ) then

! Compute pt for initial run: scaled virtual potential temperature
! defined as (virtual temp deg K)/pkz. pt will be written to restart (SJL)

!$omp parallel do private(i,j,k)
      do k = beglev, endlev
         do j = beglat, endlat
            do i = 1, plon
               pt(i,j,k) =  t3(i,k,j)*(1.+zvir*q3(i,j,k,1))/pkz(i,j,k)
            enddo
         enddo
      enddo
   endif

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))

!----------------------------------------------------------
! Allocate auxiliary variables
!----------------------------------------------------------

   allocate (dudt(plon,beglev:endlev,beglat:endlat))
   allocate (dvdt(plon,beglev:endlev,beglat:endlat))
   allocate (dummy3(plon,beglat:endlat,beglev:endlev))

!----------------------------------------------------------
! Allocate variables for secondary 2D xy decomposition
!----------------------------------------------------------

   allocate (phisxy(beglonxy:endlonxy      , beglatxy:endlatxy))
   allocate (  psxy(beglonxy:endlonxy      , beglatxy:endlatxy))
   allocate (omgaxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))    
   allocate ( u3sxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))  ! now unghosted, was N1
   allocate ( v3sxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
   allocate (delpxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
   allocate (  ptxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
   allocate (  q3xy(beglonxy:endlonxy, beglatxy:endlatxy, plev, pcnst+pnats))
   allocate ( pkzxy(beglonxy:endlonxy, beglatxy:endlatxy, plev))
   allocate (  pkxy(beglonxy:endlonxy, beglatxy:endlatxy, plevp))
   allocate (  pexy(beglonxy:endlonxy, plevp, beglatxy:endlatxy)) 
   allocate (pilnxy(beglonxy:endlonxy, plevp, beglatxy:endlatxy)) 
   allocate (dudtxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
   allocate (dvdtxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
   allocate (dummy3xy(beglonxy:endlonxy,beglatxy:endlatxy,plev))

!-----------------------------------------------------------------------
! Transpose phis if 2D decomposition
! The transposed phis (called phisxy) is needed for the geopotential
!    calculation (when using transpose option) and for remapping
! Embed in 3D array since transpose machinery cannot handle 2D arrays
!-----------------------------------------------------------------------

#if defined( SPMD )
   if (twod_decomp .eq. 1) then
!$omp parallel do private(i,j,k)
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,j,k) = phis(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(dummy3, ijk_yz_to_xy%SendDesc,              &
                      ijk_yz_to_xy%RecvDesc, dummy3xy)
      call mp_recvirr(dummy3xy, ijk_yz_to_xy%RecvDesc )
!$omp parallel do private(i,j)
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            phisxy(i,j) = dummy3xy(i,j,1)
         enddo
      enddo
   endif
#endif

!----------------------------------------------------------
! Beginning of basic time step loop
!----------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION
!
!
!     A 2D xy decomposition is used for handling the Lagrangian surface
!     remapping, the ideal physics, and (optionally) the geopotential
!     calculation.
!
!     The transpose from yz to xy decomposition takes place within dynpkg.
!     The xy decomposed variables are then transposed directly to the
!     physics decomposition within d_p_coupling.
!
!     The xy decomposed variables have names corresponding to the
!     yz decomposed variables: simply append "xy". Thus, "uxy" is the
!     xy decomposed version of "u".
!
!     To assure that the latitudinal decomposition operates
!     as efficiently as before, a separate parameter "twod_decomp" has
!     been defined; a value of 0 means the original 1D decomposition
!     (with the 2D decomposition formalism skipped), and a value of 1
!     utilizes the 2D decomposition formalism.
!
!     For questions/comments, contact Art Mirin, mirin@llnl.gov
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   do

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

!--------------------------------------------------------------------------
! Perform finite-volume dynamics -- this dynamical core contains some 
! yet to be published algorithms. Its use in the CAM is
! for software development purposes only. 
! Please contact S.-J. Lin (slin@dao.gsfc.nasa.gov)
! if you plan to use this mudule for scientific purposes. Contact S.-J. Lin
! or Will Sawyer (sawyer@dao.gsfc.nasa.gov) if you plan to modify the
! software.
!--------------------------------------------------------------------------

      call t_startf('dynpkg')

!----------------------------------------------------------
! For 2-D decomposition, phisxy is input to dynpkg, and the other
! xy variables are output. Some are computed through direct
! transposes, and others are derived.
!----------------------------------------------------------
      call dynpkg(plon,    plat,    plev,    u3s,     v3s,         &
                  ps,      delp,    pe,      pk,      pkz,         &
                  piln,    ptop,    beglat,  endlat,               &
                  beglev,  endlev,  endlevp,  phis,                &
                  omga,    cpair,   pcnst,   q3,      pdt,         &
                  omega,   rair,    rearth,  nsplit,  pcnst+pnats, &
                  pt,      full_phys,        iord,    jord,        &
                  kord,    te0,                                    &
                  u3sxy,   v3sxy,   psxy,    delpxy,  pexy,        &
                  pkxy,    pkzxy,   pilnxy,  phisxy,  omgaxy,      &
                  q3xy,    ptxy,                                   &
                  beglonxy, endlonxy, beglatxy, endlatxy )

      call t_stopf('dynpkg')

      call t_startf('physics')

      call t_startf('dp_coupling_total')
      call t_startf('d_p_coupling')
!----------------------------------------------------------
! Move data into phys_state structure.
!----------------------------------------------------------
      call d_p_coupling (ps,   u3s,  v3s,   pt,    coslon, sinlon,    &
                         q3,   omga,  phis,  pe,     piln,   pk,      &
                         pkz, phys_state, phys_tend, pbuf, full_phys, &
                         psxy,        u3sxy, v3sxy, ptxy,             &
                         q3xy, omgaxy,      phisxy,pexy,  pilnxy,     &
                         pkxy, pkzxy )
      call t_stopf('d_p_coupling')
      call t_stopf('dp_coupling_total')


!----------------------------------------------------------
! Call full physics
!----------------------------------------------------------
      if ( full_phys ) then
!
         call t_startf('physpkg')
         call physpkg (phys_state, w, dtime, phys_tend, pbuf)
         call t_stopf('physpkg')
!----------------------------------------------------------
! Otherwise, adiabatic physics to update calendar and history
!----------------------------------------------------------
      else
         call phys_adiabatic (phys_state, phys_tend)
      endif

!----------------------------------------------------------
! Update dynamics variables using phys_state & phys_tend.
! 2-D decomposition: Compute dudtxy, dvdtxy, ptxy and q3xy; for ideal
!   physics, scale ptxy by (old) pkzxy; then transpose to yz variables
! 1-D decomposition: Compute dudt, dvdt, pt and q3; for ideal physics,
!   scale pt by old pkz.
! Call uv3s_update to update u3s and v3s from dudt and dvdt.
! Call p_d_adjust to update pt, q3, pe, delp, ps, piln, pkz and pk.
! For adiabatic case, transpose to yz variables.
!----------------------------------------------------------

      call t_startf('dp_coupling_total')
      call t_startf('p_d_coupling')
      call p_d_coupling(phys_state,   phys_tend,   full_phys,      &
                        adiabatic, t3, t3old,   q3old,             &
                        q3,    pt,    dudt,    dvdt,     pkz,      &
                        q3xy,  ptxy,  dudtxy,  dvdtxy,   pkzxy,    &
                        dtime, u3s,   v3s,     u3sxy,    v3sxy,    &
                        zvir,  cappa, ptop,    pk,                 &
                        piln,  ps,                                 &
                        pe,    pexy,  delp,    delpxy     )
      call t_stopf('p_d_coupling')
      call t_stopf('dp_coupling_total')

!----------------------------------------------------------
! Monitor max/min/mean of selected fields
!
!  SEE BELOW  ****  SEE BELOW  ****  SEE BELOW

! Beware that fv_out uses both dynamics and physics instanciations.
! However, I think that they are used independently, so that the
! answers are correct. Still, this violates the notion that the
! physics state is no longer active after p_d_coupling.
!----------------------------------------------------------
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      i = pdt + ncsec      !  step complete, but nstep not incremented yet

      if ( fv_monitor .and. mod(i, freq_diag) == 0 ) then
         call fv_out( plon,    plat ,  plev,   beglat,      endlat,     &
                      ng_d,    beglev, endlev, pk,          pt  ,       &
                      ptop,    ps,     q3,     pcnst+pnats, pcnst,      &
                      delp,    pe, surface_state2d, phys_state, ncdate, &
                      i,       full_phys )
      endif

      call t_stopf('physics')

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


!----------------------------------------------------------
! History and restart logic
!----------------------------------------------------------

! Write and/or dispose history tapes if required

      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Shift time pointers (used in cloud water).  For LR, this must be 
! done after physpkg and before writing the restart file.
!
      call shift_time_indices ()

! Write restart file

      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if

! Dispose necessary files

      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if

! Increment timestep before returning to top of loop

      call advance_timestep()
      
! Write initial file if requested.

      call get_curr_date(yr, mon, day, ncsec)
      l_write_inithist = .false.
      if (inithist == '6-HOURLY'  ) then
         if (mod(ncsec, 21600) == 0) l_write_inithist = .true.
      end if
      if (ncsec == 0                               .and. inithist == 'DAILY'  ) l_write_inithist = .true.
      if (ncsec == 0 .and. day == 1                .and. inithist == 'MONTHLY') l_write_inithist = .true.
      if (ncsec == 0 .and. day == 1 .and. mon == 1 .and. inithist == 'YEARLY' ) l_write_inithist = .true.
      if (l_write_inithist) call write_inithist

! Check for end of run

      if (nlend) then
         call print_memusage ('End stepon')
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(dudt)
         deallocate(dvdt)
         deallocate(dummy3)
         deallocate(phisxy)
         deallocate(  psxy)
         deallocate(omgaxy)
         deallocate( u3sxy)
         deallocate( v3sxy)
         deallocate(delpxy)
         deallocate(  ptxy)
         deallocate(  q3xy)
         deallocate( pkzxy)
         deallocate(  pkxy)
         deallocate(  pexy)
         deallocate(pilnxy)
         deallocate(dudtxy)
         deallocate(dvdtxy)
         deallocate(dummy3xy)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

   end do  ! End of timestep loop
!EOC
end subroutine stepon
!-----------------------------------------------------------------------
