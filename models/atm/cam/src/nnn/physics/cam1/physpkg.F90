#include <misc.h>
#include <params.h>

subroutine physpkg(phys_state, gw, ztodt, phys_tend, pbuf)


!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics
! 
! Method: 
! COUP_CSM and must be checked in order to invoke the proper calling
! sequence for running the CSM model
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, plat, masterproc
   use ppgrid, only: pcols, pver
#ifdef CRM
   use buffer, only: pblht, tpert, qpert, qrs, qrl, &
                     u_crm, v_crm, w_crm, t_crm, q_crm, qn_crm, qp_crm, qrs_crm, qrl_crm, &
                     rad_buffer, qrs1, qrl1, &
		     fsds_crm,fsns_crm,fsntoa_crm,fsutoa_crm,  &
		     flwds_crm,flns_crm,flut_crm, &
		     fsdsc_crm,fsntoac_crm,flnsc_crm, flutc_crm
#else
   use buffer, only: pblht, tpert, qpert, qrs, qrl
#endif
   use check_energy, only: check_energy_gmean
   use comsrf
   use comsrfdiag
#ifdef COUP_CSM
   use ccsm_msg, only: ccsmave, dorecv, dosend, ccsmsnd, ccsmrcv
#else
   use atm_lndMod, only: atmlnd_drv
#endif
#ifdef SPMD
   use mpishorthand
#endif
   use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_allocate, pbuf_deallocate, &
                             pbuf_update_tim_idx
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physics_types,  only: physics_state, physics_tend
   use diagnostics,    only: diag_surf
   use time_manager,   only: get_nstep, is_first_step, is_first_restart_step, &
                             is_end_curr_month, get_curr_date, is_end_curr_day
   use constituents,    only: pcnst, pnats, ppcnst, qphystendnam
   use physconst, only: stebol
#if (!defined COUP_CSM)
   use ice_constants, only: TfrezK
#endif
   use history, only: outfld
#if (!defined COUP_CSM)
#if (!defined COUP_SOM)
   use sst_data, only: sstint
   use ice_data, only: iceint
#endif
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gw(plat)                    ! Gaussian weights
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,m,lat,c,lchnk                   ! indices
   integer :: lats(pcols)                       ! array of latitude indices
   integer :: lons(pcols)                       ! array of longitude indices
   integer :: ncol                              ! number of columns
   integer :: nstep                             ! current timestep number
   integer :: ncdate                            ! current date in integer format [yyyymmdd]
   integer :: ncsec                             ! current time of day [seconds]
   integer :: yr, mon, day                      ! year, month, and day components of a date
                                                
   real(r8) fsds(pcols,begchunk:endchunk)        ! Surface solar down flux
   real(r8) :: tmp(pcols,begchunk:endchunk)
   real(r8) :: tphystend(pcols,pver,begchunk:endchunk)
   real(r8) :: qphystend(pcols,pver,ppcnst,begchunk:endchunk)
   real(r8) :: dtinv ! inverse timestep
!-----------------------------------------------------------------------

   call t_startf ('physpkg_st')
   nstep = get_nstep()

   call pbuf_allocate('physpkg')

! Compute total energy of input state and previous output state
   call t_startf ('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf, ztodt, nstep)
   call t_stopf ('chk_en_gmean')

!-----------------------------------------------------------------------
! Advance time information
!-----------------------------------------------------------------------

   call advnce()
   call t_stopf ('physpkg_st')
!
!
!-----------------------------------------------------------------------
! Tendency physics before flux coupler invokation
!-----------------------------------------------------------------------
!
   call t_startf ('bc_physics')

!$OMP PARALLEL DO PRIVATE (C,NCOL)

   do c=begchunk, endchunk
      call t_startf ('tphysbc')
      ncol = get_ncols_p(c)
      tphystend(:ncol,:,c) = phys_state(c)%t(:ncol,:)
      qphystend(:ncol,:,:,c) = phys_state(c)%q(:ncol,:,:)
      call tphysbc (ztodt, pblht(1,c), tpert(1,c),             &
	              srfflx_state2d(c)%ts, srfflx_state2d(c)%sst,        &
                    qpert(1,1,c), surface_state2d(c)%precl,               &
	   	      surface_state2d(c)%precc, surface_state2d(c)%precsl,&
                      surface_state2d(c)%precsc,                          &
                    srfflx_state2d(c)%asdir, srfflx_state2d(c)%asdif,    &
                      srfflx_state2d(c)%aldir, srfflx_state2d(c)%aldif,  &
                      snowhland(1,c),                                    &
                    qrs(1,1,c), qrl(1,1,c), surface_state2d(c)%flwds,     &
                      fsns(1,c), fsnt(1,c),                               &
                    flns(1,c),    flnt(1,c), srfflx_state2d(c)%lwup,       &
                    surface_state2d(c)%srfrad, surface_state2d(c)%sols,    &
                      surface_state2d(c)%soll, surface_state2d(c)%solsd,   &
                      surface_state2d(c)%solld,                           &
                      phys_state(c), phys_tend(c),                        &
			pbuf, prcsnw(1,c), fsds(1,c), landm(1,c), landfrac(1,c), &
			ocnfrac(1,c),icefrac(1,c) &
#ifdef CRM
                      ,u_crm(1,1,1,1,c), v_crm(1,1,1,1,c), w_crm(1,1,1,1,c) &
                      ,t_crm(1,1,1,1,c), q_crm(1,1,1,1,c), qn_crm(1,1,1,1,c), qp_crm(1,1,1,1,c) &
                      ,qrs_crm(1,1,1,1,c), qrl_crm(1,1,1,1,c), rad_buffer(1,1,c) &
                      ,qrs1(1,1,c), qrl1(1,1,c) &
		      ,fsds_crm(1,1,1,c), fsns_crm(1,1,1,c) &
		      ,fsntoa_crm(1,1,1,c), fsutoa_crm(1,1,1,c)  &
		      ,flwds_crm(1,1,1,c), flns_crm(1,1,1,c), flut_crm(1,1,1,c) &
		      ,fsdsc_crm(1,1,1,c), fsntoac_crm(1,1,1,c) &
		      ,flnsc_crm(1,1,1,c), flutc_crm(1,1,1,c) &
                      ,srfflx_state2d(c)%wsx, srfflx_state2d(c)%wsy &
	              ,srfflx_state2d(c)%shf, srfflx_state2d(c)%lhf &
#endif
                      )

      call t_stopf ('tphysbc')

      if (dosw .or. dolw) then
	call output_flns_fsns_fluxes(surface_state2d(c),c)
      end if	

#if ( ! defined COUP_CSM )
!
! zero surface fluxes at beginning of each time step.  Land Ocean and Ice
! processes will will write into process specific flux variables
! at the end of the time step these separate fluxes will be combined over the
! entire grid
!
      call srfflx_state_reset (srfflx_state2d(c))
#endif

   end do
   call t_stopf ('bc_physics')

#if ( ! defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities - no flux coupler
!-----------------------------------------------------------------------
!
   if (.not. aqua_planet) then
!
! Call land model driving routine
!
#ifdef TIMING_BARRIERS
      call t_startf ('sync_tphysbc_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_tphysbc_lnd')
#endif
      call t_startf ('atmlnd_drv')
      call atmlnd_drv(nstep, iradsw, eccen, obliqr, lambm0,&
                      mvelpp,surface_state2d,srfflx_parm2d)
      call t_stopf ('atmlnd_drv')
#ifdef TIMING_BARRIERS
      call t_startf ('sync_after_lnd')
      call mpibarrier (mpicom)
      call t_stopf ('sync_after_lnd')
#endif
!
! save off albedos and longwave for som offline vars
!
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            if (landfrac(i,c) > 0.) then
               asdirlnd(i,c) = srfflx_parm2d(c)%asdir(i)
               asdiflnd(i,c) = srfflx_parm2d(c)%asdif(i)
               aldirlnd(i,c) = srfflx_parm2d(c)%aldir(i)
               aldiflnd(i,c) = srfflx_parm2d(c)%aldif(i)
               lwuplnd(i,c)  = srfflx_parm2d(c)%lwup(i)
            else
               asdirlnd(i,c) = 0. 
               asdiflnd(i,c) = 0. 
               aldirlnd(i,c) = 0. 
               aldiflnd(i,c) = 0. 
               lwuplnd(i,c)  = 0. 
            end if
         end do
!
!output shf/lhf fluxes for land model
!
         call output_shf_lhf_fluxes(srfflx_parm2d(c), c, ncol, landfrac(1,c), 'LND')
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)
      end do
   end if                    ! end of not aqua_planet if block

#if (defined COUP_SOM)
!
! Set ocean surface quantities - ocn model internal to atm
!
   if (is_end_curr_day ()) then
      call print_coverage ('icefrac', ' million km^2', icefrac, 1.d-12)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*sicthk(i,c)
         end do
      end do
      call print_coverage ('icevol ', ' 10^13m^3', tmp, 1.d-13)

      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = icefrac(i,c)*snowhice(i,c)
         end do
      end do
      call print_coverage ('snowvol', ' 10^13m^3', tmp, 1.d-13)
   end if

   call t_startf ('somint')
   call somint ()
   call t_stopf ('somint')

   call t_startf ('somoce')
   call somoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('somoce')

#else

   call t_startf ('sstint')
   call sstint ()
   call t_stopf ('sstint')
!
! iceint may change ocean fraction, so call it before camoce
!
   call t_startf ('iceint')
   call iceint ()
   call t_stopf ('iceint')

   call t_startf ('camoce')
   call camoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('camoce')
#endif
!
! Set ice surface quantities - icn model internal to atm
!
   call t_startf('camice')
   call camice (surface_state2d, srfflx_parm2d)
   call t_stopf('camice')
!
! output shf/lhf fluxes for ice/ocn/som_offline 
!
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do i=1,ncol
         if(icefrac(i,c) > 0.) then
            tsice_rad(i,c) = sqrt(sqrt(srfflx_parm2d(c)%lwup(i)/stebol))
         else
            tsice_rad(i,c) = TfrezK
         endif
      end do
      call output_shf_lhf_fluxes (srfflx_parm2d(c), c, ncol, icefrac(1,c), 'ICE')
      call output_shf_lhf_fluxes (srfflx_parm2d_ocn(c), c, ncol, ocnfrac(1,c), 'OCN')
      call output_shfoi_lhfoi_fluxes (srfflx_parm2d_ocn(c), srfflx_parm2d(c), c)

!JR SOM case: Have to wait to call update routine till after both ocean and ice have
!JR operated, since the fractions can change internal to the parameterization
      do i = 1, ncol
         srfflx_state2d(c)%sst(i) = srfflx_parm2d_ocn(c)%ts(i)
      enddo
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d_ocn(c), ocnfrac(1,c), ncol)
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), icefrac(1,c), ncol)
   end do
#endif

#if ( defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities using csm flux coupler
!-----------------------------------------------------------------------
!
! If send data to flux coupler only on radiation time steps:
!
   if (flxave) then
!
! Average the precipitation input to lsm between radiation calls.
!
      call ccsmave(iradsw, nstep, dosw)
!
! Use solar radiation flag to determine data exchange steps 
! with flux coupler. This processes are not independent since 
! instantaneous radiative fluxes are passed, valid over the 
! interval to the next radiation calculation. The same 
! considerations apply to the long and shortwave fluxes, so 
! the intervals must be the same. Data is received from the 
! coupler one step after it is sent.
!
      if (nstep == 0) then
         dorecv = .true.
         dosend = .true.
      else if (nstep == 1) then
         dorecv = .false.
         dosend = .false.
      else if ( (nstep == 2) .and. (iradsw == 1) ) then
         dorecv = .true.
         dosend = dosw
      else
         dorecv = dosend
         dosend = dosw
      end if
   endif
!
! If send data to flux coupler on every time step
!
   if (.not. flxave) then
      if (nstep /= 1) then
         dorecv = .true.
         dosend = .true.
      else 
         dorecv = .false.
         dosend = .false.
      endif
   endif
!
! Send/recv data to/from the csm flux coupler.
!
   if (dosend) call ccsmsnd ( )
   if (dorecv) call ccsmrcv ( )
#endif
!
!-----------------------------------------------------------------------
! Tendency physics after coupler 
! Not necessary at terminal timestep.
!-----------------------------------------------------------------------
!
   call t_startf ('ac_physics')

!$OMP PARALLEL DO PRIVATE (C, NCOL)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
!
! surface diagnostics for history files
!
      call diag_surf (c, ncol, srfflx_state2d(c)%shf, srfflx_state2d(c)%lhf, srfflx_state2d(c)%cflx, &
                      srfflx_state2d(c)%tref, trefmxav(1,c), trefmnav(1,c), srfflx_state2d(c)%qref, &
                      srfflx_state2d(c)%wsx, srfflx_state2d(c)%wsy, &
                      icefrac(1,c), ocnfrac(1,c), surface_state2d(c)%tssub, srfflx_state2d(c)%ts, sicthk(1,c), &
                      snowhland(1,c),snowhice(1,c), tsnam,   landfrac(1,c) , tsice(1,c), &
                      surface_state2d(c)%tbot)
      call t_startf ('tphysac')
      call tphysac (ztodt, pblht(1,c), qpert(1,1,c), tpert(1,c), srfflx_state2d(c)%shf,        &
                    srfflx_state2d(c)%wsx,srfflx_state2d(c)%wsy, srfflx_state2d(c)%cflx, sgh(1,c), srfflx_state2d(c)%lhf,        &
                    landfrac(1,c), snowhland(1,c),srfflx_state2d(c)%tref, surface_state2d(c)%precc, surface_state2d(c)%precl,    &
                    surface_state2d(c)%precsc, surface_state2d(c)%precsl, phys_state(c), phys_tend(c), pbuf, &
                    ocnfrac(1,c), fsds(1,c), icefrac(1,c), fv(1,c), ram1(1,c))
      dtinv=1./ztodt
      tphystend(:ncol,:,c) = (phys_state(c)%t(:ncol,:) - tphystend(:ncol,:,c))*dtinv
      qphystend(:ncol,:,:,c) = (phys_state(c)%q(:ncol,:,:) - qphystend(:ncol,:,:,c))*dtinv
      call outfld('TPHYSTND',tphystend(1,1,c),pcols   ,c   )
      do m=1,ppcnst
         call outfld(qphystendnam(m),qphystend(1,1,m,c),pcols   ,c   )
      enddo
      call t_stopf ('tphysac')
   end do                    ! Chunk loop

   call t_stopf('ac_physics')

   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()

end subroutine physpkg
