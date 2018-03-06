#include <misc.h>
#include <params.h>
!#define CLOUDBRAIN
#define XEONPHI ! Pritch at TACC / Stampede

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
  use crmics
  use runtime_opts, only: crminitread
   use filenames, only: ncdata
   use ioFileMod,    only: getfil

#else
   use buffer, only: pblht, tpert, qpert, qrs, qrl
#endif
   use analyses, only: analyses_int
   use runtime_opts,  only: l_analyses
   use runtime_opts, only: aqua_AndKua, aqua_uniform, aqua_uniform_sst_degC


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
#ifdef QRLDAMP
   use qrl_anncycle, only: read_qrl_anncycle, allocate_qrlanncycle,ref_qrl_anncycle_int
#endif
#ifdef FLUXDAMP
   use sfcwind_anncycle, only: read_sfcwindanncycle, allocate_sfcwindanncycle,ref_sfcwindanncycle_int, sfcwind_interference
   use phys_grid,       only: get_rlat_all_p
   use runtime_opts, only: fluxdampfac,fluxdamp_equatoronly, flux_dylat,flux_critlat_deg

#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
#include <comlun.h>
   include 'netcdf.inc'
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
   real(r8) :: uphystend(pcols,pver,begchunk:endchunk)
   real(r8) :: vphystend(pcols,pver,begchunk:endchunk)
   real(r8) :: qphystend(pcols,pver,ppcnst,begchunk:endchunk)
   real(r8) :: dtinv ! inverse timestep
   character(len=256) locfn ! local filename

#ifdef PVBUDGET
   real(r8)  :: pv1(pcols,pver,begchunk:endchunk),pv1_tmp(pcols,pver,begchunk:endchunk)
   real(r8)  :: pv2(pcols,pver,begchunk:endchunk),pv2_tmp(pcols,pver,begchunk:endchunk)
   real(r8)  :: pv3(pcols,pver,begchunk:endchunk),pv3_tmp(pcols,pver,begchunk:endchunk)
   real (r8) :: pvtot_tmp(pcols,pver,begchunk:endchunk)
#endif
  real (r8) :: aux (pcols,pver)
#ifdef FLUXDAMP
  real (r8) :: clat (pcols), dubot(pcols), dvbot(pcols), dwindbot(pcols)
  real (r8) :: lhf0 (pcols,begchunk:endchunk), shf0(pcols,begchunk:endchunk), dlhf(pcols), dshf(pcols) ! diag
  integer :: ifluxcalc
  type (surface_state) :: mysave_surface_state2d(begchunk:endchunk)
  type (srfflx_state) :: mysave_srfflx_state2d(begchunk:endchunk)
  type (srfflx_parm) :: mysave_srfflx_parm2d(begchunk:endchunk)
  type (srfflx_parm) :: mysave_srfflx_parm2d_ocn (begchunk:endchunk)
#endif

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
! If pushing model with analysis data, process data
!
   if(l_analyses) call analyses_int (ztodt)
#ifdef CRM
   if (is_first_step() .and. crminitread) then
     if (masterproc) then
       call getfil (ncdata, locfn)
       call wrap_open (locfn, NF_NOWRITE, ncid_ini)
     end if
      call read_crm_ics
!      call crmic_debug_dump ('[tphysbc] afterread')
   end if
#endif

#ifdef PVBUDGET
! Mike Pritchard, store PV before physics:
  call calculate_physics_PV (phys_state,pv1,pv2,pv3)
!$OMP PARALLEL DO PRIVATE (C,NCOL)
  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
     pvtot_tmp(:ncol,:pver,c) = phys_state(c)%pv(:ncol,:pver)
     pv1_tmp(:ncol,:pver,c) = pv1(:ncol,:pver,c)
     pv2_tmp(:ncol,:pver,c) = pv2(:ncol,:pver,c)
     pv3_tmp(:ncol,:pver,c) = pv3(:ncol,:pver,c)
  end do
#endif
!-----------------------------------------------------------------------
! Tendency physics before flux coupler invokation
!-----------------------------------------------------------------------
!

#ifdef QRLDAMP
   if (is_first_step() .or. is_first_restart_step()) then
     if (is_first_step()) then
       call allocate_qrlanncycle
     end if 
     call read_qrl_anncycle
   end if
   call ref_qrl_anncycle_int
#endif


   call t_startf ('bc_physics')


#ifndef XEONPHI
! Standard threading... here in the scope of physpkg, surrounding tphysbc

!$OMP PARALLEL DO PRIVATE (C,NCOL)
   call t_startf ('tphysbc')
   do c=begchunk, endchunk
      ncol = get_ncols_p(c)
      tphystend(:ncol,:,c) = phys_state(c)%t(:ncol,:)
      uphystend(:ncol,:,c) = phys_state(c)%u(:ncol,:)
      vphystend(:ncol,:,c) = phys_state(c)%v(:ncol,:)
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
                   ,u_crm, v_crm, w_crm, t_crm, q_crm, qn_crm, qp_crm &
                  ,qrs_crm, qrl_crm, rad_buffer, qrs1, qrl1  &
                   ,fsds_crm,fsns_crm,fsntoa_crm,fsutoa_crm  &
                   ,flwds_crm,flns_crm,flut_crm   &
                   ,fsdsc_crm,fsntoac_crm,flnsc_crm, flutc_crm & 
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
      end do ! loop over chunks
   
#else 
! For MIC / XEON PHI at TACC STAMPEDE.
! Pritch custom threading... separated chunk loops inside tphysbc.
! A minimalist chunk loop bracketing CRM work is hopefuly simple enough to 
! offload to the Many-Integrated-Core coprocessor.

   call t_startf ('tphysbc_internallythreaded')
   do c=begchunk, endchunk
      ncol = get_ncols_p(c)
      tphystend(:ncol,:,c) = phys_state(c)%t(:ncol,:)
      uphystend(:ncol,:,c) = phys_state(c)%u(:ncol,:)
      vphystend(:ncol,:,c) = phys_state(c)%v(:ncol,:)
      qphystend(:ncol,:,:,c) = phys_state(c)%q(:ncol,:,:)
      ! SR: Debug output T 1
      call outfld ('DBGT1',phys_state(c)%t(:ncol,:),pcols,c)
   end do

   call tphysbc_internallythreaded (ztodt, pblht(:,begchunk:endchunk), tpert(:,begchunk:endchunk),             &
                     srfflx_state2d(begchunk:endchunk), &
                    qpert(:,:,begchunk:endchunk),surface_state2d(begchunk:endchunk), &
                      snowhland(:,begchunk:endchunk),                                    &
                    qrs(:,:,begchunk:endchunk), qrl(:,:,begchunk:endchunk), &
                      fsns(:,begchunk:endchunk), fsnt(:,begchunk:endchunk),                               &
                    flns(:,begchunk:endchunk),    flnt(:,begchunk:endchunk), &
                      phys_state(begchunk:endchunk), phys_tend(begchunk:endchunk),  &
                       pbuf, prcsnw(:,begchunk:endchunk), fsds(:,begchunk:endchunk), landm(:,begchunk:endchunk), landfrac(:,begchunk:endchunk),&
                       ocnfrac(:,begchunk:endchunk),icefrac(:,begchunk:endchunk) &
#ifdef CRM
                   ,u_crm, v_crm, w_crm, t_crm, q_crm, qn_crm, qp_crm &
                  ,qrs_crm, qrl_crm, rad_buffer, qrs1, qrl1  &
                   ,fsds_crm,fsns_crm,fsntoa_crm,fsutoa_crm  &
                   ,flwds_crm,flns_crm,flut_crm   &
                   ,fsdsc_crm,fsntoac_crm,flnsc_crm, flutc_crm & 
#endif
)
do c=begchunk, endchunk
      ncol = get_ncols_p(c)
      ! SR: Debug output T 2
      call outfld ('DBGT2',phys_state(c)%t(:ncol,:),pcols,c)
end do                      
      call t_stopf ('tphysbc_internallythreaded')
      do c=begchunk,endchunk 
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

#endif ! XEONPHI
   call t_stopf ('bc_physics')

#if ( ! defined COUP_CSM )
!
!-----------------------------------------------------------------------
! Determine surface quantities - no flux coupler
!-----------------------------------------------------------------------

#ifdef FLUXDAMP
   if (is_first_step() .or. is_first_restart_step()) then
     if (is_first_step()) then
       call allocate_sfcwindanncycle
     end if 
     call read_sfcwindanncycle
   end if
   call ref_sfcwindanncycle_int
#endif



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
! --- pritch: this is where surface fluxes are calculated for land points by clm:
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

#ifdef FLUXDAMP
!$OMP PARALLEL DO PRIVATE (C,NCOL,I,DUBOT,DVBOT,DWINDBOT)
#else
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)
#endif
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
! pritch: this writes the LHFLXLND variable to history
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)
         ! note update_srf_fluxes from comsrf.F90 just addes 2nd arg to first
         ! first calcluation is the baseline - what woul dhave happened without
         ! interference
      end do ! do chunk
   end if                    ! end of not aqua_planet if block

#ifdef FLUXDAMP

 ! =========== pritch double iteration of ocean flux calculations =====
 ! ==== (second time using the auxiliary, interfered surface wind fields ) =====

 ! We will do this twice; second time with interfered winds; to bracket effect on energy
 do c=begchunk,endchunk
   mysave_surface_state2d(c) = surface_state2d(c)
   mysave_srfflx_state2d(c) = srfflx_state2d(c)
   mysave_srfflx_parm2d(c) = srfflx_parm2d(c)
   mysave_srfflx_parm2d_ocn(c) = srfflx_parm2d_ocn(c)
 end do

 do ifluxcalc = 1,2 ! first time just for diag
   if (ifluxcalc .eq. 2) then
      ! Second time, forget what just happened, revert to previous state:
      ! Saved right after land calcluation.
      do c=begchunk,endchunk
         surface_state2d(c) = mysave_surface_state2d(c)
         srfflx_parm2d(c) = mysave_srfflx_parm2d(c)
         srfflx_state2d(c) = mysave_srfflx_state2d(c)
         srfflx_parm2d_ocn(c) = mysave_srfflx_parm2d_ocn(c)
      end do
      ! --- pritch: this is the place to implement surface wind interference ---
      ! ASSUMPTION: modifications to surface_state2d%ubot,vbot will ONLY be felt
      ! by intended surface flux processes in camoce call below. 
       ! (THAT COULD USE SOME VALIDATION...beware unintended consequences.)
      do c=begchunk,endchunk
        ncol = get_ncols_p(c)
        call get_rlat_all_p(c, ncol, clat)
        do i=1,ncol
          call sfcwind_interference (c,ncol,surface_state2d(c)%ubot,surface_state2d(c)%vbot,clat,fluxdampfac,fluxdamp_equatoronly, flux_dylat,flux_critlat_deg,dubot,dvbot,dwindbot,ztodt)
          call outfld('DUBOT   ',dubot ,pcols,c)
          call outfld('DVBOT   ',dvbot ,pcols,c)
          call outfld('DWINDBOT',dwindbot,pcols,c)
        end do
      end do
   end if ! is second iteration? 
! note as-yet-unterminated iteration loop above inn this ifdef block
#endif

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

#ifdef FLUXDAMP
   write (6,*) 'FLUX interference not enabled for SOM.'
   call endrun
#endif


#else
! i.e. not COUPCSM (my case for FLUXDAMP)
#ifdef FLUXDAMP
   if (ifluxcalc .eq. 1) then
#endif
   call t_startf ('sstint')
   call sstint (.false.,aqua_uniform, aqua_AndKua, aqua_uniform_sst_degC)
   call t_stopf ('sstint')
!
! iceint may change ocean fraction, so call it before camoce
!
   call t_startf ('iceint')
   call iceint ()
   call t_stopf ('iceint')
#ifdef FLUXDAMP
  end if
#endif

   call t_startf ('camoce')
   call camoce (surface_state2d, srfflx_parm2d_ocn)
   call t_stopf ('camoce')

#endif
! COUPCSM if block

! Set ice surface quantities - icn model internal to atm
!
   call t_startf('camice')
   call camice (surface_state2d, srfflx_parm2d)
   call t_stopf('camice')

!
! output shf/lhf fluxes for ice/ocn/som_offline 
!
#ifndef FLUXDAMP
!$OMP PARALLEL DO PRIVATE (C,NCOL,I)
#endif
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
#ifdef FLUXDAMP
    if (ifluxcalc .eq. 2) then
      ! 2nd time, calculate the energy interference diagnostic and save to disk
      do i=1,ncol
        dlhf(i) = srfflx_state2d(c)%lhf(i) - lhf0(i,c)
        dshf(i) = srfflx_state2d(c)%shf(i) - shf0(i,c)
      end do
#ifdef MDEBUG
    if (masterproc .and. c .eq. begchunk) then
      write (6,*) 'dlhf = ',dlhf(1),', dshf = ',dshf(1)
    end if
#endif
      call outfld('DLHFLX  ',dlhf ,pcols,c)
      call outfld('DSHFLX  ',dshf ,pcols,c)

   else ! first iteration, only need to save for diags:
        do i = 1,ncol
          lhf0(i,c) = srfflx_state2d(c)%lhf(i)
          shf0(i,c) = srfflx_state2d(c)%shf(i)
        end do          
   endif ! is second iteration?
#ifdef MDEBUG
    if (masterproc .and. c .eq. begchunk) then
      write (6,*) 'iter=',ifluxcalc,', srfflx_state2d%lhf =',srfflx_state2d(c)%lhf(1)
      write (6,*) 'iter=',ifluxcalc,', srfflx_state2d%shf =',srfflx_state2d(c)%shf(1)
    end if
#endif

#endif

   end do ! End chunk loop for ocn/ice flux calculation.
#endif

#ifdef FLUXDAMP
  end do
! =========== pritch end double iteration of ocean flux calculations =====
 ! ==== (second time using the auxiliary, interfered surface wind fields ) =====
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

! MSP: Note we cannot recalculate PV after tphybc because state%t is out of date.
! (only state%s has been updated so far. The next place we can calcaulte PV is after tphysac, 
! i.e. bracketing the total physics package tendency).
!
!-----------------------------------------------------------------------
! Tendency physics after coupler 
! Not necessary at terminal timestep.
!-----------------------------------------------------------------------
!
   call t_startf ('ac_physics')
      dtinv=1./ztodt

!$OMP PARALLEL DO PRIVATE (C, NCOL, M, AUX)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
!
! surface diagnostics for history files
!
      call outfld ('DBGT3',phys_state(c)%t(:ncol,:pver),pcols,c) ! SR: Debug T 3
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
      tphystend(:ncol,:,c) = (phys_state(c)%t(:ncol,:) - tphystend(:ncol,:,c))*dtinv
      call outfld('TPHYSTND',tphystend(1,1,c),pcols   ,c   )
      uphystend(:ncol,:,c) = (phys_state(c)%u(:ncol,:) - uphystend(:ncol,:,c))*dtinv
      call outfld ('UPHYSTND',uphystend(1,1,c),pcols, c)
      vphystend(:ncol,:,c) = (phys_state(c)%v(:ncol,:) - vphystend(:ncol,:,c))*dtinv
      call outfld ('VPHYSTND',vphystend(1,1,c),pcols, c)
      qphystend(:ncol,:,:,c) = (phys_state(c)%q(:ncol,:,:) - qphystend(:ncol,:,:,c))*dtinv
      do m=1,ppcnst
         call outfld(qphystendnam(m),qphystend(1,1,m,c),pcols   ,c   )
      enddo
      call t_stopf ('tphysac')
     ! MSP added:
     ! Note that state%t IS UPDATED via tend@dtdt at the end of tphysac.
#ifdef CLOUDBRAIN
  phys_state(c)%tap = phys_state(c)%t
  phys_state(c)%qap = phys_state(c)%q(:,:,1)
  phys_state(c)%vap = phys_state(c)%v   ! SR: Add v-wind
#endif
     aux(:ncol,:pver) = phys_state(c)%t(:ncol,:pver)
     call outfld ('TAP',aux,pcols,c)
     aux(:ncol,:pver) = phys_state(c)%q(:ncol,:pver,1)
     call outfld ('QAP',aux,pcols,c)
     aux(:ncol,:pver) = phys_state(c)%q(:ncol,:pver,2)
     call outfld ('QCAP',aux,pcols,c)
     aux(:ncol,:pver) = phys_state(c)%q(:ncol,:pver,3)
     call outfld ('QIAP',aux,pcols,c)
     aux(:ncol,:pver) = phys_state(c)%u(:ncol,:pver)
     call outfld ('UAP',aux,pcols,c)
     aux(:ncol,:pver) = phys_state(c)%v(:ncol,:pver)
     call outfld ('VAP',aux,pcols,c)
   end do                    ! Chunk loop
#ifdef PVBUDGET
    ! Mike Pritchard, store PV before physics:
   call calculate_physics_PV (phys_state,pv1,pv2,pv3)
!$OMP PARALLEL DO PRIVATE (C, NCOL, AUX)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
     aux(:ncol,:pver) = phys_state(c)%pv(:ncol,:)
     call outfld ('PVAP',aux,pcols,c)
     aux(:ncol,:pver) = (phys_state(c)%pv(:ncol,:pver)-pvtot_tmp(:ncol,:pver,c) )/ztodt
     call outfld('PVPHYSTND',aux,pcols   ,c   )
     aux(:ncol,:pver) = (pv1(:ncol,:pver,c) - pv1_tmp(:ncol,:pver,c))/ztodt
     call outfld('PV1PHYSTND',aux,pcols,c)
     aux(:ncol,:pver) = (pv2(:ncol,:pver,c) - pv2_tmp(:ncol,:pver,c))/ztodt
     call outfld('PV2PHYSTND',aux,pcols,c)
     aux(:ncol,:pver) = (pv3(:ncol,:pver,c) - pv3_tmp(:ncol,:pver,c))/ztodt
     call outfld('PV3PHYSTND',aux,pcols,c)
   end do
#endif

   call t_stopf('ac_physics')

   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()

end subroutine physpkg
