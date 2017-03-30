#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: driver
!
! !INTERFACE: 
subroutine driver (doalb, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION: 
! Main CLM2 driver calling sequence:
! -> loop over subgrid hierarchy, calling for each column:
!    -> DriverInit          save of variables from previous time step
!    -> Hydrology1          canopy interception and precip on ground
!    -> Biogeophysics1      leaf temperature and surface fluxes
!    -> Biogeophysics_Lake  lake temperature and surface fluxes
!    -> Biogeophysics2      soil/snow and ground temp and update surface fluxes
!    -> Hydrology2          surface and soil hydrology
!    -> Hydrology_Lake      lake hydrology
!    -> Biogeochemistry     surface biogeochemical fluxes (LSM)
!    -> EcosystemDyn:       ecosystem dynamics: phenology, vegetation, soil carbon
!    -> SurfaceAlbedo:      albedos for next time step
!    -> BalanceCheck        check for errors in energy and water balances
!  -> write_diagnostic      output diagnostic if appropriate
!  -> Rtmriverflux          calls RTM river routing model
!  -> update_hbuf           accumulate history fields over history time interval
!  -> htapes_wrapup         write history tapes if appropriate
!  -> restart               write restart file if appropriate
!  -> inicwrt               write initial file if appropriate 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clmpoint
  use globals
  use spmdMod, only : masterproc 
  use clm_varcon, only : zlnd		
  use time_manager, only : get_step_size, get_curr_calday, get_curr_date, get_nstep
  use histFileMod, only : update_hbuf, htapes_wrapup 
  use restFileMod, only : restart
  use inicFileMod, only : inicwrt, do_inicwrite 
  use DriverInitMod, only : DriverInit
  use BalanceCheckMod, only : BalanceCheck
  use Hydrology1Mod, only : Hydrology1
  use Hydrology2Mod, only : Hydrology2
  use HydrologyLakeMod, only : HydrologyLake
  use Biogeophysics1Mod, only : Biogeophysics1
  use Biogeophysics2Mod, only : Biogeophysics2
  use BiogeophysicsLakeMod, only : BiogeophysicsLake
  use SurfaceAlbedoMod, only : SurfaceAlbedo, Snowage
  use pft2columnMod, only : pft_to_col_glob, pft_to_col_soil
#if (defined COUP_CSM)
  use clm_varctl, only : wrtdia, fsurdat, csm_doflxave
#else
  use clm_varctl, only : wrtdia, fsurdat
#endif
#if (defined DGVM)
   use EcosystemDynDGVMMod, only : EcosystemDyn
   use DGVMMod, only : lpj, lpjreset1, lpjreset2
#else
  use EcosystemDynMod, only : EcosystemDyn
  use mvegFileMod, only : interpmonthlyveg
#endif
#if (defined BGC)
  use BiogeochemistryMod, only : Biogeochemistry
  use DustMod, only : DustDryDep
#endif
  use accFldsMod, only : updateAccFlds
#if (defined RTM)
  use RtmMod, only : Rtmriverflux
#endif
#if (defined COUP_CSM)
  use clm_csmMod, only : csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
       dorecv, dosend, csmstop_now
#endif
  use lnd2atmMod, only : lnd2atm
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: doalb  !true if time for surface albedo calculation
  real(r8), intent(in) :: eccen  !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0 !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp !Earth's moving vernal equinox long. of perihelion + pi (radians)
!
! !CALLED FROM:
! program program_off (if COUP_OFFLINE cpp variable is defined)
! program program_csm (if COUP_CSM cpp variable is defined)
! subroutine atm_lnddrv in module atm_lndMod (if COUP_CAM cpp variable is defined)
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
! 2003.07.09 Peter Thornton: Reversion to old fraction snow cover algorithm
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: ci                    !indices
  integer  :: yrp1                  !year (0, ...) for nstep+1
  integer  :: monp1                 !month (1, ..., 12) for nstep+1
  integer  :: dayp1                 !day of month (1, ..., 31) for nstep+1
  integer  :: secp1                 !seconds into current date for nstep+1   
  real(r8) :: caldayp1              !calendar day for nstep+1
  logical, external :: do_restwrite !determine if time to write restart   
  type(column_type), pointer :: c   !local pointer to derived subtype
!-----------------------------------------------------------------------

  call t_startf('clm_driver')

  ! ----------------------------------------------------------------------
  ! Calendar information for current time step
  ! ----------------------------------------------------------------------

  nstep = get_nstep() 
  call get_curr_date (ctl%year, ctl%month, ctl%day, ctl%secs)

  ! ----------------------------------------------------------------------
  ! Coupler receive
  ! ----------------------------------------------------------------------

#if (defined COUP_CSM)
  ! Determine if information should be sent/received to/from flux coupler
  call csm_dosndrcv(doalb)

  ! Get atmospheric state and fluxes from flux coupler
  if (dorecv) then
     call csm_recv()
     if (csmstop_now) then
        call t_stopf('clm_driver')
        RETURN
     endif
  endif
#endif

  ! ----------------------------------------------------------------------
  ! Calendar information for next time step
  ! o caldayp1 = calendar day (1.00 -> 365.99) for cosine solar zenith angle
  !              calday is based on Greenwich time
  ! o monp1    = month (1 -> 12) for leaf area index and stem area index
  ! o dayp1    = day (1 -> 31)   for leaf area index and stem area index
  ! ----------------------------------------------------------------------

  dtime = get_step_size()
  caldayp1 = get_curr_calday(offset=int(dtime))
  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))

#if (!defined DGVM)
  ! ----------------------------------------------------------------------
  ! Determine weights for time interpolation of monthly vegetation data.
  ! This also determines whether it is time to read new monthly vegetation and
  ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
  ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
  ! weights obtained here are used in subroutine ecosystemdyn to obtain time
  ! interpolated values.
  ! ----------------------------------------------------------------------

  if (doalb) call interpMonthlyVeg (fsurdat, monp1, dayp1)
#endif

  ! ----------------------------------------------------------------------
  ! LOOP 1
  ! ----------------------------------------------------------------------

  call t_startf('clm_loop1')

!$OMP PARALLEL DO PRIVATE (ci,c)
  do ci = cols1d%beg,cols1d%end
     c => cpoint(ci)%c

     ! Initialize variables from previous time step
     call DriverInit(c)

     if (.not. c%cps%lps%lakpoi) then

        ! Determine canopy interception and precipitation onto ground surface.
        ! Determine the fraction of foliage covered by water and the fraction
        ! of foliage that is dry and transpiring. Initialize snow layer if the
        ! snow accumulation exceeds 10 mm.
        ! Hydrology1() includes a loop through all the pfts for this column
        call Hydrology1(c)

        ! Determine leaf temperature and surface fluxes based on ground
        ! temperature from previous time step.
        call Biogeophysics1(c)

     else if (c%cps%lps%lakpoi) then

        ! Determine lake temperature and surface fluxes
        call BiogeophysicsLake(c)

#if (defined BGC)
        ! Dust dry deposition routine (C. Zender's modified codes)
        call DustDryDep(c)
#endif
     endif

     if (.not. c%cps%lps%lakpoi) then

#if (defined BGC)
        ! Surface biogeochemical fluxes: 
        ! co2 respiration and plant production
        call Biogeochemistry (c)

        ! Dust dry deposition routine (C. Zender's modified codes)
        call DustDryDep(c)
#endif
        ! Ecosystem dynamics: 
        ! phenology, vegetation, soil carbon, snow fraction
        call EcosystemDyn(c, doalb, .false.)

     end if

     if (doalb) then

        ! Determine albedos for next time step
        call SurfaceAlbedo(c, caldayp1, eccen, obliqr, lambm0, mvelpp)

     endif

     if (.not. c%cps%lps%lakpoi) then

        ! Determine soil/snow temperatures including ground temperature and
        ! update surface fluxes for new ground temperature.
        call Biogeophysics2(c)

        ! Perform averaging from PFT level to column level - non lake points
        call pft_to_col_soil(c)

     endif

     ! Perform averaging from PFT level to column level - lake and non-lake
     call pft_to_col_glob(c)

  end do ! end column loop

  call t_stopf('clm_loop1')

  ! ----------------------------------------------------------------------
  ! Coupler send
  ! ----------------------------------------------------------------------

#if (defined COUP_CSM)
  ! Average fluxes over interval if appropriate
  ! Surface states sent to the flux coupler states are not time averaged
  if (csm_doflxave) call csm_flxave()

  ! Send fields to flux coupler
  ! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
  ! albedos[n+1], and ocnrof_vec[n]
  if (dosend) call csm_send()
#endif

  ! ----------------------------------------------------------------------
  ! LOOP 2
  ! ----------------------------------------------------------------------

  call t_startf('clm_loop2')
!$OMP PARALLEL DO PRIVATE (ci,c)
  do ci = cols1d%beg,cols1d%end
     c => cpoint(ci)%c

     ! Vertical (column) soil and surface hydrology
     if (.not. c%cps%lps%lakpoi) call Hydrology2 (c)

     ! Lake hydrology
     if (c%cps%lps%lakpoi) call HydrologyLake (c)

     ! Update Snow Age (needed for surface albedo calculation 
     call SnowAge (c)

     ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)  
     !c%cps%frac_sno = tanh(c%cps%snowdp / (2.5 * zlnd))
	  c%cps%frac_sno = c%cps%snowdp/(10.*zlnd + c%cps%snowdp)

     ! Check the energy and water balance
     call BalanceCheck (c)

  end do   
  call t_stopf('clm_loop2')

#if (defined OFFLINE)
  ! ----------------------------------------------------------------------
  ! Determine fields that would be sent to atm for diagnostic purposes
  ! When not in offline mode, this routine is called from clm_csmMod and
  ! lp_coupling
  ! ----------------------------------------------------------------------

  call lnd2atm()
#endif

  ! ----------------------------------------------------------------------
  ! Write global average diagnostics to standard output
  ! ----------------------------------------------------------------------

  call write_diagnostic (wrtdia, nstep)

#if (defined RTM)
  ! ----------------------------------------------------------------------
  ! Route surface and subsurface runoff into rivers
  ! ----------------------------------------------------------------------

  call t_startf('clm_rtm')
  call Rtmriverflux ()
  call t_stopf('clm_rtm')
#endif

  ! Update accumulators

  call t_startf('accum')
  call updateAccFlds()
  call t_stopf('accum')

  ! Update history buffer

  call t_startf('hbuf')
  call update_hbuf()
  call t_stopf('hbuf')

#if (defined DGVM)
  ! Call DGVM (Dynamic Global Vegetation Model) LPJ at last tstep of year
  ! Then reset vegetation distribution and some counters for the next year.
  ! NOTE: lpjreset2 is called after htape_wrapup call in order to
  ! use the old weights in December's history calculations.
  ! NOTE: lpjreset1 must be called after update_accum and update_hbuf
  ! in order to have the correct values of the accumulated variables
  ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1

  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
     if (masterproc) write(6,*) 'End of year. DGVM called now: nstep =', nstep
     call lpj()
     call lpjreset1 (caldayp1, eccen, obliqr, lambm0, mvelpp)
  end if
#endif

  call t_startf('clm_output')

  ! Create history tapes if  appropriate and  write to history tapes if appropriate

  call htapes_wrapup()

#if (defined DGVM)
  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
     call lpjreset2 ()
     if (masterproc) write(6,*) 'Annual DGVM calculations are complete'
  end if
#endif

  ! Write restart files if appropriate

  if (do_restwrite()) call restart('write')

  ! Write intial files if appropriate

  if (do_inicwrite()) call inicwrt()

  call t_stopf('clm_output')

  call t_stopf('clm_driver')

  return
end subroutine driver

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
subroutine write_diagnostic (wrtdia, nstep)
!
! !DESCRIPTION: 
! Write diagnostic output.
!
! !USES:
  use clmtype
  use clmpoint
#if (defined SPMD)
  use spmdMod, only : masterproc, gather_data_to_master, mpicom
#else
  use spmdMod, only : masterproc
#endif
  use shr_sys_mod, only : shr_sys_flush
  use system_messages
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic   
  integer, intent(in) :: nstep      !model time step
!                                   
! !REVISION HISTORY:                
! Created by Mariana Vertenstein    
!                                   
!EOP                                
!                                   
! !LOCAL VARIABLES:                 
  integer  :: gi                    !index
  integer  :: begg,endg,numg        !gridcell 1d indices
  integer  :: ier                   !error status 
  real(r8) :: tsxyav                !average ts for diagnostic output
  real(r8), pointer :: rloc(:)      !temporary buffer 
  real(r8), pointer :: rglob(:)     !temporary buffer
  type(gridcell_type), pointer :: g !local pointer to derived subtype
!------------------------------------------------------------------------

  begg = grid1d%beg
  endg = grid1d%end
  numg = grid1d%num

  if (wrtdia) then

#if (defined TIMING_BARRIERS)
     call t_startf ('sync_write_diag')
     call mpi_barrier (mpicom, ier)
     call t_stopf ('sync_write_diag')
#endif
     allocate (rloc(begg:endg), stat=ier)
     if (ier/=0) then
        write(6,*)'write_diagnostic allocation error for rloc'; call endrun
     endif
     do gi = begg,endg
        g => gpoint(gi)%g
        rloc(gi) = g%l2as%t_rad
     end do
#if (defined SPMD)
     if (masterproc) then
        allocate (rglob(numg), stat=ier)
        if (ier/=0) then
           write(6,*)'write_diagnostic allocation error for rglob'; call endrun
        endif
     endif
     call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
#else
     rglob => rloc
#endif
     if (masterproc) then
        tsxyav = 0._r8
        do gi = 1,numg
           tsxyav = tsxyav + rglob(gi)
        end do
        tsxyav = tsxyav / numg
        write (6,1000) nstep, tsxyav
        call shr_sys_flush(6)
     end if
#if (defined SPMD)
     if (masterproc) deallocate(rglob)
#endif
     deallocate (rloc)

  else

#if (!defined COUP_CAM)
     if (masterproc) then
        write(6,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(6)
     endif
#endif

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',e21.15)

end subroutine write_diagnostic
