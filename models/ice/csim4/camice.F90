#include <misc.h>
#include <params.h>

subroutine camice (srf_state, srfflx)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! CAM sea ice surface fluxes.
!
! Method: 
! 
! Author:
! 
!-----------------------------------------------------------------------
!
! $Id: camice.F90,v 1.1.4.7 2003/12/12 22:55:01 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, begchunk,endchunk
  use pmgrid,          only: masterproc
  use ice_constants,   only: Tffresh
!JR Renamed Fhnet to Focn
  use comsrf,          only: surface_state, srfflx_parm, icefrac, snowhice, sicthk, &
                             tsice, asdirice, asdifice, aldirice, aldifice, &
                             tsocn, Focn, frzmlt, aice, landfrac, ocnfrac, lwupice, &
                             update_ocnice
  use phys_grid,       only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use time_manager,    only: get_step_size, get_curr_calday, get_nstep
  use ice_dh,          only: prognostic_icesnow
  use ice_globalcalcs, only: sea_ice_energy, gmean_ice

  implicit none
!
! Input/Output arguments
!
  type(surface_state), intent(inout) :: srf_state(begchunk:endchunk)
  type(srfflx_parm), intent(inout)   :: srfflx(begchunk:endchunk)

#include <comctl.h>

!---------------------------Local variables-----------------------------
  integer :: dtime          ! timestep size [seconds]
  real(r8) rtime            ! calendar day for next timestep
  real(r8) lats(pcols)      ! 
  real(r8) lons(pcols)      ! 
  real(r8) cdaynext         ! calendar day for next timestep
  real(r8) cosznext(pcols)  ! cosine solar zenith angle next timestep

  integer ncol              ! number of columns in chunk
  integer c                 ! chunk index
  integer i                 ! temporary variables

  real(r8) :: snowfall(pcols,begchunk:endchunk)  ! total snowfall rate
  real(r8) :: fsns(pcols,begchunk:endchunk)      ! SW absorbed in ice
  real(r8) :: evap(pcols,begchunk:endchunk)      ! evaporative flux off snow and ice
  real(r8) :: aiceinit(pcols,begchunk:endchunk)
  real(r8) :: gmsie_in          ! global mean sea ice/snow internal energy on input
  real(r8) :: gmsie_out         ! global mean sea ice/snow internal energy on output
  real(r8) :: sier              ! sea ice energy rate
  real(r8) :: fluxsum           ! F_ice - F_ocn - F_frzmlt
  real(r8) :: F_ice
  real(r8) :: F_ocn
  real(r8) :: F_frzmlt
  real(r8) :: errterm           ! error term in global energy calculations
  real(r8) :: deltae(pcols,begchunk:endchunk)    ! change in energy
  real(r8) :: deltaaice(pcols,begchunk:endchunk) ! change in aice

  logical :: ice_conschk        ! whether to apply global energy conservation check

!-----------------------------------------------------------------------
!
! Calendar day for next time step
!
  call t_startf ('camice_st')
  dtime = get_step_size()
  rtime = dtime
  cdaynext = get_curr_calday(offset=dtime)
!
! set up snowfall here so it doesn't have to be private in the omp call
!
  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
     do i = 1,ncol
	if (prognostic_icesnow) then
           snowfall(i,c) = srf_state(c)%precsc(i) + srf_state(c)%precsl(i)
        else
           snowfall(i,c) = 0.	
        end if
     end do
  end do
  call t_stopf ('camice_st')

!BPB Compute initial sea ice/snow internal energy and globally average over
!BPB ocean.  gmsie_in = global mean sea ice/snow internal energy on input

#ifdef COUP_SOM
  ice_conschk = .false.
  if (ice_conschk_frq > 0) then
     ice_conschk = mod (get_nstep(), ice_conschk_frq) == 0
  end if
  if (ice_conschk) then
     gmsie_in = sea_ice_energy (srf_state, dtime, 1, deltae, deltaaice, snowfall)
  end if
#endif

!$OMP PARALLEL DO PRIVATE (C, NCOL, LATS, LONS, COSZNEXT,I)

  do c=begchunk,endchunk
     ncol = get_ncols_p(c)

! Sea ice surface fluxes and temperatures

     call seaice (c, ncol, rtime, aice(1,c), tsice(1,c), &
                  sicthk(1,c), snowhice(1,c), srf_state(c)%ubot, &
                     srf_state(c)%vbot, srf_state(c)%tbot, &
                  srf_state(c)%qbot, srf_state(c)%thbot, srf_state(c)%zbot, &
                     srf_state(c)%pbot ,srf_state(c)%flwds, &
                  srf_state(c)%sols, srf_state(c)%soll, srf_state(c)%solsd, &
                     srf_state(c)%solld, asdirice(1,c), &
                  aldirice(1,c), asdifice(1,c), aldifice(1,c), &
    	             snowfall(1,c), &
                  tsocn(1,c), frzmlt(1,c), Focn(1,c),  &
                  srf_state(c)%tssub, srfflx(c)%cflx, srfflx(c)%wsx, srfflx(c)%wsy, &
                     srfflx(c)%ts, srfflx(c)%shf, &
         	  srfflx(c)%lhf, srfflx(c)%lwup, srfflx(c)%tref, fsns(1,c), evap(1,c), &
                  aiceinit(1,c))
!
! Make all surface fractions consistent, and update physics buffer
!
     call update_ocnice (c)
!
! Albedos for next time step 
!
! Note the total albedo here that is returned to the atmosphere 
! model is based on a weighted sum of the albedo over ice and ocean
! using fractional areas from this time step. The absorbed shortwave over
! sea ice in the next step uses ice albedos that are saved at there
! present value but with a NEW fractional area that is input prior to 
! the next time through the sea ice model.  Hence
! there is a time step mismatch in the absorbed solar over sea ice. 
! CCSM would not allow such a thing, but here we are specifying sst, 
! over the ocean fraction anyway so it doesn't really matter. 

     call get_rlat_all_p(c, ncol, lats)
     call get_rlon_all_p(c, ncol, lons)
     call zenith (cdaynext, lats, lons, cosznext, ncol)
     call albice(c,ncol, &
                 srf_state(c)%tbot,snowhice(1,c),cosznext, &
                 srfflx(c)%asdir, srfflx(c)%aldir, &
                 srfflx(c)%asdif, srfflx(c)%aldif)
!
! save off ice albedos for sea ice routine per email Bitz.
! I should note that I made one change to the "physics" from John's
! original fracice implementation. John had the absorbed solar by the
! sea ice equal to the gridcell average.  This is pretty far off when
! the sea ice fraction is small. I realize that it is standard practise
! in many models, but it doesn't have to be.  Therefore I have compute a
! special srfrad over ice and I send the ice albedos to the restart
! file.
!
     do i = 1,ncol
        asdirice(i,c)=srfflx(c)%asdir(i)
        aldirice(i,c)=srfflx(c)%aldir(i)
        asdifice(i,c)=srfflx(c)%asdif(i)
        aldifice(i,c)=srfflx(c)%aldif(i)
        lwupice(i,c) =srfflx(c)%lwup(i)
     end do
  end do

#ifdef COUP_SOM

  if (ice_conschk) then

!BPB Check sea ice energy conservation: sier = sea ice energy rate

     gmsie_out = sea_ice_energy (srf_state, dtime, 2, deltae, deltaaice, snowfall)
     sier = (gmsie_out - gmsie_in) / dtime

!BPB Compute global mean F_ice, F_frzmlt and F_ocn
!BPB Then compute output sea ice/snow internal energy and globally
!BPB average over ocean.  gmsie_out = global mean sea ice/snow
!BPB internal energy on output

     call gmean_ice (srfflx, srf_state, fsns, snowfall, F_ice, &
                     F_ocn, F_frzmlt, dtime, deltae, evap, aiceinit)

     fluxsum = F_ice - F_ocn - F_frzmlt
     errterm = abs (sier/fluxsum - 1.)
     if (masterproc) then
        write(6,'(a,1p,4e12.4)') 'CAMICE: sier,fluxsum,diff,errterm=', &
                                 sier, fluxsum, fluxsum-sier, errterm
                   
        write(6,'(a,1p,5e12.4)') '        gmsie_in,gmsie_out,F_ice,F_ocn,F_frzmlt=', &
                                 gmsie_in, gmsie_out, F_ice, F_ocn, F_frzmlt
     end if
  end if

#endif

  return
end subroutine camice
