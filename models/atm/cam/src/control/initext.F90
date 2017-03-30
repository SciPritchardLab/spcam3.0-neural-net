#include <misc.h>
#include <params.h>

subroutine initext
!----------------------------------------------------------------------- 
! 
! Purpose: Initialize external models and/or boundary dataset information
! 
! Method: 
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use ppgrid, only: begchunk, endchunk
   use phys_grid, only: get_ncols_p, get_rlat_all_p, get_rlon_all_p,get_lat_all_p, get_lon_all_p
   use pspect
   use comsrf
   use runtime_opts, only: aqua_AndKua, aqua_uniform, aqua_uniform_sst_degC

   use rgrid
   use shr_orb_mod
   use ioFileMod
   use so4bnd
   use commap
#if ( ! defined COUP_CSM )
   use ice_constants, only: Tffresh, init_constants
   use shr_const_mod, only: shr_const_tkfrz
#endif
   use filenames, only: bndtvo, bndtvs
   use physconst, only: stebol
   use time_manager, only: is_first_step, get_curr_calday, get_curr_date, &
                           is_perpetual, get_perp_date
#if ( defined SPMD )
   use mpishorthand
#endif

#if (defined COUP_CSM)
   use ccsm_msg, only: ccsmini
#else
   use atm_lndMod, only: atmlnd_ini
#if ( ! defined COUP_SOM )
   use sst_data, only: sstini, sstint, sstan, sst
   use ice_data, only: iceini, iceint
   use atm_lndMod
#endif
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!--------------------------Local Variables------------------------------
!
#if (!defined COUP_CSM)
   integer  :: i,k,c        ! indices
   integer  :: ncol               ! number of columns in current chunk
   real(r8) :: tfreeze            ! freezing temp of fresh water
   real(r8) :: coszrs(pcols)      ! Cosine solar zenith angle
   real(r8) :: clat1(pcols)       ! Current latitude(radians)
   real(r8) :: clon1(pcols)       ! Current longitude(radians)
   integer  :: sghid              ! NetCDF sgh field id 
   logical  :: oro_hires          ! true => ORO came from high res topo file
   logical  :: log_print          ! Flag to print out log information or not
   integer  :: ret                ! NetCDF returned status 
   integer  :: attlen             ! NetCDF attribute length
   character(len=256) :: text     ! NetCDF attribute
#endif
   character(len=256) :: locfn    ! netcdf local filename to open 
   character*4 ncnam(5)
   integer  :: yr, mon, day, tod  ! components of a date
   real(r8) :: calday             ! current calendar day
   real(r8) :: tsice_tmp (pcols,begchunk:endchunk)
   integer  :: lchnk
   integer  :: lats(pcols)
   integer  :: lons(pcols)
   real(r8) :: tssav(pcols,begchunk:endchunk) ! cam surface temperatures
!
!-----------------------------------------------------------------------

   calday = get_curr_calday()
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Obtain datasets
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Obtain time-variant ozone and sst datatsets and do initial read of
! ozone dataset
!
   if (.not. ideal_phys) then
      if (masterproc) then
         call getfil (bndtvo, locfn)
         call wrap_open (locfn, 0, ncid_oz)
         write(6,*)'INITEXT: NCOPN returns id ',ncid_oz,' for file ',trim(locfn)
      endif
#if ( ! defined COUP_CSM )
      if (.not. aqua_planet) then
         if (masterproc) then
            call getfil(bndtvs, locfn)
            call wrap_open(locfn, 0, ncid_sst)
            write(6,*)'INITEXT: NCOPN returns id ',ncid_sst,' for file ',trim(locfn)
         endif
      endif
#endif
      call oznini ()
   endif
!
! Obtain sulfate aerosol datasets
!
   if ( doRamp_so4 ) then
      call sulfini
   end if

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Preprocessing if -- If NOT coupled to the Climate System Model (CSM)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( ! defined COUP_CSM )
!
! Determine if SGH field came from hi-res dataset
!
   if (is_first_step()) then
      if (masterproc) then
         call wrap_inq_varid (ncid_ini, 'SGH', sghid)
         ret = nf_inq_attlen (ncid_ini, sghid, 'from_hires', attlen)
         if (ret == nf_noerr .and. attlen > 256) then
            write(6,*)'INITEXT: Att length of from_hires is too long'
            call endrun
         end if
         ret = nf_get_att_text (ncid_ini, sghid, 'from_hires', text)
         if (ret == nf_noerr .and. text(1:4) == 'true')then
            oro_hires = .true.
            write(6,*)'INITEXT: attribute from_hires is true.'
            write(6,*)'         Will use tssub values to guess sea ice'
         else
            oro_hires = .false.
            write(6,*)'INITEXT: attribute from_hires is either false or not present.'
            write(6,*)'         Where sea ice exists, its initial temperature will be just below freezing'
         end if
      end if
#if ( defined SPMD )
      call mpibcast (oro_hires, 1, mpilog, 0, mpicom)
#endif
   end if
!
! Setup the characteristics of the orbit
! (Based on the namelist parameters)
!
   if (masterproc) then
      log_print = .true.
   else
      log_print = .false.
   end if
   call shr_orb_params (iyear_AD, eccen , obliq , mvelp, obliqr, &
                        lambm0, mvelpp, log_print)
!
! Initialize land model. This involves initializing land 
! albedos, surface temperature, lwup and snowh.  NOTE: On restart, 
! lwup, ts, albedos and snowh, come from the atm restart data.  
!
   if (is_first_step()) then
      do c=begchunk,endchunk
         call srfflx_state_reset (srfflx_state2d(c))
      end do
   end if
   if (.not. adiabatic .and. .not. ideal_phys .and. .not. aqua_planet) then
      call atmlnd_ini(srfflx_parm2d)

      if (is_first_step()) then
!
! save off albedos and longwave for som offline vars on restart read these from
! restart datasets
!
         do lchnk=begchunk,endchunk
            ncol = get_ncols_p(lchnk)
            do i=1,ncol
               if (landfrac(i,lchnk) > 0.) then
                  asdirlnd(i,lchnk) = srfflx_parm2d(lchnk)%asdir(i)
                  asdiflnd(i,lchnk) = srfflx_parm2d(lchnk)%asdif(i)
                  aldirlnd(i,lchnk) = srfflx_parm2d(lchnk)%aldir(i)
                  aldiflnd(i,lchnk) = srfflx_parm2d(lchnk)%aldif(i)
                  lwuplnd(i,lchnk) = srfflx_parm2d(lchnk)%lwup(i)
               else
                  asdirlnd(i,lchnk) = 0. 
                  asdiflnd(i,lchnk) =  0. 
                  aldirlnd(i,lchnk) =  0. 
                  aldiflnd(i,lchnk) =  0. 
                  lwuplnd(i,lchnk) =  0. 
               endif
            end do
         end do

      endif
   endif
!
! Save off ts here because its needed below for calculating lwup.  The 
! updated state values are summed according to the fractional area of the 
! underlying surface and only represent an actual grid box value after the 
! last surface process has been called. TS is a special case as it is 
! calculated from lwup and overwritten each time update surface
! fluxes is called.  The intermediate ts values returned from the update 
! routine are wrong until ts is calculated from the full gridbox value of lwup
! lwup is only complete after the last surface process is called or if we
! are calculating a grid point that is all land.  We will save the 
! intermediate value returned from land for our calculations below.
!
   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      tssav(:,c) = srfflx_parm2d(c)%ts(:)
      call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), landfrac(1,c), ncol)
   end do	

#if ( defined COUP_SOM )
!
! Slab ocean model: set initial surf temps for initial run. Read in 2 time slices of
! mixed layer depths and q fluxes from boundary dataset whether initial or restart
!
   call init_constants ()
   call somini ()
!
! Initialize surface and sub-surface temperatures, set new 
! new sea ice concentrations and compute longwave up over non-land
!
   if (is_first_step()) then
      write(6,*)'INITEXT: ocnfrac=',ocnfrac(1,begchunk)
      do lchnk=begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         call verify_fractions (lchnk, ncol)

         if (.not. adiabatic .and. .not. ideal_phys) then
         do i=1,ncol
            frzmlt(i,lchnk) = 0.            ! set initial freeze/melt potential to zero
            Focn(i,lchnk) = 0.
            srfflx_state2d(lchnk)%ts(i) = &
                 landfrac(i,lchnk)*srfflx_state2d(lchnk)%ts(i) + &
                 icefrac(i,lchnk)*tsice(i,lchnk) + &
                 ocnfrac(i,lchnk)*(tsocn(i,lchnk)+Tffresh)
            if (landfrac(i,lchnk).ne.1.) then
	            srfflx_state2d(lchnk)%lwup(i) = &
                         stebol*(srfflx_state2d(lchnk)%ts(i)**4)
                    srfflx_state2d(lchnk)%sst(i) = tsocn(i,lchnk) + Tffresh
            end if

            lwuplnd(i,lchnk) = 0.
            lwupice(i,lchnk) = 0.
            lwupocn(i,lchnk) = 0.
            if (landfrac(i,lchnk) > 0.) then
               lwuplnd(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
            end if
            if (icefrac(i,lchnk) > 0.) then
               lwupice(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
            end if
            if (ocnfrac(i,lchnk) > 0.) then
               lwupocn(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
            end if
         end do
         end if
      end do
   end if

#else
!
! Data ocean model: Initialize ocean/sea-ice surface datasets and determine initial sea surface 
! temperature 
! 
   if (adiabatic .or. ideal_phys) then
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         aice(:ncol,c) = 0.
         call update_ocnice (c)
      end do
   else
      call sstini
      call iceini
!
! on a restart the ice values come from the restart data sets 
!
      if (is_first_step()) then
!         call sstint(prev_timestep=.true.)
         call sstint (.true.,aqua_uniform, aqua_AndKua, aqua_uniform_sst_degC)
         tsice_tmp(:,:) = tsice(:,:)
         call iceint(prev_timestep=.true., startup=.true.)
         tsice(:,:)     = tsice_tmp(:,:)
      else
         call sstint (.false.,aqua_uniform, aqua_AndKua, aqua_uniform_sst_degC)
         do lchnk=begchunk,endchunk
            ncol = get_ncols_p(lchnk)
            do i=1,ncol
               previcefrac(i,lchnk) = icefrac(i,lchnk)
            end do
         end do
      endif
   end if
!
! Initialize surface and sub-surface temperatures, set new 
! new sea ice concentrations and compute longwave up over non-land
!
   if (is_first_step()) then
      tfreeze = shr_const_tkfrz
      do lchnk=begchunk,endchunk
         if (.not. adiabatic .and. .not. ideal_phys) then
            ncol = get_ncols_p(lchnk)
            do i=1,ncol
               srfflx_state2d(lchnk)%lwup(i) = &
                    srfflx_state2d(lchnk)%lwup(i)                      + &
                    icefrac(i,lchnk)*stebol*tsice_rad(i,lchnk)     **4 + &
                    ocnfrac(i,lchnk)*stebol* (sst(i,lchnk)+tfreeze)**4
               if (landfrac(i,lchnk).ne.1.) then
                  srfflx_state2d(lchnk)%ts(i) = &
                       sqrt(sqrt(srfflx_state2d(lchnk)%lwup(i)/stebol))
                  srfflx_state2d(lchnk)%sst(i) = sst(i,lchnk) + tfreeze
               end if
            end do
            do i=1,ncol
               lwuplnd(i,lchnk) = 0.
               lwupice(i,lchnk) = 0.
               lwupocn(i,lchnk) = 0.
               if (landfrac(i,lchnk) > 0.) then
                  lwuplnd(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
               end if
               if (icefrac(i,lchnk) > 0.) then
                  lwupice(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
               end if
               if (ocnfrac(i,lchnk) > 0.) then
                  lwupocn(i,lchnk) = srfflx_state2d(lchnk)%lwup(i)
               end if
           end do
         end if
      end do
   end if
#endif
!
! Initialize non-land albedos at NSTEP = 0.  At NSTEP = 1 and 
! beyond, albedos will be computed for the *next* timestep to 
! accomodate coupling with a single interface.
!
   if (is_first_step()) then
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         call get_rlat_all_p(c, ncol, clat1)
         call get_rlon_all_p(c, ncol, clon1)
         call zenith (calday, clat1, clon1, coszrs, ncol)
         call albocean (c, ncol, coszrs, &
                        srfflx_parm2d(c)%asdir, srfflx_parm2d(c)%aldir, &
                        srfflx_parm2d(c)%asdif, srfflx_parm2d(c)%aldif)
!
! fill in ocn albedoes for mixed layer model
!
         asdirocn(:ncol,c)= srfflx_parm2d(c)%asdir(:ncol)
         aldirocn(:ncol,c)= srfflx_parm2d(c)%aldir(:ncol)
         asdifocn(:ncol,c)= srfflx_parm2d(c)%asdif(:ncol)
         aldifocn(:ncol,c)= srfflx_parm2d(c)%aldif(:ncol)
         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), ocnfrac(1,c), ncol)
      end do

      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         call get_lat_all_p(c, ncol, lats)
         call get_lon_all_p(c, ncol, lons)
         call get_rlat_all_p(c, ncol, clat1)
         call get_rlon_all_p(c, ncol, clon1)
         call zenith (calday, clat1, clon1, coszrs, ncol)
         call albice(c,ncol, surface_state2d(c)%tbot, snowhice(1,c), coszrs, &
              srfflx_parm2d(c)%asdir, &
              srfflx_parm2d(c)%aldir, srfflx_parm2d(c)%asdif, &
              srfflx_parm2d(c)%aldif)
!
! fill in ice albedoes for therm ice model
!
         asdirice(:ncol,c)= srfflx_parm2d(c)%asdir(:ncol)
         aldirice(:ncol,c)= srfflx_parm2d(c)%aldir(:ncol)
         asdifice(:ncol,c)= srfflx_parm2d(c)%asdif(:ncol)
         aldifice(:ncol,c)= srfflx_parm2d(c)%aldif(:ncol)

         call update_srf_fluxes (srfflx_state2d(c), srfflx_parm2d(c), icefrac(1,c), ncol)
      end do
   end if
#endif
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Preprocessing if -- if coupled to (CSM)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( defined COUP_CSM )
!
! Initial communications with coupler
!
   call ccsmini
#endif
   return
end subroutine initext
