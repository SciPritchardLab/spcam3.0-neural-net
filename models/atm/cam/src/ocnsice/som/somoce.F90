#include <misc.h>
#include <params.h>

subroutine somoce (srf_state, srfflx)

!-----------------------------------------------------------------------
!
! Purpose:
! CAM Slab Ocean
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, begchunk, endchunk
  use comsrf,       only: surface_state, srfflx_parm, tsocn, Focn, frzmlt, &
                          aice,landfrac, asdirocn, aldirocn, asdifocn, aldifocn, &
                          sicthk,ocnfrac,lwupocn, icefrac
        
  use phys_grid,    only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use ocean_data,   only: qflux, mld
  use time_manager, only: get_step_size, get_curr_calday, get_nstep, is_end_curr_day
  use ice_constants,only: Tffresh, rhow, cp_ocn
  use mixed_layer_globalcalcs
  use history,      only: outfld

  implicit none
!
! Input/Output arguments
!
   type(surface_state), intent(inout) :: srf_state(begchunk:endchunk)
   type(srfflx_parm),   intent(inout) :: srfflx(begchunk:endchunk)

#include <comctl.h>

!---------------------------Local variables-----------------------------
  integer :: dtime          ! timestep size [seconds]
  real(r8) rtime            ! calendar day for next timestep
  real(r8) cdaynext         ! calendar day for next timestep
  real(r8) clat(pcols)      ! current latitudes(radians)
  real(r8) clon(pcols)      ! current longitudes(radians)
  real(r8) cosznext(pcols)  ! cosine solar zenith angle next timestep
  real(r8) qfluxgmin        ! initial qflux global mean
  real(r8) qfluxgmot        ! qflux global mean output after adjustment for cold ocean
  real(r8) qdiff            ! global mean diff qfluxgmin and qfluxgmon; also warm ocn area adjusted
  real(r8) :: qflux_adj(pcols,begchunk:endchunk) ! adjusted qflux
  real(r8) :: oie1(pcols,begchunk:endchunk)      ! ocean internal energy (timestep start)
  real(r8) :: oie2(pcols,begchunk:endchunk)      ! ocean internal energy (timestep end)
  real(r8) :: ftot(pcols,begchunk:endchunk)      ! net ocean flux
  real(r8) :: onf(pcols,begchunk:endchunk)       ! ocean net flux (sw + lw + lh + sh)
  real(r8) :: frzmltmax(pcols,begchunk:endchunk)
  real(r8) :: oierate(pcols,begchunk:endchunk)
  real(r8) :: wght
  real(r8) :: wghtw         ! global warm ocean area weight
  real(r8) :: oiegmin       ! global mean initial ocean internal energy
  real(r8) :: oiegmot       ! global mean final ocean internal energy
  real(r8) :: oierate_g     ! global change rate of ocean internal energy
  real(r8) :: ftotgmot      ! total ocean flux (rhs of energy eq. for ocean)
  real(r8) :: Focngmot
  real(r8) :: frzmltgmot
  real(r8) :: tsgmot

  integer ncol              ! number of columns in chunk
  integer c                 ! chunk index
  integer i

  logical :: som_conschk    ! Whether to apply global energy conservation check

  real(r8), parameter :: cpw  = cp_ocn*rhow  ! Specific heat of ocean water (J/m3/C)
!-----------------------------------------------------------------------
!
! Calendar day for next time step
!
  dtime = get_step_size()
  rtime = dtime
  cdaynext = get_curr_calday(offset=dtime)

  som_conschk = .false.
  if (som_conschk_frq > 0) then
     som_conschk = mod (get_nstep(), som_conschk_frq) == 0
  end if

!$OMP PARALLEL DO PRIVATE (I, C, NCOL)

  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
     do i=1,ncol
        oie1(i,c) = cpw*mld(i,c)*tsocn(i,c)
     end do

#ifdef BFBDEV45
!
! This is a roundoff change to get bfb vs. dev51
!
     do i=1,ncol
        if (landfrac(i,c) < 1.) then
           aice(i,c) = aice(i,c)*(1. - landfrac(i,c))
           aice(i,c) = aice(i,c)/(1. - landfrac(i,c))
        end if
     end do
#endif
!
! Ocean surface fluxes and temperatures
!JR Reordered frzmlt, Focn order to be what I think is correct
!
     call mixed_layer1 (c, ncol, rtime, landfrac(1,c), aice(1,c),                     &
                        srf_state(c)%ubot, srf_state(c)%vbot, srf_state(c)%tbot,      &
                           srf_state(c)%qbot, srf_state(c)%thbot,                     &
                        srf_state(c)%zbot, srf_state(c)%pbot, srf_state(c)%flwds,     &
                           srf_state(c)%sols, srf_state(c)%soll,                      &
                        srf_state(c)%solsd, srf_state(c)%solld, asdirocn(1,c),        &
                           aldirocn(1,c), asdifocn(1,c),                              &
                        aldifocn(1,c), mld(1,c), qflux(1,c), tsocn(1,c), frzmlt(1,c), &
                        Focn(1,c), srfflx(c)%cflx, srfflx(c)%wsx, srfflx(c)%wsy,      &
                           srfflx(c)%shf,                                             &
                        srfflx(c)%lhf, srfflx(c)%lwup, srfflx(c)%tref,                &
                           qflux_adj(1,c), sicthk(1,c),                               &
                        onf(1,c))
  end do

!BPB (Extrae Comment): Ensure qflux energy conservation by giving/taking
!BPB heat from warm ocean (i.e. sst > 0C) to balance that required by
!BPB changes to qflux under cold ocean (i.e. sst < 0C).
!
!BPB Ensure energy conservation of qflux:
!BPB
!BPB Compute output global mean of the adjusted qflux:
!BPB Also, compute area of sst > 0

!
! Compute global mean input qflux
!
  qfluxgmin = gmean (qflux, wght)
  qfluxgmot = gmean_warm (qflux_adj, wght, wghtw)

!BPB Qdifference:
#ifdef ZEROADJ
  qdiff = 0.
#else
  qdiff = qfluxgmin - qfluxgmot
#endif

!BPB Renormalize for area of warm ocean:
  qdiff = qdiff * wght/wghtw

!$OMP PARALLEL DO PRIVATE (I, C, NCOL, CLAT, CLON, COSZNEXT)

  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
!
! Ocean surface fluxes and temperatures
!JR Reordered frzmlt, Focn order to be what I think is correct
!
     call mixed_layer2 (c, ncol, rtime, landfrac(1,c), mld(1,c), &
                        tsocn(1,c), srfflx(c)%ts, qdiff, qflux_adj(1,c))
!
! Albedos for next time step
!
     call get_rlat_all_p(c, ncol, clat)
     call get_rlon_all_p(c, ncol, clon)
     call zenith (cdaynext, clat, clon, cosznext, ncol)
     call albocean (c,ncol, &
                    cosznext, srfflx(c)%asdir, srfflx(c)%aldir, &
                    srfflx(c)%asdif, srfflx(c)%aldif)
!
! Save off ocean albedos for next timestep
!
     do i = 1,ncol
        asdirocn(i,c) = srfflx(c)%asdir(i)
        aldirocn(i,c) = srfflx(c)%aldir(i)
        asdifocn(i,c) = srfflx(c)%asdif(i)
        aldifocn(i,c) = srfflx(c)%aldif(i)
        lwupocn(i,c) = srfflx(c)%lwup(i)
     end do
!
! Compute total internal ocean energy and total flux for conservation check.
! Note: include frzmlt in calculation only when positive because
! positive frzmlt means heating mixed layer back to freezing temperature, while
! negative frzmlt means mixed layer temperature above freezing.
!
     do i=1,ncol
        oie2(i,c)      = cpw*mld(i,c)*tsocn(i,c)
        oierate(i,c)   = (oie2(i,c) - oie1(i,c))/dtime
     end do

     call outfld ('OIE     ', oie2(1,c),   pcols, c)
     call outfld ('OIERATE ', oierate(1,c), pcols, c)

     if (som_conschk) then
        do i=1,ncol
           ftot(i,c)      = onf(i,c) + Focn(i,c) + max (frzmlt(i,c), 0._r8) - qflux_adj(i,c)
           frzmltmax(i,c) = max (frzmlt(i,c), 0._r8)
        end do
     end if
  end do

  if (som_conschk) then
     oiegmin    = gmean (oie1, wght)
     Focngmot   = gmean (Focn, wght)
     frzmltgmot = gmean (frzmltmax, wght)
     oiegmot    = gmean (oie2, wght)
     ftotgmot   = gmean (ftot, wght)

     if (is_end_curr_day ()) then
        tsgmot  = gmean (tsocn, wght)
        if (masterproc ) then
           write(6,*)'TSGMOT=',tsgmot
        end if
     end if

     oierate_g = (oiegmot-oiegmin)/dtime
     if (masterproc) then
        write(6,'(a,1p,4e16.8)') 'SOMOCE: ftotgmot,oierate,oiegmot,oiegmin=', &
                                 ftotgmot, oierate_g, oiegmot, oiegmin
        write(6,'(a,1p,2e16.8)') '        Focngmot, frzmltgmot=', Focngmot, frzmltgmot
     end if

     call check_conservation (qflux_adj, qfluxgmin, wght)
  end if

  return
end subroutine somoce
