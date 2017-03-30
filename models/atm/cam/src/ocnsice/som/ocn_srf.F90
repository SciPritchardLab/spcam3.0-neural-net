#include <misc.h>
subroutine mixed_layer1 (c, ncol, dtime, landfrac, aice,    &
                         ubot, vbot, tbot, qbot, thbot,     &
                         zbot, pbot, flwds, swvdr, swidr,   &
                         swvdf, swidf, alvdr, alidr, alvdf, &
                         alidf, depth, qflux, sst, frzmlt,  &
                         Focn, cflx, taux, tauy, shflx,     &
                         lhflx, lwup, tref, qflux_adj, hi,  &
                         onf)

!---!-------------------------------------------------------------------
!---! Calculate flux exchange over open ocean and update mixed layer sst
!---!
!---!-------------------------------------------------------------------

   use shr_kind_mod,  only: r8 => shr_kind_r8
   use ppgrid,        only: pcols
   use phys_grid,     only: get_rlat_all_p
   use constituents,  only: pcnst, pnats
   use ice_constants, only: Tffresh, tfrez, rhow, cp_ocn
   use physconst,     only: latvap
   use history,       only: outfld

   implicit none

!------------------------------Arguments--------------------------------
   integer , intent(in) :: c               ! chunk index
   integer , intent(in) :: ncol            ! number of columns this chunk

   real(r8), intent(in) :: dtime           ! model timestep
   real(r8), intent(in) :: landfrac(pcols) ! land fraction
   real(r8), intent(in) :: aice(pcols)     ! ice fraction not incl. land.
   real(r8), intent(in) :: ubot(pcols)     ! Bottom level u wind
   real(r8), intent(in) :: vbot(pcols)     ! Bottom level v wind
   real(r8), intent(in) :: tbot(pcols)     ! Bottom level temperature
   real(r8), intent(in) :: qbot(pcols)     ! Bottom level specific humidity
   real(r8), intent(in) :: thbot(pcols)    ! Bottom level potential temperature
   real(r8), intent(in) :: zbot(pcols)     ! Bottom level height above surface
   real(r8), intent(in) :: pbot(pcols)     ! Bottom level pressure
   real(r8), intent(in) :: flwds(pcols)    ! net down longwave radiation at surface
   real(r8), intent(in) :: swvdr(pcols)    ! direct beam solar radiation onto srf (sw)
   real(r8), intent(in) :: swidr(pcols)    ! direct beam solar radiation onto srf (lw)
   real(r8), intent(in) :: swvdf(pcols)    ! diffuse solar radiation onto srf (sw)
   real(r8), intent(in) :: swidf(pcols)    ! diffuse solar radiation onto srf (lw)
   real(r8), intent(in) :: alvdr(pcols)    ! ocean + ice albedo: shortwave, direct
   real(r8), intent(in) :: alidr(pcols)    ! ocean + ice albedo: longwave, direct
   real(r8), intent(in) :: alvdf(pcols)    ! ocean + ice albedo: shortwave, diffuse
   real(r8), intent(in) :: alidf(pcols)    ! ocean + ice albedo: longwave, diffuse

   real(r8), intent(in) :: depth(pcols)    ! mixed layer depth (m)
!
! Declare qflux intent(inout) because it gets set to zero over land
!
   real(r8), intent(inout) :: qflux(pcols) ! mixed layer heat flux (W/m**2, positive down)
   real(r8), intent(inout) :: sst(pcols)   ! ocean layer temp (C)
   real(r8), intent(inout) :: frzmlt(pcols)! freeze/melt potential  (W/m**2) (if >0 then freeze)
   real(r8), intent(inout) :: Focn(pcols)  ! ocean-ice heat flux for basal and lateral melt (<0)

! fluxes/quantities summed over surface types

   real(r8), intent(inout):: cflx(pcols,pcnst+pnats) ! Constituent flux (kg/m2/s)
   real(r8), intent(inout):: taux(pcols)             ! X surface stress (N/m2)
   real(r8), intent(inout):: tauy(pcols)             ! Y surface stress (N/m2)
   real(r8), intent(inout):: shflx(pcols)            ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout):: lhflx(pcols)            ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout):: lwup(pcols)             ! surface longwave up flux (W/m2)
   real(r8), intent(inout):: tref(pcols)             ! 2m reference temperature
   real(r8), intent(out):: qflux_adj(pcols)        ! adjusted qflux

   real(r8), intent(in):: hi(pcols)                ! ice thickness
   real(r8), intent(out):: onf(pcols)              ! ocean net flux

!---------------------------Local variables-----------------------------
! fluxes/quantities over sea ice only
!BPB  CCSM sign convention: positive means cooling of mixed layer; negative warming
#ifdef ZEROADJ
   real(r8), parameter :: qnh = 0.
   real(r8), parameter :: qsh = 0.
#else
   real(r8), parameter :: qnh = +15.  ! extrae under ice qflux for northern hemisphere
   real(r8), parameter :: qsh = -10.  ! extrae under ice qflux for southern hemisphere
#endif
!JR Change cpw to be consistent with what was used in defineqflux (cpocn*rhoocn)
   real(r8), parameter :: cpw  = rhow*cp_ocn  ! Specific heat of ocean water (J/m3/C)

   real(r8) :: tauxocn(pcols)       ! X surface stress (N/m2)
   real(r8) :: tauyocn(pcols)       ! Y surface stress (N/m2)
   real(r8) :: shflxocn(pcols)      ! Surface sensible heat flux (J/m2/s)
   real(r8) :: lhflxocn(pcols)      ! Surface latent   heat flux (J/m2/s)
   real(r8) :: lwupocn(pcols)
   real(r8) :: rlat(pcols)          ! latitudes for chunk (radians)

   real(r8) :: Flwdabs,Fswabs,Fswabsv,Fswabsi
   real(r8) :: ft                   ! fraction reduction of positive qflux
   real(r8) :: thickfactor          ! thickness factor for adjusting qflux
   real(r8) :: ltheat(pcols)        ! latent heat (from flxoce) latvap
   real(r8) :: sstk(pcols)          ! ocn temperature (deg K.)
   real(r8) :: qflux_ft(pcols)      ! qflux after ft-factor adjustment (history diagnostic)
   real(r8) :: qflux_thick(pcols)   ! qflux after thickness-factor adjustment (history diagnostic)

   integer :: npts                  ! number of gridcells with ocn
   integer :: indx(pcols)
   integer :: i,ii,m                ! indices

   call t_startf ('mixed_layer1')
   call get_rlat_all_p (c, ncol, rlat)

   npts = 0
   do i=1,ncol
      if (landfrac(i) < 1.) then
         npts = npts + 1
         indx(npts) = i
!JR Array depth, mixed layer depth, is array mld at higher levels.
      else
         qflux(i) = 0.     ! set qflux to zero over land
!JR Put zero in sst over total land cover.  Plus I think Tbot is in deg. K and at this
!JR point we want deg C.
         sst(i) = 0.
!JR         sst(i) = Tbot(i)
!JR Changed Fhnet to Focn
!         Fhnet(i) = 0.
         Focn(i) = 0.
         frzmlt(i) = 0.
      end if
   end do

   lwupocn(:) = 0.
   shflxocn(:) = 0.
   lhflxocn(:) = 0.
   tauxocn(:) = 0.
   tauyocn(:) = 0.
   qflux_adj(:) = 0.
   qflux_ft(:) = 0.
   qflux_thick(:) = 0.
   onf(:) = 0.

   call outfld ('QFLUX   ', qflux, pcols, c)
   call outfld ('FOCN    ', Focn, pcols, c)
   call outfld ('MLDANN  ', depth, pcols, c)

   ltheat(:ncol) = latvap

   if (npts.gt.0) then
      sstk(:ncol) = sst(:ncol) + Tffresh
      call flxoce (indx, npts, pbot, ubot, vbot, &
                   tbot, qbot, thbot, zbot, sstk, &
                   ltheat, shflx, lhflx, taux, tauy, &
                   lwup, tref)
      do ii=1,npts
         i = indx(ii)
         shflxocn(i) = -shflx(i)
         lhflxocn(i) = -lhflx(i)
         tauxocn(i)  = -taux(i)
         tauyocn(i)  = -tauy(i)
         lwupocn(i)  = -lwup(i)

! Code copied from ice_atm_flux for shortwave absorption

         Fswabsv = swvdr(i)*(1.-alvdr(i)) + swvdf(i)*(1.-alvdf(i))
         Fswabsi = swidr(i)*(1.-alidr(i)) + swidf(i)*(1.-alidf(i))
         Fswabs  = fswabsv + fswabsi

! first, compute sst change due to exchange with atm above
! over open water:

         onf(i) = (Fswabs + flwds(i) + lwupocn(i) + lhflxocn(i) + shflxocn(i))*(1.-aice(i))

         sst(i) = sst(i) + onf(i)*dtime/(cpw*depth(i))

! adjust qflux if cooling of mixed layer will occur when sst < sst0;
! Reduce qflux > 0 down to zero

         if (sst(i) < 0. .and. qflux(i) > 0. ) then
            ft  = max((1. - sst(i)/Tfrez),0._r8)
         else
            ft  = 1.
         endif
#ifdef ZEROADJ
         ft = 1.
#endif
         qflux_adj(i) = qflux(i)*ft
!
! qflux_ft: history diagnostic indicating effect of only the ft factor
!
         qflux_ft(i) = qflux_adj(i) - qflux(i)

!BPB Extrae adjustment of qflux under ice; make
!BPB adjustment additive instead of supplantive as
!BPB previously; we will see if this might need
!BPB modification.....note that the aice weighting
!BPB ensures this is an under ice only addition:

!BPB: 9/13/02: Thickness dependence added for stability.
!JR: 9/13/02:  This could probably be done more elegantly, but as is it should work.

         if (rlat(i) > 0. ) then
            if (qnh > 0.) then
               thickfactor = 1./(1. + hi(i))
            else
               thickfactor = hi(i)/(1. + hi(i))
            end if
            qflux_adj(i) = qflux_adj(i) + qnh*aice(i)*thickfactor
!
! qflux_thick: history diagnostic indicating effect of only thickfactor
!
            qflux_thick(i) = qnh*aice(i)*(thickfactor - 1.)
         else
            if (qsh > 0.) then
               thickfactor = 1./(1. + hi(i))
            else
               thickfactor = hi(i)/(1. + hi(i))
            end if
            qflux_adj(i) = qflux_adj(i) + qsh*aice(i)*thickfactor
!
! qflux_thick: history diagnostic indicating effect of only thickfactor
!
            qflux_thick(i) = qsh*aice(i)*(thickfactor - 1.)
         endif

! compute T change due to exchange with deep layers:
         sst(i) = sst(i) - qflux_adj(i)*dtime/(cpw*depth(i))

! compute T change due to heat used by ice model:
         sst(i) = sst(i) + Focn(i)*dtime/(cpw*depth(i))

! compute potential to freeze or melt ice
         frzmlt(i) = (Tfrez - sst(i))*(cpw*depth(i))/dtime

! if sst is below freezing, reset sst to Tfrez
         sst(i) = max (sst(i), Tfrez)

         cflx(i,1) = -lhflxocn(i)/ltheat(i)
      end do

!JR Code barfs under lf95 strict error checking if the following loop isnt included
! Set non-water constituent fluxes to zero

      do m=2,pcnst+pnats
         do ii=1,npts
            i = indx(ii)
            cflx(i,m) = 0.
         end do
      end do
   end if

   call outfld ('QFLUX_FT', qflux_ft, pcols, c)
   call outfld ('QFLUX_TH', qflux_thick, pcols, c)
   call outfld ('ONF     ', onf, pcols, c)
   call t_stopf ('mixed_layer1')
   return
end subroutine mixed_layer1

subroutine mixed_layer2 (c, ncol, dtime, landfrac, depth, &
                         sst, ts, qdiff, qflux_adj)
!---!-------------------------------------------------------------------
!---! Calculate flux exchange over open ocean and update mixed layer sst
!---!
!---!-------------------------------------------------------------------

   use shr_kind_mod,  only: r8 => shr_kind_r8
   use ppgrid,        only: pcols
   use constituents,  only: pcnst, pnats
   use ice_constants, only : Tffresh, rhow, cp_ocn
   use history,       only: outfld

   implicit none

!------------------------------Arguments--------------------------------
   integer , intent(in) :: c               ! chunk index
   integer , intent(in) :: ncol            ! number of columns this chunk

   real(r8), intent(in) :: dtime           ! model timestep
   real(r8), intent(in) :: landfrac(pcols) ! land fraction
   real(r8), intent(in) :: depth(pcols)    ! mixed layer depth (m)
   real(r8), intent(in) :: qdiff           ! global mean diff qfluxgmin and qfluxgmot; also warm ocn area adjusted
   real(r8), intent(inout) :: qflux_adj(pcols)

   real(r8), intent(inout) :: sst(pcols)   ! ocean layer temp (C)

! fluxes/quantities summed over surface types

   real(r8), intent(out):: ts(pcols)       ! surface temperature (K)

!---------------------------Local variables-----------------------------
! fluxes/quantities over sea ice only
   real(r8), parameter :: cpw = rhow*cp_ocn  ! Specific heat of ocean water (J/m3/C)

   integer :: npts            ! number of gridcells with ocn
   integer :: indx(pcols)
   integer :: i,ii            ! indices

   call t_startf ('mixed_layer2')
   npts = 0
   do i=1,ncol
      Ts(i) = 0.0
      if (landfrac(i) < 1.) then
         npts = npts + 1
         indx(npts) = i
      end if
   end do

   do ii=1,npts
      i = indx(ii)
      if (sst(i) > 0.) then
         sst(i) = sst(i) - qdiff*dtime/(cpw*depth(i))

!BPB Check conservation: add qdiff to qflux_adj over warm ocean,
!BPB Global mean will then be computed, and compared with the initial value

         qflux_adj(i) = qflux_adj(i) + qdiff
      end if
      Ts(i) = sst(i) + Tffresh
   end do
   call outfld ('QFLUX_A2', qflux_adj, pcols, c)
   call t_stopf ('mixed_layer2')

   return
end subroutine mixed_layer2
