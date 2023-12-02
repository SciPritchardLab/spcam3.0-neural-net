#include <misc.h>
#include <params.h>

subroutine tphysac (ztodt,   pblh,    qpert,   tpert,  shf,  &
                    taux,    tauy,    cflx,    sgh,    lhf,  &
                    landfrac,snowh,   tref,    precc,  precl,&
                    precsc, precsl,   state,   tend,    pbuf,&
                    ocnfrac, fsds, icefrac, fv, ram1 )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics after coupling to land, sea, and ice models.
! Computes the following:
!   o Radon surface flux and decay (optional)
!   o Vertical diffusion and planetary boundary layer
!   o Multiple gravity wave drag
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,             only: pcols, pver
   use chemistry,          only: trace_gas, chem_timestep_tend
   use gw_drag,            only: gw_intr
   use vertical_diffusion, only: vd_intr
   use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
                                 physics_ptend_init, physics_dme_adjust
   use phys_buffer,        only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents,       only: ppcnst, qmin
   use test_tracers,       only: trace_test1, trace_test2, trace_test3, test_tracers_timestep_tend
   use physconst,          only: zvir, gravit, rhoh2o, latvap,latice, cpair, rair
   use aerosol_intr,       only: aerosol_emis_intr, aerosol_drydep_intr, aerosol_srcsnk_intr
   use check_energy,    only: check_energy_chng
   use time_manager,    only: get_nstep

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
   real(r8), intent(in) :: landfrac(pcols)        ! Land fraction
   real(r8), intent(in) :: icefrac(pcols)         ! ice fraction
   real(r8), intent(in) :: ocnfrac(pcols)         ! Land fraction
   real(r8), intent(in) :: snowh(pcols)           ! snow depth (liquid water equivalent)
   real(r8), intent(in) :: fv(pcols)              ! for dry deposition velocities for dust from land model
   real(r8), intent(in) :: ram1(pcols)            ! for dry deposition velocities for dust from land model
   real(r8), intent(in) :: tref(pcols)            ! 2m air temperature
   real(r8), intent(in) :: precc(pcols)           ! convective precipitation
   real(r8), intent(in) :: precl(pcols)           ! large-scale precipitation
   real(r8), intent(in) :: fsds(pcols)            ! down solar flux
   real(r8), intent(in) :: precsc(pcols)           ! convective snow
   real(r8), intent(in) :: precsl(pcols)           ! large-scale snow
   real(r8), intent(out) :: pblh(pcols)           ! Planetary boundary layer height
   real(r8), intent(out) :: qpert(pcols,ppcnst)   ! Moisture/constit. perturbation (PBL)
   real(r8), intent(out) :: tpert(pcols)          ! Temperature perturbation (PBL)
   real(r8), intent(inout) :: shf(pcols)          ! Sensible heat flux (w/m^2)
   real(r8), intent(in) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(in) :: tauy(pcols)            ! Y surface stress (meridional)
   real(r8), intent(inout) :: cflx(pcols,ppcnst)  ! Surface constituent flux (kg/m^2/s)
   real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
   real(r8), intent(inout) :: lhf(pcols)          ! Latent heat flux (w/m^2)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

   integer  :: nstep                              ! current timestep number
   real(r8) :: zero(pcols)                        ! array of zeros

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns
   integer i,k,m                 ! Longitude, level indices
   integer :: yr, mon, day, tod       ! components of a date

   logical :: labort                            ! abort flag

   real(r8) tvm(pcols,pver)           ! virtual temperature
   real(r8) prect(pcols)              ! total precipitation
   real(r8) surfric(pcols)              ! surface friction velocity
   real(r8) obklen(pcols)             ! Obukhov length
   real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry

! physics buffer fields for total energy and mass adjustment
   integer itim, ifld
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
   real(r8), pointer, dimension(:,:) :: cld
!
!-----------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
!
! accumulate fluxes into net flux array
! jrm Include latent heat of fusion for snow
!
   do i=1,ncol
      tend%flx_net(i) = tend%flx_net(i) + shf(i) + (precc(i) + precl(i))*latvap*rhoh2o &
           + (precsc(i) + precsl(i))*latice*rhoh2o
   end do

! Initialize parameterization tendency structure

   call physics_ptend_init(ptend)

! emission of aerosols at surface
   call aerosol_emis_intr (state, ptend, cflx, ztodt)
   call physics_update (state, tend, ptend, ztodt)

! get nstep and zero array for energy checker
   zero = 0.
   nstep = get_nstep()

! Check if latent heat flux exceeds the total moisture content of the
! lowest model layer, thereby creating negative moisture.

   call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,          &
              state%q(1,pver,1),state%rpdel(1,pver) ,shf ,lhf ,cflx(1,1) )

!===================================================
! Source/sink terms for advected tracers.
!===================================================

   if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
      call test_tracers_timestep_tend(state, ptend, cflx, landfrac, ztodt)
      call physics_update (state, tend, ptend, ztodt)
   end if

! Advected greenhouse trace gases:

   if (trace_gas) then
      call chem_timestep_tend(state, ptend, cflx, ztodt, fh2o)
      call physics_update (state, tend, ptend, ztodt)
! Check energy integrals
      call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
   end if

!===================================================
! Vertical diffusion/pbl calculation
! Call vertical diffusion code (pbl, free atmosphere and molecular)
!===================================================

   call vd_intr (ztodt    ,state    ,taux     ,tauy     , shf    ,&
                 cflx     ,pblh     ,tpert    ,qpert    , surfric  ,&
                 obklen   ,ptend    ,cld      ,ocnfrac  , landfrac, sgh )

   call physics_update (state, tend, ptend, ztodt)
! Check energy integrals
   call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cflx(:,1), zero, zero, shf)

!  aerosol dry deposition processes
   call aerosol_drydep_intr (state, ptend, cflx(:,1), ztodt, &
       fsds, obklen, tref, surfric, prect, snowh, pblh, shf, landfrac, &
       icefrac, ocnfrac, fv, ram1)

   call physics_update (state, tend, ptend, ztodt)

!===================================================
! Gravity wave drag
!===================================================

   call gw_intr (state   ,sgh     ,pblh    ,ztodt   , ptend , landfrac)
   call physics_update (state, tend, ptend, ztodt)
! Check energy integrals
   call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)

! conversion of aerosols from one category to another by non-wet processes
   call aerosol_srcsnk_intr (state, ptend, ztodt, ocnfrac)
   call physics_update (state, tend, ptend, ztodt)

!-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   teout(:ncol) = state%te_cur(:ncol)

!*** BAB's FV heating kludge *** apply the heating as temperature tendency.
!*** BAB's FV heating kludge *** modify the temperature in the state structure
   state%t(:ncol,:pver) = tini(:ncol,:pver) + ztodt*tend%dtdt(:ncol,:pver)

! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
   call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)
!-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   if (aqua_planet) then
      labort = .false.
      do i=1,ncol
         if (ocnfrac(i) /= 1.) labort = .true.
      end do
      if (labort) then
         write(6,*) 'ERROR:  grid contains non-ocean point'
         call endrun ()
      endif
   endif

end subroutine tphysac
