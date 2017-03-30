#include <misc.h>
#include <preproc.h>

module DustMod

#if (defined BGC)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: DustMod
!
! !DESCRIPTION: 
! Routines in this module calculate Dust mobilization and dry deposition for dust.
! Simulates dust mobilization due to wind from the surface into the 
! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface dust 
! emission (kg/m**2/s) [ + = to atm].
! Calculates the turbulent component of dust dry deposition, (the turbulent deposition 
! velocity through the lowest atmospheric layer). CAM will calculate the settling 
! velocity through the whole atmospheric column. The two calculations will determine 
! the dust dry deposition flux to the surface.
!                              
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8 
  use clmtype
  use spmdMod   , only : masterproc
  use clm_varpar, only : dst_src_nbr, ndst, sz_nbr
  use clm_varcon, only : grav, istsoil
!  
! !PUBLIC TYPES
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public Dustini        ! Initialize variables used in subroutine Dust
  public DustEmission   ! Dust mobilization 
  public DustDryDep     ! Turbulent dry deposition for dust
!
! !REVISION HISTORY
! Created by Sam Levis, updated to clm2.1 by Mariana Vertenstein
! Source: C. Zender's dust model
!
!EOP
!
! Data private to this module
!
  private
  real(r8) tmp1          !Factor in saltation computation (named as in Charlie's code)
  real(r8) ovr_src_snk_mss(dst_src_nbr,ndst)  
  real(r8) dmt_vwr(ndst) ![m] Mass-weighted mean diameter resolved
  real(r8) stk_crc(ndst) ![frc] Correction to Stokes settling velocity
  real(r8) dns_aer       ![kg m-3] Aerosol density
!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DustEmission
!
! !INTERFACE:
  subroutine DustEmission (c)
!
! !DESCRIPTION: 
! Dust mobilization. This code simulates dust mobilization due to wind
! from the surface into the lowest atmospheric layer
! On output flx_mss_vrt_dst(ndst) is the surface dust emission 
! (kg/m**2/s) [ + = to atm]
! Source: C. Zender's dust model
!
! !USES
   use shr_const_mod, only : SHR_CONST_RHOFW
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout) :: c !column derived type
!
! !LOCAL VARIABLES
!
! local pointers to implicit in scalars
!
    integer , pointer :: ityplun         !landunit type
    real(r8), pointer :: tlai            !one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai            !one-sided stem area index, no burying by snow
    real(r8), pointer :: frac_sno        !fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: gwc_thr         !threshold gravimetric soil moisture based on clay content
    real(r8), pointer :: forc_rho        !density (kg/m**3)
    real(r8), pointer :: fv              !friction velocity (m/s) (for dust model)
    real(r8), pointer :: u10             !10-m wind (m/s) (created for dust model)
    real(r8), pointer :: mbl_bsn_fct     !basin factor
    real(r8), pointer :: mss_frc_cly_vld ![frc] Mass fraction clay limited to 0.20
!
! local pointers to implicit in arrays
!
    real(r8), pointer :: h2osoi_vol(:)   !volumetric soil water (0<=h2osoi_vol<=watsat)
    real(r8), pointer :: h2osoi_liq(:)   !liquid soil water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:)   !frozen soil water (kg/m2)
    real(r8), pointer :: watsat(:)       !saturated volumetric soil water

! local pointers to implicit out arrays
!
    real(r8), pointer :: flx_mss_vrt_dst(:)  !surface dust emission (kg/m**2/s) 
    real(r8), pointer :: flx_mss_vrt_dst_tot !total dust flux into atmosphere 

! !REVISION HISTORY
! Created by Sam Levis
! Migrated to new data structures by Peter Thornton and Mariana Vertenstein
! !Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: pi              !pft index
    integer  :: m,n             !indices
    real(r8) :: liqfrac         !fraction of total water that is liquid
    real(r8) :: wnd_frc_rat     ![frc] Wind friction threshold over wind friction
    real(r8) :: wnd_frc_slt_dlt ![m s-1] Friction velocity increase from saltatn
    real(r8) :: wnd_rfr_dlt     ![m s-1] Reference windspeed excess over threshld
    real(r8) :: dst_slt_flx_rat_ttl
    real(r8) :: flx_mss_hrz_slt_ttl
    real(r8) :: flx_mss_vrt_dst_ttl
    real(r8) :: frc_thr_wet_fct
    real(r8) :: frc_thr_rgh_fct
    real(r8) :: wnd_frc_thr_slt
    real(r8) :: wnd_rfr_thr_slt
    real(r8) :: wnd_frc_slt
    real(r8) :: lnd_frc_mbl
    real(r8) :: bd
    real(r8) :: gwc_sfc
!    
! constants
!
    real(r8), parameter :: cst_slt = 2.61           ![frc] Saltation constant
    real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4 ![frc] Empir. mass flx tuning eflx_lh_vegt
    real(r8), parameter :: vai_mbl_thr = 0.1        ![m2 m-2] VAI threshold quenching dust mobilization
!
! local pointers to derived subtypes
!
    type(atm2lnd_state_type), pointer :: a2ls
    type(column_pstate_type), pointer :: cps
    type(column_wstate_type), pointer :: cws
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_dflux_type)    , pointer :: pdf
!------------------------------------------------------------------------

    ! Assign local pointers to column-level derived subtypes
    a2ls => c%a2ls
    cps => c%cps
    cws => c%cws

    ! Assign local pointers to column-level implicit input arrays
    h2osoi_vol => cws%h2osoi_vol
    h2osoi_liq => cws%h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice
    watsat => cps%watsat

    ! Assign local pointers to derived type scalar members (column-level)
    ityplun => cps%lps%itype
    frac_sno => cps%frac_sno
    gwc_thr => cps%gwc_thr
    forc_rho => a2ls%forc_rho
    mbl_bsn_fct => cps%mbl_bsn_fct
    mss_frc_cly_vld => cps%mss_frc_cly_vld

    ! Begin loop through pfts
    do pi=1, cps%npfts

       ! Assign local pointers to pft-level derived subtypes
       p => c%p(pi)
       pps => p%pps
       pdf => p%pdf

       ! Assign local pointers to pft-level implicit output arrays
       flx_mss_vrt_dst => pdf%flx_mss_vrt_dst
       flx_mss_vrt_dst_tot => pdf%flx_mss_vrt_dst_tot

       ! Assign local pointers to derived type scalar members (pft-level)
       tlai => pps%tlai
       tsai => pps%tsai
       fv => pps%fv
       u10 => pps%u10

       ! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
       ! purpose: return fraction of each gridcell suitable for dust mobilization

       ! the "bare ground" fraction of the current sub-gridscale cell decreases
       ! linearly from 1 to 0 as VAI(=tlai+tsai) increases from 0 to vai_mbl_thr

       if (ityplun == istsoil) then
          if (tlai+tsai < vai_mbl_thr) then
             lnd_frc_mbl = 1.0 - (tlai+tsai)/vai_mbl_thr
          else
             lnd_frc_mbl = 0.0
          endif
          lnd_frc_mbl = lnd_frc_mbl * (1.0 - frac_sno)
       else          !if ice sheet, wetland, or lake, no dust allowed
          lnd_frc_mbl = 0.0
       end if

       if (lnd_frc_mbl>1.0 .or. lnd_frc_mbl<0.0) then
          write (6,*)'Error dstmbl: pft index1d= ',p%pps%index1d, &
               ' lnd_frc_mbl= ',lnd_frc_mbl
          call endrun 
       end if

       ! reset history output variables before next if-statement to avoid output = inf

       flx_mss_vrt_dst_tot = 0.0
       do n = 1, ndst
          flx_mss_vrt_dst(n) = 0.0
       end do

       ! only perform the following calculations if lnd_frc_mbl is non-zero 

       if (lnd_frc_mbl /= 0.0) then

          ! the following comes from subr. frc_thr_rgh_fct_get
          ! purpose: compute factor by which surface roughness increases threshold
          !          friction velocity (currently a constant)
          
          frc_thr_rgh_fct = 1.0
          
          ! the following comes from subr. frc_thr_wet_fct_get
          ! purpose: compute factor by which soil moisture increases threshold friction
          !          velocity
          ! adjust threshold velocity for inhibition by moisture
          ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
          ! water content
          
          if (lnd_frc_mbl > 0.0) then
             bd = (1.-watsat(1))*2.7e3 ![kg m-3] Bulk density of dry surface soil
             gwc_sfc=h2osoi_vol(1)*SHR_CONST_RHOFW/bd ![kg kg-1] Gravimetric H2O cont
             if (gwc_sfc > gwc_thr) then
                frc_thr_wet_fct = sqrt(1.0 + 1.21 * (100.0*(gwc_sfc-gwc_thr))**0.68)
             else
                frc_thr_wet_fct = 1.0
             end if
             ! slevis: adding liqfrac here, because related to effects from soil water
             liqfrac = max( 0.0, min( 1.0, h2osoi_liq(1) / &
                  (h2osoi_ice(1)+h2osoi_liq(1)+1.0e-6) ) )
          end if
          
          ! the following lines come from subr. dst_mbl
          ! purpose: adjust threshold friction velocity to acct for moisture and
          !          roughness. The ratio tmp1 / sqrt(forc_rho) comes from
          !          subr. wnd_frc_thr_slt_get which computes dry threshold
          !          friction velocity for saltation
          
          wnd_frc_thr_slt = tmp1 / sqrt(forc_rho) * frc_thr_wet_fct * frc_thr_rgh_fct
          
          ! reset these variables which will be updated in the following if-block
          
          wnd_frc_slt = fv
          flx_mss_hrz_slt_ttl = 0.0
          flx_mss_vrt_dst_ttl = 0.0
          
          if (lnd_frc_mbl > 0.0) then
             
             ! the following line comes from subr. dst_mbl
             ! purpose: threshold saltation wind speed
             
             wnd_rfr_thr_slt = u10 * wnd_frc_thr_slt / fv
             
             ! the following if-block comes from subr. wnd_frc_slt_get 
             ! purpose: compute the saltating friction velocity
             ! theory: saltation roughens the boundary layer, AKA "Owen's effect"
             
             if (u10 >= wnd_rfr_thr_slt) then
                wnd_rfr_dlt = u10 - wnd_rfr_thr_slt
                wnd_frc_slt_dlt = 0.003 * wnd_rfr_dlt * wnd_rfr_dlt
                wnd_frc_slt = fv + wnd_frc_slt_dlt
             endif
             
             ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
             ! purpose: compute vertically integrated streamwise mass flux of particles
             
             if (wnd_frc_slt > wnd_frc_thr_slt) then
                wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
                flx_mss_hrz_slt_ttl = cst_slt * forc_rho * (wnd_frc_slt**3.0) * &
                     (1.0 - wnd_frc_rat) * (1.0 + wnd_frc_rat) * (1.0 + wnd_frc_rat) / grav
                
                ! the following loop originates from subr. dst_mbl
                ! purpose: apply land sfc and veg limitations and global tuning factor
                ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect 
                ! of frozen soil
                
                flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl * mbl_bsn_fct * &
                     flx_mss_fdg_fct * liqfrac
             endif

             ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
             ! purpose: diagnose total vertical mass flux of dust from vertically
             !          integrated streamwise mass flux
             
             dst_slt_flx_rat_ttl = 100.0 * exp( log(10.0) * (13.4*mss_frc_cly_vld - 6.0) )
             flx_mss_vrt_dst_ttl = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl
             
          endif !lnd_frc_mbl > 0.0
          
          ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
          ! purpose: partition total vertical mass flux of dust into transport bins
          
          do n = 1, ndst
             if (lnd_frc_mbl > 0.0) then
                do m = 1, dst_src_nbr
                   flx_mss_vrt_dst(n) = flx_mss_vrt_dst(n) + &
                        ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl
                end do
                flx_mss_vrt_dst_tot = flx_mss_vrt_dst_tot + flx_mss_vrt_dst(n)
             end if
          end do

       endif   ! end of lnd_frc_mbl not zero if-block  

    end do ! end of pft loop   

    return
  end subroutine DustEmission

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine DustDryDep(c)
!
! !INTERFACE:
!
  subroutine DustDryDep(c)
!
! !DESCRIPTION: 
!
    ! Determine Turbulent dry deposition for dust. Calculate the turbulent 
    ! component of dust dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the dust dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, dustini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
!
! !USES
!
    use shr_const_mod, ONLY: SHR_CONST_PI, SHR_CONST_RDAIR, SHR_CONST_BOLTZ
!
! !ARGUMENTS:
!
    implicit none
    type (column_type),target,intent(inout) :: c !column derived type
!
! !LOCAL VARIABLES
! The following only list local variables that have implicit input/output effects
!
! local pointers to implicit in scalars
    real(r8), pointer :: forc_t      !atm temperature (K)
    real(r8), pointer :: forc_pbot   !atm pressure (Pa)
    real(r8), pointer :: forc_rho    !atm density (kg/m**3)
    real(r8), pointer :: fv          !friction velocity (m/s)
    real(r8), pointer :: ram1        !aerodynamical resistance (s/m)
    real(r8), pointer :: vlc_trb(:)  !Turbulent deposn velocity (m/s)
    real(r8), pointer :: vlc_trb_1   !Turbulent deposition velocity 1
    real(r8), pointer :: vlc_trb_2   !Turbulent deposition velocity 2
    real(r8), pointer :: vlc_trb_3   !Turbulent deposition velocity 3
    real(r8), pointer :: vlc_trb_4   !Turbulent deposition velocity 4
!
! !REVISION HISTORY
! Created by Sam Levis
!
!EOP
!------------------------------------------------------------------------

!------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,pi          !indices
    real(r8) :: vsc_dyn_atm   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(ndst) ![frc] Slip correction factor
    real(r8) :: vlc_grv(ndst) ![m s-1] Settling velocity
    real(r8) :: rss_lmn(ndst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp           !temporary 

    ! constants
    real(r8),parameter::shm_nbr_xpn_lnd=-2./3. ![frc] shm_nbr_xpn over land

    ! local pointers to derived subtypes
    type(atm2lnd_state_type), pointer :: a2ls
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(column_pstate_type), pointer :: cps
    type(pft_dflux_type)    , pointer :: pdf
!------------------------------------------------------------------------

    ! assign local pointers to derived subtypes
    a2ls => c%a2ls
    cps => c%cps

    ! assign local pointers to derived type scalar members
    forc_pbot => a2ls%forc_pbot
    forc_rho => a2ls%forc_rho
    forc_t => a2ls%forc_t

    do pi = 1,cps%npfts !begin pft loop

       ! Assign local pointers to pft-level derived subtypes
       p => c%p(pi)
       pps => p%pps
       pdf => p%pdf

       ! Assign local pointers to derived type scalar members (pft-level)
       fv => pps%fv
       ram1 => pps%ram1
       vlc_trb => pdf%vlc_trb
       vlc_trb_1 => pdf%vlc_trb_1
       vlc_trb_2 => pdf%vlc_trb_2
       vlc_trb_3 => pdf%vlc_trb_3
       vlc_trb_4 => pdf%vlc_trb_4

       ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
       ! when code asks to use midlayer density, pressure, temperature,
       ! I use the data coming in from the atmosphere, ie forc_t, forc_pbot, forc_rho

       ! Quasi-laminar layer resistance: call rss_lmn_get
       ! Size-independent thermokinetic properties
       vsc_dyn_atm = 1.72e-5_r8 * ((forc_t/273.0_r8)**1.5_r8) * 393.0_r8 / &
            (forc_t+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
       mfp_atm = 2.0_r8 * vsc_dyn_atm / &   ![m] SeP97 p. 455
            (forc_pbot*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*forc_t)))
       vsc_knm_atm = vsc_dyn_atm / forc_rho ![m2 s-1] Kinematic viscosity of air

       do m = 1, ndst
          slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm * &
               (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
               dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
               grav * slp_crc(m) / vsc_dyn_atm ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(m) = vlc_grv(m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
       end do

       do m = 1, ndst
          stk_nbr = vlc_grv(m) * fv * fv / (grav*vsc_knm_atm)    ![frc] SeP97 p.965
          dff_aer = SHR_CONST_BOLTZ * forc_t * slp_crc(m) / &    ![m2 s-1]
               (3.0_r8*SHR_CONST_PI*vsc_dyn_atm*dmt_vwr(m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm / dff_aer                        ![frc] SeP97 p.972
          shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv) ![s m-1] SeP97 p.972,965
       end do

       ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)
       do m = 1, ndst
          rss_trb = ram1 + rss_lmn(m) + ram1*rss_lmn(m)*vlc_grv(m) ![s m-1]
          vlc_trb(m) = 1.0_r8 / rss_trb                            ![m s-1]
       end do
       vlc_trb_1 = vlc_trb(1)
       vlc_trb_2 = vlc_trb(2)
       vlc_trb_3 = vlc_trb(3)
       vlc_trb_4 = vlc_trb(4)

    end do !end pft loop   

    return
  end subroutine DustDryDep

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Dustini()
!
! !INTERFACE:
!
  subroutine Dustini()
!
! !DESCRIPTION: 
!
! Compute source efficiency factor from topography
! Initialize other variables used in subroutine Dust:
! ovr_src_snk_mss(m,n) and tmp1.
! Define particle diameter and density needed by atm model
! as well as by dry dep model
! Source: Paul Ginoux (for source efficiency factor)
! Modifications by C. Zender and later by S. Levis
! Rest of subroutine from C. Zender's dust model
!
! !USES
!
    use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_RDAIR
    use clmpoint, only : cpoint
!
! !ARGUMENTS:
!
    implicit none
!
! !REVISION HISTORY
! Created by Samual Levis
!
!EOP
!------------------------------------------------------------------------

!------------------------------------------------------------------------
    !Local Variables
    integer  :: ci,m,n                  !indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ![frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ![frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ![frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 !Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 !Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ![m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ![m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ![m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ![m] Width of size bin
    real(r8) :: slp_crc(ndst)           ![frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ![m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ![m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ![m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ![frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ![frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     !temporary 
    real(r8) :: ln_gsd                  ![frc] ln(gsd)
    real(r8) :: gsd_anl                 ![frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ![m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ![m] Number median particle diameter
    real(r8) :: lgn_dst                 !Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ![frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ![frc] Current relative accuracy
    real(r8) :: itr_idx                 ![idx] Counting index
    real(r8) :: dns_mdp                 ![kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ![m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ![kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ![m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            !Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      !Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ![m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ![m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ![m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ![m] Size Bin widths
    
    ! constants
    real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ![m] Mass median diameter
         (/ 0.832e-6 , 4.82e-6 , 19.38e-6 /)        !BSM96 p. 73 Table 2
    real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ![frc] Geometric std deviation
         (/ 2.10     ,  1.90   , 1.60     /)        !BSM96 p. 73 Table 2
    real(r8) :: mss_frc_src(dst_src_nbr) =    &     ![frc] Mass fraction 
         (/ 0.036, 0.957, 0.007 /)                  !BSM96 p. 73 Table 2
    real(r8) :: dmt_grd(5) =                  &     ![m] Particle diameter grid
         (/ 0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6    ![m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0         ![kg m-3] Density of optimal saltation particles

    ! local pointers to derived subtypes
    type(column_type), pointer :: c !local pointer to derived subtype

    ! declare erf intrinsic function
    real(r8) :: dum     !dummy variable for erf test
#if (defined AIX) 
#define ERF erf
#else
#define ERF derf
    real(r8) derf
#endif
!------------------------------------------------------------------------

    ! Sanity check on erf: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

    dum = 1.0
    if (abs(0.8427-ERF(dum))/0.8427>0.001) then
       write (6,*) 'erf(1.0) = ',ERF(dum)
       write (6,*) 'Dustini: Error function error'
       call endrun
    endif
    dum = 0.0
    if (ERF(dum) /= 0.0) then
       write (6,*) 'erf(0.0) = ',ERF(dum)
       write (6,*) 'Dustini: Error function error'
       call endrun
    endif

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
       sqrt2lngsdi = sqrt(2.0) * log(gsd_anl_src(m))
       do n = 1, ndst
          lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
          lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
          ovr_src_snk_frc = 0.5 * (ERF(lndmaxjovrdmdni/sqrt2lngsdi) - &
                                   ERF(lndminjovrdmdni/sqrt2lngsdi))
          ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
       enddo
    enddo

    ! The following code from subroutine wnd_frc_thr_slt_get was placed 
    ! here because tmp1 needs to be defined just once

    ryn_nbr_frc_thr_prx_opt = 0.38 + 1331.0 * (100.0*dmt_slt_opt)**1.56

    if (ryn_nbr_frc_thr_prx_opt < 0.03) then
       write (6,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
       call endrun
    else if (ryn_nbr_frc_thr_prx_opt < 10.0) then
       ryn_nbr_frc_thr_opt_fnc = -1.0 + 1.928 * (ryn_nbr_frc_thr_prx_opt**0.0922)
       ryn_nbr_frc_thr_opt_fnc = 0.1291 * 0.1291 / ryn_nbr_frc_thr_opt_fnc
    else
       ryn_nbr_frc_thr_opt_fnc = 1.0 - 0.0858 * exp(-0.0617*(ryn_nbr_frc_thr_prx_opt-10.0))
       ryn_nbr_frc_thr_opt_fnc = 0.120 * 0.120 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
    endif

    icf_fct = 1.0 + 6.0e-07 / (dns_slt * grav * (dmt_slt_opt**2.5))
    dns_fct = dns_slt * grav * dmt_slt_opt
    tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

    ! Set basin factor to 1 for now

    do ci = cols1d%beg,cols1d%end
       c => cpoint(ci)%c
       c%cps%mbl_bsn_fct = 1.0
    end do

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)
    
    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here ndst = dst_nbr
    
    ! Override automatic grid with preset grid if available
    if (ndst == 4) then
       do n = 1, ndst
          dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
          dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
          dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
          dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
       end do
    else
       write (6,*) 'Dustini error: ndst must equal to 4 with current code'
       call endrun                                      !see more comments above
    endif                                               !endif ndst == 4

    ! Bin physical properties
    gsd_anl = 2.0      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e+3   ! [kg m-3] Aerosol density
    ! Set a fundamental statistic for each bin
    dmt_vma = 2.524e-6 ! [m] Mass median diameter analytic She84 p.75 Table1
    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)
    dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]

    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl
    do n = 1, ndst
       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do

       lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*SHR_CONST_PI))
       dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved
       do m = 1, sz_nbr
          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)
          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
       end do
       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do

    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)

    eps_max = 1.0e-4_r8
    dns_mdp = 100000. / (295.0*SHR_CONST_RDAIR) ![kg m-3] const prs_mdp & tpt_vrt
    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_r8 * ((295.0/273.0_r8)**1.5_r8) * 393.0_r8 / &
         (295.0+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000.*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*295.0)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, ndst
       slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
            (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
            grav * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, ndst

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]
       do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_r8) then
             cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                  (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                  9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                  (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_r8) then
             cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
          else
             write (6,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
             write(6,*)'Dustini error: Reynolds number too large in stk_crc_get()'
             call endrun 
          endif

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)
          vlc_grv(m) = sqrt(4.0_r8 * grav * dmt_vwr(m) * slp_crc(m) * dns_aer / &
               (3.0_r8*cff_drg_grv(m)*dns_mdp))   
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          endif
          if (itr_idx > 20) then
             write (6,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
             goto 100                                        !to next iteration
          endif
          itr_idx = itr_idx + 1

       end do                                                !end while
100    continue   !Label to jump to when iteration does not converge
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities
    do m = 1, ndst
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

    return
  end subroutine Dustini

#endif

end module DustMod






