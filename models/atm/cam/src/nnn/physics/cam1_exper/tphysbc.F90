#include <misc.h>
#include <params.h>
#define PCWDETRAIN
#define RADTIME 900.
#define SP_DIR_NS
#define MOMTRANS

subroutine tphysbc (ztodt,   pblht,   tpert,   ts,      sst,     &
                    qpert,   precl,   precc,   precsl,  precsc,  &
                    asdir,   asdif,   aldir,   aldif,   snowh,   &
                    qrs,     qrl,     flwds,   fsns,    fsnt,    &
                    flns,    flnt,    lwup,    srfrad,  sols,    &
                    soll,    solsd,   solld,   state,   tend,    &
                    pbuf,    prcsnw,  fsds ,   landm,   landfrac,&
		    ocnfrac, icefrac  &
#ifdef CRM
                   ,u_crm, v_crm, w_crm, t_crm, q_crm, qn_crm, qp_crm &
                   ,qrs_crm, qrl_crm, rad_buffer, qrs1, qrl1  &
		   ,fsds_crm,fsns_crm,fsntoa_crm,fsutoa_crm  &
		   ,flwds_crm,flns_crm,flut_crm   &
		   ,fsdsc_crm,fsntoac_crm,flnsc_crm, flutc_crm &
	           ,taux, tauy, shf, lhf &
#endif
                      )
!-----------------------------------------------------------------------
!
! Purpose:
! Tendency physics BEFORE coupling to land, sea, and ice models.
!
! Method:
! Call physics subroutines and compute the following:
!     o cloud calculations (cloud fraction, emissivity, etc.)
!     o radiation calculations
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use cldcond,         only: cldcond_tend, cldcond_zmconv_detrain, cldcond_sediment
   use param_cldoptics, only: param_cldoptics_calc
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use diagnostics,     only: diag_dynvar
   use history,         only: outfld, fillvalue
   use physconst,       only: gravit, latvap, cpair, tmelt, cappa, zvir, rair, rga
   use radheat,         only: radheat_net
   use constituents,    only: pcnst, pnats, ppcnst, qmin
   use constituents,    only: dcconnam, cnst_get_ind
   use zm_conv,         only: zm_conv_evap, zm_convr
   use time_manager,    only: is_first_step, get_nstep, get_curr_calday
   use moistconvection, only: cmfmca
   use check_energy,    only: check_energy_chng, check_energy_fix
   use dycore,          only: dycore_is
   use cloudsimulatorparms, only: doisccp
   use cloudsimulator, only: ccm_isccp
   use aerosol_intr, only: aerosol_wet_intr
#ifdef CRM
   use history,         only: outfldcol
   use crmdims,       only: crm_nx, crm_ny, crm_nz
   use buffer,        only: nrad_buffer
   use pkg_cldoptics, only: cldefr, cldems, cldovrlap
   use check_energy,  only: check_energy_timestep_init
   use cloudsimulatorparms
   use cloudsimulator, only: crm_isccp
#endif

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: ts(pcols)                      ! surface temperature
   real(r8), intent(in) :: sst(pcols)                     ! sea surface temperature
   real(r8), intent(inout) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(inout) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess
   real(r8), intent(in) :: asdir(pcols)                  ! Albedo: shortwave, direct
   real(r8), intent(in) :: asdif(pcols)                  ! Albedo: shortwave, diffuse
   real(r8), intent(in) :: aldir(pcols)                  ! Albedo: longwave, direct
   real(r8), intent(in) :: aldif(pcols)                  ! Albedo: longwave, diffuse
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: qrs(pcols,pver)            ! Shortwave heating rate
   real(r8), intent(inout) :: qrl(pcols,pver)            ! Longwave  heating rate
   real(r8), intent(inout) :: flwds(pcols)               ! Surface longwave down flux
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: lwup(pcols)                    ! Surface longwave up flux
   real(r8), intent(out) :: srfrad(pcols)                 ! Net surface radiative flux (watts/m**2)
   real(r8), intent(inout) :: sols(pcols)                   ! Direct beam solar rad. onto srf (sw)
   real(r8), intent(inout) :: soll(pcols)                   ! Direct beam solar rad. onto srf (lw)
   real(r8), intent(inout) :: solsd(pcols)                  ! Diffuse solar radiation onto srf (sw)
   real(r8), intent(inout) :: solld(pcols)                  ! Diffuse solar radiation onto srf (lw)
   real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
   real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
   real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
   real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
   real(r8), intent(out) :: prcsnw(pcols)                 ! snowfall rate (precsl + precsc)
   real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction
   real(r8), intent(in) :: ocnfrac(pcols)                ! land fraction
   real(r8), intent(in) :: icefrac(pcols)                ! land fraction
#ifdef CRM
   real(r8), intent(inout) :: u_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: v_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: w_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: t_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: q_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: qn_crm  (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: qp_crm (pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: qrs_crm(pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: qrl_crm(pcols, crm_nx, crm_ny, crm_nz)
   real(r8), intent(inout) :: rad_buffer(pcols,nrad_buffer)
   real(r8), intent(inout) :: qrs1(pcols,pver)
   real(r8), intent(inout) :: qrl1(pcols,pver)
   real(r8), intent(inout) :: fsds_crm(pcols,crm_nx,crm_ny)   ! Flux Shortwave Downwelling Surface
   real(r8), intent(inout) :: fsns_crm(pcols,crm_nx,crm_ny)   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsntoa_crm(pcols,crm_nx,crm_ny) ! Net column abs solar flux at model top
   real(r8), intent(inout) :: fsutoa_crm(pcols,crm_nx,crm_ny) ! Flux Shortwave Upwelling TOA
   real(r8), intent(inout) :: fsntoac_crm(pcols,crm_nx,crm_ny)! Clear sky total column abs solar flux
   real(r8), intent(inout) :: fsdsc_crm(pcols,crm_nx,crm_ny)  ! Clear sky downard solar flux surface 
   real(r8), intent(inout) :: flwds_crm(pcols,crm_nx,crm_ny)  ! Surface longwave down flux
   real(r8), intent(inout) :: flns_crm(pcols,crm_nx,crm_ny)   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flut_crm(pcols,crm_nx,crm_ny)   ! Outgoing lw flux at model top
   real(r8), intent(inout) :: flutc_crm(pcols,crm_nx,crm_ny)  ! Clear sky outgoing lw flux at model top
   real(r8), intent(inout) :: flnsc_crm(pcols,crm_nx,crm_ny)  ! Clear sky lw flux at srf (up-down)
   real(r8), intent(in)    :: taux(pcols)   ! surface stress (zonal) (N/m2)
   real(r8), intent(in)    :: tauy(pcols)   ! surface stress (zonal) (N/m2)
   real(r8), intent(in)    :: shf(pcols)   ! surface sensible heat flux (W/m2)
   real(r8), intent(in)    :: lhf(pcols)   ! surface latent heat flux (W/m2) 
#endif

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: rhdfda(pcols,pver)            ! dRh/dcloud, old 
   real(r8) :: rhu00 (pcols,pver)            ! Rh threshold for cloud, old

   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: zmrprd(pcols,pver)            ! rain production in ZM convection
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c
   real(r8) :: cmfsl(pcols,pver)             ! Moist convection lw stat energy flux
   real(r8) :: cmflq(pcols,pver)             ! Moist convection total water flux
   real(r8) :: dtcond(pcols,pver)            ! dT/dt due to moist processes
   real(r8) :: dqcond(pcols,pver,ppcnst)     ! dq/dt due to moist processes

   real(r8) cldst(pcols,pver)
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                      !       "     low  cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfdqr2(pcols,pver)               ! dq/dt due to moist convective rainout
   real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
   real(r8) cmfsl2(pcols,pver)                ! Moist convection lw stat energy flux
   real(r8) cmflq2(pcols,pver)                ! Moist convection total water flux
   real(r8) cnt(pcols)                        ! Top level of convective activity
   real(r8) cnb(pcols)                        ! Lowest level of convective activity
   real(r8) cnt2(pcols)                       ! Top level of convective activity
   real(r8) cnb2(pcols)                       ! Bottom level of convective activity
   real(r8) concld(pcols,pver)             
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) prect(pcols)                      ! total (conv+large scale) precip rate
   real(r8) dlf2(pcols,pver)                   ! dq/dt due to rainout terms
   real(r8) qpert2(pcols,ppcnst)              ! Perturbation q
   real(r8) rtdt                              ! 1./ztodt
   real(r8) tpert2(pcols)                     ! Perturbation T
   real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
!                                             !    maximally overlapped region.
!                                             !    0->pmxrgn(i,1) is range of pressure for
!                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                             !    2nd region, etc
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
   integer  i,k,m                             ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
                                           
!  real(r8) engt                              ! Thermal   energy integral
!  real(r8) engk                              ! Kinetic   energy integral
!  real(r8) engp                              ! Potential energy integral
   real(r8) rel(pcols,pver)                   ! Liquid cloud particle effective radius
   real(r8) rei(pcols,pver)                   ! Ice effective drop size (microns)
   real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
   real(r8) clc(pcols)                        ! Total convective cloud (cloud scheme)
   real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
!
   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for 

! physics buffer fields to compute tendencies for cloud condensation package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: qcwat, tcwat, lcwat, cld

! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
!                                          
! Used for OUTFLD only                     
!                                          
   real(r8) icwmr1(pcols,pver)                ! in cloud water mixing ration for zhang scheme
   real(r8) icwmr2(pcols,pver)                ! in cloud water mixing ration for hack scheme
   real(r8) fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble
   real(r8) timestep(pcols)
!
!     Variables for doing deep convective transport outside of zm_convr
!
   real(r8) mu2(pcols,pver)
   real(r8) eu2(pcols,pver)
   real(r8) du2(pcols,pver)
   real(r8) md2(pcols,pver)
   real(r8) ed2(pcols,pver)
   real(r8) dp(pcols,pver)
   real(r8) dsubcld(pcols)
   real(r8) conicw(pcols,pver)
   real(r8) cmfdqrt(pcols,pver)               ! dq/dt due to moist convective rainout

! stratiform precipitation variables
   real(r8) :: prec_pcw(pcols)                ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)                ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)                ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)                ! snow from cloud ice sedimentation

! convective precipitation variables
   real(r8) :: prec_zmc(pcols)                ! total precipitation from ZM convection
   real(r8) :: snow_zmc(pcols)                ! snow from ZM convection
   real(r8) :: prec_cmf(pcols)                ! total precipitation from Hack convection
   real(r8) :: snow_cmf(pcols)                ! snow from Hack convection

! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: cape(pcols)                   ! cape
   real(r8) :: flx_cnd(pcols)
   real(r8) :: flx_heat(pcols)
   logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps

   integer jt(pcols)
   integer maxg(pcols)
   integer ideep(pcols)
   integer lengath
   real(r8) cldc(pcols,pver)
   real(r8) nevapr(pcols,pver)
   real(r8) qme(pcols,pver)
   real(r8) prain(pcols,pver)
   real(r8) cflx(pcols,ppcnst)

#ifdef FILLDEBUG
   real(r8) fakefld(pcols,pver)
   integer lons(pcols)
#endif
#ifdef CRM

   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsdsc(pcols)         ! Clear sky flux Shortwave Downwelling Surface
   real(r8) :: fsutoa(pcols)        ! Flux Shortwave Upwelling TOA
   real(r8) :: fsutoac(pcols)       ! Clear sky Flux Shortwave Upwelling TOA
   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: flut(pcols)          ! Upward flux at top of model
   real(r8) :: flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) :: flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) :: lwcf(pcols)          ! longwave cloud forcing
   real(r8) :: swcf(pcols)          ! shortwave cloud forcing
   real(r8) :: solin(pcols)         ! Solar incident flux
   real(r8) :: flwdsc(pcols)        ! Clear-sky Surface longwave down flux

   type(physics_state) :: state_save
   type(physics_tend ) :: tend_save
   real(r8) tmp_crm(pcols, crm_nx, crm_ny, crm_nz)
   real(r8) qc_crm (pcols, crm_nx, crm_ny, crm_nz)
   real(r8) qi_crm (pcols, crm_nx, crm_ny, crm_nz)
   real(r8) qpc_crm(pcols, crm_nx, crm_ny, crm_nz)
   real(r8) qpi_crm(pcols, crm_nx, crm_ny, crm_nz)
   real(r8) prec_crm(pcols, crm_nx, crm_ny)
   real(r8) mctot(pcols,pver)          ! total cloud mass flux
   real(r8) mcup(pcols,pver)           ! cloud updraft mass flux
   real(r8) mcdn(pcols,pver)           ! cloud downdraft mass flux
   real(r8) mcuup(pcols,pver)          ! unsaturated updraft mass flux
   real(r8) mcudn(pcols,pver)          ! unsaturated downdraft mass flux
   real(r8) crm_qc(pcols,pver)         ! cloud water
   real(r8) crm_qi(pcols,pver)         ! cloud ice
   real(r8) crm_qs(pcols,pver)         ! snow
   real(r8) crm_qg(pcols,pver)         ! graupel
   real(r8) crm_qr(pcols,pver)         ! rain
   real(r8) flux_qt(pcols,pver)        ! nonprecipitating water flux
   real(r8) flux_qp(pcols,pver)        ! precipitating water flux
   real(r8) precflux(pcols,pver)       ! precipitation flux
   real(r8) qt_ls(pcols,pver)        ! water tendency due to large-scale
   real(r8) qt_trans(pcols,pver)     ! nonprecip water tendency due to transport
   real(r8) qp_trans(pcols,pver)     ! precip water tendency due to transport
   real(r8) qp_fall(pcols,pver)      ! precip water tendency due to fall-out
   real(r8) qp_evp(pcols,pver)       ! precip water tendency due to evap
   real(r8) qp_src(pcols,pver)       ! precip water tendency due to conversion
   real(r8) t_ls(pcols,pver)        ! tendency of crm's liwse due to large-scale
   real(r8) t_rad (crm_nx, crm_ny, crm_nz) ! rad temperuture
   real(r8) qv_rad(crm_nx, crm_ny, crm_nz) ! rad vapor
   real(r8) qc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud water
   real(r8) qi_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice
   real(r8) trad(pcols,pver)
   real(r8) qvrad(pcols,pver,ppcnst)
   real(r8) fice(pcols,pver)                  ! Ice fraction from ice and liquid mixing ratios
   real(r8) cldn(pcols,pver) ! cloud top pdf
   real(r8) cldr(pcols,pver) ! cloud fraction based on -30dBZ radar reflectivity
   real(r8) cldtop(pcols,pver)
   real(r8) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
   real(r8) :: gicewp(pcols,pver)      ! grid-box cloud ice water path
   real(r8) :: gliqwp(pcols,pver)      ! grid-box cloud liquid water path
   real(r8) :: gwp   (pcols,pver)      ! grid-box cloud (total) water path
   real(r8) :: tgicewp(pcols)          ! Vertically integrated ice water path
   real(r8) :: tgliqwp(pcols)          ! Vertically integrated liquid water path
   real(r8) :: tgwp   (pcols)          ! Vertically integrated (total) cloud water path
   real(r8) stat_buffer(pcols,19) ! One-column CRM statistics for the ARM diagnostics
   real(r8) cld_crm(crm_nx, crm_ny, crm_nz)
   real(r8) cliqwp_crm(crm_nx, crm_ny, crm_nz)
   real(r8) cicewp_crm(crm_nx, crm_ny, crm_nz)
   real(r8) rel_crm(crm_nx, crm_ny, crm_nz)
   real(r8) rei_crm(crm_nx, crm_ny, crm_nz)
   real(r8) emis_crm(crm_nx, crm_ny, crm_nz)
   real(r8) fq_isccp_s1(pcols,ntau*npres) !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types
   real(r8) totalcldarea(pcols) !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  This should
					!  equal the sum over all entries of fq_isccp
   real(r8) lowcldarea(pcols), midcldarea(pcols), hghcldarea(pcols)
   real(r8) meantaucld(pcols) !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
   real(r8) meanttop(pcols) !  mean cloud top temp (k) - linear averaging
   real(r8) meanptop(pcols) !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
   real(r8) cloudy(pcols)
   integer  nbreak
   logical  doinitial
   real(r8) qtot
   real(r8) coef
   real(r8) kabs                   ! longwave absorption coeff (m**2/g)
! CRM column radiation stuff:
   real(r8) qrs_tmp(pcols,pver)            ! Shortwave heating rate
   real(r8) qrl_tmp(pcols,pver)            ! Longwave  heating rate
   real(r8) flwdsc_crm(pcols,crm_nx,crm_ny)  ! Clear-sky surface longwave down flux
   real(r8) sols_crm(pcols,crm_nx,crm_ny)    ! Direct beam solar rad. onto srf (s w)
   real(r8) soll_crm(pcols,crm_nx,crm_ny)    ! Direct beam solar rad. onto srf (l w)
   real(r8) solsd_crm(pcols,crm_nx,crm_ny)   ! Diffuse solar radiation onto srf ( sw)
   real(r8) solld_crm(pcols,crm_nx,crm_ny)   ! Diffuse solar radiation onto srf (
   real(r8) fsnt_crm(pcols,crm_nx,crm_ny)    ! Net downward solar flux at model top
   real(r8) fsntc_crm(pcols,crm_nx,crm_ny)   ! Clear sky net downward solar flux at model top
   real(r8) fsnsc_crm(pcols,crm_nx,crm_ny)   ! Clear sky net downward Shortwave flux Surface
   real(r8) fsutoac_crm(pcols,crm_nx,crm_ny) ! Clear sky Shortwave Flux Upwelling TOA
   real(r8) flnt_crm(pcols,crm_nx,crm_ny)    ! Net Upward lw flux at top of model
   real(r8) flntc_crm(pcols,crm_nx,crm_ny)   ! Net Upward Clear Sky flux at top of model
   real(r8) solin_crm(pcols,crm_nx,crm_ny)   ! Solar incident flux
   real(r8) prectend(pcols) ! tendency in precipitating water and ice
   real(r8) precstend(pcols) ! tendency in precipitating ice
   real(r8) wtricesink(pcols) ! sink of water vapor + cloud water + cloud ice
   real(r8) icesink(pcols) ! sink of
   real(r8) tau00  ! surface stress
   real(r8) wnd  ! surface wnd
   real(r8) bflx   ! surface buoyancy flux (Km/s)
   real(r8) taux_crm(pcols)  ! zonal CRM surface stress perturbation
   real(r8) tauy_crm(pcols)  ! merid CRM surface stress perturbation
   real(r8) z0m(pcols)  ! surface momentum roughness length

   real timing_factor(pcols) ! factor for crm cpu-usage: 1 means no subcycling

   integer ii, jj, mm
   integer iii,lll
#endif
!
!-----------------------------------------------------------------------
   zero = 0.
!
   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1./ztodt

   nstep = get_nstep()
   calday = get_curr_calday()
!
! Output NSTEP for debugging
!
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)

#ifdef FILLDEBUG
   call get_lon_all_p(lchnk, ncol, lons)
   do k=1,pver
      do i=1,ncol
         if (nstep > lons(i)) then
            fakefld(i,k) = 1.d0/nstep
         else
            fakefld(i,k) = fillvalue
         end if
      end do
   end do
   call outfld('FAKEFLD ',fakefld, pcols,lchnk)
#endif

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QCWAT')
   qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('TCWAT')
   tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('LCWAT')
   lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.

   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure
!
! Make sure that input tracers are all positive (probably unnecessary)
!
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              ppcnst,qmin  ,state%q )
!
! Setup q and t accumulation fields
!
   dqcond(:ncol,:,:) = state%q(:ncol,:,:)
   dtcond(:ncol,:)   = state%s(:ncol,:)

   fracis (:ncol,:,1:ppcnst) = 1.

!===================================================
! Global mean total energy fixer
!===================================================
   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
   qini(:ncol,:pver) = state%q(:ncol,:pver,1)

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )

!
!===================================================
! Dry adjustment
!===================================================

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

   call t_startf ('dadadj')
   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call t_stopf ('dadadj')
   call physics_update (state, tend, ptend, ztodt)

#ifdef CRM
! Save the state and tend variables to overwrite conventional physics effects
! leter before calling the superparameterization. Conventional moist
! physics is allowed to compute tendencies due to conventional
! moist physics for diagnostics purposes. -Marat

    state_save = state
    tend_save = tend

#endif
!
!===================================================
! Moist convection
!===================================================
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   qpert(:ncol,2:ppcnst) = 0.0
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')
   call zm_convr( lchnk,    ncol, &
                  state%t,   state%q,    prec_zmc,   cnt,     cnb,      &
                  pblht,   state%zm, state%phis,    state%zi,   ptend%q(:,:,1),     &
                  ptend%s, state%pmid,   state%pint,  state%pdel,       &
                  .5*ztodt,cmfmc,    cmfcme,             &
                  tpert,   dlf,      pflx,    zdu,     zmrprd,   &
                  mu2,      md2,     du2,     eu2,     ed2,      &
                  dp,       dsubcld, jt,      maxg,    ideep,    &
                  lengath, icwmr1,   rliq, cape    )

!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) * 100./gravit

!  Momentum transport by deep convection, Marat Khairoutdinov, May 2005
#ifdef MOMTRANS
   ptend%u(:,:) = 0.
   ptend%v(:,:) = 0.
   do i=1,ncol 
     if(prec_zmc(i)*86400.*1000..gt.20.) then   
!       print*, 'prec_zmc=',prec_zmc(i)*86400.*1000., int(cnt(i)), state%zm(i,int(cnt(i)))
       call convmom(state%zm(i,:),state%zi(i,:),state%pdel(i,:),state%u(i,:),cape, cmfmc(i,:),&
                    state%zm(i,int(cnt(i))), pver, ptend%u(i,:))  
       call convmom(state%zm(i,:),state%zi(i,:),state%pdel(i,:),state%v(i,:),cape, cmfmc(i,:),&
                    state%zm(i,int(cnt(i))), pver, ptend%v(i,:))  
     end if
   end do
   ptend%lu    = .TRUE.
   ptend%lv    = .TRUE.
   call outfld('UCONVMOM',ptend%u,pcols   ,lchnk   )
   call outfld('VCONVMOM',ptend%v,pcols   ,lchnk   )
#endif


!   if(lchnk.eq.478.and.nstep.le.96) then
!        i = 10
!        write(22,*) &
!            state%t(i,:),   state%q(i,:,:ppcnst), pblht(i),   &
!            state%zm(i,:), state%phis(i),    state%zi(i,:), tpert(i), qpert(i,:), &
!            state%pmid(i,:),   state%pint(i,:),  state%pdel(i,:), ztodt, &
!            state%omega(i,:),tcwat(i,:), qcwat(i,:), lcwat(i,:), cld(i,:),&
!            landfrac(i),snowh(i),icefrac(i),landm(i),ts(i), sst(i)
!
!   end if

   ptend%name  = 'zm_convr'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   cmfsl (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.
   cmflq (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.

   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf('zm_convr')

   call physics_update(state, tend, ptend, ztodt)
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
   call zm_conv_evap(state, ptend, zmrprd, cld, ztodt, prec_zmc, snow_zmc, .false.)
   call physics_update(state, tend, ptend, ztodt)
! Check energy integrals, including "reserved liquid"
   flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
   call check_energy_chng(state, tend, "zm_evap", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)

! Transport cloud water and ice only
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ptend%name = 'convtran1'
   ptend%lq(ixcldice) = .true.
   ptend%lq(ixcldliq) = .true.
   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                  du2,     eu2,     ed2,     dp,      dsubcld,  &
                  jt,      maxg,    ideep,   1,       lengath,  &
                  nstep,   fracis,  ptend%q   )
   call t_stopf ('convtran1')
   call physics_update (state, tend, ptend, ztodt)

!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
!
   call t_startf('cmfmca')
   tpert2(:ncol  ) =0.
   qpert2(:ncol,:) = qpert(:ncol,:)  ! BAB Why is this not zero, if tpert2=0???
   call cmfmca (lchnk,   ncol, &
                nstep,   ztodt,   state%pmid,  state%pdel,   &
                state%rpdel,   state%zm,      tpert2,  qpert2,  state%phis,     &
                pblht,   state%t,   state%q,   ptend%s,   ptend%q,      &
                cmfmc2,  cmfdqr2, cmfsl2,  cmflq2,  prec_cmf,   &
                dlf2,     cnt2,    cnb2,    icwmr2   , rliq2)
   ptend%name  = 'cmfmca'
   ptend%ls    = .TRUE.
   ptend%lq(:) = .TRUE.


! Add shallow cloud water detrainment to cloud water detrained from ZM
   dlf(:ncol,:pver) = dlf(:ncol,:pver) + dlf2(:ncol,:pver)
   rliq(:ncol) = rliq(:ncol) + rliq2(:ncol)
   
   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
   call outfld('CMFDQ   ',ptend%q(1,1,1),pcols   ,lchnk   )
   call t_stopf('cmfmca')
   call physics_update (state, tend, ptend, ztodt)
!
! Determine the phase of the precipitation produced and add latent heat of fusion
   call zm_conv_evap(state, ptend, cmfdqr2, cld, ztodt, prec_cmf, snow_cmf, .true.)
   call physics_update(state, tend, ptend, ztodt)
   flx_cnd(:ncol) = prec_cmf(:ncol) + rliq2(:ncol)
   call check_energy_chng(state, tend, "hk_evap", nstep, ztodt, zero, flx_cnd, snow_cmf, zero)
!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
   do i=1,ncol
      if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
   end do
!
   cmfmc (:ncol,:pver) = cmfmc (:ncol,:pver) + cmfmc2 (:ncol,:pver)
   cmfsl (:ncol,:pver) = cmfsl (:ncol,:pver) + cmfsl2 (:ncol,:pver)
   cmflq (:ncol,:pver) = cmflq (:ncol,:pver) + cmflq2 (:ncol,:pver)
   call outfld('CMFMC' , cmfmc  , pcols, lchnk)
!  output new partition of cloud condensate variables, as well as precipitation 
   call outfld('QC      ',dlf2           ,pcols   ,lchnk   )
   call outfld('PRECSH  ',prec_cmf      ,pcols   ,lchnk       )
   call outfld('CMFDQR', cmfdqr2, pcols, lchnk)
   call outfld('CMFSL' , cmfsl  , pcols, lchnk)
   call outfld('CMFLQ' , cmflq  , pcols, lchnk)
   call outfld('DQP'   , dlf2    , pcols, lchnk)

! Allow the cloud liquid drops and ice particles to sediment
! Occurs before adding convectively detrained cloud water, because the phase of the
! of the detrained water is unknown.
   call t_startf('cldwat_sediment')
   call cldcond_sediment(state, ptend, ztodt,cld, icefrac, landfrac, ocnfrac, prec_sed, &
                         snow_sed, landm, snowh)

   call physics_update(state, tend, ptend, ztodt)
   call t_stopf('cldwat_sediment')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_sediment", nstep, ztodt, zero, prec_sed, snow_sed, zero)

! Put the detraining cloud water from convection into the cloud and environment. 
   call t_startf('cldwat_detrain')
   call cldcond_zmconv_detrain(dlf, cld, state, ptend)
   call physics_update(state, tend, ptend, ztodt)
   call t_stopf('cldwat_detrain')

! check energy integrals, reserved liquid has now been used
   flx_cnd(:ncol) = -rliq(:ncol)
   call check_energy_chng(state, tend, "cldwat_detrain", nstep, ztodt, zero, flx_cnd, zero, zero)
!
! cloud fraction after transport and convection,
! derive the relationship between rh and cld from 
! the employed cloud scheme
!
   call t_startf('cldnrh')
!   if(lchnk.eq.478.and.nstep.lt.96) then
!        i = 10
!        write(22,*) &
!        state%t(i,:),   state%q(i,:,:), cnt(i),     cnb(i),cmfmc(i,:),   cmfmc2(i,:),zdu(i,:),dlf(i,:)
!   end if
   call cldnrh(lchnk,   ncol,                                &
               state%pmid,    state%t,   state%q(1,1,1),   state%omega, &
               cnt,     cnb,     cld,    clc,     state%pdel,   &
               cmfmc,   cmfmc2, landfrac,snowh,   concld,  cldst,    &
               ts,      sst, state%pint(1,pverp),       zdu,  ocnfrac, &
               rhdfda,   rhu00 , state%phis)
   call t_stopf('cldnrh')
#ifdef CRM
   call outfld('_CONCLD  ',concld, pcols,lchnk)
   call outfld('_CLDST   ',cldst,  pcols,lchnk)
   call outfld('_CNVCLD  ',clc,    pcols,lchnk)
#else
   call outfld('CONCLD  ',concld, pcols,lchnk)
   call outfld('CLDST   ',cldst,  pcols,lchnk)
   call outfld('CNVCLD  ',clc,    pcols,lchnk)
#endif
! cloud water and ice parameterizations
   call t_startf('cldwat_tend')
   call cldcond_tend(state, ptend, ztodt, &
       tcwat, qcwat, lcwat, prec_pcw, snow_pcw, icefrac, rhdfda, rhu00, cld, nevapr, prain, qme, snowh)
   call physics_update (state, tend, ptend, ztodt)
!   if(lchnk.eq.478.and.nstep.le.96) then
!        i = 10
!	write(22,*) state%t(i,:),   state%q(i,:,:), cld(i,:), prec_zmc(i), prec_cmf(i), prec_sed(i), &
!          prec_pcw(i)
!
!	if(nstep.eq.96) close(22)
!   end if
   call t_stopf('cldwat_tend')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_pcw, snow_pcw, zero)

! Save off q and t after cloud water
   do k=1,pver
      qcwat(:ncol,k) = state%q(:ncol,k,1)
      tcwat(:ncol,k) = state%t(:ncol,k)
      lcwat(:ncol,k) = state%q(:ncol,k,ixcldice) + state%q(:ncol,k,ixcldliq)
   end do
!
!  aerosol wet chemistry determines scavenging fractions, and transformations
   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
   call get_rlat_all_p(lchnk, ncol, clat)
   conicw(:ncol,:) = icwmr1(:ncol,:) + icwmr2(:ncol,:)
   cmfdqrt(:ncol,:) = zmrprd(:ncol,:) + cmfdqr2(:ncol,:)
   call aerosol_wet_intr (state, ptend, cflx, nstep, ztodt, lat, clat, lon,&
        qme, prain, &
       nevapr, cldc, cld, fracis, calday, cmfdqrt, conicw)
   call physics_update (state, tend, ptend, ztodt)

!
!     Convective transport of all trace species except water vapor and
!     cloud liquid and ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   ptend%name  = 'convtran2'
   ptend%lq(:) = .true.
   ptend%lq(ixcldice) = .false.
   ptend%lq(ixcldliq) = .false.
   call t_startf ('convtran2')
   call convtran (lchnk,                                           &
                  ptend%lq,state%q, ppcnst,     mu2,     md2,      &
                  du2,     eu2,     ed2,        dp,      dsubcld,  &
                  jt,      maxg,    ideep,      1,       lengath,  &
                  nstep,   fracis,  ptend%q)
   call t_stopf ('convtran2')

   call physics_update (state, tend, ptend, ztodt)
!
! Compute rates of temperature and constituent change due to moist processes
!
   dtcond(:ncol,:) = (state%s(:ncol,:) - dtcond(:ncol,:))*rtdt / cpair
   dqcond(:ncol,:,:) = (state%q(:ncol,:,:) - dqcond(:ncol,:,:))*rtdt
   call outfld('DTCOND  ',dtcond, pcols   ,lchnk   )
   do m=1,ppcnst
      call outfld(dcconnam(m),dqcond(1,1,m),pcols   ,lchnk )
      call outfld('DQCOND  ',dqcond(1,1,m),pcols   ,lchnk   )
   end do

! Compute total convective and stratiform precipitation and snow rates
   do i=1,ncol
      precc (i) = prec_zmc(i) + prec_cmf(i)
      precl (i) = prec_sed(i) + prec_pcw(i)
      precsc(i) = snow_zmc(i) + snow_cmf(i)
      precsl(i) = snow_sed(i) + snow_pcw(i)
! jrm These checks should not be necessary if they exist in the parameterizations
      if(precc(i).lt.0.) precc(i)=0.
      if(precl(i).lt.0.) precl(i)=0.
      if(precsc(i).lt.0.) precsc(i)=0.
      if(precsl(i).lt.0.) precsl(i)=0.
      if(precsc(i).gt.precc(i)) precsc(i)=precc(i)
      if(precsl(i).gt.precl(i)) precsl(i)=precl(i)
! end jrm
   end do
   prcsnw(:ncol) = precsc(:ncol) + precsl(:ncol)   ! total snowfall rate: needed by slab ocean model
!
!===================================================
! Moist physical parameteriztions complete: 
! send dynamical variables, and derived variables to history file
!===================================================
!
#ifndef CRM
   call diag_dynvar (lchnk, ncol, state)
#endif
!
!===================================================
! Radiation computations
!===================================================
!
! Cosine solar zenith angle for current time step
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call zenith (calday, clat, clon, coszrs, ncol)

   if (dosw .or. dolw) then

! Compute cloud water/ice paths and optical properties for input to radiation
      call t_startf('cldoptics')
      call param_cldoptics_calc(1, ncol, state, cld, landfrac, landm,icefrac, &
                                cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh)
      call t_stopf('cldoptics')
!
! Complete radiation calculations
!
      call t_startf ('radctl')
      call radctl (lchnk, 1, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cld, cicewp, cliqwp, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, state%zm, state, fsds &
#ifdef CRM
                  ,fsntoa ,fsntoac ,fsdsc   ,flwdsc  ,fsntc   ,fsnsc   , &
                   fsutoa ,fsutoac ,flut    ,flutc   ,flntc   ,flnsc   ,solin   , &
                   .true., doabsems, dosw, dolw &
#endif
             )
      call t_stopf ('radctl')



!
! Cloud cover diagnostics
! radctl can change pmxrgn and nmxrgn so cldsav needs to follow
! radctl.
!
      call cldsav (lchnk, 1, ncol, cld, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
!
! Dump cloud field information to history tape buffer (diagnostics)
!
#ifdef CRM
      call outfld('_CLOUD  ',cld,  pcols,lchnk)
      call outfld('_CLDTOT ',cltot  ,pcols,lchnk)
      call outfld('_CLDLOW ',cllow  ,pcols,lchnk)
      call outfld('_CLDMED ',clmed  ,pcols,lchnk)
      call outfld('_CLDHGH ',clhgh  ,pcols,lchnk)
#else
      call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
      call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
      call outfld('CLDMED  ',clmed  ,pcols,lchnk)
      call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
      call outfld('CLOUD   ',cld    ,pcols,lchnk)
#endif
      if (doisccp) then
         call ccm_isccp(lchnk, ncol, state%pmid, state%pint, state%q, &
                        state%t, ts, concld,   &
                        cld, cliqwp, cicewp, rel, rei, emis, cltot, coszrs)
      end if
   else

! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k =1 , pver
            do i = 1, ncol
               qrs(i,k) = qrs(i,k)/state%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state%pdel(i,k)
            end do
         end do
      end if

   end if


#ifdef CRM
!========================================================
!========================================================
!=   Warning! You're entering a no-spin zone ! ==========
!========================================================
!========================================================
!========================================================
!========================================================
!  CRM (Superparameterization).
! Author: Marat Khairoutdinov marat@atmos.colostate.edu
!========================================================

   call t_startf ('crm')

! Forget all changes to the state due to conventional physics above:

   state = state_save
   tend  = tend_save

! Initialize stuff:

   if(is_first_step()) then

      call check_energy_timestep_init(state, tend, pbuf)

      do k=1,crm_nz
        do i=1,ncol
             m = pver-k+1
#ifdef SP_DIR_NS
             u_crm  (i,:,:,k) = state%v(i,m)
             v_crm  (i,:,:,k) = state%u(i,m)
#else
             u_crm  (i,:,:,k) = state%u(i,m)  
             v_crm  (i,:,:,k) = state%v(i,m)  
#endif
             w_crm  (i,:,:,k) = 0.
             t_crm  (i,:,:,k) = state%t(i,m)
             q_crm  (i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
             qn_crm  (i,:,:,k) = state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
             qp_crm (i,:,:,k) = 0.
             qc_crm (i,:,:,k) = 0.
             qi_crm (i,:,:,k) = 0.
             qpc_crm (i,:,:,k) = 0.
             qpi_crm (i,:,:,k) = 0.
             qrs_crm(i,:,:,k) = (qrs(i,m)) /cpair
             qrl_crm(i,:,:,k) = (qrl(i,m)) /cpair
          end do
      end do
      rad_buffer(:,:) = 0.
! use radiation from grid-cell mean radctl on first time step
      coef = dble(crm_nx*crm_ny)
      do i=1,ncol
         prec_crm (i,:,:) = 0.
         qrs1(i,:pver) = qrs(i,:pver) * coef * state%pdel(i,:pver)
         qrl1(i,:pver) = qrl(i,:pver) * coef * state%pdel(i,:pver)
         rad_buffer(i,1) =  solin(i) * coef
         rad_buffer(i,2) =  fsds(i) * coef
         rad_buffer(i,3) =  fsdsc(i) * coef
         rad_buffer(i,4) =  fsnt(i) * coef
         rad_buffer(i,5) =  fsns(i) * coef
         rad_buffer(i,6) =  fsntc(i) * coef
         rad_buffer(i,7) =  fsnsc(i) * coef
         rad_buffer(i,8) =  fsntoa(i) * coef
         rad_buffer(i,9) =  fsntoac(i) * coef
         rad_buffer(i,10) = sols(i) * coef
         rad_buffer(i,11) = soll(i) * coef
         rad_buffer(i,12) = solsd(i) * coef
         rad_buffer(i,13) = solld(i) * coef
         rad_buffer(i,14) = flnt(i) * coef
         rad_buffer(i,15) = flut(i) * coef
         rad_buffer(i,16) = flntc(i) * coef
         rad_buffer(i,17) = flutc(i) * coef
         rad_buffer(i,18) = flns(i) * coef
         rad_buffer(i,19) = flnsc(i) * coef
         rad_buffer(i,20) = flwds(i) * coef
         rad_buffer(i,21) = fsutoa(i) * coef
         rad_buffer(i,22) = fsutoac(i) * coef
         rad_buffer(i,23) = flwdsc(i) * coef
         fsds_crm(i,:,:) = fsds(i)
         fsns_crm(i,:,:) = fsns(i)
         fsdsc_crm(i,:,:) = fsdsc(i)
         fsntoa_crm(i,:,:) = fsntoa(i)
         fsntoac_crm(i,:,:) = fsntoac(i)
         fsutoa_crm(i,:,:) = fsutoa(i)
         flwds_crm(i,:,:) = flwds(i)
         flns_crm(i,:,:) = flns(i)
         flnsc_crm(i,:,:) = flnsc(i)
         flnt_crm(i,:,:) = flnt(i)
         flntc_crm(i,:,:) = flntc(i)
         flut_crm(i,:,:) = flut(i)
         flutc_crm(i,:,:) = flutc(i)
         lwcf(i) = flutc(i) - flut(i)
         swcf(i) = fsntoa(i) - fsntoac(i)
      enddo
      lwcf(:) = flutc(:) - flut(:)
      swcf(:) = fsntoa(:) - fsntoac(:)
      ptend%q(:,:,1) = 0.
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.
      ptend%s(:,:) = 0.
      precc(:) = 0.
      precl(:) = 0.
      precsc(:) = 0.
      precsl(:) = 0.
      cltot(:) = 0.
      clhgh(:) = 0.
      clmed(:) = 0.
      cllow(:) = 0.
      cld(:,:) = 0.
      cldr(:,:) = 0.
      cldtop(:,:) = 0.
      gicewp(:,:)=0
      gliqwp(:,:)=0
      stat_buffer(:,:) = 0.
      mctot(:,:) = 0.
      mcup(:,:) = 0.
      mcdn(:,:) = 0.
      mcuup(:,:) = 0.
      mcudn(:,:) = 0.
      crm_qc(:,:) = 0.
      crm_qi(:,:) = 0.
      crm_qs(:,:) = 0.
      crm_qg(:,:) = 0.
      crm_qr(:,:) = 0.
      flux_qt(:,:) = 0.
      flux_qp(:,:) = 0.
      precflux(:,:) = 0.
      qt_ls(:,:) = 0.
      qt_trans(:,:) = 0.
      qp_trans(:,:) = 0.
      qp_fall(:,:) = 0.
      qp_evp(:,:) = 0.
      qp_src(:,:) = 0.
      z0m(:) = 0.
      taux_crm(:) = 0.
      tauy_crm(:) = 0.
      t_ls(:,:) = 0.


   else

      ptend%q(:,:,1) = 0.  ! necessary?
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.
      ptend%s(:,:) = 0. ! necessary?
     trad(:ncol,:)  = state%t(:ncol,:)
     qvrad(:ncol,:,1) = state%q(:ncol,:,1)
     cwp   = 0.
     fice  = 0.
     cldn  = 0.
     emis  = 0.
     gicewp = 0.
     gliqwp = 0.
     cicewp = 0.
     cliqwp = 0.
     stat_buffer(:ncol,:) = 0.
     cltot = 0.
     clhgh = 0.
     clmed = 0.
     cllow = 0.

! Recall previous step's rad statistics for correct averaging on a current time step::

     qrs = qrs1
     qrl = qrl1
     solin  (:)   = rad_buffer(:,1)
     fsds   (:)   = rad_buffer(:,2)
     fsdsc  (:)   = rad_buffer(:,3)
     fsnt   (:)   = rad_buffer(:,4)
     fsns   (:)   = rad_buffer(:,5)
     fsntc  (:)   = rad_buffer(:,6)
     fsnsc  (:)   = rad_buffer(:,7)
     fsntoa (:)   = rad_buffer(:,8)
     fsntoac(:)   = rad_buffer(:,9)
     sols   (:)   = rad_buffer(:,10)
     soll   (:)   = rad_buffer(:,11)
     solsd  (:)   = rad_buffer(:,12)
     solld  (:)   = rad_buffer(:,13)
     flnt   (:)   = rad_buffer(:,14)
     flut   (:)   = rad_buffer(:,15)
     flntc  (:)   = rad_buffer(:,16)
     flutc  (:)   = rad_buffer(:,17)
     flns   (:)   = rad_buffer(:,18)
     flnsc  (:)   = rad_buffer(:,19)
     flwds  (:)   = rad_buffer(:,20)
     fsutoa (:)   = rad_buffer(:,21)
     fsutoac(:)   = rad_buffer(:,22)
     flwdsc (:)   = rad_buffer(:,23)

! Initialize save-buffers

     rad_buffer(:,:) = 0.
     qrs1(:,:) = 0.
     qrl1(:,:) = 0.
!
! superparameterization radiation cycling starts here:
!

     do i = 1,ncol

        doinitial = .true.
        nbreak = nint(ztodt/RADTIME)

        do m=1,crm_nz
           k = pver-m+1
           qrs_crm(i,:,:,m) = qrs_crm(i,:,:,m) / state%pdel(i,k) ! for energy conservation
           qrl_crm(i,:,:,m) = qrl_crm(i,:,:,m) / state%pdel(i,k) ! for energy conservation
        end do

	tau00 = sqrt(taux(i)**2 + tauy(i)**2)
	wnd = sqrt(state%u(i,pver)**2 + state%v(i,pver)**2)
        bflx = shf(i)/cpair + 0.61*state%t(i,pver)*lhf(i)/latvap

        do  mm=1,nbreak
	  
   call t_startf ('crm_call')
	  call crm (lchnk, i, &
               state%t(i,:), state%q(i,:,1), state%q(i,:,ixcldliq), &
#ifdef SP_DIR_NS
               state%q(i,:,ixcldice), state%v(i,:), state%u(i,:),&
#else
               state%q(i,:,ixcldice), state%u(i,:), state%v(i,:),&
#endif
               state%ps(i), state%pmid(i,:), state%pdel(i,:), state%phis(i), &
               state%zm(i,:), state%zi(i,:), ztodt, pver, &
               ptend%q(i,:,1), ptend%q(i,:,ixcldliq), ptend%q(i,:,ixcldice), ptend%s(i,:), &
               u_crm(i,:,:,:), v_crm(i,:,:,:), w_crm(i,:,:,:), &
               t_crm(i,:,:,:), q_crm(i,:,:,:), qn_crm(i,:,:,:), qp_crm(i,:,:,:) ,&
               qc_crm(i,:,:,:), qi_crm(i,:,:,:), qpc_crm(i,:,:,:), qpi_crm(i,:,:,:),  &
               prec_crm(i,:,:), qrs_crm(i,:,:,:), qrl_crm(i,:,:,:), &
	       fsds_crm(i,:,:), fsns_crm(i,:,:), fsntoa_crm(i,:,:), fsutoa_crm(i,:,:),  &
	       flwds_crm(i,:,:), flns_crm(i,:,:), flut_crm(i,:,:),    &
	       fsntoac_crm(i,:,:), fsdsc_crm(i,:,:), flutc_crm(i,:,:), flnsc_crm(i,:,:), &
               t_rad, qv_rad(:,:,:), qc_rad(:,:,:), qi_rad(:,:,:), &
               precc(i), precl(i), precsc(i), precsl(i), &
               cltot(i), clhgh(i), clmed(i), cllow(i), &
               stat_buffer(i,:), cld(i,:), cldr(i,:), cldtop(i,:) , &
	       gicewp(i,:), gliqwp(i,:),     &
               mctot(i,:), mcup(i,:), mcdn(i,:), mcuup(i,:), mcudn(i,:), &
               crm_qc(i,:),crm_qi(i,:),crm_qs(i,:),crm_qg(i,:),crm_qr(i,:), &
               flux_qt(i,:), flux_qp(i,:), precflux(i,:), qt_ls(i,:), &
               qt_trans(i,:), qp_trans(i,:), qp_fall(i,:), &
               qp_evp(i,:), qp_src(i,:), t_ls(i,:), prectend(i), precstend(i), &
               nbreak, doinitial,  ocnfrac(i), wnd, tau00, bflx, taux_crm(i), tauy_crm(i), z0m(i), &
               timing_factor(i) &
               )

   call t_stopf('crm_call')
!
! CRM column-by-column radiation calculations
!
!
! Loop over CRM columns:

          do jj=1,crm_ny
           do ii=1,crm_nx

!
!          Compute liquid and ice water paths for a given CRM column
            do m=1,crm_nz
               k = pver-m+1
               trad(i,k) = t_rad(ii,jj,m)
               qvrad(i,k,1) = qv_rad(ii,jj,m)
               qtot = qc_rad(ii,jj,m) + qi_rad(ii,jj,m)
               if(qtot.gt.1.e-9) then
                 fice(i,k) = qi_rad(ii,jj,m)/qtot
                 cldn(i,k) = 0.999_r8
		 cld_crm(ii,jj,m)=0.999_r8
                 cicewp(i,k) = qi_rad(ii,jj,m)*state%pdel(i,k)/gravit*1000.0 &
		           / max(0.01_r8,cldn(i,k)) ! In-cloud ice water path.
                 cliqwp(i,k) = qc_rad(ii,jj,m)*state%pdel(i,k)/gravit*1000.0 &
		           / max(0.01_r8,cldn(i,k)) ! In-cloud liquid water path.
               else
                 fice(i,k)=0.
                 cldn(i,k)=0.
		 cld_crm(ii,jj,m)=0.
                 cicewp(i,k) = 0.           ! In-cloud ice water path.
                 cliqwp(i,k) = 0.           ! In-cloud liquid water path.
               end if
               cwp(i,k) = cicewp(i,k) + cliqwp(i,k)
	       cliqwp_crm(ii,jj,m)=cliqwp(i,k)
	       cicewp_crm(ii,jj,m)=cicewp(i,k)
            end do

!           Cloud water and ice particle sizes
            call cldefr(lchnk, i, i, landfrac, state%t, rel, rei, state%ps, state%pmid, landm, icefrac, snowh)

!           Cloud emissivity.
            call cldems(lchnk, i, i, cwp, fice, rei, emis)

            if(doisccp)then
             do m=1,crm_nz
               k = pver-m+1
	       rel_crm(ii,jj,m)=rel(i,k)
	       rei_crm(ii,jj,m)=rei(i,k)
	       emis_crm(ii,jj,m)=emis(i,k)
             end do
	    endif

            call cldovrlap(lchnk,i,i,state%pint,cldn,nmxrgn,pmxrgn)

!
!   Compute radiation:
   call t_startf ('crmrad_call')

            call radctl (lchnk, i, i, lwup, emis, state%pmid, &
                state%pint, state%lnpmid, state%lnpint, trad, qvrad, &
                cldn, cicewp, cliqwp, coszrs, asdir, asdif, &
                aldir, aldif, pmxrgn, nmxrgn, &
                fsns_crm(1,ii,jj), fsnt_crm(1,ii,jj), flns_crm(1,ii,jj), flnt_crm(1,ii,jj), &
                qrs_tmp, qrl_tmp, flwds_crm(1,ii,jj), rel, rei,  &
                sols_crm(1,ii,jj), soll_crm(1,ii,jj), solsd_crm(1,ii,jj), solld_crm(1,ii,jj),  &
                landfrac, state%zm, state, fsds_crm(1,ii,jj), &
                fsntoa_crm(1,ii,jj)  ,fsntoac_crm(1,ii,jj) ,fsdsc_crm(1,ii,jj),   &
		flwdsc_crm(1,ii,jj),fsntc_crm(1,ii,jj) ,fsnsc_crm(1,ii,jj), &
                fsutoa_crm(1,ii,jj) ,fsutoac_crm(1,ii,jj) ,flut_crm(1,ii,jj) , &
		flutc_crm(1,ii,jj) ,flntc_crm(1,ii,jj) ,flnsc_crm(1,ii,jj) ,solin_crm(1,ii,jj), &
                .false., .false.,.true.,.true.)
   call t_stopf('crmrad_call')
		
            do m=1,crm_nz
               k = pver-m+1
               qrs_crm(i,ii,jj,m) = (qrs_tmp(i,k)) / cpair
               qrl_crm(i,ii,jj,m) = (qrl_tmp(i,k)) / cpair
            end do

            if(ii.eq.1.and.jj.eq.1) then
              stat_buffer(i,10) = stat_buffer(i,10) + fsds_crm(i,ii,jj)
              stat_buffer(i,11) = stat_buffer(i,11) + flwds_crm(i,ii,jj)
              stat_buffer(i,12) = stat_buffer(i,12) + fsutoa_crm(i,ii,jj)
              stat_buffer(i,13) = stat_buffer(i,13) + flut_crm(i,ii,jj)
              stat_buffer(i,14) = stat_buffer(i,14) + fsdsc_crm(i,ii,jj)
              stat_buffer(i,15) = stat_buffer(i,15) + fsutoac_crm(i,ii,jj)
            end if

            if(mm.ne.nbreak) then

             qrs(i,:) = qrs(i,:) + qrs_tmp(i,:)  * state%pdel(i,:)
             qrl(i,:) = qrl(i,:) + qrl_tmp(i,:)  * state%pdel(i,:)
             solin(i) = solin(i) + solin_crm(i,ii,jj)
             fsnt(i) = fsnt(i) + fsnt_crm(i,ii,jj)
             fsns(i) = fsns(i) + fsns_crm(i,ii,jj)
             fsds(i) = fsds(i) + fsds_crm(i,ii,jj)
             fsdsc(i) = fsdsc(i) + fsdsc_crm(i,ii,jj)
             fsutoa(i) = fsutoa(i) + fsutoa_crm(i,ii,jj)
             fsutoac(i) = fsutoac(i) + fsutoac_crm(i,ii,jj)
             fsntc(i) = fsntc(i) + fsntc_crm(i,ii,jj)
             fsnsc(i) = fsnsc(i) + fsnsc_crm(i,ii,jj)
             fsntoa(i) = fsntoa(i) + fsntoa_crm(i,ii,jj)
             fsntoac(i) = fsntoac(i) + fsntoac_crm(i,ii,jj)
             sols(i) = sols(i) + sols_crm(i,ii,jj)
             soll(i) = soll(i) + soll_crm(i,ii,jj)
             solsd(i) = solsd(i) + solsd_crm(i,ii,jj)
             solld(i) = solld(i) + solld_crm(i,ii,jj)
             flut(i) = flut(i) + flut_crm(i,ii,jj)
             flutc(i) = flutc(i) + flutc_crm(i,ii,jj)
             flnt(i) = flnt(i) + flnt_crm(i,ii,jj)
             flns(i) = flns(i) + flns_crm(i,ii,jj)
             flntc(i) = flntc(i) + flntc_crm(i,ii,jj)
             flnsc(i) = flnsc(i) + flnsc_crm(i,ii,jj)
             flwds(i) = flwds(i) + flwds_crm(i,ii,jj)
             flwdsc(i) = flwdsc(i) + flwdsc_crm(i,ii,jj)

!
!   Save to the buffer for use in the next step for averaging:
!
            else

             qrs1(i,:) = qrs1(i,:) + qrs_tmp(i,:)
             qrl1(i,:) = qrl1(i,:) + qrl_tmp(i,:)
             rad_buffer(i,1) =  rad_buffer(i,1) + solin_crm(i,ii,jj)
             rad_buffer(i,2) =  rad_buffer(i,2) + fsds_crm(i,ii,jj)
             rad_buffer(i,3) =  rad_buffer(i,3) + fsdsc_crm(i,ii,jj)
             rad_buffer(i,4) =  rad_buffer(i,4) + fsnt_crm(i,ii,jj)
             rad_buffer(i,5) =  rad_buffer(i,5) + fsns_crm(i,ii,jj)
             rad_buffer(i,6) =  rad_buffer(i,6) + fsntc_crm(i,ii,jj)
             rad_buffer(i,7) =  rad_buffer(i,7) + fsnsc_crm(i,ii,jj)
             rad_buffer(i,8) =  rad_buffer(i,8) + fsntoa_crm(i,ii,jj)
             rad_buffer(i,9) =  rad_buffer(i,9) + fsntoac_crm(i,ii,jj)
             rad_buffer(i,10) = rad_buffer(i,10) + sols_crm(i,ii,jj)
             rad_buffer(i,11) = rad_buffer(i,11) + soll_crm(i,ii,jj)
             rad_buffer(i,12) = rad_buffer(i,12) + solsd_crm(i,ii,jj)
             rad_buffer(i,13) = rad_buffer(i,13) + solld_crm(i,ii,jj)
             rad_buffer(i,14) = rad_buffer(i,14) + flnt_crm(i,ii,jj)
             rad_buffer(i,15) = rad_buffer(i,15) + flut_crm(i,ii,jj)
             rad_buffer(i,16) = rad_buffer(i,16) + flntc_crm(i,ii,jj)
             rad_buffer(i,17) = rad_buffer(i,17) + flutc_crm(i,ii,jj)
             rad_buffer(i,18) = rad_buffer(i,18) + flns_crm(i,ii,jj)
             rad_buffer(i,19) = rad_buffer(i,19) + flnsc_crm(i,ii,jj)
             rad_buffer(i,20) = rad_buffer(i,20) + flwds_crm(i,ii,jj)
             rad_buffer(i,21) = rad_buffer(i,21) + fsutoa_crm(i,ii,jj)
             rad_buffer(i,22) = rad_buffer(i,22) + fsutoac_crm(i,ii,jj)
             rad_buffer(i,23) = rad_buffer(i,23) + flwdsc_crm(i,ii,jj)

           end if

          end do ! ii
         end do ! jj

        end do ! mm

        if (doisccp) then
            call crm_isccp (state%pmid(i,:), state%pint(i,:), qv_rad, t_rad, ts(i),    &
                cld_crm, cliqwp_crm, cicewp_crm, rel_crm, rei_crm, emis_crm, coszrs(i),    &
	        fq_isccp_s1(i,:), totalcldarea(i), lowcldarea(i), midcldarea(i), hghcldarea(i), &
                meantaucld(i), meanptop(i), meanttop(i), cloudy(i))
        end if

        qrs1(i,:) = qrs1(i,:) * state%pdel(i,:)
        qrl1(i,:) = qrl1(i,:) * state%pdel(i,:)
        coef = 1._r8/dble(nbreak*crm_nx*crm_ny)
        qrs(i,:) = qrs(i,:) * coef / state%pdel(i,:)
        qrl(i,:) = qrl(i,:) * coef / state%pdel(i,:)
        solin(i) = solin(i) * coef
        fsnt(i) = fsnt(i) * coef
        fsns(i) = fsns(i) * coef
        fsntoa(i) = fsntoa(i) * coef
        fsntoac(i) = fsntoac(i) * coef
        fsds(i) = fsds(i) * coef
        fsdsc(i) = fsdsc(i) * coef
        fsutoa(i) = fsutoa(i) * coef
        fsutoac(i) = fsutoac(i) * coef
        fsntc(i) = fsntc(i) * coef
        fsnsc(i) = fsnsc(i) * coef
        sols(i) = sols(i) * coef
        soll(i) = soll(i) * coef
        solsd(i) = solsd(i) * coef
        solld(i) = solld(i) * coef
        flut(i) = flut(i) * coef
        flutc(i) = flutc(i) * coef
        flntc(i) = flntc(i) * coef
        flnsc(i) = flnsc(i) * coef
        flnt(i) = flnt(i) * coef
        flns(i) = flns(i) * coef
        flwds(i) = flwds(i) * coef
        flwdsc(i) = flwdsc(i) * coef
        lwcf(i) = flutc(i) - flut(i)
        swcf(i) = fsntoa(i) - fsntoac(i)
        stat_buffer(i,10:15) = stat_buffer(i,10:15)/dble(nbreak)

     end do ! i

     if (doisccp) then
        if (any(coszrs(:ncol) > 0.)) then
           call outfld('FISCCP1 ',fq_isccp_s1, pcols,lchnk)
           call outfld('TCLDAREA',totalcldarea,pcols,lchnk)
           call outfld('LCLDAREA',lowcldarea,pcols,lchnk)
           call outfld('MCLDAREA',midcldarea,pcols,lchnk)
           call outfld('HCLDAREA',hghcldarea,pcols,lchnk)
           call outfld('MEANPTOP',meanptop    ,pcols,lchnk)
           call outfld('MEANTAU ',meantaucld  ,pcols,lchnk)
           call outfld('MEANTTOP',meanttop    ,pcols,lchnk)
           call outfld('CLOUDY  ',cloudy      ,pcols,lchnk)
        end if
     end if
!
!  subtract radiative heating tendency from the CRM tendency:
!  it will be added later:

     ptend%s(:ncol,pver-crm_nz+1:pver) = ptend%s(:ncol,pver-crm_nz+1:pver) -  &
        (qrs(:ncol,pver-crm_nz+1:pver) + qrl(:ncol,pver-crm_nz+1:pver))

     ptend%name  = 'crm'
     ptend%ls    = .TRUE.
     ptend%lq(1) = .TRUE.
     ptend%lq(ixcldliq) = .TRUE.
     ptend%lq(ixcldice) = .TRUE.



   call outfld('PRES    ',state%pmid ,pcols   ,lchnk   )
   call outfld('DPRES   ',state%pdel ,pcols   ,lchnk   )
   call outfld('HEIGHT  ',state%zm   ,pcols   ,lchnk   )
   call outfld('CRM_U   ',u_crm          ,pcols   ,lchnk   )
   call outfld('CRM_V   ',v_crm          ,pcols   ,lchnk   )
   call outfld('CRM_W   ',w_crm          ,pcols   ,lchnk   )
   call outfld('CRM_TABS',t_crm          ,pcols   ,lchnk   )
   call outfld('CRM_QV  ',(q_crm-qc_crm-qi_crm)*1000.,pcols   ,lchnk   )
   call outfld('CRM_QC  ',qc_crm*1000.   ,pcols   ,lchnk   )
   call outfld('CRM_QI  ',qi_crm*1000.   ,pcols   ,lchnk   )
   call outfld('CRM_QPC ',qpc_crm*1000.  ,pcols   ,lchnk   )
   call outfld('CRM_QPI ',qpi_crm*1000.  ,pcols   ,lchnk   )
   call outfld('CRM_PREC',prec_crm       ,pcols   ,lchnk   )
   call outfld('CRM_QRS ',qrs_crm        ,pcols   ,lchnk   )
   call outfld('CRM_QRL ',qrl_crm        ,pcols   ,lchnk   )
   call outfld('CRM_FSNT',fsntoa_crm     ,pcols   ,lchnk   )
   call outfld('CRMFSNTC',fsntoac_crm    ,pcols   ,lchnk   )
   call outfld('CRM_FSUT',fsutoa_crm     ,pcols   ,lchnk   )
   call outfld('CRM_FSNS',fsns_crm       ,pcols   ,lchnk   )
   call outfld('CRMFSDSC',fsdsc_crm      ,pcols   ,lchnk   )
   call outfld('CRM_FSDS',fsds_crm       ,pcols   ,lchnk   )
   call outfld('CRM_FLUT',flut_crm       ,pcols   ,lchnk   )
   call outfld('CRMFLUTC',flutc_crm      ,pcols   ,lchnk   )
   call outfld('CRM_FLNS',flns_crm       ,pcols   ,lchnk   )
   call outfld('CRMFLNSC',flnsc_crm      ,pcols   ,lchnk   )
   call outfld('CRM_FLDS',flwds_crm      ,pcols   ,lchnk   )


   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('SPDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('SPDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('SPDQC   ',ptend%q(1,1,ixcldliq) ,pcols   ,lchnk   )
   call outfld('SPDQI   ',ptend%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('SPMC    ',mctot          ,pcols   ,lchnk   )
   call outfld('SPMCUP  ',mcup           ,pcols   ,lchnk   )
   call outfld('SPMCDN  ',mcdn           ,pcols   ,lchnk   )
   call outfld('SPMCUUP ',mcuup          ,pcols   ,lchnk   )
   call outfld('SPMCUDN ',mcudn          ,pcols   ,lchnk   )
   call outfld('SPQC    ',crm_qc         ,pcols   ,lchnk   )
   call outfld('SPQI    ',crm_qi         ,pcols   ,lchnk   )
   call outfld('SPQS    ',crm_qs         ,pcols   ,lchnk   )
   call outfld('SPQG    ',crm_qg         ,pcols   ,lchnk   )
   call outfld('SPQR    ',crm_qr         ,pcols   ,lchnk   )
   call outfld('SPQTFLX ',flux_qt        ,pcols   ,lchnk   )
   call outfld('SPQPFLX ',flux_qp        ,pcols   ,lchnk   )
   call outfld('SPPFLX  ',precflux       ,pcols   ,lchnk   )
   call outfld('SPQTLS  ',qt_ls          ,pcols   ,lchnk   )
   call outfld('SPQTTR  ',qt_trans       ,pcols   ,lchnk   )
   call outfld('SPQPTR  ',qp_trans       ,pcols   ,lchnk   )
   call outfld('SPQPEVP ',qp_evp         ,pcols   ,lchnk   )
   call outfld('SPQPFALL',qp_fall        ,pcols   ,lchnk   )
   call outfld('SPQPSRC ',qp_src         ,pcols   ,lchnk   )
   call outfld('SPTLS   ',t_ls           ,pcols   ,lchnk   )

   call outfld('CLOUD   ',cld,  pcols,lchnk)
   call outfld('CLOUDR  ',cldr,  pcols,lchnk)
   call outfld('CLOUDTOP',cldtop, pcols,lchnk)
   call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
   call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
   call outfld('CLDMED  ',clmed  ,pcols,lchnk)
   call outfld('CLDLOW  ',cllow  ,pcols,lchnk)

   call outfld('CLDLOWR ',stat_buffer(:,1),pcols,lchnk)
   call outfld('CLDMEDR ',stat_buffer(:,2),pcols,lchnk)
   call outfld('CLDHGHR ',stat_buffer(:,3),pcols,lchnk)
   call outfld('CLDTOTR ',stat_buffer(:,4),pcols,lchnk)
   call outfld('LWPR    ',stat_buffer(:,5),pcols,lchnk)
   call outfld('IWPR    ',stat_buffer(:,6),pcols,lchnk)
   call outfld('LOBR    ',stat_buffer(:,7),pcols,lchnk)
   call outfld('HIBR    ',stat_buffer(:,8),pcols,lchnk)
   call outfld('PRECR   ',stat_buffer(:,9),pcols,lchnk)
   call outfld('FSDSR   ',stat_buffer(:,10),pcols,lchnk)
   call outfld('FLDSR   ',stat_buffer(:,11),pcols,lchnk)
   call outfld('FSUTR   ',stat_buffer(:,12),pcols,lchnk)
   call outfld('FLUTR   ',stat_buffer(:,13),pcols,lchnk)
   call outfld('FSDSCR  ',stat_buffer(:,14),pcols,lchnk)
   call outfld('FSUTCR  ',stat_buffer(:,15),pcols,lchnk)
   call outfld('CLDLOWRD',stat_buffer(:,16),pcols,lchnk)
   call outfld('CLDMEDRD',stat_buffer(:,17),pcols,lchnk)
   call outfld('CLDHGHRD',stat_buffer(:,18),pcols,lchnk)
   call outfld('CLDTOTRD',stat_buffer(:,19),pcols,lchnk)

   call outfld('QRS     ',qrs/cpair  ,pcols,lchnk)
   call outfld('SOLIN   ',solin ,pcols,lchnk)
   call outfld('FSDS    ',fsds  ,pcols,lchnk)
   call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
   call outfld('FSUTOA  ',fsutoa,pcols,lchnk)
   call outfld('FSUTOAC ',fsutoac,pcols,lchnk)
   call outfld('FSNT    ',fsnt  ,pcols,lchnk)
   call outfld('FSNS    ',fsns  ,pcols,lchnk)
   call outfld('FSNTC   ',fsntc ,pcols,lchnk)
   call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
   call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
   call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
   call outfld('SOLS    ',sols  ,pcols,lchnk)
   call outfld('SOLL    ',soll  ,pcols,lchnk)
   call outfld('SOLSD   ',solsd ,pcols,lchnk)
   call outfld('SOLLD   ',solld ,pcols,lchnk)
   call outfld('QRL     ',qrl/cpair ,pcols,lchnk)
   call outfld('FLNT    ',flnt  ,pcols,lchnk)
   call outfld('FLUT    ',flut  ,pcols,lchnk)
   call outfld('FLUTC   ',flutc ,pcols,lchnk)
   call outfld('FLNTC   ',flntc ,pcols,lchnk)
   call outfld('FLNS    ',flns  ,pcols,lchnk)
   call outfld('FLNSC   ',flnsc ,pcols,lchnk)
   call outfld('FLWDS   ',flwds ,pcols,lchnk)
   call outfld('FLWDSC  ',flwdsc,pcols,lchnk)
   call outfld('LWCF    ',lwcf  ,pcols,lchnk)
   call outfld('SWCF    ',swcf  ,pcols,lchnk)

   call outfld('Z0M     ',z0m  ,pcols,lchnk)
   call outfld('TAUX_CRM',taux_crm  ,pcols,lchnk)
   call outfld('TAUY_CRM',tauy_crm  ,pcols,lchnk)

   call outfld('TIMINGF ',timing_factor  ,pcols,lchnk)

! Compute energy and water integrals of input state

     call check_energy_timestep_init(state, tend, pbuf)
     call physics_update (state, tend, ptend, ztodt)
!    check energy integrals
     wtricesink(:ncol) = precc(:ncol) + precl(:ncol) + prectend(:ncol)*1.e-3 ! include precip storage term
     icesink(:ncol) = precsc(:ncol) + precsl(:ncol) + precstend(:ncol)*1.e-3   ! conversion of ice to snow
!     write(6,'(a,12e10.3)')'prect=',(prect(i),i=1,12)
     call check_energy_chng(state, tend, "crm", nstep, ztodt, zero, wtricesink, icesink, zero)

!
! Compute liquid water paths (for diagnostics only)
!
    tgicewp(:ncol) = 0.
    tgliqwp(:ncol) = 0.
    do k=1,pver
       do i = 1,ncol
          cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cld(i,k)) ! In-cloud ice water path.
          cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cld(i,k)) ! In-cloud liquid water path.
          tgicewp(i)  = tgicewp(i) + gicewp(i,k) ! grid cell mean ice water path.
          tgliqwp(i)  = tgliqwp(i) + gliqwp(i,k) ! grid cell mean ice water path.
       end do
    end do
    tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
    gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
    cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

    call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
    call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
    call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
    call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
    call outfld('ICLDLWP' ,cwp    , pcols,lchnk)


   end if ! (is_first_step())

   do m=1,crm_nz
      k = pver-m+1
      do i=1,ncol
         qrs_crm(i,:,:,m) = qrs_crm(i,:,:,m) * state%pdel(i,k) ! for energy conservation
         qrl_crm(i,:,:,m) = qrl_crm(i,:,:,m) * state%pdel(i,k) ! for energy conservation
   end do
   end do

   call t_stopf('crm')


   call diag_dynvar (lchnk, ncol, state)

!========================================================
!========================================================
!========================================================
! End of superparameterization zone.
#endif
!
! Compute net flux
! Since fsns, fsnt, flns, and flnt are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   do i=1,ncol
      tend%flx_net(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do
!
! Compute net radiative heating
!
   call radheat_net (state, ptend, qrl, qrs)
!
! Add radiation tendencies to cummulative model tendencies and update profiles
!
   call physics_update(state, tend, ptend, ztodt)

! check energy integrals
   call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, tend%flx_net)
!
! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   srfrad(:ncol) = fsns(:ncol) + flwds(:ncol)
   call outfld('SRFRAD  ',srfrad,pcols,lchnk)
!
! Save atmospheric fields to force surface models
!
   call srfxfer (lchnk, ncol, state%ps, state%u(1,pver), state%v(1,pver),    &
                 state%t(1,pver), state%q(1,pver,1), state%exner(1,pver), state%zm(1,pver), &
                    state%pmid,      &
                 state%rpdel(1,pver))

!---------------------------------------------------------------------------------------
! Save history variables. These should move to the appropriate parameterization interface
!---------------------------------------------------------------------------------------

   call outfld('PRECL   ',precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',precc   ,pcols   ,lchnk       )
   call outfld('PRECSL  ',precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',precsc  ,pcols   ,lchnk       )
   
   prect(:ncol) = precc(:ncol) + precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )

#if ( defined COUP_CSM )
   call outfld('PRECLav ',precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',precc   ,pcols   ,lchnk   )
#endif
!     
! Compute heating rate for dtheta/dt
!
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )

! convert radiative heating rates to Q*dp for energy conservation
   if (conserve_energy) then
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if
   

   return
end subroutine tphysbc
