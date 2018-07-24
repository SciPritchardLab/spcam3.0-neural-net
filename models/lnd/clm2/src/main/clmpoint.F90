#include <misc.h>
#include <preproc.h>

module clmpoint

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmpoint
!
! !DESCRIPTION: 
! Module containing one-dimensional real and integer pointers for clmtype 
! entries. Use of these "pointer" arrays significantly speeds up model I/O.
!
! !USES:
  use clmtype
!
! !PUBLIC TYPES:
  implicit none
  save
!
! 1d real and integer pointers for clmtype entries
!
  type value_pointer_is
     integer , pointer :: p
  end type value_pointer_is

  type value_pointer_ia
     integer , dimension(:), pointer:: p
  end type value_pointer_ia

  type value_pointer_rs
     real(r8), pointer :: p
  end type value_pointer_rs

  type value_pointer_ra
     real(r8), dimension(:), pointer:: p
  end type value_pointer_ra

  type clmpoint_is
     type (value_pointer_is), dimension(:), pointer :: val
  end type clmpoint_is

  type clmpoint_ia
     type (value_pointer_ia), dimension(:), pointer :: val
  end type clmpoint_ia

  type clmpoint_rs
     type (value_pointer_rs), dimension(:), pointer :: val
  end type clmpoint_rs

  type clmpoint_ra
     type (value_pointer_ra), dimension(:), pointer :: val
  end type clmpoint_ra

  integer, parameter :: max_mapflds = 500
  type (clmpoint_rs) :: clmptr_rs(max_mapflds)
  type (clmpoint_ra) :: clmptr_ra(max_mapflds)

! 1d gridcell pointer
!
  type gridcell_single_pointer
     type(gridcell_type), pointer :: g
  end type gridcell_single_pointer
  type (gridcell_single_pointer), dimension(:), pointer :: gpoint
!
! 1d landunit pointer
!
  type landunit_single_pointer
     type(landunit_type), pointer :: l
  end type landunit_single_pointer
  type (landunit_single_pointer), dimension(:), pointer :: lpoint
!
! 1d column pointer
!
  type column_single_pointer
     type(column_type), pointer :: c
  end type column_single_pointer
  type (column_single_pointer), dimension(:), pointer :: cpoint
!
! 1d pft pointer
!
  type pft_single_pointer
     type(pft_type), pointer :: p
  end type pft_single_pointer
  type (pft_single_pointer), dimension(:), pointer :: ppoint
!
! pft indices
!
  integer, parameter :: not_set = -999
  integer :: lastindex = 1 

  integer :: ip_pps_rsw = not_set
  integer :: ip_pps_emv = not_set
  integer :: ip_pps_z0mv = not_set
  integer :: ip_pps_z0hv = not_set
  integer :: ip_pps_z0qv = not_set
  integer :: ip_pps_dewmx = not_set
  integer :: ip_pps_rssun = not_set
  integer :: ip_pps_rssha = not_set
  integer :: ip_pps_laisun = not_set
  integer :: ip_pps_laisha = not_set
  integer :: ip_pps_btran = not_set
  integer :: ip_pps_fsun = not_set
  integer :: ip_pps_tlai = not_set
  integer :: ip_pps_tsai = not_set
  integer :: ip_pps_elai = not_set
  integer :: ip_pps_esai = not_set
  integer :: ip_pps_igs = not_set
  integer :: ip_pps_stembio = not_set
  integer :: ip_pps_rootbio = not_set
  integer :: ip_pps_fwet = not_set
  integer :: ip_pps_fdry = not_set
  integer :: ip_pps_dt_veg = not_set
  integer :: ip_pps_htop = not_set
  integer :: ip_pps_hbot = not_set
  integer :: ip_pps_z0m = not_set
  integer :: ip_pps_displa = not_set
  integer :: ip_pps_ndvi = not_set
  integer :: ip_pps_u10 = not_set
  integer :: ip_pps_fv = not_set
  integer :: ip_pps_area = not_set
  integer :: ip_pps_wt = not_set
  integer :: ip_pps_wtxy = not_set
  integer :: ip_pps_rootfr = not_set
  integer :: ip_pps_rootr = not_set
  integer :: ip_pps_albd = not_set
  integer :: ip_pps_albi = not_set
  integer :: ip_pps_fabd = not_set
  integer :: ip_pps_fabi = not_set
  integer :: ip_pps_ftdd = not_set
  integer :: ip_pps_ftid = not_set
  integer :: ip_pps_ftii = not_set
  integer :: ip_pes_t_ref2m = not_set
  integer :: ip_pes_t_ref2m_min = not_set
  integer :: ip_pes_t_ref2m_max = not_set
  integer :: ip_pes_q_ref2m = not_set
  integer :: ip_pes_t_veg = not_set
  integer :: ip_pes_t_rad = not_set
  integer :: ip_pws_h2ocan = not_set
  integer :: ip_pef_sabg = not_set
  integer :: ip_pef_sabv = not_set
  integer :: ip_pef_fsa = not_set
  integer :: ip_pef_fsr = not_set
  integer :: ip_pef_fsds_vis_d = not_set
  integer :: ip_pef_fsds_nir_d = not_set
  integer :: ip_pef_fsds_vis_i = not_set
  integer :: ip_pef_fsds_nir_i = not_set
  integer :: ip_pef_fsds_vis_d_ln = not_set
  integer :: ip_pef_fsds_nir_d_ln = not_set
  integer :: ip_pef_fsr_vis_d = not_set
  integer :: ip_pef_fsr_nir_d = not_set
  integer :: ip_pef_fsr_vis_i = not_set
  integer :: ip_pef_fsr_nir_i = not_set
  integer :: ip_pef_fsr_vis_d_ln = not_set
  integer :: ip_pef_fsr_nir_d_ln = not_set
  integer :: ip_pef_parsun = not_set
  integer :: ip_pef_parsha = not_set
  integer :: ip_pef_dlrad = not_set
  integer :: ip_pef_ulrad = not_set
  integer :: ip_pef_eflx_lh_tot = not_set
  integer :: ip_pef_eflx_lh_grnd = not_set
  integer :: ip_pef_eflx_soil_grnd = not_set
  integer :: ip_pef_eflx_sh_tot = not_set
  integer :: ip_pef_eflx_sh_grnd = not_set
  integer :: ip_pef_eflx_sh_veg = not_set
  integer :: ip_pef_eflx_lh_vege = not_set
  integer :: ip_pef_eflx_lh_vegt = not_set
  integer :: ip_pef_cgrnd = not_set
  integer :: ip_pef_cgrndl = not_set
  integer :: ip_pef_cgrnds = not_set
  integer :: ip_pef_eflx_gnet = not_set
  integer :: ip_pef_dgnetdT = not_set
  integer :: ip_pef_eflx_lwrad_out = not_set
  integer :: ip_pef_eflx_lwrad_net = not_set
  integer :: ip_pmf_taux = not_set
  integer :: ip_pmf_tauy = not_set
  integer :: ip_pwf_qflx_prec_intr = not_set
  integer :: ip_pwf_qflx_prec_grnd = not_set
  integer :: ip_pwf_qflx_rain_grnd = not_set
  integer :: ip_pwf_qflx_snow_grnd = not_set
  integer :: ip_pwf_qflx_snowcap = not_set
  integer :: ip_pwf_qflx_evap_veg = not_set
  integer :: ip_pwf_qflx_tran_veg = not_set
  integer :: ip_pwf_qflx_evap_soi = not_set
  integer :: ip_pwf_qflx_evap_tot = not_set
  integer :: ip_pwf_qflx_evap_grnd = not_set
  integer :: ip_pwf_qflx_dew_grnd = not_set
  integer :: ip_pwf_qflx_sub_snow = not_set
  integer :: ip_pwf_qflx_dew_snow = not_set
  integer :: ip_pwf_qflx_evap_can = not_set
  integer :: ip_pcf_psnsun = not_set
  integer :: ip_pcf_psnsha = not_set
  integer :: ip_pcf_fpsn = not_set
  integer :: ip_pcf_frm = not_set
  integer :: ip_pcf_frmf = not_set
  integer :: ip_pcf_frms = not_set
  integer :: ip_pcf_frmr = not_set
  integer :: ip_pcf_frg = not_set
  integer :: ip_pcf_dmi = not_set
  integer :: ip_pcf_fco2 = not_set
  integer :: ip_pcf_fmicr = not_set
  integer :: ip_pvf_vocflx_tot = not_set
  integer :: ip_pvf_vocflx_1 = not_set
  integer :: ip_pvf_vocflx_2 = not_set
  integer :: ip_pvf_vocflx_3 = not_set
  integer :: ip_pvf_vocflx_4 = not_set
  integer :: ip_pvf_vocflx_5 = not_set
  integer :: ip_pdf_flx_mss_vrt_dst_tot = not_set
  integer :: ip_pdf_vlc_trb_1 = not_set
  integer :: ip_pdf_vlc_trb_2 = not_set
  integer :: ip_pdf_vlc_trb_3 = not_set
  integer :: ip_pdf_vlc_trb_4 = not_set
  integer :: ip_pdgvs_t_mo = not_set
  integer :: ip_pdgvs_t_mo_min = not_set
  integer :: ip_pdgvs_t10 = not_set
  integer :: ip_pdgvs_fnpsn10 = not_set
  integer :: ip_pdgvs_prec365 = not_set
  integer :: ip_pdgvs_agdd0 = not_set
  integer :: ip_pdgvs_agdd5 = not_set
  integer :: ip_pdgvs_agddtw = not_set
  integer :: ip_pdgvs_agdd = not_set
!
! column indices
!
  integer :: ic_cps_pps_a_rsw = not_set
  integer :: ic_cps_pps_a_emv = not_set
  integer :: ic_cps_pps_a_z0mv = not_set
  integer :: ic_cps_pps_a_z0hv = not_set
  integer :: ic_cps_pps_a_z0qv = not_set
  integer :: ic_cps_pps_a_dewmx = not_set
  integer :: ic_cps_pps_a_rssun = not_set
  integer :: ic_cps_pps_a_rssha = not_set
  integer :: ic_cps_pps_a_laisun = not_set
  integer :: ic_cps_pps_a_laisha = not_set
  integer :: ic_cps_pps_a_btran = not_set
  integer :: ic_cps_pps_a_fsun = not_set
  integer :: ic_cps_pps_a_tlai = not_set
  integer :: ic_cps_pps_a_tsai = not_set
  integer :: ic_cps_pps_a_elai = not_set
  integer :: ic_cps_pps_a_esai = not_set
  integer :: ic_cps_pps_a_igs = not_set
  integer :: ic_cps_pps_a_stembio = not_set
  integer :: ic_cps_pps_a_rootbio = not_set
  integer :: ic_cps_pps_a_fwet = not_set
  integer :: ic_cps_pps_a_fdry = not_set
  integer :: ic_cps_pps_a_dt_veg = not_set
  integer :: ic_cps_pps_a_htop = not_set
  integer :: ic_cps_pps_a_hbot = not_set
  integer :: ic_cps_pps_a_z0m = not_set
  integer :: ic_cps_pps_a_displa = not_set
  integer :: ic_cps_pps_a_ndvi = not_set
  integer :: ic_cps_pps_a_u10 = not_set
  integer :: ic_cps_pps_a_fv = not_set
  integer :: ic_cps_pps_a_area = not_set
  integer :: ic_cps_pps_a_wt = not_set
  integer :: ic_cps_pps_a_wtxy = not_set
  integer :: ic_cps_vwc_thr = not_set
  integer :: ic_cps_mss_frc_cly_vld = not_set
  integer :: ic_cps_mbl_bsn_fct = not_set
  integer :: ic_cps_snowdp = not_set
  integer :: ic_cps_snowage = not_set
  integer :: ic_cps_frac_sno = not_set
  integer :: ic_cps_sfact = not_set
  integer :: ic_cps_sfactmax = not_set
  integer :: ic_cps_emg = not_set
  integer :: ic_cps_z0mg = not_set
  integer :: ic_cps_z0hg = not_set
  integer :: ic_cps_z0qg = not_set
  integer :: ic_cps_htvp = not_set
  integer :: ic_cps_beta = not_set
  integer :: ic_cps_zii = not_set
  integer :: ic_cps_ur = not_set
  integer :: ic_cps_wf = not_set
  integer :: ic_cps_area = not_set
  integer :: ic_cps_wt = not_set
  integer :: ic_cps_wtxy = not_set
  integer :: ic_cps_pps_a_rootfr = not_set
  integer :: ic_cps_pps_a_rootr = not_set
  integer :: ic_cps_frac_iceold = not_set
  integer :: ic_cps_eff_porosity = not_set
  integer :: ic_cps_rootr_column = not_set
  integer :: ic_cps_albd = not_set
  integer :: ic_cps_albi = not_set
  integer :: ic_cps_fabd = not_set
  integer :: ic_cps_fabi = not_set
  integer :: ic_cps_ftdd = not_set
  integer :: ic_cps_ftid = not_set
  integer :: ic_cps_ftii = not_set
  integer :: ic_cps_albgrd = not_set
  integer :: ic_cps_albgri = not_set
  integer :: ic_ces_pes_a_t_ref2m = not_set
  integer :: ic_ces_pes_a_q_ref2m = not_set
  integer :: ic_ces_pes_a_t_veg = not_set
  integer :: ic_ces_t_grnd = not_set
  integer :: ic_ces_dt_grnd = not_set
  integer :: ic_ces_t_rad = not_set
  integer :: ic_ces_t_snow = not_set
  integer :: ic_ces_thv = not_set
  integer :: ic_ces_thm = not_set
  integer :: ic_ces_t_lake = not_set
  integer :: ic_ces_t_soisno = not_set
  integer :: ic_ces_tssbef = not_set
  integer :: ic_cws_pws_a_h2ocan = not_set
  integer :: ic_cws_h2osno = not_set
  integer :: ic_cws_qg = not_set
  integer :: ic_cws_dqgdT = not_set
  integer :: ic_cws_snowice = not_set
  integer :: ic_cws_snowliq = not_set
  integer :: ic_cws_h2osoi_liq = not_set
  integer :: ic_cws_h2osoi_ice = not_set
  integer :: ic_cws_h2osoi_vol = not_set
  integer :: ic_ccs_soilc = not_set
  integer :: ic_cef_pef_a_sabg = not_set
  integer :: ic_cef_pef_a_sabv = not_set
  integer :: ic_cef_pef_a_fsds_vis_d = not_set
  integer :: ic_cef_pef_a_fsds_nir_d = not_set
  integer :: ic_cef_pef_a_fsds_vis_i = not_set
  integer :: ic_cef_pef_a_fsds_nir_i = not_set
  integer :: ic_cef_pef_a_fsds_vis_d_ln = not_set
  integer :: ic_cef_pef_a_fsds_nir_d_ln = not_set
  integer :: ic_cef_pef_a_fsr_vis_d = not_set
  integer :: ic_cef_pef_a_fsr_nir_d = not_set
  integer :: ic_cef_pef_a_fsr_vis_i = not_set
  integer :: ic_cef_pef_a_fsr_nir_i = not_set
  integer :: ic_cef_pef_a_fsr_vis_d_ln= not_set
  integer :: ic_cef_pef_a_fsr_nir_d_ln= not_set
  integer :: ic_cef_pef_a_fsa = not_set
  integer :: ic_cef_pef_a_fsr = not_set
  integer :: ic_cef_pef_a_parsun = not_set
  integer :: ic_cef_pef_a_parsha = not_set
  integer :: ic_cef_pef_a_dlrad = not_set
  integer :: ic_cef_pef_a_ulrad = not_set
  integer :: ic_cef_pef_a_eflx_lh_tot = not_set
  integer :: ic_cef_pef_a_eflx_lh_grnd = not_set
  integer :: ic_cef_pef_a_eflx_soil_grnd = not_set
  integer :: ic_cef_pef_a_eflx_sh_tot = not_set
  integer :: ic_cef_pef_a_eflx_sh_grnd = not_set
  integer :: ic_cef_pef_a_eflx_sh_veg = not_set
  integer :: ic_cef_pef_a_eflx_lh_vege = not_set
  integer :: ic_cef_pef_a_eflx_lh_vegt = not_set
  integer :: ic_cef_pef_a_cgrnd = not_set
  integer :: ic_cef_pef_a_cgrndl = not_set
  integer :: ic_cef_pef_a_cgrnds = not_set
  integer :: ic_cef_pef_a_eflx_gnet = not_set
  integer :: ic_cef_pef_a_dgnetdT = not_set
  integer :: ic_cef_pef_a_eflx_lwrad_out = not_set
  integer :: ic_cef_pef_a_eflx_lwrad_net = not_set
  integer :: ic_cef_eflx_snomelt = not_set
  integer :: ic_cef_eflx_impsoil = not_set
  integer :: ic_cmf_pmf_a_taux = not_set
  integer :: ic_cmf_pmf_a_tauy = not_set
  integer :: ic_cwf_pwf_a_qflx_prec_intr = not_set
  integer :: ic_cwf_pwf_a_qflx_prec_grnd = not_set
  integer :: ic_cwf_pwf_a_qflx_rain_grnd = not_set
  integer :: ic_cwf_pwf_a_qflx_snow_grnd = not_set
  integer :: ic_cwf_pwf_a_qflx_snowcap = not_set
  integer :: ic_cwf_pwf_a_qflx_evap_veg = not_set
  integer :: ic_cwf_pwf_a_qflx_tran_veg = not_set
  integer :: ic_cwf_pwf_a_qflx_evap_soi = not_set
  integer :: ic_cwf_pwf_a_qflx_evap_tot = not_set
  integer :: ic_cwf_pwf_a_qflx_evap_grnd = not_set
  integer :: ic_cwf_pwf_a_qflx_dew_grnd = not_set
  integer :: ic_cwf_pwf_a_qflx_sub_snow = not_set
  integer :: ic_cwf_pwf_a_qflx_dew_snow = not_set
  integer :: ic_cwf_pwf_a_qflx_evap_can = not_set
  integer :: ic_cwf_qflx_infl = not_set
  integer :: ic_cwf_qflx_surf = not_set
  integer :: ic_cwf_qflx_drain = not_set
  integer :: ic_cwf_qflx_top_soil = not_set
  integer :: ic_cwf_qflx_snomelt = not_set
  integer :: ic_cwf_qflx_qrgwl = not_set
  integer :: ic_cwf_qmelt = not_set
  integer :: ic_ccf_pcf_a_psnsun = not_set
  integer :: ic_ccf_pcf_a_psnsha = not_set
  integer :: ic_ccf_pcf_a_fpsn = not_set
  integer :: ic_ccf_pcf_a_frm = not_set
  integer :: ic_ccf_pcf_a_frmf = not_set
  integer :: ic_ccf_pcf_a_frms = not_set
  integer :: ic_ccf_pcf_a_frmr = not_set
  integer :: ic_ccf_pcf_a_frg = not_set
  integer :: ic_ccf_pcf_a_dmi = not_set
  integer :: ic_ccf_pcf_a_fco2 = not_set
  integer :: ic_ccf_pcf_a_fmicr = not_set
  integer :: ic_cebal_errsoi = not_set
  integer :: ic_cebal_errseb = not_set
  integer :: ic_cebal_errsol = not_set
  integer :: ic_cwbal_errh2o = not_set
!
! grid cell indices
!
  integer :: ig_gwf_qchan2 = not_set
  integer :: ig_gwf_qchocn2 = not_set
  integer :: ig_a2ls_forc_t = not_set
  integer :: ig_a2ls_forc_u = not_set
  integer :: ig_a2ls_forc_v = not_set
  integer :: ig_a2ls_forc_wind = not_set
  integer :: ig_a2ls_forc_q = not_set
  integer :: ig_a2ls_forc_hgt = not_set
  integer :: ig_a2ls_forc_hgt_u = not_set
  integer :: ig_a2ls_forc_hgt_t = not_set
  integer :: ig_a2ls_forc_hgt_q = not_set
  integer :: ig_a2ls_forc_pbot = not_set
  integer :: ig_a2ls_forc_th = not_set
  integer :: ig_a2ls_forc_vp = not_set
  integer :: ig_a2ls_forc_rho = not_set
  integer :: ig_a2ls_forc_co2 = not_set
  integer :: ig_a2ls_forc_o2 = not_set
  integer :: ig_a2ls_forc_psrf = not_set
  integer :: ig_a2lf_forc_lwrad = not_set
  integer :: ig_a2lf_forc_rain = not_set
  integer :: ig_a2lf_forc_snow = not_set
  integer :: ig_a2lf_forc_solad = not_set
  integer :: ig_a2lf_forc_solai = not_set
  integer :: ig_a2lf_forc_solar = not_set
  integer :: ig_l2as_t_rad = not_set
  integer :: ig_l2as_t_ref2m = not_set
  integer :: ig_l2as_q_ref2m = not_set
  integer :: ig_l2as_h2osno = not_set
  integer :: ig_l2as_albd = not_set
  integer :: ig_l2as_albi = not_set
  integer :: ig_l2af_taux = not_set
  integer :: ig_l2af_tauy = not_set
  integer :: ig_l2af_eflx_lh_tot = not_set
  integer :: ig_l2af_eflx_sh_tot = not_set
  integer :: ig_l2af_eflx_lwrad_out = not_set
  integer :: ig_l2af_qflx_evap_tot = not_set
  integer :: ig_l2af_fsa = not_set
!
! !PUBLIC MEMBER FUNCTIONS:
  public clmpoint_init  ! Initialize pointer arrays and allocate necessary memory
  public pointer_index  ! Set value for above pointer values 
!                              
! !REVISION HISTORY:
! Created by Mariana Vertenstein 
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS
!----------------------------------------------------------------------- 

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clmpoint_init
!
! !INTERFACE:
  subroutine clmpoint_init()
!
! !DESCRIPTION: 
! Initialize pointer arrays and allocate necessary memory
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: gi,li,ci,pi               !indices
    integer :: gindex                    !gridcell index
    integer :: lindex                    !land unit index
    integer :: cindex	 		 !column index
    integer :: pindex			 !pft index
    integer :: begg, endg                !per-process beginning and ending gridcell indices
    integer :: begl, endl                !per-process beginning and ending land unit indices
    integer :: begc, endc                !per-process beginning and ending column indices 
    integer :: begp, endp                !per-process beginning and ending pft indices
    integer :: ier                       !error status 
    type(gridcell_type), pointer :: g 	 !local pointer to derived subtype
    type(landunit_type), pointer :: l 	 !local pointer to derived subtype
    type(column_type)  , pointer :: c 	 !local pointer to derived subtype
    type(pft_type)     , pointer :: p 	 !local pointer to derived subtype
!------------------------------------------------------------------------

    ! Determine per-process beginning and ending indices

    begp = pfts1d%beg
    endp = pfts1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begl = land1d%beg
    endl = land1d%end
    begg = grid1d%beg
    endg = grid1d%end

    ! Allocate gridcell, landunit, column and pft pointers
    
    allocate (gpoint(begg:endg), &
              lpoint(begl:endl), &
              cpoint(begc:endc), &
              ppoint(begp:endp), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clmpoint_init() : gpoint,lpoint,cpoint,ppoint allocation error'
       call endrun
    end if
    
    ! Set up gridcell, landunit, column and pft pointers
    
    do gi = 1,clm%mps%ngridcells
       g => clm%g(gi)
       gindex = g%gps%index1d
       gpoint(gindex)%g => g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lindex = l%lps%index1d
          lpoint(lindex)%l => l
          do ci = 1,l%lps%ncolumns
             c => l%c(ci)
             cindex = c%cps%index1d
             cpoint(cindex)%c => c
             do pi = 1,c%cps%npfts 
                p => c%p(pi)
                pindex = p%pps%index1d
                ppoint(pindex)%p => p
             end do
          end do
       end do
    end do

  end subroutine clmpoint_init

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pointer_index
!
! !INTERFACE:
  integer function pointer_index ()
!
! !DESCRIPTION: 
! Set the current pointer index and increment the value of the index.
!
! !ARGUMENTS:
   implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

   pointer_index = lastindex
   lastindex = lastindex + 1
   if (lastindex > max_mapflds) then
      write(6,*)' pointer_index error: ',&
           ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
      call endrun
   endif

  end function pointer_index

end module clmpoint
