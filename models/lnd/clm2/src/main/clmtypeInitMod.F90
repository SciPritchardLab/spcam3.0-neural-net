module clmtypeInitMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Initialize clmtype components to signaling nan 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use nanMod
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: clmtypeInit
!                              
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_pstate_type
  private :: init_pft_ecophys_constants
  private :: init_pft_DGVMecophys_constants
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_pdgvstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_landunit_pstate_type
  private :: init_gridcell_pstate_type
  private :: init_atm2lnd_state_type
  private :: init_lnd2atm_state_type
  private :: init_atm2lnd_flux_type
  private :: init_lnd2atm_flux_type
  private :: init_gridcell_wflux_type
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clmtypeInit
!
! !INTERFACE:
  subroutine clmtypeInit()
!
! !DESCRIPTION: 
! Initialize clmtype components to signaling nan 
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
! *%area, *%wt, *%wtxy, *%ixy, *%jxy, *%mxy, *%index1d, *%ifspecial
! *%ityplun, *%itypveg, *%ngridcells, *%nlandunits, *%ncolumns, *%npfts
!
! !USES:
  use clmpoint
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
    integer :: gi,li,ci,pi            !indices
    integer :: begg, endg             !per-process beginning and ending gridcell indices
    integer :: begl, endl             !per-process beginning and ending land unit indices
    integer :: begc, endc             !per-process beginning and ending column indices 
    integer :: begp, endp             !per-process beginning and ending pft indices
    type(gridcell_type), pointer :: g !pointer to derived subtype
    type(landunit_type), pointer :: l !pointer to derived subtype
    type(column_type)  , pointer :: c !pointer to derived subtype
    type(pft_type)     , pointer :: p !pointer to derived subtype
!------------------------------------------------------------------------

    ! Set beginning and ending indices

    begg = grid1d%beg
    endg = grid1d%end
    begl = land1d%beg
    endl = land1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begp = pfts1d%beg
    endp = pfts1d%end

    ! energy balance structures (all levels)
    
    call init_energy_balance_type(clm%mebal)

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_energy_balance_type(g%gebal)
    end do

    do li = begl,endl
       l => lpoint(li)%l
       call init_energy_balance_type(l%lebal)
    end do

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_energy_balance_type(c%cebal)
    end do

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_energy_balance_type(p%pebal)
    end do

    ! water balance structures (all levels)

    call init_water_balance_type(clm%mwbal) 

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_water_balance_type(g%gwbal)
    end do

    do li = begl,endl
       l => lpoint(li)%l
       call init_water_balance_type(l%lwbal)
    end do

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_water_balance_type(c%cwbal)
    end do

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_water_balance_type(p%pwbal)
    end do

    ! pft physical state variables at pft level 
    ! and averaged to the column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_pstate_type(p%pps)
    end do

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_pft_pstate_type(c%cps%pps_a)
    end do

    ! pft ecophysiological constants 

    call init_pft_ecophys_constants()

    ! pft DGVM-specific ecophysiological constants 
    
    call init_pft_DGVMecophys_constants()

    ! pft energy state variables at the pft level and averaged to the column

    do pi = begp,endp  
       p => ppoint(pi)%p
       call init_pft_estate_type(p%pes)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_estate_type(c%ces%pes_a)
    end do

    ! pft water state variables at the pft level and averaged to the column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_wstate_type(p%pws)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_wstate_type(c%cws%pws_a)
    end do

    ! pft DGVM state variables at pft level and averaged to column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_pdgvstate_type(p%pdgvs)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_pdgvstate_type(c%cdgvs%pdgvs_a)
    end do

    ! pft energy flux variables at pft level and averaged to column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_eflux_type(p%pef)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_eflux_type(c%cef%pef_a)
    end do

    ! pft momentum flux variables at pft level and averaged to the column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_mflux_type(p%pmf)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_mflux_type(c%cmf%pmf_a)
    end do

    ! pft water flux variables 

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_wflux_type(p%pwf)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_wflux_type(c%cwf%pwf_a)
    end do

    ! pft carbon flux variables at pft level and averaged to column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_cflux_type(p%pcf)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_cflux_type(c%ccf%pcf_a)
    end do

    ! pft VOC flux variables at pft level and averaged to column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_vflux_type(p%pvf)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_vflux_type(c%cvf%pvf_a)
    end do

    ! pft dust flux variables at pft level and averaged to column

    do pi = begp,endp
       p => ppoint(pi)%p
       call init_pft_dflux_type(p%pdf)
    end do

    do ci = begc,endc  
       c => cpoint(ci)%c
       call init_pft_dflux_type(c%cdf%pdf_a)
    end do

    ! column physical state variables at column level and averaged to
    ! landunit, gridcell and model level

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_pstate_type(c%cps)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_pstate_type(l%lps%cps_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_pstate_type(g%gps%cps_a)
    end do

    call init_column_pstate_type(clm%mps%cps_a)

    ! column energy state variables at column level and averaged to
    ! landunit, gridcell and model level
    
    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_estate_type(c%ces)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_estate_type(l%les%ces_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_estate_type(g%ges%ces_a)
    end do

    call init_column_estate_type(clm%mes%ces_a)

    ! column water state variables at column level and averaged to
    ! landunit, gridcell and model level

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_wstate_type(c%cws)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_wstate_type(l%lws%cws_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_wstate_type(g%gws%cws_a)
    end do

    call init_column_wstate_type(clm%mws%cws_a)

    ! column carbon state variables at column level and averaged to
    ! landunit, gridcell and model level

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_cstate_type(c%ccs)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_cstate_type(l%lcs%ccs_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_cstate_type(g%gcs%ccs_a)
    end do

    call init_column_cstate_type(clm%mcs%ccs_a)

    ! column energy flux variables at column level and averaged to
    ! landunit, gridcell and model level

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_eflux_type(c%cef)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_eflux_type(l%lef%cef_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_eflux_type(g%gef%cef_a)
    end do

    call init_column_eflux_type(clm%mef%cef_a)

    ! column water flux variables at column level and averaged to
    ! landunit, gridcell and model level

    do ci = begc,endc
       c => cpoint(ci)%c
       call init_column_wflux_type(c%cwf)
    end do

    do li= begl,endl
       l => lpoint(li)%l
       call init_column_wflux_type(l%lwf%cwf_a)
    end do

    do gi= begg,endg
       g => gpoint(gi)%g
       call init_column_wflux_type(g%gwf%cwf_a)
    end do

    call init_column_wflux_type(clm%mwf%cwf_a)

    ! land unit physical state variables 

    do li= begl,endl
       l => lpoint(li)%l
       call init_landunit_pstate_type(l%lps)
    end do

    ! gridcell physical state variables 

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_gridcell_pstate_type(g%gps)
    end do

    ! gridcell atmosphere->land state variables 

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_atm2lnd_state_type(g%a2ls)
    end do

    ! gridcell land->atmosphere state variables 
    
    do gi = begg,endg
       g => gpoint(gi)%g
       call init_lnd2atm_state_type(g%l2as)
    end do

    ! gridcell atmosphere->land flux variables 

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_atm2lnd_flux_type(g%a2lf)
    end do

    ! gridcell land->atmosphere flux variables  

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_lnd2atm_flux_type(g%l2af)
    end do

    ! gridcell: water flux variables

    do gi = begg,endg
       g => gpoint(gi)%g
       call init_gridcell_wflux_type(g%gwf)
    end do

  end subroutine clmtypeInit

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(ebal)
!
! !DESCRIPTION: 
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    type(energy_balance_type), target, intent(inout):: ebal 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ebal%errsoi = nan       
    ebal%errseb = nan       
    ebal%errsol = nan       
    ebal%errlon = nan       
    ebal%acc_errsoi = 0._r8 
    ebal%acc_errseb = 0._r8 
    ebal%acc_errsol = 0._r8 
  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(wbal)
!
! !DESCRIPTION: 
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    type(water_balance_type), target, intent(inout):: wbal 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    wbal%begwb = nan        
    wbal%endwb = nan        
    wbal%errh2o = nan       
    wbal%acc_errh2o = 0._r8 
  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(pps)
!
! !DESCRIPTION: 
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
    type (pft_pstate_type), target, intent(inout):: pps 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pps%frac_veg_nosno = bigint	   
    pps%frac_veg_nosno_alb = bigint  
    pps%rsw = nan			
    pps%emv = nan			
    pps%z0mv = nan		
    pps%z0hv = nan		
    pps%z0qv = nan		
    pps%rootfr(:nlevsoi) = nan
    pps%rootr (:nlevsoi) = nan
    pps%dewmx = nan		
    pps%rssun = nan		
    pps%rssha = nan		
    pps%laisun = nan		
    pps%laisha = nan		
    pps%btran = nan		
    pps%fsun = nan		
    pps%tlai = nan		
    pps%tsai = nan		
    pps%elai = nan		
    pps%esai = nan		
    pps%igs = nan			
    pps%stembio = nan		
    pps%rootbio = nan		
    pps%fwet = nan		
    pps%fdry = nan		
    pps%dt_veg = nan		
    pps%htop = nan		
    pps%hbot = nan		
    pps%z0m = nan			
    pps%displa = nan		
    pps%albd(:numrad) = nan		
    pps%albi(:numrad) = nan		
    pps%fabd(:numrad) = nan		
    pps%fabi(:numrad) = nan		
    pps%ftdd(:numrad) = nan		
    pps%ftid(:numrad) = nan		
    pps%ftii(:numrad) = nan		
    pps%ndvi = nan		
    pps%u10 = nan			
    pps%fv = nan			
    pps%ram1 = nan			
  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION: 
! Initialize pft physical state
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
    integer :: ivt     !vegetation type index
!------------------------------------------------------------------------
    do ivt = 0,numpft 
       pftcon(ivt)%ncorn = bigint		
       pftcon(ivt)%nwheat = bigint		
       pftcon(ivt)%noveg = bigint	        
       pftcon(ivt)%ntree = bigint		
       pftcon(ivt)%smpmax = -1.5e5           
       pftcon(ivt)%foln = nan	        
       pftcon(ivt)%dleaf = nan		
       pftcon(ivt)%c3psn = nan		
       pftcon(ivt)%vcmx25 = nan		
       pftcon(ivt)%mp = nan			
       pftcon(ivt)%qe25 = nan		
       pftcon(ivt)%xl = nan			
       pftcon(ivt)%rhol(:numrad) = nan        
       pftcon(ivt)%rhos(:numrad) = nan        
       pftcon(ivt)%taul(:numrad) = nan        
       pftcon(ivt)%taus(:numrad) = nan        
       pftcon(ivt)%z0mr = nan		
       pftcon(ivt)%displar = nan		
       pftcon(ivt)%roota_par = nan		
       pftcon(ivt)%rootb_par = nan		
       pftcon(ivt)%sla = nan			
    end do
  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: iinit_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION: 
! Initialize pft physical state
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
    integer :: ivt     !vegetation type index
!------------------------------------------------------------------------
    do ivt = 0,numpft 
       dgv_pftcon(ivt)%respcoeff = nan      
       dgv_pftcon(ivt)%flam = nan           
       dgv_pftcon(ivt)%resist = nan         
       dgv_pftcon(ivt)%l_turn = nan         
       dgv_pftcon(ivt)%l_long = nan         
       dgv_pftcon(ivt)%s_turn = nan         
       dgv_pftcon(ivt)%r_turn = nan         
       dgv_pftcon(ivt)%l_cton = nan         
       dgv_pftcon(ivt)%s_cton = nan         
       dgv_pftcon(ivt)%r_cton = nan         
       dgv_pftcon(ivt)%l_morph = nan        
       dgv_pftcon(ivt)%l_phen = nan         
       dgv_pftcon(ivt)%lmtorm  = nan        
       dgv_pftcon(ivt)%crownarea_max = nan  
       dgv_pftcon(ivt)%init_lai = nan       
       dgv_pftcon(ivt)%x  = nan             
       dgv_pftcon(ivt)%tcmin = nan          
       dgv_pftcon(ivt)%tcmax = nan          
       dgv_pftcon(ivt)%gddmin = nan         
       dgv_pftcon(ivt)%twmax = nan          
       dgv_pftcon(ivt)%lm_sapl = nan
       dgv_pftcon(ivt)%sm_sapl = nan
       dgv_pftcon(ivt)%hm_sapl = nan
       dgv_pftcon(ivt)%rm_sapl = nan
       dgv_pftcon(ivt)%tree = .false.
       dgv_pftcon(ivt)%summergreen = .false.
       dgv_pftcon(ivt)%raingreen = .false.
       dgv_pftcon(ivt)%reinickerp = nan
       dgv_pftcon(ivt)%wooddens = nan	 
       dgv_pftcon(ivt)%latosa = nan	    
       dgv_pftcon(ivt)%allom1 = nan	    
       dgv_pftcon(ivt)%allom2 = nan     
       dgv_pftcon(ivt)%allom3 = nan	    
    end do
  end subroutine init_pft_DGVMecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_energy_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(pes)
!
! !DESCRIPTION: 
! Initialize pft energy state
!
! !ARGUMENTS:
    implicit none
    type (pft_estate_type), target, intent(inout):: pes 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pes%t_ref2m = nan         
    pes%t_ref2m_min = nan
    pes%t_ref2m_max = nan
    pes%t_ref2m_min_inst = nan
    pes%t_ref2m_max_inst = nan
    pes%q_ref2m = nan         
    pes%t_veg = nan           
    pes%t_rad = nan	
  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(pws)
!
! !DESCRIPTION: 
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    type (pft_wstate_type), target, intent(inout):: pws !pft water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pws%h2ocan = nan  
  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(pdgvs)
!
! !DESCRIPTION: 
! Initialize pft DGVM state variables
!
! !ARGUMENTS:
    implicit none
    type (pft_dgvstate_type), target, intent(inout):: pdgvs 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pdgvs%agdd0 = nan              
    pdgvs%agdd5  = nan             
    pdgvs%agddtw = nan             
    pdgvs%agdd = nan               
    pdgvs%t10 = nan                
    pdgvs%t_mo = nan               
    pdgvs%t_mo_min = nan           
    pdgvs%fnpsn10 = nan            
    pdgvs%prec365 = nan            
    pdgvs%agdd20 = nan             
    pdgvs%tmomin20 = nan           
    pdgvs%t10min = nan             
    pdgvs%tsoi25 = nan             
    pdgvs%annpsn = nan             
    pdgvs%annpsnpot = nan          
    pdgvs%present = .false.            
    pdgvs%dphen = nan              
    pdgvs%leafon = nan             
    pdgvs%leafof = nan             
    pdgvs%nind = nan               
    pdgvs%lm_ind = nan             
    pdgvs%sm_ind = nan             
    pdgvs%hm_ind = nan             
    pdgvs%rm_ind = nan             
    pdgvs%lai_ind = nan            
    pdgvs%fpcinc = nan             
    pdgvs%fpcgrid = nan            
    pdgvs%crownarea = nan          
    pdgvs%bm_inc = nan             
    pdgvs%afmicr = nan             
    pdgvs%firelength  = nan        
    pdgvs%litterag = nan           
    pdgvs%litterbg = nan           
    pdgvs%cpool_fast = nan         
    pdgvs%cpool_slow = nan         
    pdgvs%k_fast_ave = nan         
    pdgvs%k_slow_ave = nan         
    pdgvs%litter_decom_ave = nan   
    pdgvs%turnover_ind = nan       
  end subroutine init_pft_pdgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(pef)
!
! !DESCRIPTION: 
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_eflux_type), target, intent(inout):: pef 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pef%sabg = nan		
    pef%sabv = nan		
    pef%fsa = nan			
    pef%fsr = nan			
    pef%parsun = nan      
    pef%parsha = nan      
    pef%dlrad = nan    
    pef%ulrad = nan    
    pef%eflx_lh_tot = nan		
    pef%eflx_lh_grnd = nan	
    pef%eflx_soil_grnd = nan	
    pef%eflx_sh_tot = nan 
    pef%eflx_sh_grnd = nan   
    pef%eflx_sh_veg = nan 
    pef%eflx_lh_vege = nan   
    pef%eflx_lh_vegt = nan   
    pef%cgrnd = nan    
    pef%cgrndl = nan      
    pef%cgrnds = nan      
    pef%eflx_gnet = nan    
    pef%dgnetdT = nan      
    pef%eflx_lwrad_out = nan	
    pef%eflx_lwrad_net = nan	
    pef%fsds_vis_d = nan
    pef%fsds_nir_d = nan
    pef%fsds_vis_i = nan
    pef%fsds_nir_i = nan
    pef%fsr_vis_d = nan
    pef%fsr_nir_d = nan
    pef%fsr_vis_i = nan
    pef%fsr_nir_i = nan
    pef%fsds_vis_d_ln = nan
    pef%fsds_nir_d_ln = nan
    pef%fsr_vis_d_ln = nan
    pef%fsr_nir_d_ln = nan
  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(pmf)
!
! !DESCRIPTION: 
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_mflux_type), target, intent(inout):: pmf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pmf%taux = nan  
    pmf%tauy = nan  
  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(pwf)
!
! !DESCRIPTION: 
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_wflux_type), target, intent(inout):: pwf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pwf%qflx_prec_intr = nan	
    pwf%qflx_prec_grnd = nan	
    pwf%qflx_rain_grnd = nan	
    pwf%qflx_snow_grnd = nan	
    pwf%qflx_snowcap = nan	
    pwf%qflx_evap_veg = nan	
    pwf%qflx_tran_veg = nan	
    pwf%qflx_evap_can = nan       
    pwf%qflx_evap_soi = nan	
    pwf%qflx_evap_tot = nan	
    pwf%qflx_evap_grnd = nan	
    pwf%qflx_dew_grnd = nan	
    pwf%qflx_sub_snow = nan	
    pwf%qflx_dew_snow = nan	
  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(pcf)
!
! !DESCRIPTION: 
! Initialize pft carbon flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_cflux_type), target, intent(inout):: pcf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pcf%psnsun = nan		
    pcf%psnsha = nan		
    pcf%fpsn = nan		
    pcf%frm = nan		
    pcf%frmf = nan		
    pcf%frms = nan		
    pcf%frmr = nan		
    pcf%frg = nan		
    pcf%dmi = nan		
    pcf%fco2 = nan		
    pcf%fmicr = nan		
  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(pvf)
!
! !DESCRIPTION: 
! Initialize pft VOC flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_vflux_type), target, intent(inout):: pvf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pvf%vocflx_tot = nan   	
    pvf%vocflx(:nvoc) = nan	
    pvf%vocflx_1 = nan      	
    pvf%vocflx_2 = nan      	
    pvf%vocflx_3 = nan      	
    pvf%vocflx_4 = nan      	
    pvf%vocflx_5 = nan      	
  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(pdf)
!
! !DESCRIPTION: 
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    type (pft_dflux_type), target, intent(inout):: pdf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    pdf%flx_mss_vrt_dst(:ndst) = nan 
    pdf%flx_mss_vrt_dst_tot = nan   
    pdf%vlc_trb(:ndst) = nan
    pdf%vlc_trb_1 = nan
    pdf%vlc_trb_2 = nan
    pdf%vlc_trb_3 = nan
    pdf%vlc_trb_4 = nan
  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(cps)
!
! !DESCRIPTION: 
! Initialize column physical state variables
!
! !ARGUMENTS:
    implicit none
    type (column_pstate_type), target, intent(inout):: cps 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    cps%isoicol = bigint		     
    cps%bsw   (:nlevsoi) = nan          
    cps%watsat(:nlevsoi) = nan          
    cps%hksat (:nlevsoi) = nan          
    cps%sucsat(:nlevsoi) = nan          
    cps%csol  (:nlevsoi) = nan          
    cps%tkmg  (:nlevsoi) = nan          
    cps%tkdry (:nlevsoi) = nan          
    cps%tksatu(:nlevsoi) = nan          
    cps%smpmin = nan		     
    cps%gwc_thr = nan		
    cps%mss_frc_cly_vld = nan	
    cps%mbl_bsn_fct = nan		
    cps%do_capsnow = .false.
    cps%snowdp = nan		
    cps%snowage = nan		
    cps%frac_sno = nan		
    cps%zi(-nlevsno+0:nlevsoi) = nan          
    cps%dz(-nlevsno+1:nlevsoi) = nan          
    cps%z (-nlevsno+1:nlevsoi) = nan          
    cps%frac_iceold(-nlevsno+1:nlevsoi) = nan 
    cps%imelt(-nlevsno+1:nlevsoi) = bigint    
    cps%eff_porosity(:nlevsoi) = nan           
    cps%sfact = nan		
    cps%sfactmax = nan		
    cps%emg = nan			
    cps%z0mg = nan		
    cps%z0hg = nan		
    cps%z0qg = nan		
    cps%htvp = nan		
    cps%beta = nan		
    cps%zii = nan			
    cps%ur = nan			
    cps%albgrd(:numrad) = nan 
    cps%albgri(:numrad) = nan 
    cps%rootr_column(:nlevsoi) = nan		
    cps%wf = nan                  
  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(ces)
!
! !DESCRIPTION: 
! Initialize column energy state variables
!
! !ARGUMENTS:
    implicit none
    type (column_estate_type), target, intent(inout):: ces 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ces%t_grnd = nan 			
    ces%dt_grnd	 = nan		        
    ces%t_rad = nan	
    ces%t_soisno(-nlevsno+1:nlevsoi) = nan
    ces%t_lake(1:nlevlak)= nan            
    ces%tssbef(-nlevsno+1:nlevsoi) = nan    
    ces%t_snow = nan			
    ces%thv = nan       			
    ces%thm = nan			        
  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(cws)
!
! !DESCRIPTION: 
! Initialize column water state variables
!
! !ARGUMENTS:
    implicit none
    type (column_wstate_type), target, intent(inout):: cws !column water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    cws%h2osno = nan			
    cws%h2osoi_liq(-nlevsno+1:nlevsoi) = nan 
    cws%h2osoi_ice(-nlevsno+1:nlevsoi) = nan 
    cws%h2osoi_vol(:nlevsoi) = nan
    cws%h2osno_old = nan			
    cws%qg = nan				
    cws%dqgdT = nan			
    cws%snowice = nan		
    cws%snowliq = nan		
  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(ccs)
!
! !DESCRIPTION: 
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    type (column_cstate_type), target, intent(inout):: ccs 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    ccs%soilc = nan			
  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(cef)
!
! !DESCRIPTION: 
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    type (column_eflux_type), target, intent(inout):: cef 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    cef%eflx_snomelt = nan	
    cef%eflx_impsoil = nan	
  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(cwf)
!
! !DESCRIPTION: 
! Initialize column water flux variables
!
! !ARGUMENTS:
    implicit none
    type (column_wflux_type), target, intent(inout):: cwf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    cwf%qflx_infl = nan		
    cwf%qflx_surf = nan		
    cwf%qflx_drain = nan		
    cwf%qflx_top_soil = nan	
    cwf%qflx_snomelt = nan	
    cwf%qflx_qrgwl = nan		
    cwf%qmelt = nan		
  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(lps)
!
! !DESCRIPTION: 
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    type (landunit_pstate_type), target, intent(inout):: lps 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    lps%lakpoi = .false.		     
  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(gps)
!
! !DESCRIPTION: 
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    type (gridcell_pstate_type), target, intent(inout):: gps 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    gps%lat = nan		    
    gps%lon = nan		    
    gps%latdeg = nan	    
    gps%londeg = nan	    
    gps%wtfact = nan	    
    gps%landfrac = nan
  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_state_type
!
! !INTERFACE:
  subroutine init_atm2lnd_state_type(a2ls)
!
! !DESCRIPTION: 
! Initialize atmospheric state variables required by the land
!
! !ARGUMENTS:
    implicit none
    type (atm2lnd_state_type), target, intent(inout):: a2ls 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
#if (defined OFFLINE)
    a2ls%flfall = nan
#endif
    a2ls%forc_t = nan		
    a2ls%forc_u = nan		
    a2ls%forc_v = nan		
    a2ls%forc_wind = nan       
    a2ls%forc_q = nan		
    a2ls%forc_hgt = nan	
    a2ls%forc_hgt_u = nan	
    a2ls%forc_hgt_t = nan	
    a2ls%forc_hgt_q = nan	
    a2ls%forc_pbot = nan	
    a2ls%forc_th = nan		
    a2ls%forc_vp = nan		
    a2ls%forc_rho = nan	
    a2ls%forc_co2 = nan	
    a2ls%forc_o2 = nan		
    a2ls%forc_psrf = nan	
  end subroutine init_atm2lnd_state_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_state_type
!
! !INTERFACE:
  subroutine init_lnd2atm_state_type(l2as)
!
! !DESCRIPTION: 
! Initialize land state variables required by the atmosphere
!
! !ARGUMENTS:
    implicit none
    type (lnd2atm_state_type), target, intent(inout):: l2as 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    l2as%t_rad = nan	
    l2as%t_ref2m = nan 	
    l2as%q_ref2m = nan 	
    l2as%h2osno = nan 	
    l2as%albd(:numrad) = nan
    l2as%albi(numrad) = nan
#if (defined DUST)
    l2as%fv = nan
    l2as%ram1 = nan
#endif
  end subroutine init_lnd2atm_state_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_flux_type
!
! !INTERFACE:
  subroutine init_atm2lnd_flux_type(a2lf)
!
! !DESCRIPTION: 
! Initialize atmospheric fluxes required by the land
!
! !ARGUMENTS:
    implicit none
    type (atm2lnd_flux_type), target, intent(inout):: a2lf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    a2lf%forc_lwrad = nan
    a2lf%forc_solad(numrad) = nan
    a2lf%forc_solai(numrad) = nan
    a2lf%forc_solar = nan        
    a2lf%forc_rain = nan		
    a2lf%forc_snow = nan		
  end subroutine init_atm2lnd_flux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_flux_type
!
! !INTERFACE:
  subroutine init_lnd2atm_flux_type(l2af)
!
! !DESCRIPTION: 
! Initialize land fluxes required by the atmosphere
!
! !ARGUMENTS:
    implicit none
    type (lnd2atm_flux_type), target, intent(inout):: l2af 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    l2af%taux = nan		
    l2af%tauy = nan		
    l2af%eflx_lwrad_out = nan	
    l2af%eflx_sh_tot = nan	
    l2af%eflx_lh_tot = nan	
    l2af%qflx_evap_tot = nan
    l2af%fsa = nan
#if (defined DUST)
    l2af%flxdst(:) = nan
#endif
  end subroutine init_lnd2atm_flux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(gwf)
!
! !DESCRIPTION: 
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    type (gridcell_wflux_type), target, intent(inout):: gwf 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    gwf%qchan2 = nan              
    gwf%qchocn2 = nan             
  end subroutine init_gridcell_wflux_type

end module clmtypeInitMod


