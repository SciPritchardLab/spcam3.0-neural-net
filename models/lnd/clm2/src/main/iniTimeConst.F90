#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: iniTimeConst
!
! !INTERFACE:
subroutine iniTimeConst 
!
! !DESCRIPTION: 
! Initialize time invariant clm variables
! 1) removed references to shallow lake - since it is not used
! 2) ***Make c%z, c%zi and c%dz allocatable depending on if you
!    have lake or soil
! 3) rootfr only initialized for soil points
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clmtype
  use globals
  use clm_varpar, only : nlevsoi, nlevlak, lsmlon, lsmlat, numpft
  use clm_varcon, only : istice, istdlak, istwet, isturb, &
       zlak, dzlak, zsoi, dzsoi, zisoi, spval 
  use pftvarcon, only : ncorn, nwheat, noveg, ntree, roota_par, rootb_par,  &
       z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
       qe25, vcmx25, mp, c3psn, &
       pftpar , tree   , summergreen, raingreen  , sla     , &
       lm_sapl, sm_sapl, hm_sapl    , rm_sapl    , latosa  , &
       allom1 , allom2 , allom3     , reinickerp , wooddens
  use clm_varsur, only : soic2d, sand3d, clay3d, latixy, longxy, landfrac
  use clm_varctl, only : nsrest
  use time_manager, only : get_step_size
  use shr_const_mod, only : SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod.
!
! !REVISION HISTORY:
! Created by Gordon Bonan. 
! Updated to clm2.1 data structrues by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: i,j,ib,lev          !indices
  integer  :: gi,li,ci,pi         !indices 
  integer  :: ltype               !land unit type index
  integer  :: ivt                 !vegetation type index
  real(r8) :: bd                  !bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm                 !mineral conductivity
  real(r8) :: xksat               !maximum hydraulic conductivity of soil [mm/s]
  real(r8) :: scalez = 0.025      !Soil layer thickness discretization (m)     
  real(r8) :: hkdepth = 0.5       !Length scale for Ksat decrease (m)
  real(r8) :: watsat,clay,sand,zi,tkmg       !temporaries
  type(gridcell_type)       , pointer :: g   !local pointer to derived subtype
  type(landunit_type)       , pointer :: l   !local pointer to derived subtype
  type(column_type)         , pointer :: c   !local pointer to derived subtype
  type(pft_type)            , pointer :: p   !local pointer to derived subtype
  type(gridcell_pstate_type), pointer :: gps !local pointer to derived subtype
  type(column_pstate_type)  , pointer :: cps !local pointer to derived subtype
  type(landunit_pstate_type), pointer :: lps !local pointer to derived subtype
  type(pft_pstate_type)     , pointer :: pps !local pointer to derived subtype
  type(pft_epc_type)   , pointer :: pepc     !pointer to array of ecophysiological constants 
  type(pft_dgvepc_type), pointer :: pdgvepc  !pointer to array of DGVM ecophysiological constants 
!-----------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! Initialize time constant arrays of ecophysiological constants and
  ! arrays of dgvm ecophysiological constants 
  ! --------------------------------------------------------------------

   do ivt = 0,numpft 
      pftcon(ivt)%ncorn = ncorn
      pftcon(ivt)%nwheat = nwheat
      pftcon(ivt)%noveg = noveg
      pftcon(ivt)%ntree = ntree
      pftcon(ivt)%z0mr = z0mr(ivt)
      pftcon(ivt)%displar = displar(ivt)
      pftcon(ivt)%dleaf = dleaf(ivt)
      pftcon(ivt)%xl = xl(ivt)
      do ib = 1,numrad 
         pftcon(ivt)%rhol(ib) = rhol(ivt,ib)
         pftcon(ivt)%rhos(ib) = rhos(ivt,ib)
         pftcon(ivt)%taul(ib) = taul(ivt,ib)
         pftcon(ivt)%taus(ib) = taus(ivt,ib)
      end do
      pftcon(ivt)%qe25 = qe25(ivt)      
      pftcon(ivt)%vcmx25 = vcmx25(ivt)  
      pftcon(ivt)%mp = mp(ivt)          
      pftcon(ivt)%c3psn = c3psn(ivt)    
      pftcon(ivt)%sla = sla(ivt)
   end do

   do ivt = 0,numpft 
      dgv_pftcon(ivt)%respcoeff = pftpar(ivt,5)
      dgv_pftcon(ivt)%flam = pftpar(ivt,6)
      dgv_pftcon(ivt)%resist = pftpar(ivt,8)
      dgv_pftcon(ivt)%l_turn = pftpar(ivt,9)
      dgv_pftcon(ivt)%l_long = pftpar(ivt,10)
      dgv_pftcon(ivt)%s_turn = pftpar(ivt,11)
      dgv_pftcon(ivt)%r_turn = pftpar(ivt,12)
      dgv_pftcon(ivt)%l_cton = pftpar(ivt,13)
      dgv_pftcon(ivt)%s_cton = pftpar(ivt,14)
      dgv_pftcon(ivt)%r_cton = pftpar(ivt,15)
      dgv_pftcon(ivt)%l_morph = pftpar(ivt,16)
      dgv_pftcon(ivt)%l_phen = pftpar(ivt,17)
      dgv_pftcon(ivt)%lmtorm = pftpar(ivt,18)
      dgv_pftcon(ivt)%crownarea_max = pftpar(ivt,20)
      dgv_pftcon(ivt)%init_lai = pftpar(ivt,21)
      dgv_pftcon(ivt)%x = pftpar(ivt,22)
      dgv_pftcon(ivt)%tcmin = pftpar(ivt,28)
      dgv_pftcon(ivt)%tcmax = pftpar(ivt,29)
      dgv_pftcon(ivt)%gddmin = pftpar(ivt,30)
      dgv_pftcon(ivt)%twmax = pftpar(ivt,31)
      dgv_pftcon(ivt)%lm_sapl = lm_sapl(ivt)
      dgv_pftcon(ivt)%sm_sapl = sm_sapl(ivt)
      dgv_pftcon(ivt)%hm_sapl = hm_sapl(ivt)
      dgv_pftcon(ivt)%rm_sapl = rm_sapl(ivt)
      dgv_pftcon(ivt)%tree = tree(ivt)
      dgv_pftcon(ivt)%summergreen = summergreen(ivt)
      dgv_pftcon(ivt)%raingreen = raingreen(ivt)
      dgv_pftcon(ivt)%reinickerp = reinickerp
      dgv_pftcon(ivt)%wooddens = wooddens
      dgv_pftcon(ivt)%latosa = latosa
      dgv_pftcon(ivt)%allom1 = allom1
      dgv_pftcon(ivt)%allom2 = allom2
      dgv_pftcon(ivt)%allom3 = allom3
   end do

   ! --------------------------------------------------------------------
   ! Define layer structure for soil and lakes 
   ! Vertical profile of snow is initialized in routine iniTimeVar
   ! --------------------------------------------------------------------

   ! check that lake and soil levels are the same for now

   if (nlevlak /= nlevsoi) then
      write(6,*)'number of soil levels and number of lake levels must be the same'
      write(6,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
      call endrun
   endif

   ! Lake layers (assumed same for all lake patches)
   
   dzlak(1) = 0.1
   dzlak(2) = 1.
   dzlak(3) = 2.
   dzlak(4) = 3.
   dzlak(5) = 4.
   dzlak(6) = 5.
   dzlak(7) = 7.
   dzlak(8) = 7.
   dzlak(9) = 10.45
   dzlak(10)= 10.45
   
   zlak(1) =  0.05
   zlak(2) =  0.6
   zlak(3) =  2.1
   zlak(4) =  4.6
   zlak(5) =  8.1
   zlak(6) = 12.6
   zlak(7) = 18.6
   zlak(8) = 25.6
   zlak(9) = 34.325
   zlak(10)= 44.775

   ! Soil layers and interfaces (assumed same for all non-lake patches)
   ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil

   do j = 1, nlevsoi
      zsoi(j) = scalez*(exp(0.5*(j-0.5))-1.)    !node depths
   enddo
   
   dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
   do j = 2,nlevsoi-1
      dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1)) 
   enddo
   dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)
   
   zisoi(0) = 0.
   do j = 1, nlevsoi-1
      zisoi(j) = 0.5*(zsoi(j)+zsoi(j+1))         !interface depths
   enddo
   zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5*dzsoi(nlevsoi)
   
   ! --------------------------------------------------------------------
   ! Initialize soil and lake levels
   ! Initialize soil color, thermal and hydraulic properties 
   ! --------------------------------------------------------------------

   ! Grid level initialization
   do gi = 1,clm%mps%ngridcells

      ! Assign local pointer for simpler referencing
      g => clm%g(gi)					
      gps => g%gps
      
      ! Set xy indices
      i = gps%ixy
      j = gps%jxy
      
      ! Set latitudes and longitudes
      gps%lat    = latixy(i,j) * SHR_CONST_PI/180.  
      gps%lon    = longxy(i,j) * SHR_CONST_PI/180.
      gps%latdeg = latixy(i,j) 
      gps%londeg = longxy(i,j) 
      gps%landfrac = landfrac(i,j)
      
      ! Initialize fraction of model area with high water table 
      gps%wtfact = 0.3    
      
      ! Landunit level initialization
      do li = 1,gps%nlandunits
         
         ! Assign local pointer for simpler referencing
         l => g%l(li)
         lps => l%lps
         
         ! Set landunit type type
         ltype = lps%itype
         
         ! Define lake type
         if (ltype == istdlak) then
            lps%lakpoi  = .true.
         else
            lps%lakpoi  = .false.
         endif
         
         ! Column level initialization
         do ci = 1,lps%ncolumns
            
            ! Assign local pointer for simpler referencing
            c => l%c(ci)
            cps => c%cps
            
            ! Initialize restriction for min of soil potential (mm) 
            cps%smpmin = -1.e8   

            ! Soil color
            cps%isoicol = soic2d(i,j)
         
            ! Soil hydraulic and thermal properties
            if (ltype==istdlak .or. ltype==istwet .or. &
                 ltype==istice .or. ltype==isturb ) then 
               do lev = 1,nlevsoi
                  cps%bsw(lev) = spval
                  cps%watsat(lev) = spval
                  cps%hksat(lev) = spval
                  cps%sucsat(lev) = spval
                  cps%tkmg(lev) = spval
                  cps%tksatu(lev) = spval
                  cps%tkdry(lev) = spval
                  cps%csol(lev) = spval
               end do
            else 
               do lev = 1,nlevsoi
                  clay = clay3d(i,j,lev)
                  sand = sand3d(i,j,lev)
                  watsat = 0.489 - 0.00126*sand
                  bd = (1.-watsat)*2.7e3
                  zi = zisoi(lev)  
                  xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s
                  tkm = (8.80*sand+2.92*clay)/(sand+clay)          ! W/(m K)
                  
                  cps%bsw(lev) = 2.91 + 0.159*clay
                  cps%watsat(lev) = watsat
                  cps%hksat(lev) = xksat * exp(-zi/hkdepth)
                  cps%sucsat(lev) = 10. * ( 10.**(1.88-0.0131*sand) )
                  cps%tkmg(lev) = tkm ** (1.- watsat)           
                  cps%tksatu(lev) = cps%tkmg(lev)*0.57**watsat
                  cps%tkdry(lev) = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)  
                  cps%csol(lev) = (2.128*sand+2.385*clay) / (sand+clay)*1.e6  ! J/(m3 K)
               end do
            endif
        
            ! Define lake or non-lake levels layers
            if (ltype == istdlak) then
               cps%z(1:nlevlak) = zlak(1:nlevlak)
               cps%dz(1:nlevlak) = dzlak(1:nlevlak)
            else
               cps%z(1:nlevsoi) = zsoi(1:nlevsoi)
               cps%zi(0:nlevsoi) = zisoi(0:nlevsoi)
               cps%dz(1:nlevsoi) = dzsoi(1:nlevsoi)
            end if

            ! Initialize terms needed for dust model
            clay = clay3d(i,j,1)
            cps%gwc_thr = 0.17 + 0.14*clay*0.01    
            cps%mss_frc_cly_vld = min(clay*0.01_r8, 0.20_r8) 
            
            ! PFT level initialization
            do pi = 1,cps%npfts 
               
               ! Assign local pointer for simpler referencing
               p => c%p(pi)
               pps => p%pps
               
               ! Vegetation type
               ivt = pps%itype
               
               ! Initialize maximum allowed dew
               pps%dewmx  = 0.1     
               
               ! Set up pointer to appropriate array element of ecophysiological constants 
               p%pepc => pftcon(ivt)
               p%pdgvepc => dgv_pftcon(ivt)
               
               ! Initialize root fraction (computing from surface, d is depth in meter):
               ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
               ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with 
               ! beta & d_obs given in Zeng et al. (1998).
               if (ivt /= noveg) then
                  do lev = 1, nlevsoi-1
                     pps%rootfr(lev) = .5*( exp(-roota_par(ivt)*cps%zi(lev-1))  &
                                          + exp(-rootb_par(ivt)*cps%zi(lev-1))  &
                                          - exp(-roota_par(ivt)*cps%zi(lev  ))  &
                                          - exp(-rootb_par(ivt)*cps%zi(lev  )) )
                 end do
                 pps%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*cps%zi(nlevsoi-1))  &
                                          + exp(-rootb_par(ivt)*cps%zi(nlevsoi-1)) )
              else
                 pps%rootfr(1:nlevsoi) = spval
              endif

           end do ! end pft level initialization
        end do ! end column level initialization
     end do ! end landunit level initialization
   end do ! end grid level initialization  

   ! --------------------------------------------------------------------
   ! Initialize other misc derived type components - note that for a 
   ! restart run, dtime is set in routine restrd()
   ! --------------------------------------------------------------------

   if (nsrest == 0) dtime = get_step_size()
     
end subroutine iniTimeConst 
