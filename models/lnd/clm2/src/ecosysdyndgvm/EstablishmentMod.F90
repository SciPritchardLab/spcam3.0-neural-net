#include <misc.h>
#include <preproc.h>

module EstablishmentMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: EstablishmentMod
! 
! !DESCRIPTION: 
! Calculates establishment of new pfts
! Called once per year
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Establishment
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Establishment
!
! !INTERFACE:
  subroutine Establishment ()
!
! !DESCRIPTION: 
! Calculates establishment of new pfts
! Called once per year
!
! !USES:
    use clmtype
    use clm_varpar, only : numpft
    use clm_varcon, only : istsoil
    use pftvarcon, only :  noveg
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_PI, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. establishment)
! 3/4/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit scalars
!
    logical , pointer :: present   !true=> PFT present in patch
    logical , pointer :: tree      !true=> tree is present 
    integer , pointer :: itypveg   !vegetation type for this pft
    integer , pointer :: ityplun   !landunit type
    real(r8), pointer :: wtxy      !pft weight relative to grid cell
    real(r8), pointer :: tmomin20  !20-yr running mean of tmomin
    real(r8), pointer :: agdd20    !20-yr running mean of agdd
    real(r8), pointer :: agddtw    !accumulated growing degree days above twmax
    real(r8), pointer :: nind      !number of individuals (#/m**2)
    real(r8), pointer :: litterag  !above ground litter
    real(r8), pointer :: litterbg  !below ground litter
    real(r8), pointer :: lm_ind    !individual leaf mass
    real(r8), pointer :: sm_ind    !individual sapwood mass
    real(r8), pointer :: hm_ind    !individual heartwood mass
    real(r8), pointer :: rm_ind    !individual root mass
    real(r8), pointer :: fpcgrid   !foliar projective cover on gridcell (fraction)
    real(r8), pointer :: prec365   !365-day running mean of tot. precipitation
    real(r8), pointer :: crownarea !area that each individual tree takes up (m^2)
    real(r8), pointer :: htop      !canopy top (m)
    real(r8), pointer :: sla       !sp. leaf area [m2 leaf g-1 carbon]
    real(r8), pointer :: lai_ind   !LAI per individual
    real(r8), pointer :: crownarea_max !tree maximum crown area [m2]
    real(r8), pointer :: lm_sapl   !leaf mass of sapling
    real(r8), pointer :: sm_sapl   !stem mass of sapling 
    real(r8), pointer :: hm_sapl   !heartwood mass of sapling
    real(r8), pointer :: rm_sapl   !root mass of saping
    real(r8), pointer :: reinickerp!parameter in allometric equation
    real(r8), pointer :: wooddens  !wood density (gC/m3)
    real(r8), pointer :: latosa    !ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b)
    real(r8), pointer :: allom1    !parameter in allometric
    real(r8), pointer :: allom2    !parameter in allometric
    real(r8), pointer :: allom3    !parameter in allometric
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: gi,li,ci,pi   ! grid and subgrid type indices
    integer  :: m                  !index
    integer  :: ngrass             !counter 
    integer  :: npft_estab         !counter
    logical  :: grid_present(0:numpft)  !true=>pft is present in gridcell
    logical  :: grid_survive(0:numpft)  !true=>pft survives in gridcell
    logical  :: grid_estab(0:numpft)    !true=>pft is established in grirdcell
    real(r8) :: fpc_tree_total     !total fractional cover of trees in vegetated portion of gridcell
    real(r8) :: fpc_grass_total    !total fractional cover of grass in vegetated portion of gridcell
    real(r8) :: fpc_total          !old-total fractional vegetated portion of gridcell (without bare ground) 
    real(r8) :: fpc_total_new      !new-total fractional vegetated portion of gridcell (without bare ground) 
    real(r8) :: grid_tmomin20      !20 year running mean of minimum monthly temperature
    real(r8) :: grid_agdd20        !20 year running mean of growing degree days 
    real(r8) :: grid_agddtw        !growing degree base tw
    real(r8) :: fpc_ind            !individual foliage projective cover
    real(r8) :: fpcgridtemp        !temporary  
    real(r8) :: estab_rate         !establishment rate
    real(r8) :: estab_grid         !establishment rate on grid cell
    real(r8) :: bare_max           !maximum bare soil 
    real(r8) :: bare               !fractional cover of bare soil 
    real(r8) :: nind_old           !old number of individuals 
    real(r8) :: sm_ind_temp        !temporary    
    real(r8) :: stemdiam           !stem diameter
    real(r8) :: tcmin              !PFT-specific minimum coldest-month temp.
    real(r8) :: tcmax              !PFT-specific maximum coldest-month temp.
    real(r8) :: gddmin             !PFT-specific minimum GDD
!
! minimum individual density for persistence of PFT (indiv/m2)
    real(r8), parameter :: nind_min = 1.0e-10 
!
! minimum precip. for establishment (mm/s)
    real(r8), parameter :: prec_min_estab = 100./(365.*SHR_CONST_CDAY) 
!
! maximum sapling establishment rate (indiv/m2)
    real(r8), parameter :: estab_max = 0.24   
!
! local pointers to derived subtypes
    type(gridcell_type) , pointer :: g
    type(landunit_type) , pointer :: l
    type(column_type)   , pointer :: c
    type(pft_type)      , pointer :: p
    type(model_pstate_type)   , pointer :: mps
    type(gridcell_pstate_type), pointer :: gps
    type(landunit_pstate_type), pointer :: lps   
    type(column_pstate_type)  , pointer :: cps
    type(pft_pstate_type)     , pointer :: pps 
    type(pft_epc_type)        , pointer :: pepc  
    type(pft_dgvstate_type)   , pointer :: pdgvs 
    type(pft_dgvepc_type)     , pointer :: pdgvepc 
!-----------------------------------------------------------------------


    ! original local variables
    !-----------------------------------------------------------------------

    ! **********************************************************************
    ! My version of LPJ's subr. bioclim

    ! Limits based on 20-year running averages of coldest-month mean
    ! temperature and growing degree days (5 degree base).
    ! For SURVIVAL, coldest month temperature and GDD should be
    ! at least as high as PFT-specific limits.
    ! For REGENERATION, PFT must be able to survive AND coldest month
    ! temperature should be no higher than a PFT-specific limit.

    ! assign local pointers to derived subtypes
    mps => clm%mps

    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)
    do gi = 1,mps%ngridcells
       ! assign local gridcell pointers
       g => clm%g(gi)
       gps => g%gps

       ! initialize gridcell-level metrics
       grid_tmomin20 = 0.
       grid_agdd20 = 0.
       grid_agddtw = 0.

       do m = 0,numpft
          grid_present(m) = .false.
          grid_survive(m) = .false.
          grid_estab(m)   = .false.
       end do

       ngrass = 0
       npft_estab = 0
       fpc_tree_total = 0.0
       fpc_grass_total = 0.0
       fpc_total = 0.0
       fpc_total_new = 0.0

       ! Begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! Assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! Assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pdgvs => p%pdgvs

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                wtxy => pps%wtxy
                tmomin20 => pdgvs%tmomin20
                agdd20 => pdgvs%agdd20
                agddtw => pdgvs%agddtw 
                present => pdgvs%present

                ! Calculate the grid-average bioclimate variables for 
                ! survival and establishment
                ! note: wtxy is the pft weight relative to the parent
                ! gridcell, as opposed to wt, which is relative to the
                ! parent column
                grid_tmomin20 = grid_tmomin20 + wtxy *tmomin20
                grid_agdd20 = grid_agdd20 + wtxy * agdd20
                grid_agddtw = grid_agddtw + wtxy * agddtw

                ! set the presence of pft for this gridcell
                ! Note: modifies the pft-level itypveg if present is false
                ! Refresh the pointers to ecophysiological constant structures 
                ! when the vegetation type changes
                grid_present(itypveg) = present
                if (.not. present) then
                   itypveg = 0
                   p%pdgvepc => dgv_pftcon(itypveg)
                   p%pepc => pftcon(itypveg)
                endif

             enddo ! end pft loop
          enddo ! end column loop
       enddo ! end landunit loop

       ! Must go thru all 16 pfts and decide which can/cannot establish or survive
       ! Determine present, survive, estab
       ! note - if tmomin20 > tcmax then  crops, shrubs and 2nd boreal
       ! summergreen tree cannot exist yet (see EcosystemDynini) because 
       ! they often coexist using up all pft patches and allowing for no bare
       ! ground. Thus they make fpc_grid_total < 1. Solve by allowing one 
       ! one more pft patch per grid cell than the number of pfts.

       do m = 1, numpft
          tcmin  = dgv_pftcon(m)%tcmin + SHR_CONST_TKFRZ !PFT-specific minimum coldest-month temp.
          tcmax  = dgv_pftcon(m)%tcmax + SHR_CONST_TKFRZ !PFT-specific maximum coldest-month temp.
          gddmin = dgv_pftcon(m)%gddmin                  !PFT-specific minimum GDD
          if (grid_tmomin20 >= tcmin) then
             grid_survive(m) = .true.
             if (grid_tmomin20 <= tcmax &
                  .and. grid_agdd20 >= gddmin &
                  .and. nint(grid_agddtw) == 0) then
                grid_estab(m) = .true.
             endif
          endif
       end do

       ! Go back down to the pft level...
       ! Begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! Assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Assign local pointers to derived type scalar members
          ityplun => lps%itype

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                present => pdgvs%present
                nind => pdgvs%nind
                litterag => pdgvs%litterag
                litterbg => pdgvs%litterbg
                lm_ind => pdgvs%lm_ind
                sm_ind => pdgvs%sm_ind
                hm_ind => pdgvs%hm_ind
                rm_ind => pdgvs%rm_ind
                fpcgrid => pdgvs%fpcgrid
                prec365 => pdgvs%prec365
                crownarea => pdgvs%crownarea
                tree => pdgvepc%tree

                ! Case 1 -- pft ceases to exist -kill PFTs not adapted to current climate 
                ! Refresh the pointers to ecophysiological constant structures 
                ! when the vegetation type changes
                if (present .and. &
                     (.not. grid_survive(itypveg) .or. nind<nind_min)) then
                   present = .false.
                   grid_present(itypveg) = .false.

                   ! Add killed biomass to litter
                   if (tree) then
                      litterag = litterag + nind * (lm_ind + sm_ind + hm_ind)
                   else              !if grass
                      litterag = litterag + nind * lm_ind
                   endif
                   litterbg = litterbg + nind * rm_ind

                   fpcgrid = 0.0
                   itypveg = 0    
                   p%pdgvepc => dgv_pftcon(itypveg)
                   p%pepc => pftcon(itypveg)
                   pdgvepc => p%pdgvepc
                   tree => pdgvepc%tree 
                end if

                ! Case 2 -- pft begins to exist - introduce newly "adapted" PFTs
                ! Refresh the pointers to ecophysiological constant structures 
                ! when the vegetation type changes
                if (ityplun == istsoil) then
                   if (.not. present .and. prec365 >= prec_min_estab) then
                      do m = 1, numpft
                         if (itypveg /= m .and. .not. present) then
                            if (.not. grid_present(m) .and. grid_estab(m)) then
                               present = .true.
                               grid_present(m) = .true.
                               itypveg = m
                               p%pdgvepc => dgv_pftcon(itypveg)
                               p%pepc => pftcon(itypveg)
                               pdgvepc => p%pdgvepc
                               tree => pdgvepc%tree 

                               if (tree) then
                                  nind = 0.0
                               else
                                  nind = 1.0 !each grass PFT = 1 "individual"
                               endif

                               lm_ind = 0.0
                               sm_ind = 0.0
                               rm_ind = 0.0
                               hm_ind = 0.0
                               fpcgrid = 0.0

                               if (.not. tree) crownarea = 1.0
                            end if   !conditions suitable for establishment
                         end if   !no pft present and pft 'm' was not here until now
                      end do   !numpft
                   end if   !more conditions for establishment
                end if   !if soil

                ! Case 3 -- some pfts continue to exist (no change) and
                ! some pfts continue to not exist (no change)
                ! Do nothing for this case

                ! Sapling and grass establishment
                ! Calculate total woody FPC and number of woody PFTs present and
                ! able to establish
                if (present) then
                   if (tree) then
                      fpc_tree_total = fpc_tree_total + fpcgrid
                      if (grid_estab(itypveg)) npft_estab = npft_estab + 1
                   else if (.not.tree .and. itypveg > 0) then !grass
                      ngrass = ngrass + 1
                      fpc_grass_total = fpc_grass_total + fpcgrid
                   endif
                   fpc_total = fpc_total + fpcgrid
                endif

                ! These establishment counters at the grid level are
                ! required for the next steps, so close the subgrid loops

             end do ! pft loop
          end do ! column loop
       end do ! landunit loop

       ! Now that ngrass, npft_estab, fpc_tree_total, fpc_grass_total,
       ! and fpc_total are complete for this gridcell, go back down 
       ! to the pft level...
       ! Begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! Assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pepc => p%pepc
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                htop => pps%htop
                sla => pepc%sla
                present => pdgvs%present
                nind => pdgvs%nind
                litterag => pdgvs%litterag
                litterbg => pdgvs%litterbg
                lm_ind => pdgvs%lm_ind
                sm_ind => pdgvs%sm_ind
                hm_ind => pdgvs%hm_ind
                rm_ind => pdgvs%rm_ind
                fpcgrid => pdgvs%fpcgrid
                prec365 => pdgvs%prec365
                crownarea => pdgvs%crownarea
                lai_ind => pdgvs%lai_ind
                tree => pdgvepc%tree
                crownarea_max => pdgvepc%crownarea_max
                lm_sapl => pdgvepc%lm_sapl
                sm_sapl => pdgvepc%sm_sapl
                hm_sapl => pdgvepc%hm_sapl
                rm_sapl => pdgvepc%rm_sapl
                reinickerp => pdgvepc%reinickerp
                wooddens => pdgvepc%wooddens 
                latosa => pdgvepc%latosa
                allom1 => pdgvepc%allom1
                allom2 => pdgvepc%allom2
                allom3 => pdgvepc%allom3

                ! Prohibit establishment under extreme temperature or water stress.
                if (prec365 >= prec_min_estab .and. npft_estab > 0) then

                   ! Calculate establishment rate over available space, per tree PFT
                   ! Maximum establishment rate reduced by shading as tree FPC approaches 1
                   ! Total establishment rate partitioned equally among regenerating woody PFTs
                   estab_rate = estab_max * (1.0-exp(5.0*(fpc_tree_total-1.0))) / &
                        real(npft_estab)

                   ! Calculate grid-level establishment rate per woody PFT
                   ! Space available for woody PFT establishment is proportion of grid cell
                   ! not currently occupied by woody PFTs
                   estab_grid = estab_rate * (1.0-fpc_tree_total)

                else !if unsuitable climate for establishment

                   estab_grid = 0.0

                endif

                if (present .and. tree .and. grid_estab(itypveg)) then

                   ! Add new saplings to current population
                   nind_old = nind
                   nind = nind_old + estab_grid

                   lm_ind = (lm_ind * nind_old + lm_sapl * estab_grid) / nind
                   sm_ind_temp = (sm_ind * nind_old + sm_sapl * estab_grid) / nind
                   hm_ind = (hm_ind * nind_old + hm_sapl * estab_grid) / nind
                   rm_ind = (rm_ind * nind_old + rm_sapl * estab_grid) / nind

                   ! Calculate height, diameter and crown area for new average
                   ! individual such that the basic allometric relationships (A-C below)
                   ! are satisfied.
                   ! (A) (leaf area) = latosa * (sapwood xs area)
                   !        (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
                   ! (B) (leaf mass) = lmtorm * (root mass)
                   ! (C) height = allom2 * (stem diameter)**allom3  (source?)
                   ! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
                   !                                 crownarea_max)
                   ! From (A),
                   !  (1) sap_xsa = lm_ind * sla / latosa
                   !  (2) wooddens = (sm_ind + hm_ind) / stemvolume
                   !  (3) stemvolume = stem_xsa * height
                   ! From (1), (2) & (3),
                   !  (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
                   !  (5) stem_xsa = SHR_CONST_PI * (stemdiam**2) / 4
                   ! From (5),
                   !  (6) stemdiam = ( 4 * stem_xsa / SHR_CONST_PI )**0.5
                   ! From (4) & (6),
                   !  (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
                   !        SHR_CONST_PI )**0.5
                   ! From (C) & (7),
                   !  (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
                   !        ( allom2 * stemdiam**allom3 ) / SHR_CONST_PI )**0.5
                   ! From (8),
                   !  (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / SHR_CONST_PI /
                   !        allom2 )**( 1 / (2 + allom3) )

                   stemdiam = (4.0*(sm_ind_temp+hm_ind)/wooddens/SHR_CONST_PI/allom2)** &
                        (1.0/(2.0+allom3))              !Eqn 9
                   htop = allom2 * stemdiam**allom3  !Eqn C
                   crownarea = min(crownarea_max, allom1*stemdiam**reinickerp)!Eqn D

                   ! Recalculate sapwood mass, transferring excess sapwood to heartwood
                   ! compartment, if necessary to satisfy Eqn A
                   sm_ind = lm_ind * htop * wooddens * sla / latosa
                   hm_ind = hm_ind + (sm_ind_temp-sm_ind)

                   ! Update LAI and FPC
                   if (crownarea > 0.0) then
                      lai_ind = lm_ind * sla / crownarea
                   else
                      lai_ind = 0.0
                   endif

                   fpc_ind = 1.0 - exp(-0.5*lai_ind)
                   fpcgrid = crownarea *nind * fpc_ind

                   fpc_total_new = fpc_total_new + fpcgrid

                endif ! add new saplings block

                ! close the subgrid loops to update fpc_total_new

             enddo ! pft loop
          enddo ! column loop
       enddo ! landunit loop

       ! Go back down to the pft level...
       ! Begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! Assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! Assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pepc => p%pepc
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                sla => pepc%sla
                present => pdgvs%present
                nind => pdgvs%nind
                fpcgrid => pdgvs%fpcgrid
                litterag => pdgvs%litterag
                litterbg => pdgvs%litterbg
                lm_ind => pdgvs%lm_ind
                sm_ind => pdgvs%sm_ind
                hm_ind => pdgvs%hm_ind
                rm_ind => pdgvs%rm_ind
                crownarea => pdgvs%crownarea
                tree => pdgvepc%tree
                lm_sapl => pdgvepc%lm_sapl
                rm_sapl => pdgvepc%rm_sapl

                if (fpc_total_new > 0.95) then
                   if (tree .and. present) then
                      nind_old = nind
                      nind = nind / (fpc_total_new/0.95)
                      fpcgrid = fpcgrid / (fpc_total_new/0.95)
                      litterag = litterag + (nind_old-nind) * (lm_ind+sm_ind+hm_ind)
                      litterbg = litterbg + (nind_old-nind) * rm_ind
                   endif
                   fpc_total = 0.95
                endif

                ! Section for grasses
                if (present .and. .not. tree) then
                   if (grid_estab(itypveg)) then
                      ! Grasses can establish in non-vegetated areas
                      if (ngrass > 0) then
                         bare = (1.0-fpc_total) / real(ngrass)
                      else
                         bare = 0.0
                      endif

                      bare_max = (-2.0 * crownarea *                     &
                           log(max(1.0_r8-bare-fpcgrid, 0.000001_r8)) / sla - lm_ind) /  &
                           lm_sapl

                      bare = max(0.0_r8, min(bare, bare_max))

                      lm_ind = lm_ind + bare * lm_sapl
                      rm_ind = rm_ind + bare * rm_sapl
                   endif

                   if (lm_ind <= 0.0) then
                      present = .false.
                      litterbg = litterbg + rm_ind * nind
                   endif
                endif

             enddo ! pft loop
          enddo ! column loop
       enddo ! landunit loop

       ! Recalculate fpc's and do error check
       fpc_total = 0.0
       fpc_total_new = 0.0

       ! Go back down to the pft level...
       ! Begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! Assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pepc => p%pepc
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                sla => pepc%sla
                present => pdgvs%present
                crownarea => pdgvs%crownarea
                lai_ind => pdgvs%lai_ind
                lm_ind => pdgvs%lm_ind
                fpcgrid => pdgvs%fpcgrid
                nind => pdgvs%nind

                if (present) then
                   if (crownarea > 0.0) then
                      lai_ind = lm_ind * sla / crownarea
                   else
                      lai_ind = 0.0
                   endif
                   fpc_ind = 1.0 - exp(-0.5*lai_ind)
                   fpcgrid = crownarea * nind * fpc_ind
                   fpc_total = fpc_total + fpcgrid
                else
                   fpcgrid = 0.0
                endif

             enddo ! pft loop
          enddo ! column loop
       enddo ! landunit loop

       if (fpc_total < 0.0) then
          write(6,*) 'Error in Establishment: fpc_total is',fpc_total,gi
          call endrun
       end if

       ! Adjustment b/c fpc_total > 1. This can happen, because grasses are allowed
       ! to grow on the bare soil as determined before updating the trees
       ! This may imply that LPJ allows a brief coexistence of some
       ! trees and grasses. LSM cannot allow such coexistence.
       if (fpc_total > 1.0) then
          if (fpc_total > 1.15) then
             write(6,*) 'Error in Establishment: fpc_total is',fpc_total,gi
             call endrun
          end if

          ! Go back down to the pft level...
          ! Begin loop through landunits on this gridcell
          do li = 1, gps%nlandunits
             ! Assign local landunit pointers
             l => g%l(li)
             lps => l%lps

             ! Begin loop through columns on this landunit
             do ci = 1, lps%ncolumns
                ! Assign local column pointers
                c => l%c(ci)
                cps => c%cps

                do pi=1, cps%npfts
                   ! Assign local pointers to pft-level derived subtypes
                   p => c%p(pi)
                   pps => p%pps
                   pdgvs => p%pdgvs

                   ! Assign local pointers to derived type scalar members
                   itypveg => pps%itype
                   fpcgrid => pdgvs%fpcgrid
                   litterag => pdgvs%litterag
                   litterbg => pdgvs%litterbg
                   lm_ind => pdgvs%lm_ind
                   rm_ind => pdgvs%rm_ind
                   nind => pdgvs%nind

                   if (itypveg >= 12 .and. fpcgrid > 0.0) then
                      fpcgridtemp = fpcgrid
                      fpcgrid = max(0.0_r8, fpcgrid-(fpc_total-1.0_r8))
                      if (fpcgrid == 0.0) then
                         present = .false.
                         litterag = litterag + nind * lm_ind
                         litterbg = litterbg + nind * rm_ind
                      else !recalc lai_ind for consistency with fpcgrid
                         lai_ind = -2.0 * log(1.0 - fpcgrid/(crownarea*nind))
                      end if
                      fpc_total = fpc_total - fpcgridtemp + fpcgrid
                   end if

                enddo ! pft loop
             enddo ! column loop
          enddo ! landunit loop
       end if

       ! Next lines avoid fpcgrid=0 with present=.T. and
       ! itypveg=0 with lai>0. Should C pools be reset, to balance carbon?
       ! go back down to the pft level...
       ! begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! Assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! Begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! Assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pps => p%pps
                pdgvs => p%pdgvs

                ! Assign local pointers to derived type scalar members
                itypveg => pps%itype
                fpcgrid => pdgvs%fpcgrid
                present => pdgvs%present
                lai_ind => pdgvs%lai_ind

                if (fpcgrid == 0.) present = .false.
                if (.not. present) lai_ind = 0.0
                if (fpc_total < 1.0) then
                   if (itypveg == noveg) then
                      fpcgrid = 1.0 - fpc_total
                      fpc_total = 1.0
                   end if
                end if
                fpc_total_new = fpc_total_new + fpcgrid

             enddo ! pft loop
          enddo ! column loop
       enddo ! landunit loop

       if (abs(fpc_total_new - 1.0) > 1.0e-6) then
          write(6,*) 'Error in Establishment: fpc_total_new =',fpc_total_new,gi
          call endrun
       end if

    end do ! gridcell loop

    return
  end subroutine Establishment

#endif 

end module EstablishmentMod
