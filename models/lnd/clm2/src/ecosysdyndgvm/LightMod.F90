#include <misc.h>
#include <preproc.h>

module LightMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: LightMod
! 
! !DESCRIPTION: 
! Calculate light competition
! Update fpc (for establishment routine)
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
  public :: Light
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
! !IROUTINE: Light
!
! !INTERFACE:
  subroutine Light ()
!
! !DESCRIPTION: 
! Calculate light competition
! Update fpc (for establishment routine)
! Called once per year
!
! !USES:
    use clmtype
    use clmpoint, only : gpoint
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine light)
! 3/4/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    logical , pointer:: present   !whether PFT present in patch
    real(r8), pointer:: fpcinc    !foliar projective cover increment (fraction) 
    real(r8), pointer:: sm_ind    !individual stem mass
    real(r8), pointer:: hm_ind    !individual heartwood mass
    real(r8), pointer:: crownarea !area that each individual tree takes up (m^2)
    real(r8), pointer:: sla       !sp. leaf area [m2 leaf g-1 carbon]
    logical , pointer:: tree      !whether this pft is a tree type
!
! local pointers to implicit inout scalars
!
    real(r8), pointer:: fpcgrid   !foliar projective cover on gridcell (fraction)
    real(r8), pointer:: nind      !number of individuals
    real(r8), pointer:: litterag  !above ground litter
    real(r8), pointer:: litterbg  !below ground litter
    real(r8), pointer:: lm_ind    !individual leaf mass
    real(r8), pointer:: rm_ind    !individual root mass
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: gi,li,ci,pi       !indices
    integer  :: begg,endg         !gridcell 1d indices
    real(r8) :: fpc_tree_total
    real(r8) :: fpc_inc_tree
    real(r8) :: fpc_grass_total
    integer  :: ntree
    real(r8) :: excess
    real(r8) :: nind_kill
    real(r8) :: lm_old
    real(r8) :: lm_kill
    real(r8) :: rm_kill
    real(r8) :: lai_ind
    real(r8) :: fpc_ind
    real(r8), parameter :: fpc_tree_max = 0.95  !maximum total tree FPC

! local pointers to derived subtypes
    type(gridcell_type), pointer :: g 
    type(landunit_type), pointer :: l 
    type(column_type)  , pointer :: c 
    type(pft_type)     , pointer :: p 
    type(gridcell_pstate_type), pointer :: gps   
    type(landunit_pstate_type), pointer :: lps   
    type(column_pstate_type)  , pointer :: cps   
    type(pft_epc_type)        , pointer :: pepc  
    type(pft_dgvstate_type)   , pointer :: pdgvs 
    type(pft_dgvepc_type)     , pointer :: pdgvepc 
!-----------------------------------------------------------------------

    ! Set 1d indices

    begg = grid1d%beg
    endg = grid1d%end

    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)
    
    do gi = begg,endg
       ! set up pointers to gridcell 
       g => gpoint(gi)%g
       gps => g%gps

       ! initialize gridcell-level metrics
       fpc_tree_total = 0.
       fpc_inc_tree = 0.
       fpc_grass_total = 0.
       ntree = 0

       ! begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! assign local pointers to derived type scalar members
                present => pdgvs%present
                fpcgrid => pdgvs%fpcgrid
                fpcinc => pdgvs%fpcinc
                tree => pdgvepc%tree

                if (present) then
                   if (tree) then
                      ntree = ntree + 1
                      fpc_tree_total = fpc_tree_total + fpcgrid
                      fpc_inc_tree = fpc_inc_tree + fpcinc
                   else ! if grass
                      fpc_grass_total = fpc_grass_total + fpcgrid
                   endif
                endif
             enddo ! end pft loop

          enddo ! end column loop

       enddo ! end landunit loop

       ! the gridcell level metrics are now in place, continue...
       ! begin loop through landunits on this gridcell
       do li = 1, gps%nlandunits
          ! assign local landunit pointers
          l => g%l(li)
          lps => l%lps

          ! begin loop through columns on this landunit
          do ci = 1, lps%ncolumns
             ! assign local column pointers
             c => l%c(ci)
             cps => c%cps

             do pi=1, cps%npfts
                ! Assign local pointers to pft-level derived subtypes
                p => c%p(pi)
                pepc => p%pepc
                pdgvs => p%pdgvs
                pdgvepc => p%pdgvepc

                ! assign local pointers to derived type scalar members
                present => pdgvs%present
                nind => pdgvs%nind
                fpcgrid => pdgvs%fpcgrid
                fpcinc => pdgvs%fpcinc
                litterag => pdgvs%litterag
                litterbg => pdgvs%litterbg
                lm_ind => pdgvs%lm_ind
                sm_ind => pdgvs%sm_ind
                hm_ind => pdgvs%hm_ind
                rm_ind => pdgvs%rm_ind
                crownarea => pdgvs%crownarea
                tree => pdgvepc%tree
                sla => pepc%sla

                ! LIGHT COMPETITION
                if (present) then

                   if (tree) then
                      if (fpc_tree_total > fpc_tree_max) then    ! case (1)

                         if (fpc_inc_tree > 0.0) then
                            excess = (fpc_tree_total - fpc_tree_max) * &
                                 fpcinc / fpc_inc_tree
                         else
                            excess = (fpc_tree_total - fpc_tree_max) / &
                                 real(ntree)
                         endif

                         ! Reduce individual density (and thereby gridcell-level biomass)
                         ! so that total tree FPC reduced to 'fpc_tree_max'

                         nind_kill = nind * excess / fpcgrid
                         nind = nind - nind_kill

                         ! Transfer lost biomass to litter

                         litterag = litterag + nind_kill * &
                              (lm_ind + sm_ind + hm_ind)
                         litterbg = litterbg + nind_kill * rm_ind

                      endif

                   else !if grass

                      if (fpc_grass_total > &
                           (1.0-min(fpc_tree_total, fpc_tree_max))) then

                         ! grass competes with itself if total fpc exceeds 1
                         ! NEEDS COMMENTS !!!!!!????!!!!

                         excess = (min(fpc_tree_total, fpc_tree_max) + &
                              fpc_grass_total - 1.0) * fpcgrid / fpc_grass_total
                         lm_old = lm_ind
                         lm_ind = -2.0 * log(1.0-(fpcgrid - excess)) / sla
                         lm_kill = lm_old - lm_ind
                         rm_kill = rm_ind * lm_kill/lm_old
                         rm_ind = rm_ind - rm_kill

                         ! Transfer lost biomass to litter

                         litterag = litterag + lm_kill
                         litterbg = litterbg + rm_kill
                      endif

                   endif

                   ! update fpc (for establishment routine)
                   ! slevis: lai_ind is local here

                   if (crownarea > 0.0) then
                      lai_ind = lm_ind * sla / crownarea
                   else
                      lai_ind = 0.0
                   endif

                   fpc_ind = 1.0 - exp(-0.5*lai_ind)
                   fpcgrid = crownarea * nind * fpc_ind

                endif ! present

             enddo ! pft loop

          enddo ! column loop

       enddo ! landunit loop

    enddo ! gridcell loop

  end subroutine Light

#endif

end module LightMod





