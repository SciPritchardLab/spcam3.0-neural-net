#include <misc.h>
#include <preproc.h>

module MortalityMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: MortalityMod
! 
! !DESCRIPTION: 
! Tree background and stress mortality
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
  public :: Mortality
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
! !IROUTINE: Mortality
!
! !INTERFACE:
  subroutine Mortality (bm_inc   , nind    , turnover_ind, &
                        lm_ind   , sm_ind  , hm_ind      , &
                        rm_ind   , sla     , litter_ag   , &
                        litter_bg, present , tree        , &
                        agddtw   , mort_max)
!
! !DESCRIPTION: 
!
! !ARGUMENTS:
    implicit none
    logical , intent(in)    :: tree
    real(r8), intent(in)    :: bm_inc
    real(r8), intent(in)    :: lm_ind
    real(r8), intent(in)    :: sm_ind
    real(r8), intent(in)    :: hm_ind
    real(r8), intent(in)    :: rm_ind
    real(r8), intent(in)    :: turnover_ind
    real(r8), intent(in)    :: sla
    real(r8), intent(in)    :: agddtw
    real(r8), intent(in)    :: mort_max !asymptotic maximum mortality rate (/yr)
    logical , intent(inout) :: present
    real(r8), intent(inout) :: nind
    real(r8), intent(inout) :: litter_ag
    real(r8), intent(inout) :: litter_bg
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. mortality)
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), parameter :: k_mort = 0.3 !coefficient of growth efficiency in mortality equation
    real(r8), parameter :: ramp_agddtw = 300.0
    real(r8) :: greffic
    real(r8) :: bm_delta  !net individual living biomass increment 
                          !(incorporating loss through leaf, root and sapwood turnover) (gC)
    real(r8) :: mort      !tree mortality rate
    real(r8) :: nind_kill !reduction in individual density due to mortality (indiv/m2)
    real(r8) :: heatstress
!-----------------------------------------------------------------------

    if (present .and. tree) then

       heatstress = min(1.0, agddtw / ramp_agddtw)

       ! Calculate net individual living biomass increment
       bm_delta = max(0.0, bm_inc / nind - turnover_ind)

       ! Calculate growth efficiency (net biomass increment per unit leaf area)
       greffic = bm_delta / lm_ind / sla

       ! Mortality rate inversely related to growth efficiency (Prentice et al 1993)
       mort = mort_max / (1.0 + k_mort * greffic)

       ! Reduce individual density (=> gridcell-level biomass) by mortality rate
       mort = min(1.0, mort + heatstress)
       nind_kill = nind * mort
       nind = nind - nind_kill

       ! Transfer lost biomass to litter
       litter_ag = litter_ag + nind_kill * (lm_ind + sm_ind + hm_ind)
       litter_bg = litter_bg + nind_kill * rm_ind

    endif

    if (nind == 0.0) present =.false.

    return
  end subroutine Mortality

#endif

end module MortalityMod
