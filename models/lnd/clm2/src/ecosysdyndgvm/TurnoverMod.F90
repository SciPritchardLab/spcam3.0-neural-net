#include <misc.h>
#include <preproc.h>

module TurnoverMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: TurnoverMod
! 
! !DESCRIPTION: 
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
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
  public :: Turnover
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
! !IROUTINE: Turnover
!
! !INTERFACE:
  subroutine Turnover (l_turn   , s_turn   , r_turn, &
                       litter_ag, litter_bg, lm_ind, &
                       sm_ind   , hm_ind   , rm_ind, &
                       nind     , present  , turnover_ind)
!
! !DESCRIPTION: 
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
! Called once per year
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(out)   :: turnover_ind
    real(r8), intent(inout) :: litter_ag
    real(r8), intent(inout) :: litter_bg
    real(r8), intent(inout) :: lm_ind
    real(r8), intent(inout) :: sm_ind
    real(r8), intent(inout) :: hm_ind
    real(r8), intent(inout) :: rm_ind
    real(r8), intent(in)    :: nind
    real(r8), intent(in)    :: l_turn
    real(r8), intent(in)    :: s_turn
    real(r8), intent(in)    :: r_turn
    logical , intent(in)    :: present
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. turnover)
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) :: l_torate
    real(r8) :: s_torate
    real(r8) :: r_torate
    real(r8) :: lm_turn
    real(r8) :: sm_turn
    real(r8) :: rm_turn
!-----------------------------------------------------------------------

    ! Turnover rates are reciprocals of tissue longevity

    l_torate = 1.0 / l_turn
    s_torate = 1.0 / s_turn
    r_torate = 1.0 / r_turn

    if (present) then

       ! Calculate the biomass turnover in this year

       lm_turn = lm_ind * l_torate
       sm_turn = sm_ind * s_torate
       rm_turn = rm_ind * r_torate

       ! Update the pools

       lm_ind = lm_ind - lm_turn
       sm_ind = sm_ind - sm_turn
       rm_ind = rm_ind - rm_turn

       ! Convert sapwood to heartwood

       hm_ind = hm_ind + sm_turn

       ! Transfer to litter pools

       litter_ag = litter_ag + lm_turn * nind
       litter_bg = litter_bg + rm_turn * nind

       ! Record total turnover

       turnover_ind = lm_turn + sm_turn + rm_turn

    endif

    return
  end subroutine Turnover

#endif

end module TurnoverMod

