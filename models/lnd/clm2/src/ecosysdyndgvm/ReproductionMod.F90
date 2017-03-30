#include <misc.h>
#include <preproc.h>

module ReproductionMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: ReproductionMod
! 
! !DESCRIPTION: 
! Deduction of reproduction costs from annual biomass increment
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
  public :: Reproduction
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
! !IROUTINE: Reproduction
!
! !INTERFACE:
  subroutine Reproduction (bm_inc, litter_ag, present)
!
! !DESCRIPTION: 
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: litter_ag
    real(r8), intent(inout) :: bm_inc
    logical , intent(in)    :: present
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. reproduction)
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) :: reprod
    real(r8), parameter :: reprod_cost = 0.1 !proportion of NPP lost to reproduction (Harper 1977)
!-----------------------------------------------------------------------

    if (present) then

       ! Calculate allocation to reproduction
       ! Reproduction costs taken simply as a constant fraction of annual NPP

       reprod = max(bm_inc * reprod_cost, 0.0)

       ! assume the costs go to reproductive structures which will
       ! eventually enter the litter pool

       litter_ag = litter_ag + reprod

       ! Reduce biomass increment by reproductive cost

       bm_inc = bm_inc - reprod

    endif

    return
  end subroutine Reproduction

#endif

end module ReproductionMod

