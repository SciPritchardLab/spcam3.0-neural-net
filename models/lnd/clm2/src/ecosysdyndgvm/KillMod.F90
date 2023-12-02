#include <misc.h>
#include <preproc.h>

module KillMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: KillMod
! 
! !DESCRIPTION: 
! Removal of PFTs with negative annual C increment
! NB: PFTs newly beyond their bioclimatic limits are removed in
! subroutine establishment
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
  public :: Kill
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
! !IROUTINE: Kill
!
! !INTERFACE:
  subroutine Kill (bm_inc, litter_ag, litter_bg, &
                   lm_ind, sm_ind   , hm_ind   , &
                   rm_ind, nind     , present  , &
                   tree)
!
! !DESCRIPTION: 
! Removal of PFTs with negative annual C increment
! NB: PFTs newly beyond their bioclimatic limits are removed in
! subroutine establishment
! Called once per year
!
! !ARGUMENTS:
    implicit none
    logical , intent(inout) :: present
    real(r8), intent(inout) :: litter_ag
    real(r8), intent(inout) :: litter_bg
    real(r8), intent(in)    :: bm_inc
    real(r8), intent(in)    :: lm_ind
    real(r8), intent(in)    :: sm_ind
    real(r8), intent(in)    :: hm_ind
    real(r8), intent(in)    :: rm_ind
    real(r8), intent(in)    :: nind
    logical , intent(in)    :: tree
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. kill)
!
!EOP
!-----------------------------------------------------------------------

    if (present) then

       if (bm_inc < 0.0) then !negative C increment this year

          present = .false.   !remove PFT

          ! Transfer killed biomass to litter

          if (tree) then !redundant if block? (slevis)
             litter_ag = litter_ag + (lm_ind + sm_ind + hm_ind) * nind
          else !if grass
             litter_ag = litter_ag + lm_ind * nind
          endif

          litter_bg = litter_bg + rm_ind * nind

       endif

    endif

    return
  end subroutine Kill

#endif

end module KillMod
