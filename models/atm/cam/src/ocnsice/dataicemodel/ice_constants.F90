  module ice_constants
!----------------------------------------------------------------------- 
! 
! Purpose: ice_constants module. Mimic required behavior in data-ice-model
!		as for the CCSM sea-ice model CSIM.
!
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

      real (r8), parameter :: Tffresh = 273.15_r8    ! Freezing temp of fresh ice (K)
      real (r8), parameter :: TfrezK = Tffresh - 1.8_r8
     contains
!
! Dummy call to get dataice model working
!
     subroutine init_constants
       return
     end subroutine init_constants
  end module ice_constants
