!----------------------------------------------------------------------- 
! 
! Purpose: ice data handling module. Mimic required data in data-ice-model
!	as needed from the CCSM sea-ice model CSIM.
!
! Author: Erik Kluzek
! 
!-----------------------------------------------------------------------
      module ice_dh

      implicit none
!
! Add flag if prognostic snow on ice, and flag to reset csim iceprops.
!
      logical :: prognostic_icesnow, reset_csim_iceprops

      contains

      logical function icemodel_is( name )
!
! Input arguments
!
      character(len=*) :: name
      
      if (name == 'ccmdata' .or. name == 'data' .or. name == 'CCMDATA') then
         icemodel_is = .true.
      else
         icemodel_is = .false.
      end if
      
      return

      end function icemodel_is

      end module ice_dh
