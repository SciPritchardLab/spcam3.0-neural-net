#include <misc.h>
#include <preproc.h>

module DriverInitMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: DriverInitMod
! 
! !DESCRIPTION: 
! Initialization of driver variables needed from previous timestep
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: DriverInit
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DriverInit
!
! !INTERFACE:
  subroutine DriverInit (c)
!
! !DESCRIPTION: 
! Initialization of driver variables needed from previous timestep
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    real(r8), pointer :: endwb         !water mass end of the time step
    real(r8), pointer :: begwb         !water mass begining of the time step
    integer , pointer :: snl           !number of snow layers
    real(r8), pointer :: h2osno        !snow water (mm H2O)
    real(r8), pointer :: h2ocan        !total canopy water (mm H2O)
!
! local pointers to original implicit out scalars
!
    logical , pointer :: do_capsnow    !true => do snow capping
    real(r8), pointer :: h2osno_old    !snow water (mm H2O) at previous time step 
!
! local pointers to original implicit in arrays
!
    real(r8), dimension(:), pointer :: h2osoi_ice   !ice lens (kg/m2)
    real(r8), dimension(:), pointer :: h2osoi_liq   !liquid water (kg/m2)
!
! local pointers to original implicit out arrays
!
    real(r8), dimension(:), pointer :: frac_iceold  !fraction of ice relative to the tot water
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: ci,pi,j,n         !indices
    type(pft_type)          , pointer :: p    !local pointers 
    type(pft_pstate_type)   , pointer :: pps  !local pointers 
    type(column_pstate_type), pointer :: cps  !local pointers 
    type(column_wstate_type), pointer :: cws  !local pointers 
    type(water_balance_type), pointer :: cwbal!local pointers to derived subtypes
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes
    cps => c%cps
    cws => c%cws
    cwbal => c%cwbal
    snl => cps%snl

    ! Assign local pointers to derived type members
    h2osno => cws%h2osno         
    h2osno_old => cws%h2osno_old  
    h2ocan => cws%pws_a%h2ocan
    do_capsnow => cps%do_capsnow  
    frac_iceold => cps%frac_iceold
    h2osoi_ice => cws%h2osoi_ice  
    h2osoi_liq => cws%h2osoi_liq  
    endwb => cwbal%endwb
    begwb => cwbal%begwb

    ! Save snow mass at previous time step
    h2osno_old = h2osno

    ! Initialize fraction of vegetation not covered by snow
    cps%pps_a%frac_veg_nosno = cps%pps_a%frac_veg_nosno_alb

    do pi = 1, cps%npfts
       p => c%p(pi)
       pps => p%pps
       pps%frac_veg_nosno = pps%frac_veg_nosno_alb
    end do
    
    ! Decide whether to cap snow
    if (h2osno > 1000.) then
       do_capsnow = .true.
    else
       do_capsnow = .false.
    end if

    ! Set ice-fraction and water balance for non-lake points
    if (.not. cps%lps%lakpoi) then
       
       ! Initial set of previous time-step variables
       ! ice fraction of snow at previous time step
       n = snl+1
       do j=n,0
          frac_iceold(j) = h2osoi_ice(j)/(h2osoi_liq(j)+h2osoi_ice(j))
       end do

       ! Determine beginning water balance (at previous time step)
       begwb = h2ocan + h2osno
       do j=1, nlevsoi
          begwb = begwb + h2osoi_ice(j) + h2osoi_liq(j)
       end do

    endif

  end subroutine DriverInit

end module DriverInitMod
