#include <misc.h>
#include <preproc.h>

module BiogeochemistryMod

#if (defined BGC)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: BiogeochemistryMod
! 
! !DESCRIPTION: 
! Calculates surface biogeochemical fluxes.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeochemistry
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
! !IROUTINE: Biogeochemistry
!
! !INTERFACE:
  subroutine Biogeochemistry (c)
!
! !DESCRIPTION: 
! Calculates surface biogeochemical fluxes
!
! !USES:
    use clmtype
    use globals
#if (defined DGVM)
    use shr_const_mod, only: SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_varpar, only : nlevsoi
#endif
    use DustMod, only : DustEmission
    use VOCEmissionMod, only : VOCEmission
!
! !ARGUMENTS:
    implicit none
    type (column_type), target, intent(inout) :: c !column derived type
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan and Sam Levis
! 1/31/02, PET: migrated to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: itypveg      !vegetation type for current pft 
    real(r8), pointer :: psnsun       !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha       !shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: t_veg        !vegetation temperature (Kelvin)
#if (defined DGVM)
    real(r8), pointer :: fpcgrid      !foliar projective cover on gridcell (fraction)
    real(r8), pointer :: nind         !number of individuals
    real(r8), pointer :: dphen        !phenology [0 to 1]
    real(r8), pointer :: lm_ind       !individual leaf mass
    real(r8), pointer :: sm_ind       !individual stem mass
    real(r8), pointer :: rm_ind       !individual root mass
    real(r8), pointer :: respcoeff    !respiration coefficient (LPJ)
    real(r8), pointer :: l_cton       !c/n for leaves (LPJ)
    real(r8), pointer :: s_cton       !c/n for stems (LPJ)
    real(r8), pointer :: r_cton       !c/n for roots (LPJ)
#endif
!
! local pointers to implicit inout scalars
!
    real(r8), pointer :: fmicr        !microbial respiration (umol CO2 /m**2 /s)
#if (defined DGVM)
    real(r8), pointer :: bm_inc       !biomass increment
    real(r8), pointer :: afmicr       !annual microbial respiration
#endif
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: frmf         !leaf maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frms         !stem maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frmr         !root maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frm          !total maintenance respiration (umol CO2 /m**2/s)
    real(r8), pointer :: frg          !growth respiration (umol CO2 /m**2 /s)
    real(r8), pointer :: dmi          !total dry matter production (ug /m**2 /s)
    real(r8), pointer :: fco2         !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
#if (defined DGVM)
    real(r8), pointer :: tsoi25       !soil temperature to 0.25 m (Kelvin)
#endif
!
! local pointers to implicit in arrays
!
    real(r8), pointer :: fpsn               !photosynthesis (umol CO2 /m**2 /s)
    real(r8), dimension(:), pointer :: z    !layer thickness (m)
    real(r8), dimension(:), pointer :: dz   !layer depth (m)
    real(r8), dimension(:), pointer :: t_soisno   !soil temperature (Kelvin)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  i,j           !indices
    real(r8) tf            !temperature factor
    real(r8) dmcf          !co2-to-biomass conversion (ug biomass / umol CO2)
#if (defined DGVM)
    real(r8), parameter :: k = 0.0548 / SHR_CONST_CDAY !from [/day] to [/second]
    real(r8) tsoi,dep
#endif
    integer:: pi           !pft index
!
! local pointers to derived subtypes
!
    type(column_pstate_type), pointer :: cps
    type(column_estate_type), pointer :: ces
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_estate_type)   , pointer :: pes
    type(pft_dgvstate_type) , pointer :: pdgvs
    type(pft_dgvepc_type)   , pointer :: pdgvepc
    type(pft_cflux_type)    , pointer :: pcf
!-----------------------------------------------------------------------

    ! Assign local pointers to column-level derived subtypes
    cps => c%cps
    ces => c%ces

    ! Assign local pointers to column-level implicit input arrays
    z           => cps%z
    dz          => cps%dz
    t_soisno    => ces%t_soisno


    ! Begin loop through pfts
    do pi=1, cps%npfts

       ! Assign local pointers to pft-level derived subtypes
       p => c%p(pi)
       pps => p%pps
       pes => p%pes
       pcf => p%pcf
#if (defined DGVM)
       pdgvs => p%pdgvs
       pdgvepc => p%pdgvepc
#endif

       ! assign local pointers to derived type scalar members (pft-level)
       itypveg => pps%itype
       t_veg => pes%t_veg
       psnsun => pcf%psnsun
       psnsha => pcf%psnsha
       fmicr => pcf%fmicr
       fpsn => pcf%fpsn
       frmf => pcf%frmf
       frms => pcf%frms
       frmr => pcf%frmr
       frm => pcf%frm
       frg => pcf%frg
       dmi => pcf%dmi
       fco2 => pcf%fco2
#if (defined DGVM)
       fpcgrid=> pdgvs%fpcgrid
       nind => pdgvs%nind
       dphen => pdgvs%dphen
       lm_ind => pdgvs%lm_ind
       sm_ind => pdgvs%sm_ind
       rm_ind => pdgvs%rm_ind
       respcoeff => pdgvepc%respcoeff
       l_cton => pdgvepc%l_cton
       s_cton => pdgvepc%s_cton
       r_cton => pdgvepc%r_cton
       bm_inc => pdgvs%bm_inc
       afmicr => pdgvs%afmicr
       tsoi25 => pdgvs%tsoi25
#endif

       ! begin calculations for this pft
       dmcf = 28.5

       ! determine vegetation type
       i = itypveg

#if (!defined DGVM)

       frmf = 0.
       frms = 0.
       frmr = 0.

#else

       ! maintenance respiration: LPJ equations w/ units of [gC m-2 gridcell s-1]
       ! converted to LSM units of [umol CO2 m-2 patch s-1]

       ! Soil temperature to a depth of 0.25 m.
       tsoi = 0.
       dep = 0.
       do j = 1, nlevsoi
          if (z(j)+0.5*dz(j) <= 0.25) then
             tsoi = tsoi + t_soisno(j)*dz(j)
             dep = dep + dz(j)
          end if
       end do
       if (dep /= 0.) then
          tsoi25 = tsoi/dep
       else
          tsoi25 = t_soisno(1)
       end if

       if (i > 0 .and. fpcgrid > 0.0) then
          if (t_veg >= SHR_CONST_TKFRZ-40.) then
             tf = exp(308.56 * (1.0/56.02 - 1.0/(t_veg-227.13)))
          else
             tf = 0.0
          end if

          frmf = respcoeff * k * lm_ind * nind / l_cton * tf * dphen * 2.0e6 / dmcf / fpcgrid
          frms = respcoeff * k * sm_ind * nind / s_cton * tf * 2.0e6 / dmcf / fpcgrid

          if (tsoi25 >= SHR_CONST_TKFRZ-40.) then
             tf = exp(308.56 * (1.0/56.02 - 1.0/(tsoi25-227.13)))
          else
             tf = 0.0
          end if

          frmr = respcoeff * k * rm_ind * nind / r_cton * tf * dphen * 2.0e6 / dmcf / fpcgrid
       else
          frmf = 0.0
          frms = 0.0
          frmr = 0.0
       end if
#endif

       frm  = frmf + frms + frmr          

       ! growth respiration and production
#if (defined DGVM)
       frg = 0.25 * max(fpsn - frm, 0.0)      !changed to match LPJ
       dmi = (fpsn - frm - frg) * dmcf
#else
       frg = 0.                     
       dmi = 0.     
#endif

#if (defined DGVM)
       !bm_inc=[gC/m2 patch area] from dmi=[ug dry matter/m2 patch area/s]
       bm_inc = bm_inc + dmi * dtime * 0.5e-6 

       ! microbial respiration

       ! DGVM calculates clm%fmicr in LitterSOM; problem with units in relation
       ! to history grid averages; {fmicr}=[gC/m2 gridcell vegetated area] in
       ! LPJ calculation => sum(fmicr) over a gridcell would give the total for
       ! the gridcell; it would be wrong to convert to [/m2 patch area] because
       ! soil carbon continues to exist even when the plants above it die and
       ! fpcgrid goes to 0; how to reconcile with history calculation which
       ! will take the following values and weight them by the corresponding
       ! weights of their patches? Could chg. soil carbon and plant litter to
       ! gridcell level pools; this will affect fco2, as well; could treat this
       ! issue simultaneously with soil water, which we also want converted to
       ! the gridcell level for plant competition. (slevis)

       afmicr = afmicr + fmicr ![gC/m2 gridcell vegetated area]     
#else
       fmicr = 0.
#endif

       ! net CO2 flux
       fco2 = -fpsn + frm + frg + fmicr

    end do ! end of pft loop

    ! The following two routines (dust and VOCs) could just as well
    ! be called from driver...

    ! dust mobilization routine (C. Zender's modified codes)
    call DustEmission(c)

    ! VOC emission routine (A. Guenther's model)
    call VOCEmission(c)

    return
  end subroutine Biogeochemistry

#endif

end module BiogeochemistryMod
