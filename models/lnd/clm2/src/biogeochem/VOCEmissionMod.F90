#include <misc.h>
#include <preproc.h>

module VOCEmissionMod

#if (defined BGC)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: VOCEmissionMod
! 
! !DESCRIPTION: 
! Volatile organic compound emission
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
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
! !IROUTINE: VOCEmission
!
! !INTERFACE:
  subroutine VOCEmission(c)
!
! !DESCRIPTION: 
! Volatile organic compound emission
! This code simulates volatile organic compound emissions
! following the algorithm presented in Guenther, A., 1999: Modeling
! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
! Reactive Hydrocarbons in the Atmosphere, Ch. 3
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    E = epsilon * gamma * density * delta
! The factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. We neglect this factor at the present time.
! This factor is discussed in Guenther (1997).
! Subroutine written to operate at the patch level.
! IN FINAL IMPLEMENTATION, REMEMBER:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! Output: vocflx(nvoc) !VOC flux [ug C m-2 h-1]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use shr_const_mod, only : SHR_CONST_RGAS
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c		!column derived type
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
! 2/1/02, Peter Thornton: migration to new data structure
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: itypveg      !vegetation type for current pft 
    real(r8), pointer :: t_veg        !vegetation temperature (Kelvin)
    real(r8), pointer :: fsun         !sunlit fraction of canopy
    real(r8), pointer :: elai         !one-sided leaf area index with burying by snow
!
! local pointers to implicit in arrays
!
    real(r8), dimension(:), pointer :: forc_solad  !direct beam radiation (visible only)
    real(r8), dimension(:), pointer :: forc_solai  !diffuse radiation     (visible only)
!
! local pointers to original implicit out arrays
!
    real(r8), dimension(:), pointer :: vocflx      !VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_tot                !VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_1                  !VOC flux(1) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_2                  !VOC flux(2) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_3                  !VOC flux(3) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_4                  !VOC flux(4) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_5                  !VOC flux(5) [ug C m-2 h-1]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: n       !loop index
    real(r8) :: epsilon !emission factor [ug g-1 h-1]
    real(r8) :: gamma   !activity factor (instantaneous light and temp. condns)
    real(r8) :: density !source density factor [g dry wgt foliar mass/m2 ground]
    real(r8) :: cl
    real(r8) :: ct
    real(r8) :: par
    real(r8) :: reciprod
    integer  :: pi        !pft loop index

    ! constants
    real(r8), parameter :: alpha = 0.0027 !empirical coefficient
    real(r8), parameter :: cl1 = 1.066    !empirical coefficient
    real(r8), parameter :: ct1 = 95000.0  !empirical coefficient [J mol-1]
    real(r8), parameter :: ct2 = 230000.0 !empirical coefficient [J mol-1]
    real(r8), parameter :: ct3 = 0.961    !empirical coefficient
    real(r8), parameter :: tm  = 314.0    !empirical coefficient [K]
    real(r8), parameter :: R   = SHR_CONST_RGAS*0.001 !univ. gas constant [J K-1 mol-1]
    real(r8), parameter :: tstd = 303.0   !std temperature [K]
    real(r8), parameter :: bet = 0.09     !beta empirical coefficient [K-1]

    ! Specific leaf areas [m2 leaf g-1 carbon]
    ! These are the values from my version of genesis-ibis / 1000.
    ! With DGVM defined, use LPJ's sla [m2 leaf g-1 carbon]
    ! Divide by 2 in the equation to get dry weight foliar mass from grams carbon
    real(r8) :: hardwire_sla(0:16)
    real(r8) :: sla

    ! local pointers to derived subtypes
    type(atm2lnd_flux_type) , pointer :: a2lf
    type(column_pstate_type), pointer :: cps
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_estate_type)   , pointer :: pes
    type(pft_vflux_type)    , pointer :: pvf
!-----------------------------------------------------------------------

#if (!defined DGVM)
    hardwire_sla( 0) = 0.
    hardwire_sla( 1) = 0.0125 !needleleaf
    hardwire_sla( 2) = 0.0125 !Gordon Bonan suggests NET = 0.0076
    hardwire_sla( 3) = 0.0125 !Gordon Bonan suggests NDT = 0.0200
    hardwire_sla( 4) = 0.0250 !broadleaf
    hardwire_sla( 5) = 0.0250 !Gordon Bonan suggests BET = 0.0178
    hardwire_sla( 6) = 0.0250 !Gordon Bonan suggests BDT = 0.0274
    hardwire_sla( 7) = 0.0250
    hardwire_sla( 8) = 0.0250
    hardwire_sla( 9) = 0.0250
    hardwire_sla(10) = 0.0250
    hardwire_sla(11) = 0.0250
    hardwire_sla(12) = 0.0200 !grass
    hardwire_sla(13) = 0.0200
    hardwire_sla(14) = 0.0200
    hardwire_sla(15) = 0.0200
    hardwire_sla(16) = 0.0200 !numpft = 16
#endif

    ! Assign local pointers to column-level derived subtypes
    a2lf => c%a2lf
    cps => c%cps

    ! Assign local pointers to derived subtypes components (column-level)
    forc_solad => a2lf%forc_solad
    forc_solai => a2lf%forc_solai

    ! Begin loop through pfts
    do pi=1,cps%npfts

       ! Assign local pointers to pft-level derived subtypes
       p => c%p(pi)
       pps => p%pps
       pes => p%pes
       pvf => p%pvf

       ! Assign local pointers to original implicit out arrays
       vocflx => pvf%vocflx
       vocflx_tot => pvf%vocflx_tot
       vocflx_1 => pvf%vocflx_1
       vocflx_2 => pvf%vocflx_2
       vocflx_3 => pvf%vocflx_3
       vocflx_4 => pvf%vocflx_4
       vocflx_5 => pvf%vocflx_5

       ! Assign local pointers to derived subtypes components (pft-level)
       itypveg => pps%itype
       t_veg => pes%t_veg
       fsun => pps%fsun
       elai => pps%elai

       ! Determine specific leaf array 
#if (!defined DGVM)
       sla = hardwire_sla(itypveg)
#else
       sla = p%pepc%sla
#endif

       ! Begin loop through voc species
       do n = 1, nvoc

          ! epsilon: use values from table 3 in Guenther (1997) which originate in
          ! -------  Guenther et al. (1995). In the comments below, I mention the pft
          !          category as described in table 3. Some values were taken directly
          !          from Guenther et al. (1995). Units: [ug g-1 h-1]
          !          Values were updated on 1/2002 (Guenther, personal communication)

          ! isoprenes:

          if (n == 1) then
             if (itypveg == 1) then       !needleleaf evergreen temperate
                epsilon = 2.
             else if (itypveg == 2) then  !needleleaf evergreen boreal
                epsilon = 4.
             else if (itypveg == 3) then  !needleleaf deciduous
                epsilon = 0.
             else if (itypveg == 4) then  !broadleaf evergreen tropical
                epsilon = 24.
             else if (itypveg >= 5 .and. itypveg <= 11) then !other woody veg
                epsilon = 24.
             else if (itypveg >= 12 .and. itypveg <= 17) then !grass & crop
                epsilon = 0.
             else
                epsilon = 0.
             end if

             ! monoterpenes:
          else if (n == 2) then
             if (itypveg >= 1 .and. itypveg <= 2) then !needleleaf evergreen
                epsilon = 2.0
             else if (itypveg == 3) then                   !needleleaf deciduous
                epsilon = 1.6
             else if (itypveg == 4) then                   !broadleaf everg trop
                epsilon = 0.4
             else if (itypveg >= 5 .and. itypveg <= 11) then !other woody veg
                epsilon = 0.8
             else if (itypveg >= 12 .and. itypveg <= 17) then !grass & crop
                epsilon = 0.1
             else
                epsilon = 0.0
             end if

             ! other VOCs (OVOCs)
          else if (n == 3) then

             epsilon = 1.0                 !Guenther (personal communication)

             ! other reactive VOCs (ORVOCs)
          else if (n == 4) then

             epsilon = 1.0                 !Guenther (personal communication)

             ! CO
          else if (n == 5) then

             epsilon = 0.3                 !Guenther (personal communication)

          end if

          ! gamma: Activity factor. Units [dimensionless]

          ! isoprenes:
          if (n == 1) then

             ! scale total incident par by fraction of sunlit leaves (added on 1/2002)
             ! multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)
             ! got this value from subr. Stomata

             reciprod = 1. / (R * t_veg * tstd)
             ct = exp(ct1 * (t_veg - tstd) * reciprod) / &
                  (ct3 + exp(ct2 * (t_veg - tm) * reciprod))

             par = (forc_solad(1) + fsun * forc_solai(1)) * 4.6 
             cl = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)
             gamma = cl * ct !gamma = 1 under std temp & light condns

             par = ((1. - fsun) * forc_solai(1)) * 4.6
             cl = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)
             gamma = gamma + cl * ct !gamma(sun) + gamma(sha)

             ! monoterpenes, OVOCs, and ORVOCs (Guenther, 1999 and 1995):
          else

             gamma = exp(bet * (t_veg - tstd)) 

          end if

          ! density: Source density factor [g dry weight foliar mass m-2 ground]
          if (itypveg > 0) then
             density = elai / (sla * 0.5)
          else
             density = 0.
          end if

          ! calculate the voc flux
          vocflx(n) = epsilon * gamma * density

       end do ! voc species loop

       ! calculate total voc flux and individual components for history output

       vocflx_tot = 0._r8
       do n = 1, nvoc
          vocflx_tot = vocflx_tot + vocflx(n)
       end do
       vocflx_1 = vocflx(1)
       vocflx_2 = vocflx(2)
       vocflx_3 = vocflx(3)
       vocflx_4 = vocflx(4)
       vocflx_5 = vocflx(5)

    end do ! pft loop

    return
  end subroutine VOCEmission

#endif

end module VOCEmissionMod
