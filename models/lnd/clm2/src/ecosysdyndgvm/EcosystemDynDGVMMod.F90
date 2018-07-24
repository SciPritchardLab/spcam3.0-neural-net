#include <misc.h>
#include <preproc.h>

module EcosystemDynDGVMMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: EcosystemDynDGVMMod
! 
! !DESCRIPTION: 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: EcosystemDynini ! LPJ and other DGVM related initializations
  public  :: EcosystemDyn    ! Ecosystem dynamics: phenology, vegetation
!
! !PUBLIC MEMBER FUNCTIONS:
  private :: Phenology     ! Compute summer and drought phenology
  private :: FireSeason    ! Calculate length of fire season in a year
  private :: LitterSOM     ! Calculate litter and soil decomposition
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
! !IROUTINE: EcosystemDynini
!
! !INTERFACE:
  subroutine EcosystemDynini ()
!
! !DESCRIPTION: 
! LPJ and other DGVM related initializations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use nanMod
    use clmtype
    use clmpoint, only : ppoint
    use clm_varpar, only: numpft, npftpar, maxpatch_pft
    use pftvarcon, only: pftpar , tree   , summergreen, raingreen  , sla     , &
         lm_sapl, sm_sapl, hm_sapl    , rm_sapl    , latosa  , &
         allom1 , allom2 , allom3     , reinickerp , wooddens, &
         noveg
    use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from LPJ initialization subroutines)
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: pi,n                   !indices
    integer :: begp,endp              !beginning and ending pft indices
    type(gridcell_type), pointer :: g !local pointer to derived subtype
    type(landunit_type), pointer :: l !local pointer to derived subtype
    type(column_type)  , pointer :: c !local pointer to derived subtype
    type(pft_type)     , pointer :: p !local pointer to derived subtype
    type(model_pstate_type)   , pointer :: mps   !local pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps   !local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps   !local pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps   !local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps   !local pointer to derived subtype
    type(pft_dgvstate_type)   , pointer :: pdgvs !local pointer to derived subtype
    integer  :: pft
    real(r8) :: x
    real(r8) :: stemdiam
    real(r8) :: height_sapl
    real(r8) :: table(0:numpft,1:npftpar)
!-----------------------------------------------------------------------

    ! PFT PARAMETERS (follows LPJ subroutine pftparameters)

    !  1  fraction of roots in upper soil layer
    !  2  plants with C4 (1) or C3 (0) photosynthetic pathway
    !  3  water scalar value at which leaves shed by drought deciduous PFT
    !  4  canopy conductance component (gmin, mm/s) not associated with
    !     photosynthesis (Haxeltine & Prentice 1996, Table 4)
    !  5  maintenance respiration coefficient
    !  6  flammability threshold
    !  7  maximum foliar N content (mg/g)
    !     (Haxeltine & Prentice 1996a, Fig 4)
    !  8  fire resistance index
    !  9  leaf turnover period (years)
    ! 10  leaf longevity (years)
    ! 11  sapwood turnover period (sapwood converted to heartwood) (years)
    ! 12  root turnover period (years)
    ! 13  leaf C:N mass ratio
    ! 14  sapwood C:N mass ratio
    ! 15  root C:N mass ratio
    ! 16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
    ! 17  phenology type: evergreen (1), summergreen (2), raingreen (3),
    !     any type (4)
    ! 18  leaf to root ratio under non-water stressed conditions
    ! 19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
    ! 20  tree maximum crown area (m2)
    ! 21  sapling (or grass on initialisation) LAI
    ! 22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    ! 23  boreal pft (1), non-boreal pft (0)
    ! 24  low temperature limit for CO2 uptake
    ! 25  lower range of temperature optimum for photosynthesis
    ! 26  upper range of temperature optimum for photosynthesis
    ! 27  high temperature limit for CO2 unptake

    !BIOCLIMATIC LIMITS

    ! 28 minimum coldest monthly mean temperature
    ! 29 maximum coldest monthly mean temperature
    ! 30 minimum growing degree days (at or above 5 deg C)
    ! 31 upper limit of temperature of the warmest month (twmax)
    ! 32 lower limit of growth efficiency (g/m2)


    ! ---------------------------------------------------------------------
    !      1      2      3      4      5      6      7      8          PFT
    ! ---------------------------------------------------------------------

    data ((table(pft,n),n=1,8),pft=0,numpft) /               &
          inf,   inf,   inf,   inf,   inf,  0.15,   inf,   inf, &      !  0
         0.70,   0.0,  0.00,   0.3,  1.20,  0.15, 100.0,  0.12, &      !  1
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  2 was 1.20
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  3 was 1.20
         0.85,   0.0,  0.00,   0.5,  0.50,  0.15, 100.0,  0.12, &      !  4 was 0.20
         0.70,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.50, &      !  5
         0.70,   0.0,  0.35,   0.5,  0.50,  0.15, 100.0,  0.50, &      !  6 was 0.20
         0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      !  7
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      !  8 was 1.20
         0.85,   0.0,  0.00,   0.5,  1.20,  0.15, 100.0,  0.12, &      !  9 (=4or5)
         0.80,   0.0,  0.00,   0.5,  1.20,  0.15, 120.0,  0.12, &      ! 10
         0.90,   0.0,  0.00,   0.3,  0.60,  0.15, 100.0,  0.12, &      ! 11 was 1.20
         0.90,   0.0,  0.35,   0.5,  0.60,  0.15, 100.0,  1.00, &      ! 12 was 1.20
         0.90,   0.0,  0.35,   0.5,  0.60,  0.15, 100.0,  1.00, &      ! 13 was 1.20
         0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! 14
         0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00, &      ! 15 (.not.
         0.90,   1.0,  0.35,   0.5,  1.20,  0.15, 100.0,  1.00/        ! 16 present)


    ! ---------------------------------------------------------------------
    !      9     10     11     12     13     14     15     16     17   PFT
    ! ---------------------------------------------------------------------

    data ((table(pft,n),n=9,17),pft=0,numpft) /                     &
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 0
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! 1
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, & ! 2
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! 3
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 4
         1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 5
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, & ! 6
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & ! 7
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & ! 8
         2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, & ! 9
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, & !10
         1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, & !11
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !12
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !13
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !14
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0, & !15
         1.0,  1.00,   1.0,   2.0,  29.0, 330.0,  29.0,   3.0,   4.0/   !16


    ! ------------------------------------------------------
    !       18      19     20      21    22     23     PFT
    ! ------------------------------------------------------

    data ((table(pft,n),n=18,23),pft=0,numpft) /  &
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 0
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 1
         1.0, 1000.0,  15.0,  1.500,  1.2,   1.0, &  ! 2
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! 3
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 4
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 5
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 6
         1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  ! 7
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  ! 8
         1.0, 1000.0,  15.0,  1.500,  1.2,   0.0, &  ! 9
         1.0,  200.0,  15.0,  1.500,  1.2,   0.0, &  !10
         1.0,  200.0,  15.0,  1.500,  1.2,   1.0, &  !11
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  !12
        0.75,  100.0,   0.0,  0.001,  1.2,   1.0, &  !13
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  !14
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0, &  !15
        0.75,  100.0,   0.0,  0.001,  1.2,   0.0/    !16

    ! -------------------------------------
    !      24     25     26      27    PFT
    ! -------------------------------------
    data ((table(pft,n),n=24,27),pft=0,numpft) / &
          inf,   inf,   inf,    inf, & ! 0
         -4.0,  20.0,  30.0,   42.0, & ! 1
         -4.0,  15.0,  25.0,   38.0, & ! 2
         -4.0,  15.0,  25.0,   38.0, & ! 3
          2.0,  25.0,  30.0,   55.0, & ! 4
         -4.0,  20.0,  30.0,   42.0, & ! 5
          2.0,  25.0,  30.0,   55.0, & ! 6
         -4.0,  20.0,  25.0,   38.0, & ! 7
         -4.0,  15.0,  25.0,   38.0, & ! 8
          2.0,  25.0,  30.0,   55.0, & ! 9
         -4.0,  20.0,  25.0,   38.0, & !10
         -4.0,  15.0,  25.0,   38.0, & !11
         -4.0,  10.0,  30.0,   45.0, & !12
         -4.0,  10.0,  30.0,   45.0, & !13
          6.0,  20.0,  45.0,   55.0, & !14
          6.0,  20.0,  45.0,   55.0, & !15
          6.0,  20.0,  45.0,   55.0/   !16

    ! --------------------------------------------------------
    !      28       29      30       31      32     PFT
    ! --------------------------------------------------------
    data ((table(pft,n),n=28,npftpar),pft=0,numpft) / &
         inf,     inf,    inf,  1000.0,    inf, & !  0
         -2.0,    22.0,  900.0,  1000.0,    0.0, & !  1
        -32.5,    -2.0,  600.0,    23.0,    0.0, & !  2
       9999.9,    -2.0,  350.0,    23.0,    0.0, & !  3 (was -1000.0)
         15.5,  1000.0,    0.0,  1000.0,    0.0, & !  4
          3.0,    18.8, 1200.0,  1000.0,    0.0, & !  5
         15.5,  1000.0,    0.0,  1000.0,    0.0, & !  6
        -17.0,    15.5, 1200.0,  1000.0,    0.0, & !  7
      -1000.0,    -2.0,  350.0,    23.0,    0.0, & !  8
       9999.9,  1000.0,    0.0,  1000.0,    0.0, & !  9 (was    15.5)
       9999.9,    15.5, 1200.0,  1000.0,    0.0, & ! 10 (was   -17.0)
       9999.9,    -2.0,  350.0,    23.0,    0.0, & ! 11 (was -1000.0)
      -1000.0,   -17.0,    0.0,  1000.0,    0.0, & ! 12 (an LSM type)
        -17.0,    15.5,    0.0,  1000.0,    0.0, & ! 13 (was -1000.0)
         15.5,  1000.0,    0.0,  1000.0,    0.0, & ! 14
       9999.9,  1000.0,    0.0,  1000.0,    0.0, & ! 15 (was    15.5)
       9999.9,  1000.0,    0.0,  1000.0,    0.0/   ! 16 (was    15.5)
    !----------------------------------------------------------------------------

    pftpar(noveg,:) = table(noveg,:)
    sla(noveg)      = inf
    lm_sapl(noveg)  = inf
    rm_sapl(noveg)  = inf
    sm_sapl(noveg)  = inf
    hm_sapl(noveg)  = inf
    tree(noveg)        = .false.
    summergreen(noveg) = .false.
    raingreen(noveg)   = .false.

    do pft = 1,numpft

       ! Transfer parameter values to array pftpar

       do n = 1, npftpar
          pftpar(pft,n) = table(pft,n)
       enddo

       ! Assign leaf and phenology logicals

       if (pftpar(pft,16) <= 2.0) then     !woody vegetation: trees, shrubs
          tree(pft) = .true.
       else                                !non woody vegetation: grasses
          tree(pft) = .false.
       endif

       if     (pftpar(pft,17) == 1.0) then !evergreen
          summergreen(pft) = .false.
          raingreen(pft)   = .false.
       elseif (pftpar(pft,17) == 2.0) then !summergreen
          summergreen(pft) = .true.
          raingreen(pft)   = .false.
       elseif (pftpar(pft,17) == 3.0) then !raingreen
          summergreen(pft) = .false.
          raingreen(pft)   = .true.
       else                                !any of the above
          summergreen(pft) = .true.
          raingreen(pft)   = .true.
       endif

       ! Calculate specific leaf area (SLA) for each PFT from leaf longevity
       ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
       ! Equation based on Reich et al 1997, Fig 1f:

       ! SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

       ! SLA in m2/gC, leaf_longevity in years

       sla(pft) = 2.0e-4 * exp(6.15 - 0.46*log(pftpar(pft,10)*12.0))

       ! Define initial mass structure

       if (tree(pft)) then !woody PFTs

          ! Calculate leafmass for a sapling individual
          !  (1) lai = leafmass * sla / (crown area)
          !  (2) (leaf area) = latosa * (sapwood xs area)
          !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
          !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
          !         (Reinickes theory)
          ! From (1),
          !  (4) leafmass = lai * (crown area) / sla
          ! From (1) & (3),
          !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
          ! From (2),
          !  (6) leafmass = latosa * (sapwood xs area) / sla
          !  (7) (sapwood xs area) = SHR_CONST_PI * (sapwood diameter)**2 / 4
          ! From (6) and (7),
          !  (8) leafmass = latosa * SHR_CONST_PI * (sapwood diameter)**2 / 4 / sla
          ! From (8),
          !  (9) (sapwood diameter) = [ 4 * leafmass * sla / SHR_CONST_PI / latosa ]**0.5
          ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
          ! Define x,
          ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
          !          (sapwood diameter)
          ! From (10) & (11),
          ! (12) (stem diameter) = x * (sapwood diameter)
          ! From (5), (9) & (12),
          ! (13) leafmass = lai * allom1 * x**reinickerp *
          !               (4*leafmass*sla/SHR_CONST_PI/latosa)**(reinickerp*0.5) / sla
          ! From (13),
          ! (14) leafmass = [ lai * allom1 * x**reinickerp *
          !      (4*sla/SHR_CONST_PI/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

          x = pftpar(pft,22)

          lm_sapl(pft) = (pftpar(pft,21) * allom1 * x**reinickerp *          &
               (4.0 * sla(pft) / SHR_CONST_PI / latosa)**(reinickerp * 0.5) / sla(pft))** &
               (2.0/(2.0-reinickerp)) !eqn 14

          ! Calculate sapling stem diameter
          ! From (9) & (12),
          ! (15) (stem diameter) = x * [ 4 * leafmass * sla / SHR_CONST_PI / latosa ]**0.5

          stemdiam = x * (4.0*lm_sapl(pft)*sla(pft)/SHR_CONST_PI/latosa)**0.5 !Eqn 15

          ! Calculate sapling height
          ! (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl = allom2 * stemdiam**allom3 !Eqn 16

          ! Calculate sapling sapwood mass
          ! (17) (sapwood volume) = height * (sapwood xs area)
          ! (18) (sapwood xs area) = leafmass * sla / latosa
          ! From (17) & (18),
          ! (19) (sapwood volume) = height * leafmass * sla / latosa
          ! (20) (sapwood mass) = (wood density) * (sapwood volume)
          ! From (19) & (20),
          ! (21) (sapwood mass) = (wood density) * height * leafmass * sla / latosa

          sm_sapl(pft)=wooddens*height_sapl*lm_sapl(pft)*sla(pft)/latosa !Eqn 21

          ! Calculate sapling heartwood mass
          ! From (11),
          ! (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft) = (x-1.0) * sm_sapl(pft) !Eqn 22

       else !grass PFTs

          lm_sapl(pft) = pftpar(pft,21) / sla(pft)

       endif

       ! Calculate sapling or initial grass rootmass
       ! (23) lmtorm = (leafmass) / (rootmass)
       ! where lmtorm=pftpar(pft,18)

       rm_sapl(pft) = lm_sapl(pft) / pftpar(pft,18) !From Eqn 23

    enddo ! pft loop

    ! ---------------------------------------------------------------
    ! Some of the following comes from LPJ subroutine initgrid
    ! ---------------------------------------------------------------

    ! Set beginning and ending pft indices

    begp = pfts1d%beg
    endp = pfts1d%end

    ! Loop over pfts

    do pi = begp,endp
       p => ppoint(pi)%p
       pps => p%pps
       pdgvs => p%pdgvs
       
       ! used in Phenology
       pdgvs%t10min = 1.0e+36
       pdgvs%lai_ind = 0.0
       
       ! updated in Phenology
       pdgvs%dphen = 0.0
       pdgvs%leafon = 0.0
       pdgvs%leafof = 0.0
       
       ! accumulated in FireSeason (must reset at end of every year)
       pdgvs%firelength = 0.0       
       
       ! used in FireSeason; updated in LitterSOM; updated in annual portion of LPJ
       pdgvs%litterag = 0.0
       pdgvs%litterbg = 0.0
       
       ! updated in LitterSOM
       pdgvs%cpool_fast = 0.0
       pdgvs%cpool_slow = 0.0
       pdgvs%k_fast_ave = 0.0
       pdgvs%k_slow_ave = 0.0
       pdgvs%litter_decom_ave = 0.0
       p%pcf%fmicr = 0.0 !initialize b/c use in Biogeochemistry before LitterSOM
       
       ! used and updated in annual portion of LPJ
       pdgvs%present  = .false.
       pdgvs%nind     = 0.
       pdgvs%lm_ind   = 0.
       pdgvs%sm_ind   = 0.
       pdgvs%hm_ind   = 0.
       pdgvs%rm_ind   = 0.
       pdgvs%tmomin20 = SHR_CONST_TKFRZ - 5. !initialize this way for Phenology code
       pdgvs%agdd20   = 0.
       pdgvs%t_mo_min = 1.0e+36
       
       ! already a variable in LSM but now updated in annual portion of LPJ
       ! note - fpcgrid is not relevant where ist/=istsoil (consistent with subr. surfrd)
       pps%htop = 0.
       pps%tsai = 0.
       pdgvs%fpcgrid = 1.0 / maxpatch_pft !
       !
       ! accumulated in Biogeochemistry and used/reset in annual portion of LPJ
       pdgvs%bm_inc = 0.0
       pdgvs%afmicr = 0.0
       
       ! accumulated in Stomata and used/reset in annual portion of LPJ
       pdgvs%annpsn = 0.0
       pdgvs%annpsnpot = 0.0
       
    end do   !end of pft loop

  end subroutine EcosystemDynini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EcosystemDyn
!
! !INTERFACE:
  subroutine EcosystemDyn (c, doalb, endofyr)
!
! !DESCRIPTION: 
! Ecosystem dynamics: phenology, vegetation
! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
! height (htop)
!
! !USES:
    use clmtype
    use clmpoint
    use globals
    use shr_const_mod, only: SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    type (column_type),target,intent(inout):: c	!column derived type
    logical, intent(in):: doalb   !true = surface albedo calculation time step
    logical, intent(in):: endofyr
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrated to new data structure.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: snowdp   !snow height (m)
    integer , pointer :: itypveg  !vegetation type for this pft
    real(r8), pointer :: lai_ind  !LAI per individual
    real(r8), pointer :: dphen    !phenology [0 to 1]
    real(r8), pointer :: fpcgrid  !foliar projective cover on gridcell (fraction)
!
! local pointers to implicit in/out scalars
!
    real(r8), pointer :: htop     !canopy top (m)
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: tlai     !one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai     !one-sided stem area index, no burying by snow
    real(r8), pointer :: hbot     !canopy bottom (m)
    real(r8), pointer :: elai     !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai     !one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb !frac of vegetation not covered by snow [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8) :: ol                !thickness of canopy layer covered by snow (m)
    real(r8) :: fb                !fraction of canopy layer covered by snow
    integer  :: pi                !pft index
    integer  :: begp,endp         !pft beginning and ending indices

! local pointers to derived subtypes
    type(column_pstate_type), pointer :: cps
    type(pft_type)          , pointer :: p
    type(pft_pstate_type)   , pointer :: pps
    type(pft_dgvstate_type) , pointer :: pdgvs
!-----------------------------------------------------------------------

    ! Assign local pointers to column-level derived subtypes
    cps => c%cps

    ! Assign local pointers to derived type scalar members
    snowdp => cps%snowdp

    ! Begin loop through pfts
    do pi=1, cps%npfts

       ! Assign local pointers to pft-level derived subtypes
       p => c%p(pi)
       pps => p%pps
       pdgvs => p%pdgvs

       ! today's phenology: returns dphen
       ! with small chg in subr, could be called every tstep if preferred

       if (nstep==0 .or. &
           mod(nstep-1, nint(SHR_CONST_CDAY/dtime)) == 0. .or. endofyr) then
          call Phenology(p)
       end if

       ! fire season; returns firelength for use at end of yr in subr. Fire
       ! litter and soil decomposition; returns litterag,bg and fmicr

       if (nstep > 1) then
          call FireSeason(p)
          call LitterSOM(p, ctl%year)
       end if

       ! obtain vegetation structure here (tlai+sai,hbot+top)

       if (doalb .or. endofyr) then

          ! Assign local pointers to derived type scalar members (pft-level)
          tlai => pps%tlai 
          tsai => pps%tsai 
          elai => pps%elai 
          esai => pps%esai 
          htop => pps%htop 
          hbot => pps%hbot 
          frac_veg_nosno_alb => pps%frac_veg_nosno_alb 
          itypveg => pps%itype
          htop    => pps%htop
          lai_ind => pdgvs%lai_ind
          dphen   => pdgvs%dphen
          fpcgrid => pdgvs%fpcgrid

          if (itypveg > 0 .and. fpcgrid > 0.0) then
             if (itypveg <= 11) then          !woody vegetation (following IBIS)
                tlai = lai_ind * dphen
                tsai = 0.250 * lai_ind
                hbot = max(0.0, min(3.00, htop-1.00))
             else                             !grasses (following IBIS)
                tlai = lai_ind * dphen / fpcgrid
                tsai = 0.050 * tlai
                htop = max(0.25, tlai * 0.25)
                hbot = max(0.0, min(0.05, htop-0.20))
             end if
          else
             tlai = 0.0
             tsai = 0.0
             htop = 0.0
             hbot = 0.0
          end if

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05 but non zero, set to 0.05 to prevent numerical
          ! problems associated with very small lai and sai.

          ol = min( max(snowdp-hbot,0._r8), htop-hbot)
          fb = 1. - ol / max(1.e-06,htop-hbot)
          elai = max(tlai*fb,0.0_r8)
          esai = max(tsai*fb,0.0_r8)
          if (elai > 0.0 .and. elai < 0.05) elai = 0.05
          if (esai > 0.0 .and. esai < 0.05) esai = 0.05

          ! Fraction of vegetation free of snow
          if ((elai + esai) >= 0.05) then
             frac_veg_nosno_alb = 1
          else
             frac_veg_nosno_alb = 0
          endif

       end if !end of if-doalb block

    end do ! end of pft loop

  end subroutine EcosystemDyn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Phenology
!
! !INTERFACE:
  subroutine Phenology (p)
!
! !DESCRIPTION: 
! Summer and drought phenology. Called once per day.
!
! !USES:
    use clmtype
    use shr_const_mod, only : SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    type (pft_type),target,intent(inout):: p		!pft derived type
!
! !CALLED FROM:
! subroutine Ecosysdyn in this module
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Jon Foley's IBIS subroutine pheno)
! 2/1/02, Peter Thornton: Migrated to new data structure
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: itypveg   !vegetation type for this pft
    real(r8), pointer :: t10       !10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: agdd0     !accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5     !accumulated growing degree days above -5
    real(r8), pointer :: fnpsn10   !10-day running mean net photosynthesis
    real(r8), pointer :: l_long    !leaf longevity [years]
    real(r8), pointer :: tmomin20  !20 year running mean of monthly minimum	
    logical , pointer :: tree
    logical , pointer :: raingreen
    logical , pointer :: summergreen
!
! local pointers to implicit in/out scalars
!
    real(r8), pointer :: t10min    !annual minimum of 10-day running mean (K)
    real(r8), pointer :: leafon    !leafon days
    real(r8), pointer :: leafof    !leafoff days
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: dphen     !phenology [0 to 1]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8), parameter :: ddfacu = 1.0 / 15.0 !'drop day factor' causes pheno-
    real(r8), parameter :: ddfacl = 1.0 /  5.0 !logy to switch in 15 or 5 days
    real(r8) :: tthreshold

! local pointers to derived subtypes
    type(pft_pstate_type)  ,pointer :: pps
    type(pft_dgvstate_type),pointer :: pdgvs
    type(pft_dgvepc_type)  ,pointer :: pdgvepc
!-----------------------------------------------------------------------

    ! Assign local pointers to pft-level derived subtypes
    pps => p%pps
    pdgvepc => p%pdgvepc
    pdgvs => p%pdgvs

    ! Assign local pointers to derived type scalar members (pft-level)
    itypveg => pps%itype
    t10 => pdgvs%t10
    agdd0 => pdgvs%agdd0
    agdd5 => pdgvs%agdd5
    fnpsn10 => pdgvs%fnpsn10
    l_long => pdgvepc%l_long
    tree => pdgvepc%tree
    raingreen => pdgvepc%raingreen
    summergreen => pdgvepc%summergreen
    t10min => pdgvs%t10min
    leafon => pdgvs%leafon
    leafof => pdgvs%leafof
    dphen => pdgvs%dphen
    tmomin20 => pdgvs%tmomin20 	

    t10min = min (t10min, t10) !reset to 1.0e+36 once per year

    if (tree .and. .not.summergreen .and. .not.raingreen) then

       dphen = 1.0

    else if (tree .and. summergreen) then

       ! ---------------------------------------------------------------------
       ! * * * upper canopy winter phenology * * *
       ! ---------------------------------------------------------------------
       ! temperature threshold for budburst and senescence
       ! temperature threshold is assumed to be SHR_CONST_TKFRZ
       ! or 5 degrees warmer than the coldest monthly temperature
       ! tmomin20 is initialized to tfrz-5.0
       tthreshold = max (SHR_CONST_TKFRZ, tmomin20 + 5.0)

       ! determine if growing degree days are initiated
       ! slevis:*t10 = 10-day running mean air temperature (K)
       ! determine leaf display
       if (t10 < tthreshold) then
          dphen = max (0.0, dphen - ddfacu)
       else
          dphen = min (1.0, max (0.0, agdd0 - 100.0) / 50.0)
       endif

    else if (itypveg > 0 .and. .not.tree) then !NB: grass has no specific phenology

       ! ---------------------------------------------------------------------
       ! * * * lower canopy phenology * * *
       ! ---------------------------------------------------------------------
       ! temperature threshold for budburst and senescence
       ! temperature threshold is assumed to be SHR_CONST_TKFRZ
       tthreshold = SHR_CONST_TKFRZ

       ! determine leaf display
       if (t10 < tthreshold) then             !cold phenology for grasses
          dphen = max (0.0, dphen - ddfacl)   !slevis: made ddfacl=1/5
       else if (fnpsn10 < 0.0) then           !drought phenology for grasses
          dphen = max (0.1, dphen - ddfacl)   !slevis: made ddfacl=1/5
       else
          !        dphen = min (1.0, max (0.0, agdd5 - 150.0) / 50.0)
          dphen = min (1.0, dphen + ddfacl)   !slevis: try this line instead
       endif

    end if

    if (tree .and. raingreen) then

       ! ---------------------------------------------------------------------
       ! * * * upper canopy drought phenology * * *
       ! ---------------------------------------------------------------------

       if (fnpsn10 <  0.0) dphen = max (0.1, dphen - ddfacu)
       if (fnpsn10 >= 0.0) dphen = min (1.0, dphen + ddfacu)

       ! trying out enforced drought phenology

       if (dphen > 0.95) leafon = leafon + 1.0
       if (leafon >= 365.0*l_long) then
          dphen = 0.1
          leafof = leafof + 1.0
          if (leafof >= 365.0*l_long) then
             dphen = 1.0
             leafof = 0.0
             leafon = 1.0
          end if
       end if

    end if

  end subroutine Phenology

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FireSeason
!
! !INTERFACE: 
  subroutine FireSeason (p)
!
! !DESCRIPTION: 
! Calculate length of fire season in a year
! Orig. code was called once per day.
! slevis adapted to call every tstep.
! Orig. code operated on a grid cell basis.
! slevis adapted to operate on a patch basis.
!
! !USES:
    use clmtype
    use globals
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_CDAY, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    type (pft_type),target,intent(inout):: p		!pft derived type
!
! !CALLED FROM:
! subroutine Ecosysdyn in this module
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine fire)
!
! !LOCAL VARIABLES:
!
! local pointers toimplicit in arguments
!
    integer , pointer :: itypveg      !vegetation type for this pft
    real(r8), pointer :: t_ref2m      !2 m height surface air temperature (Kelvin)
    real(r8), pointer :: litterag     !above ground litter
    real(r8), pointer :: wf           !soil water as frac. of whc for top 0.5 m
!
! local pointers toimplicit in/out arguments
!
    real(r8), pointer :: firelength   !fire season in days
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8):: flam
    real(r8):: fire_prob

! local pointers to derived subtypes
    type(pft_pstate_type)  , pointer :: pps
    type(pft_estate_type)  , pointer :: pes
    type(pft_dgvstate_type), pointer :: pdgvs
    type(pft_dgvepc_type)  , pointer :: pdgvepc
!-----------------------------------------------------------------------

    ! Assign local pointers to pft-level derived subtypes
    pps => p%pps
    pes => p%pes
    pdgvs => p%pdgvs
    pdgvepc => p%pdgvepc

    ! Assign local pointers to derived type scalar members (pft-level)
    itypveg => pps%itype
    t_ref2m => pes%t_ref2m
    litterag => pdgvs%litterag
    wf => pps%cps%wf
    firelength => pdgvs%firelength

    ! Copy local variables
    flam = pdgvepc%flam

    ! Calculate the length of the fire season (in days)
    ! Calculate today's fire probability, fire_prob
    ! Assume fire is only possible when temperature is above SHR_CONST_TKFRZ
    ! slevis: *wf is top 0.5 m soil water as a fraction of the whc
    !         *divide fire_prob (days) by tsteps/day to get fire_prob (tsteps)
    !         *else need daily avg t_ref2m and wf to calc. fire_prob

    if (t_ref2m > SHR_CONST_TKFRZ .and. litterag > 0.0) then
       fire_prob=EXP((-SHR_CONST_PI/4.0) * (max(0.0,wf)/flam)**2) * dtime / SHR_CONST_CDAY
    else
       fire_prob=0.0
    endif

    firelength = firelength + fire_prob !reset 1/yr in subroutine lpj

  end subroutine FireSeason

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: LitterSOM
!
! !INTERFACE:
  subroutine LitterSOM (p, kyr)
!
! !DESCRIPTION: 
! Litter and soil decomposition
! Incorporates analytical solution for soil pool sizes
! once litter inputs are (assumed to be) at equilibrium,
! reducing spin-up time for carbon fluxes due to soil respiration.
!
! !USES:
    use clmtype
    use globals
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    type (pft_type), target, intent(inout) :: p !pft derived type
    integer, intent(in) :: kyr                  !year (0, ...) for nstep+1
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine littersom)
! 2/1/02, Peter Thornton: Migrate to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to  implicit in arguments
!
    real(r8), pointer :: wf           !soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: tsoi25       !soil temperature to 0.25 m (Kelvin)
!
! local pointers to  implicit in/out arguments
!
    real(r8), pointer :: litterag     !above ground litter
    real(r8), pointer :: litterbg     !below ground litter
    real(r8), pointer :: cpool_fast   !fast carbon pool 
    real(r8), pointer :: cpool_slow   !slow carbon pool
    real(r8), pointer :: k_fast_ave   !decomposition rate
    real(r8), pointer :: k_slow_ave   !decomposition rate
    real(r8), pointer :: litter_decom_ave   !decomposition rate
!
! local pointers to  implicit out arguments
!
    real(r8), pointer :: fmicr        !microbial respiration (umol CO2 /m**2 /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer , parameter :: soil_equil_year = 400     !number of years until pool sizes for soil decomposition solved analytically
    real(r8), parameter :: k_litter10 = 0.5          !litter decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_fast10 = 0.03      !fast pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_slow10 = 0.001     !slow pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: atmfrac = 0.7             !fraction of litter decomp. going directly into the atmosphere
    real(r8), parameter :: soilfrac = 1.0 - atmfrac  !fraction of litter decomp. going to soil C pools
    real(r8), parameter :: fastfrac = 0.985          !fraction of litter entering fast soil decomposition pool
    real(r8), parameter :: slowfrac = 1.0 - fastfrac !fraction of litter entering slow soil decomposition pool
!
    real(r8) :: temp_resp          !temperature response of decomposition
    real(r8) :: moist_resp         !moisture response of decomposition
    real(r8) :: k_litter           !litter decomposition rate (/tstep)
    real(r8) :: k_fast             !fast pool decomposition rate (/tstep)
    real(r8) :: k_slow             !slow pool decomposition rate (/tstep)
    real(r8) :: litter_decom       !litter decomposition
    real(r8) :: litter_decom_ag    !above-ground litter decomposition
    real(r8) :: litter_decom_bg    !below-ground litter decomposition
    real(r8) :: cflux_litter_soil  !litter decomposition flux to soil
    real(r8) :: cflux_litter_atmos !litter decomposition flux to atmosphere
    real(r8) :: cflux_fast_atmos   !soil fast pool decomposition flux to atmos.
    real(r8) :: cflux_slow_atmos   !soil slow pool decomposition flux to atmos.
!
! local pointers to derived subtypes
    type(pft_pstate_type)  , pointer :: pps   
    type(pft_dgvstate_type), pointer :: pdgvs
    type(pft_cflux_type)   , pointer :: pcf
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes
    pps   => p%pps
    pdgvs => p%pdgvs
    pcf   => p%pcf

    ! Assign local pointers to derived type scalar members
    wf => pps%cps%wf
    tsoi25 => pdgvs%tsoi25
    litterag => pdgvs%litterag
    litterbg => pdgvs%litterbg
    cpool_fast => pdgvs%cpool_fast
    cpool_slow => pdgvs%cpool_slow
    k_fast_ave => pdgvs%k_fast_ave
    k_slow_ave => pdgvs%k_slow_ave
    litter_decom_ave => pdgvs%litter_decom_ave
    fmicr => pcf%fmicr

    ! Temperature response function is a modified Q10 relationship
    ! (Lloyd & Taylor 1994)
    ! slevis: Original code used monthly avg soil temp (K); I use tstep value

    if (tsoi25 <= SHR_CONST_TKFRZ - 40.0) then !avoid division by zero
       temp_resp=0.0
    else                            !Lloyd & Taylor 1994
       temp_resp=exp(308.56*((1.0/56.02)-(1.0/(tsoi25-227.13))))
    endif

    ! Moisture response based on soil layer 1 moisture content (Foley 1995)
    ! slevis: Orig. code used monthly soil water in upper 0.5 m (fraction of whc)
    !         I use the tstep value
    moist_resp = 0.25 + 0.75 * wf

    ! Original divided by 12 to get monthly decomposition rates (k, /month)
    ! as a function of temperature and moisture
    ! slevis: make rates /tstep by dividing by the number of tsteps per year
    k_litter = k_litter10    * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365.)
    k_fast   = k_soil_fast10 * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365.)
    k_slow   = k_soil_slow10 * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365.)

    ! Calculate monthly litter decomposition using equation
    !   (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
    ! from (1),
    !   (2) c = c0*exp(-kt) where c0=initial pool size
    ! from (2), decomposition in any month given by
    !   (3) delta_c = c0 - c0*exp(-k)
    ! from (3)
    !   (4) delta_c = c0*(1.0-exp(-k))
    litter_decom_ag = litterag * (1.0-exp(-k_litter))  !eqn 4
    litter_decom_bg = litterbg * (1.0-exp(-k_litter))
    litter_decom    = litter_decom_ag + litter_decom_bg

    ! Update the litter pools
    litterag = litterag - litter_decom_ag
    litterbg = litterbg - litter_decom_bg

    ! Calculate carbon flux to atmosphere and soil
    cflux_litter_atmos = atmfrac  * litter_decom
    cflux_litter_soil  = soilfrac * litter_decom

    ! Further subdivide soil fraction between fast and slow soil pools
    cpool_fast = cpool_fast + fastfrac * cflux_litter_soil
    cpool_slow = cpool_slow + slowfrac * cflux_litter_soil

    ! Calculate monthly soil decomposition to the atmosphere
    cflux_fast_atmos = cpool_fast * (1.0-exp(-k_fast))  !eqn 4
    cflux_slow_atmos = cpool_slow * (1.0-exp(-k_slow))  !eqn 4

    ! Update the soil pools
    cpool_fast = cpool_fast - cflux_fast_atmos
    cpool_slow = cpool_slow - cflux_slow_atmos

    ! Calculate heterotrophic respiration (in LSM referred to as microbial)
    fmicr = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos

    ! Empty soil pools below a minimum threshold
    if (cpool_fast < 1.0e-5) cpool_fast = 0.0
    if (cpool_slow < 1.0e-5) cpool_slow = 0.0

    if (kyr <= soil_equil_year) then

       ! Update running average respiration rates and litter input
       ! slevis: had to multiply the denominator to chg units from years to tsteps
       k_fast_ave       = k_fast_ave       + k_fast / &
            (real(soil_equil_year) * 365. * SHR_CONST_CDAY / dtime)
       k_slow_ave       = k_slow_ave       + k_slow / &
            (real(soil_equil_year) * 365. * SHR_CONST_CDAY / dtime)
       litter_decom_ave = litter_decom_ave + litter_decom / &
            (real(soil_equil_year) * 365. * SHR_CONST_CDAY / dtime)

    else if (kyr == soil_equil_year+1) then

       ! SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
       ! Analytical solution of differential flux equations for fast and slow
       ! soil carbon pools.  Implemented after (soil_equil_year) simulation
       ! years, when annual litter inputs should be close to equilibrium.  Assumes
       ! average climate (temperature and soil moisture) from all years up to
       ! soil_equil_year.
       ! slevis: next could be done once

       ! Analytically calculate pool sizes this year only
       ! Rate of change of soil pool size = litter input - decomposition
       !   (5) dc/dt = litter_decom - kc
       ! At equilibrium,
       !   (6) dc/dt = 0
       ! From (5) & (6),
       !   (7) c = litter_decom / k
       cpool_fast=soilfrac*fastfrac*litter_decom_ave/k_fast_ave !eqn 7
       cpool_slow=soilfrac*slowfrac*litter_decom_ave/k_slow_ave !eqn 7

    endif

  end subroutine LitterSOM

#endif

end module EcosystemDynDGVMMod




