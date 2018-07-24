#include "misc.h"

module dynamics_vars
!BOP
!
! !MODULE: dynamics_vars --- Lin-Rood specific variables and methods
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! WS: ultimately we would like to uncouple fvcore from 
!     pmgrid since the latter has dependencies on params.h
!     All of the follow variables could go in spmd_dyn
!
   use pmgrid, only : twod_decomp, strip3dxyz, strip3dxzy,      &
                      strip3dxzyp, strip3dxyzp, strip3kxyz,     &
                      strip3kxyzp, strip3kxzy, strip3kxzyp

! !PUBLIC MEMBER FUNCTIONS:
   public dynamics_init, dynamics_clean

! !PUBLIC DATA MEMBERS:

!--------------------------------------------------------
! Local variables specific to the Lin-Rood dynamical core
!--------------------------------------------------------

! Temporary variables to be removed at a later time
   real(r8), allocatable, save :: yzt(:,:,:)
   real(r8), allocatable, save :: yzt2(:,:,:)
   real(r8), allocatable, save :: q3t(:,:,:,:)

! Variables set by dynamics_init
   real(r8) dtime        ! large time step
   integer  iord         ! order of LR scheme in X
   integer  jord         ! order of LR scheme in Y
   integer  im, jm, km   ! global dimensions
   integer  nq           ! number of tracers (adv. and non-adv)
   integer  ifirstxy, ilastxy !
   integer  jfirstxy, jlastxy !
   integer  jfirstyz, jlastyz !
   integer  kfirstyz, klastyz !


! General variables set by setrig
   real (r8) pi
   real(r8) dl           ! Radians per (finite-volume) longitude
   real(r8) dp           ! Radians per (finite-volume) latitude

! Geometric arrays set by setrig
   real(r8), allocatable, save :: cosp(:)   ! Cos of latitude angle - volume mean
   real(r8), allocatable, save :: sinp(:)   ! Sin of latitude angle - volume mean
   real(r8), allocatable, save :: cose(:)   ! Cos of finite-volume edges
   real(r8), allocatable, save :: sine(:)   ! Sin of finite-volume edges

! Variables set by set_eta
   real(r8) ptop         ! pressure at top of atmosphere
   real(r8) pint         ! pressure at sig-p interface
   real(r8), allocatable, save :: ak(:)   ! A's of the ETA-coordinate
   real(r8), allocatable, save :: bk(:)   ! B's of the ETA-coordinate
   integer  ks           ! Total number of pure-P layers

! Scalars set in dynpkg_init

   integer icd, jcd

   integer ng_c       ! ghost zone needed by the c-gird dynamics
   integer ng_d       ! ghost zone needed by the d-gird dynamics
   integer ng_s       ! for certain arrays, max(ng_c+1,ng_d)

   real (r8) acap     ! scaled polar cap area
   real (r8) rcap     ! inverse of scaled polar cap area

! Arrays initialized by dynpkg_init

   real(r8), allocatable, save :: coslon(:) ! Cos of longitudes - volume center
   real(r8), allocatable, save :: sinlon(:) ! Sin of longitudes - volume center

   real(r8), allocatable, save :: cosl5(:)  ! Cos of longitudes - volume center
   real(r8), allocatable, save :: sinl5(:)  ! Sin of longitudes - volume center

   real(r8), allocatable, save :: acosp(:)

! Scalars initialized by d_split
   integer   ns       ! total number of splits for Lagrangian dynamics

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in 
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        lr\_init    &  Initialize the Lin-Rood variables  \\ \hline
!        lr\_clean   &  Deallocate all internal data structures \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.06.06   Sawyer     Consolidated from various code snippets
!   01.07.12   Sawyer     Removed CCM common blocks comtim.h and commap.h
!   03.06.25   Sawyer     Cleaned up, used ParPatternCopy (Create)
!   03.07.23   Sawyer     Removed dependencies on params.h, constituents
!   03.08.05   Sawyer     Removed rayf_init and hswf_init, related vars
!   03.08.13   Sawyer     Removed xyt (no longer used)
!
! !BUGS:
!   o Where does the value of ns0 come from??
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_init --- initialize the lin-rood dynamical core
!
! !INTERFACE: 
   subroutine dynamics_init( dtime_in, iord_in, jord_in, nsplit_in,  &
                             im_in, jm_in, km_in, nq_in,             &
                             ifirstxy_in, ilastxy_in,                &
                             jfirstxy_in, jlastxy_in,                &
                             jfirstyz_in, jlastyz_in,                &
                             kfirstyz_in, klastyz_in )

! !USES:
      implicit none

! !INPUT PARAMETERS:
      real (r8), intent(in) :: dtime_in   !  Large time step
      integer, intent(in)   :: iord_in    !  Order of LR scheme in X
      integer, intent(in)   :: jord_in    !  Order of LR scheme in Y
      integer, intent(in)   :: nsplit_in  !  Order of LR scheme in Y
      integer, intent(in)   :: im_in, jm_in, km_in      !  Global dims
      integer, intent(in)   :: nq_in                    !  No. tracers
      integer, intent(in)   :: ifirstxy_in, ilastxy_in  !  Interval
      integer, intent(in)   :: jfirstxy_in, jlastxy_in  !  Interval
      integer, intent(in)   :: jfirstyz_in, jlastyz_in  !  Interval
      integer, intent(in)   :: kfirstyz_in, klastyz_in  !  Interval

! !DESCRIPTION:
!
!   Initialize Lin-Rood specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Create
!   03.07.31   Sawyer     Added the 'layout' arguments
!   03.08.05   Sawyer     Removed hswf_init and rayf_init
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
      dtime = dtime_in
      iord  = iord_in
      jord  = jord_in
      im    = im_in
      jm    = jm_in
      km    = km_in
      nq    = nq_in
      ifirstxy = ifirstxy_in
      ilastxy  = ilastxy_in
      jfirstxy = jfirstxy_in
      jlastxy  = jlastxy_in
      jfirstyz = jfirstyz_in
      jlastyz  = jlastyz_in
      kfirstyz = kfirstyz_in
      klastyz  = klastyz_in

      call setrig
      call set_eta
      call dynpkg_init
      call d_split(nsplit_in)
#if defined( SPMD )
      call spmd_vars_init()
#endif
      return
!EOC
   end subroutine dynamics_init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_clean -- clean up Lin-Rood-specific variables
!
! !INTERFACE: 
   subroutine dynamics_clean

      implicit none

! !DESCRIPTION:
!
! Clean up (deallocate) Lin-Rood-specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! Temporary data structures
#if defined(SPMD)
   call spmd_vars_clean()
#endif
   return
!EOC
   end subroutine dynamics_clean
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  setrig --- Specify the grid attributes
!
! !INTERFACE:
      subroutine setrig

! !USES:
      implicit none

!
! !DESCRIPTION:
!
!   Specify the grid attributes, such as the spacing between
!   grid points in latitude and longitude, the sines and cosines of
!   latitude at cell midpoints and edges.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j
      real*8 pi, ph5      ! This is to ensure 64-bit for any choice of r8

      pi  = 4.0 * atan(1.)
      dl  = (pi+pi)/im
      dp  = pi/(jm-1)

      allocate(cosp(jm))
      allocate(sinp(jm))
      allocate(cose(jm))
      allocate(sine(jm))

      do j=2,jm
         ph5  = -0.5d0*pi + ((j-1)-0.5d0)*(pi/(jm-1))
         sine(j) = sin(ph5)
      enddo

      cosp( 1) =  0.
      cosp(jm) =  0.

      do j=2,jm-1
         cosp(j) = (sine(j+1)-sine(j)) / dp
      enddo

! Define cosine at edges..

      do j=2,jm
         cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      enddo
         cose(1) = cose(2)

         sinp( 1) = -1.
         sinp(jm) =  1.

      do j=2,jm-1
         sinp(j) = 0.5 * (sine(j) + sine(j+1))
      enddo

      return
!EOC
      end subroutine setrig
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  set_eta --- Define vertical coordinate
!
! !INTERFACE:
      subroutine set_eta

! !USES:
      implicit none

!
! !DESCRIPTION:
!
!   Specify the vertical coordinate system.  Currently this is a 
!   dual pressure - sigma coordinate system (???), which transitions at
!   level ks, but it could be just about anything reasonable.
!
!   Choices for vertical resolutions are as follows:
!   \begin{tabular}{l}
!     NCAR: 18, 26, and 30 \\
!     NASA DAO: smoothed version of CAM's 30-level, 32, 55, 64, and 96 \\
!     New 32-layer setup with top at 0.4 mb for high horizontal
!     resolution runs. pint = 176.93 mb \\
!     Revised 55-level eta with pint at 176.93 mb  SJL: 2000-03-20
!   \end{tabular}
! 
! !REVISION HISTORY: 
!   98.01.15    Lin        Creation
!   ongoing     Lin        Fine tuning
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! NCAR specific
      real(r8) a18(19),b18(19)              ! CCM3
      real(r8) a26(27),b26(27)              ! CAM
      real(r8) a30(31),b30(31)              ! CCM3

! NASA only
      real(r8) a30m(31),b30m(31)            ! smoothed CAM 30-L
      real(r8) a32(33),b32(33)
      real(r8) a32old(33),b32old(33)
      real(r8) a55(56),b55(56)
      real(r8) a55old(56),b55old(56)
      real(r8) a64(65),b64(65)
      real(r8) a96(97),b96(97)

      integer k

! *** NCAR settings ***

      data a18 /291.70,  792.92,  2155.39,  4918.34,  8314.25, &
               7993.08, 7577.38,  7057.52,  6429.63,  5698.38, &
               4879.13, 3998.95,  3096.31,  2219.02,  1420.39, &
               754.13,  268.38,   0.0000,   0.0000 /

      data b18 /0.0000,    0.0000,    0.0000,   0.0000,   0.0000, &
                0.0380541, 0.0873088, 0.1489307, 0.2232996,       &
                0.3099406, 0.4070096, 0.5112977, 0.6182465,       &
                0.7221927, 0.8168173, 0.8957590, 0.9533137,       &
                0.9851122, 1.0  /
     
      data a26 /219.4067,  489.5209,   988.2418,   1805.201,      &
                2983.724,  4462.334,   6160.587,   7851.243,      &
                7731.271,  7590.131,   7424.086,   7228.744,      &
                6998.933,  6728.574,   6410.509,   6036.322,      &
                5596.111,  5078.225,   4468.96,    3752.191,      &
                2908.949,  2084.739,   1334.443,   708.499,       &
                252.136,   0.,         0. /

      data b26 /0.,         0.,         0.,         0.,           &
                0.,         0.,         0.,         0.,           &
                0.01505309, 0.03276228, 0.05359622, 0.07810627,   &
                0.1069411,  0.14086370, 0.180772,   0.227722,     &
                0.2829562,  0.3479364,  0.4243822,  0.5143168,    &
                0.6201202,  0.7235355,  0.8176768,  0.8962153,    &
                0.9534761,  0.9851122,  1.        /

      data a30 /225.523952394724, 503.169186413288, 1015.79474285245,  &
               1855.53170740604, 3066.91229343414,  4586.74766123295,  &
               6332.34828710556, 8070.14182209969,  9494.10423636436,  &
              11169.321089983,  13140.1270627975,  15458.6806893349,   &
              18186.3352656364, 17459.799349308,   16605.0657629967,   &
              15599.5160341263, 14416.541159153,   13024.8308181763,   &
              11387.5567913055,  9461.38575673103,  7534.44507718086,  &
               5765.89405536652, 4273.46378564835,  3164.26791250706,  &
               2522.12174236774, 1919.67375576496,  1361.80268600583,  &
                853.108894079924, 397.881818935275,    0.,             &
                  0.  /

      data b30 /0.,                 0.,                                 &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.03935482725501,&
                0.085653759539127,  0.140122056007385, 0.20420117676258,&
                0.279586911201477,  0.368274360895157, 0.47261056303978,&
                0.576988518238068,  0.672786951065063, 0.75362843275070,&
                0.813710987567902,  0.848494648933411, 0.88112789392471,&
                0.911346435546875,  0.938901245594025, 0.96355980634689,&
                0.985112190246582,  1.   /

! *** NASA DAO settings ***

! Smoothed CAM's 30-Level setup
      data a30m / 300.00000,     725.00000,    1500.00000,     &
             2600.00000,    3800.00000,    5050.00000,         &
             6350.00000,    7750.00000,    9300.00000,         &
            11100.00000,   13140.00000,   15458.00000,         &
            18186.33580,   20676.23761,   22275.23783,         &
            23025.65071,   22947.33569,   22038.21991,         &
            20274.24578,   17684.31619,   14540.98138,         &
            11389.69990,    8795.97971,    6962.67963,         &
             5554.86684,    4376.83633,    3305.84967,         &
             2322.63910,    1437.78398,     660.76994,         &
                0.00000 /

      data b30m / 0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00719,       0.02895,       &
                  0.06586,       0.11889,       0.18945,       &
                  0.27941,       0.38816,       0.50692,       &
                  0.61910,       0.70840,       0.77037,       &
                  0.81745,       0.85656,       0.89191,       &
                  0.92421,       0.95316,       0.97850,       &
                  1.00000 /

      data a32/40.00000,     100.00000,     200.00000,         &
            370.00000,     630.00000,    1000.00000,           &
           1510.00000,    2160.00000,    2900.00000,           &
           3680.00000,    4535.00000,    5505.00000,           &
           6607.26750,    7851.22980,    9236.56610,           &
          10866.34270,   12783.70000,   15039.30000,           &
          17693.00000,   20119.20876,   21686.49129,           &
          22436.28749,   22388.46844,   21541.75227,           &
          19873.78342,   17340.31831,   13874.44006,           &
          10167.16551,    6609.84274,    3546.59643,           &
           1270.49390,       0.00000,       0.00000   /

      data b32/0.00000,       0.00000,       0.00000,          &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00696,       0.02801,            &
             0.06372,       0.11503,       0.18330,            &
             0.27033,       0.37844,       0.51046,            &
             0.64271,       0.76492,       0.86783,            &
             0.94329,       0.98511,       1.00000   /

      data a32old /300.0000,  454.1491,  652.5746,  891.9637, 1159.7102, &
             1492.8248, 1902.5026, 2400.4835, 2998.6740, 3708.6584,      &
             4541.1041, 5505.0739, 6607.2675, 7851.2298, 9236.5661,      &
            10866.3427, 12420.400, 13576.500, 14365.400, 14807.800,      &
             14915.500, 14691.400, 14129.400, 13214.800, 11923.200,      &
             10220.700,  8062.000,  5849.500,  3777.000,  2017.200,      &
               720.600,     0.000,     0.000 /

      data b32old /0.00, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,     &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,     &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,     &
              0.948792 , 0.9851119, 1.0000000 /

      data a55/ 1.00000,       2.00000,       3.27000,              &
              4.75850,       6.60000,       8.93450,                &
             11.97030,      15.94950,      21.13490,                &
             27.85260,      36.50410,      47.58060,                &
             61.67790,      79.51340,     101.94420,                &
            130.05080,     165.07920,     208.49720,                &
            262.02120,     327.64330,     407.65670,                &
            504.68050,     621.68000,     761.98390,                &
            929.29430,    1127.68880,    1364.33920,                &
           1645.70720,    1979.15540,    2373.03610,                &
           2836.78160,    3380.99550,    4017.54170,                &
           4764.39320,    5638.79380,    6660.33770,                &
           7851.22980,    9236.56610,   10866.34270,                &
          12783.70000,   15039.30000,   17693.00000,                &
          20119.20876,   21686.49129,   22436.28749,                &
          22388.46844,   21541.75227,   19873.78342,                &
          17340.31831,   13874.44006,   10167.16551,                &
           6609.84274,    3546.59643,    1270.49390,                &
              0.00000,       0.00000   /

      data b55 / 0.00000,       0.00000,       0.00000,         &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00696,       0.02801,       0.06372,           &
               0.11503,       0.18330,       0.27033,           &
               0.37844,       0.51046,       0.64271,           &
               0.76492,       0.86783,       0.94329,           &
               0.98511,       1.00000  /

      data a55old /1.0000,    2.0000,    3.2700,    4.7585,     6.6000, &
              8.9345,   11.9703,   15.9495,   21.1349,    27.8526,      &
             36.5041,   47.5806,   61.6779,   79.5134,   101.9442,      &
            130.0508,  165.0792,  208.4972,  262.0212,   327.6433,      &
            407.6567,  504.6805,  621.6800,  761.9839,   929.2943,      &
           1127.6888, 1364.3392, 1645.7072, 1979.1554,  2373.0361,      &
           2836.7816, 3380.9955, 4017.5417, 4764.3932,  5638.7938,      &
           6660.3377, 7851.2298, 9236.5661,10866.3427, 12420.400 ,      &
          13576.500 , 14365.400, 14807.800, 14915.500 , 14691.400,      &
          14129.400 , 13214.800, 11923.200, 10220.700 ,  8062.000,      &
           5849.500 ,  3777.000,  2017.200,   720.600,      0.000,      &
              0.000 /

      data b55old /   0.0000000, 0.0000000, 0.0000000,                  &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,    &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,    &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,    &
              0.948792 , 0.9851119, 1.0000000 /

      data a64/1.00000,       1.80   ,       2.80086,     &
             3.93309,       5.20139,       6.77626,       &
             8.69654,      10.99483,      13.81736,       &
            17.26058,      21.43286,      26.45448,       &
            32.45730,      39.58402,      47.98678,       &
            57.82525,      69.26401,      82.46925,       &
            97.60468,     114.82686,     135.08787,       &
           158.92390,     186.96575,     219.95555,       &
           258.76633,     304.42522,     358.14053,       &
           421.33383,     495.67748,     583.13893,       &
           686.03282,     807.08215,     949.49044,       &
          1117.02644,    1314.12387,    1545.99882,       &
          1818.78771,    2139.70974,    2517.25793,       &
          2961.42386,    3483.96212,    4098.70138,       &
          4821.91034,    5672.72831,    6673.67169,       &
          7851.22983,    9236.56613,   10866.34270,       &
         12783.69059,   15039.35130,   17693.01955,       &
         20814.92310,   23609.16397,   25271.17281,       &
         25844.93368,   25345.63084,   23760.05052,       &
         21046.23129,   17132.35351,   12832.00555,       &
          8646.27815,    5012.23907,    2299.34286,       &
           781.15294,       0.00000  /

      data b64/0.00000,       0.00000,       0.00000,     &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00879,       0.03537,       &
             0.08047,       0.14526,       0.23147,       &
             0.34138,       0.47789,       0.61606,       &
             0.74456,       0.85318,       0.93300,       &
             0.97730,       1.00000  /

      data a96/  1.00000,       2.32782,       3.34990,   &
               4.49484,       5.62336,       6.93048,     &
               8.41428,      10.06365,      11.97630,     &
              14.18138,      16.70870,      19.58824,     &
              22.84950,      26.52080,      30.62845,     &
              35.19588,      40.24273,      45.78375,     &
              51.82793,      58.43583,      65.62319,     &
              73.40038,      81.77154,      90.73373,     &
             100.27628,     110.82243,     122.47773,     &
             135.35883,     149.59464,     165.32764,     &
             182.71530,     201.93164,     223.16899,     &
             246.63988,     272.57922,     301.24661,     &
             332.92902,     367.94348,     406.64044,     &
             449.40720,     496.67181,     548.90723,     &
             606.63629,     670.43683,     740.94727,     &
             818.87329,     904.99493,    1000.17395,     &
            1105.36304,    1221.61499,    1350.09326,     &
            1492.08362,    1649.00745,    1822.43469,     &
            2014.10168,    2225.92627,    2460.02905,     &
            2718.75195,    3004.68530,    3320.69092,     &
            3669.93066,    4055.90015,    4482.46240,     &
            4953.88672,    5474.89111,    6050.68994,     &
            6687.04492,    7390.32715,    8167.57373,     &
            9026.56445,    9975.89648,   11025.06934,     &
           12184.58398,   13466.04785,   14882.28320,     &
           16447.46289,   18177.25781,   20088.97461,     &
           21886.89453,   23274.16602,   24264.66602,     &
           24868.31641,   25091.15430,   24935.41016,     &
           24399.52148,   23478.13281,   22162.01758,     &
           20438.00586,   18288.83984,   15693.01172,     &
           12624.54199,    9584.35352,    6736.55713,     &
            4231.34326,    2199.57910,     747.11890,     &
              0.00000 /

      data b96/0.00000,       0.00000,       0.00000,     &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00315,       0.01263,       0.02853,      &
              0.05101,       0.08030,       0.11669,      &
              0.16055,       0.21231,       0.27249,      &
              0.34169,       0.42062,       0.51005,      &
              0.61088,       0.70748,       0.79593,      &
              0.87253,       0.93400,       0.97764,      &
              1.00000 /

      allocate(ak(km+1))
      allocate(bk(km+1))

      select case (km)

! *** Original CCM3 18-Level setup ***
        case (18)
          ks = 4
          do k=1,km+1
            ak(k) = a18(k)
            bk(k) = b18(k)
          enddo

        case (26)
! CAM 26-Level setup ***
          ks = 7
          do k=1,km+1
            ak(k) = a26(k)
            bk(k) = b26(k)
          enddo

        case (30)
! CAM 30-Level setup ***
          ks = 12
          do k=1,km+1
            ak(k) = a30(k)
            bk(k) = b30(k)
          enddo

! *** Revised 32-L setup with ptop at 0.4 mb ***
        case (32)
          ks = 18
          do k=1,km+1
            ak(k) = a32(k)
            bk(k) = b32(k)
          enddo

! *** Revised 55-L setup with ptop at 0.01 mb ***
        case (55)
          ks = 41
          do k=1,km+1
            ak(k) = a55(k)
            bk(k) = b55(k)
          enddo

! *** Others ***
        case (64)
          ks = 51
          do k=1,km+1
            ak(k) = a64(k)
            bk(k) = b64(k)
          enddo

        case (96)
          ks = 77
          do k=1,km+1
            ak(k) = a96(k)
            bk(k) = b96(k)
          enddo

      end select

          ptop = ak(1)
          pint = ak(ks+1) 

      return
!EOC
      end subroutine set_eta
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  dynpkg_init --- Initialization for dynamics package
!
! !INTERFACE:
subroutine dynpkg_init

! !USES:
   implicit none

! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the Rayleigh friction
! 
! !REVISION HISTORY: 
!   00.01.10    Grant        Creation using code from SJ Lin
!   01.03.26    Sawyer       Added ProTeX documentation
!   01.06.06    Sawyer       Modified for dynamics_vars
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, imh
      real (r8) zam5, zamda

      allocate(coslon(im))
      allocate(sinlon(im))
      allocate(cosl5(im))
      allocate(sinl5(im))
      allocate(acosp(jm))

      if( iord <= 2 ) then
         icd =  1
      else
         icd = -2
      endif
 
      if( jord <= 2 ) then
         jcd =  1
      else
         jcd =  -2
      endif

#if defined( SPMD )
!
! Calculate the ghost region sizes for the SPMD version (tricky stuff)
!
      ng_c = min(abs(jcd ), 2)
      ng_d = min(abs(jord), 3)    ! SJL: number of max ghost latitudes
      ng_d = max(ng_d, 2)
      ng_s = max( ng_c+1, ng_d )
#else
      ng_c = 0
      ng_d = 0                   ! No ghosting necessary for pure SMP runs
      ng_s = 0
#endif

!
! Pole cap area and inverse
      acap = im*(1.+sine(2)) / dp
      rcap = 1.d0 / acap
 
      imh = im/2
      if(im .ne. 2*imh) then
         write(6,*) 'im must be an even integer'
         stop
      endif
 
! Define logitude at the center of the volume
! i=1, Zamda = -pi
 
      do i=1,imh
         zam5          = ((i-1)-0.5d0) * dl
         cosl5(i)      =  cos(zam5)
         cosl5(i+imh)  = -cosl5(i)
         sinl5(i)      =  sin(zam5)
         sinl5(i+imh)  = -sinl5(i)
         zamda         = (i-1)*dl
         coslon(i)     =  cos(zamda)
         coslon(i+imh) = -coslon(i)
         sinlon(i)     =  sin(zamda)
         sinlon(i+imh) = -sinlon(i)
      enddo

      do j=2,jm-1
         acosp(j) = 1.d0 / cosp(j)
      enddo
      acosp( 1) = rcap * im
      acosp(jm) = rcap * im
      return
!EOC
end subroutine dynpkg_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  d_split --- find proper value for nsplit if not specified
!
! !INTERFACE:
      subroutine d_split(nsplit_in)
!
! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in)   :: nsplit_in  !  Small time steps in dtime

! !DESCRIPTION:
!
!    If nsplit=0 (module variable) then determine a good value 
!    for ns (used in dynpkg) based on resolution and the large-time-step 
!    (pdt). The user may have to set this manually if instability occurs.
!
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real (r8)   pdt                       ! Time-step in seconds
                                            ! Negative dt (backward in time
                                            ! integration) is allowed
      real (r8)   dim
      real (r8)   dim0                      ! base dimension
      real (r8)   dt0                       ! base time step
      real (r8)   ns0                       ! base nsplit for base dimension

      parameter ( dim0 = 180.  )
      parameter ( dt0  = 1800. )
      parameter ( ns0  = 4.    )

      if ( nsplit_in == 0 ) then
          pdt = int(dtime)   ! dtime is a variable internal to this module
          dim    = max ( im, 2*(jm-1) )
          ns = int ( ns0*abs(pdt)*dim/(dt0*dim0) + 0.75 )
          ns = max ( 1, ns )   ! for cases in which dt or dim is too small
      else
          ns = nsplit_in
      endif
      write(6,*) 'Lagrangian time splits (NSPLIT) =', ns

      return
!EOC
      end subroutine d_split
!---------------------------------------------------------------------

#if defined(SPMD)
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  spmd_vars_init --- Initialization of SPMD-related variables
!
! !INTERFACE:
subroutine spmd_vars_init

! !USES:
   use spmd_dyn, only : ghostu_yz, ghostv_yz, ghostq_yz,                &
                        ghostpe_xy, ghostpe_yz, ghostp_yz,              &
                        ghostpe1_yz, m_ttrans, q_ttrans, r_ttrans,      &
                        u_to_uxy, uxy_to_u, v_to_vxy, vxy_to_v,         &
                        ijk_yz_to_xy, ijk_xy_to_yz, ikj_yz_to_xy,       &
                        ikj_xy_to_yz, q_to_qxy, qxy_to_q,               &
                        pe_to_pexy, pexy_to_pe, pt_to_ptxy, ptxy_to_pt, &
                        pkxy_to_pkc, q3_to_qxy3, qxy3_to_q3, r_to_rxy,  &
                        rxy_to_r
   use ghostmodule, only : ghostcreate
   use parutilitiesmodule, only : gid, commglobal, parpatterncreate
   implicit none

!------------------------------Commons----------------------------------


! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the SPMD related variables.
!   This has to be done in this module since certain variables
!   (in particular the ghost sizes {\tt ng\_d, ng\_s} are first
!   defined here.
! 
! !REVISION HISTORY: 
!   02.11.08    Sawyer       Creation
!   03.05.07    Sawyer       Use ParPatternCopy for q_to_qxy, etc.
!   03.07.23    Sawyer       Removed dependency on constituents module
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! Temporary data structures
      allocate(   yzt(im,jfirstyz:jlastyz,kfirstyz:klastyz) )
      allocate(   yzt2(im,jfirstyz:jlastyz,kfirstyz:klastyz) )
      allocate(   q3t(im,jfirstyz:jlastyz,kfirstyz:klastyz,nq) )

      if ( twod_decomp == 1 ) then
! Initialize ghost regions
!
        call t_startf('ghost_creation')
        call ghostcreate( strip3dxyz, gid, im, 1, im, .true., &
                          jm, jfirstyz-ng_d, jlastyz+ng_s, .false., &
                          km, kfirstyz, klastyz, .false., ghostu_yz )
        call ghostcreate( strip3dxyz, gid, im, 1, im, .true., &
                          jm, jfirstyz-ng_s, jlastyz+ng_d, .false., &
                          km, kfirstyz, klastyz, .false., ghostv_yz )
        call ghostcreate( strip3dxzyp, gid, im, 1, im, .true., &
                          km+1, kfirstyz, klastyz+1, .false., &
                          jm, jfirstyz, jlastyz, .false., ghostpe_yz)
        call ghostcreate( strip3kxzyp, gid, im, ifirstxy, ilastxy, .true.,  &
                          km+1, 1, km+1, .false.,   &
                          jm, jfirstyz, jlastyz, .false., ghostpe_xy)
        call ghostcreate( strip3dxyz, gid, im, 1, im, .true., &
                          jm, jfirstyz-ng_d, jlastyz+ng_d, .false., &
                          km, kfirstyz, klastyz, .false., ghostp_yz )
        call ghostcreate( strip3dxyzp, gid, im, 1, im, .true., &
                          jm, jfirstyz, jlastyz, .false.,       &
                          km+1, kfirstyz, klastyz+1, .false., ghostpe1_yz)
        call t_stopf('ghost_creation')

! Initialize transposes
!
        call t_startf('transpose_creation')
!       call parpatterncreate(commglobal, ghostu_yz, strip3kxyz, u_to_uxy)
!       call parpatterncreate(commglobal, strip3kxyz, ghostu_yz, uxy_to_u)
!       call parpatterncreate(commglobal, ghostv_yz, strip3kxyz, v_to_vxy)
!       call parpatterncreate(commglobal, strip3kxyz, ghostv_yz, vxy_to_v)
        call parpatterncreate(commglobal, strip3dxyz, strip3kxyz, ijk_yz_to_xy)
        call parpatterncreate(commglobal, strip3kxyz, strip3dxyz, ijk_xy_to_yz)
!       call parpatterncreate(commglobal, strip3dxzy, strip3kxzy, ikj_yz_to_xy)
        call parpatterncreate(commglobal, strip3kxzy, strip3dxzy, ikj_xy_to_yz)
        call t_startf('q_to_qxy_creation')
        if (q_ttrans .ne. 0) then
          call parpatterncreate( ijk_yz_to_xy, q_to_qxy, m_ttrans )
          call parpatterncreate( ijk_xy_to_yz, qxy_to_q, m_ttrans )
        endif
        if (r_ttrans .ne. 0) then
          call parpatterncreate( ijk_yz_to_xy, r_to_rxy, r_ttrans )
          call parpatterncreate( ijk_xy_to_yz, rxy_to_r, r_ttrans )
        endif
        call t_stopf('q_to_qxy_creation')
!       call parpatterncreate(commglobal, ghostpe_yz, strip3kxzyp, pe_to_pexy)
        call parpatterncreate(commglobal, strip3kxzyp, ghostpe_yz, pexy_to_pe)
!       call parpatterncreate(commglobal, ghostp_yz, strip3kxyz, pt_to_ptxy)
!       call parpatterncreate(commglobal, strip3kxyz, ghostp_yz, ptxy_to_pt)
        call parpatterncreate(commglobal, strip3kxyzp, ghostpe1_yz, pkxy_to_pkc)
        call t_stopf('transpose_creation')
      endif
   return
!EOC
end subroutine spmd_vars_init
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  spmd_vars_clean --- Clean the SPMD-related variables
!
! !INTERFACE:
subroutine spmd_vars_clean

! !USES:
   use spmd_dyn, only : ghostu_yz, ghostv_yz,                           &
                        ghostq_yz, ghostpe_yz, ghostpe_xy,              &
                        ghostp_yz, ghostpe1_yz,                         &
                        m_ttrans, q_ttrans, r_ttrans,                   &
                        u_to_uxy, uxy_to_u, v_to_vxy, vxy_to_v,         &
                        ijk_yz_to_xy, ijk_xy_to_yz, ikj_yz_to_xy,       &
                        ikj_xy_to_yz, q_to_qxy, qxy_to_q,               &
                        pe_to_pexy, pexy_to_pe, pt_to_ptxy, ptxy_to_pt, &
                        pkxy_to_pkc, q3_to_qxy3, qxy3_to_q3,            &
                        r_to_rxy, rxy_to_r
   use ghostmodule, only : ghostfree
   use parutilitiesmodule, only : commglobal, parpatternfree
   implicit none

!------------------------------Commons----------------------------------


! !DESCRIPTION:
! 
!   {\bf Purpose:} Clean the SPMD related variables.
! 
! !REVISION HISTORY: 
!   02.11.08    Sawyer       Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! Temporary data structures
   deallocate( yzt )
   deallocate( yzt2 )
   deallocate( q3t )

      if ( twod_decomp == 1 ) then
! Clean the ghost regions
!
        call ghostfree( ghostu_yz )
        call ghostfree( ghostv_yz )
!       call ghostfree( ghostq_yz )
        call ghostfree( ghostpe_yz )
        call ghostfree( ghostpe_xy )
        call ghostfree( ghostp_yz )
        call ghostfree( ghostpe1_yz )
! Clean transposes
!
!       call parpatternfree(commglobal,u_to_uxy)
!       call parpatternfree(commglobal,uxy_to_u)
!       call parpatternfree(commglobal,v_to_vxy)
!       call parpatternfree(commglobal,vxy_to_v)
        call parpatternfree(commglobal,ijk_yz_to_xy)
        call parpatternfree(commglobal,ijk_xy_to_yz)
        call parpatternfree(commglobal,ikj_xy_to_yz)
!       call parpatternfree(commglobal,ikj_yz_to_xy)
!       call parpatternfree(commglobal,q3_to_qxy3)
!       call parpatternfree(commglobal,qxy3_to_q3)
        if (q_ttrans .ne. 0) then
          call parpatternfree(commglobal,q_to_qxy)
          call parpatternfree(commglobal,qxy_to_q)
        endif
        if (r_ttrans .ne. 0) then
          call parpatternfree(commglobal,r_to_rxy)
          call parpatternfree(commglobal,rxy_to_r)
        endif
!       call parpatternfree(commglobal,pe_to_pexy)
        call parpatternfree(commglobal,pexy_to_pe)
!       call parpatternfree(commglobal,pt_to_ptxy)
!       call parpatternfree(commglobal,ptxy_to_pt)
        call parpatternfree(commglobal,pkxy_to_pkc)
      endif
   return
!EOC
end subroutine spmd_vars_clean
!-----------------------------------------------------------------------
#endif
end module dynamics_vars

