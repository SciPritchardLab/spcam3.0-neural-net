
module sulchem

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Contains the sulfur cycle chemistry codes.
   ! 
   ! Author: 
   ! Chemistry code is from Mary Barth.
   ! Module coded by B. Eaton.
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,      only: pcols, pver
   use pmgrid

   implicit none
   save
   private
   public :: &
      inisulchem,   &! Initialize sulfur cycle chemistry.
      chemwdepdr     ! Chemistry and wet deposition driver.

   ! private module data

   integer, parameter :: DBLKIND = selected_real_kind(12)

   integer, parameter ::&
      nz=18,              &
      nza=8,              &
      nalb=2,             &
      nz1=nz-1,           &
      nza1=nza-1,         &
      nz2=51,             &
      nsiz =nz *nza*nalb, &
      nsiz2=nz2*nza*nalb

   real(r8) ::&
      jper2(nsiz2),          &! j values interpolated to new ht grid
      delz(nz1),             &! Difference between tabl. z's
      delang(nza1),          &! Difference between tabl. sec(zen.ang)
      delalb,                &! Difference in tabl. albedo = .5-.2
      pie,                   &! 3.14159...
      dayspy                  ! Number of days per 1 year

   real(r8), parameter :: mwa2so2 = 28.966/64. ! ratio molecular wt of air to so2
   real(r8), parameter :: mwa2h2o2 = 28.966/34. ! ratio molecular wt of air to h2o2


   real(r8), parameter, dimension(nsiz) :: jper =    &! Tabulated j values
      (/ -11.63010025,-11.53059959,-11.47570038,-11.44130039,-11.40859985,-11.40690041,-11.42399979,-11.46020031, &
         -11.48509979,-11.49100018,-11.46759987,-11.41129971,-11.32050037,-11.20450020,-10.96660042,-10.65079975, &
         -10.31389999,-10.15680027,-11.93789959,-11.81239986,-11.73900032,-11.68900013,-11.62989998,-11.60620022, &
         -11.60490036,-11.61820030,-11.62650013,-11.62040043,-11.58790016,-11.52509975,-11.43029976,-11.31110001, &
         -11.06599998,-10.75179958,-10.39190006,-10.20240021,-12.20409966,-12.06019974,-11.97229958,-11.90939999, &
         -11.82750034,-11.78409958,-11.76560020,-11.75730038,-11.74950027,-11.73139954,-11.68959999,-11.62010002, &
         -11.52070045,-11.39859962,-11.14780045,-10.83170033,-10.46100044,-10.24470043,-12.51459980,-12.35309982, &
         -12.24969959,-12.17240047,-12.06379986,-11.99619961,-11.95629978,-11.92070007,-11.89220047,-11.85890007, &
         -11.80539989,-11.72710037,-11.62139988,-11.49520016,-11.23849964,-10.91720009,-10.54150009,-10.29699993, &
         -13.14570045,-12.95919991,-12.83010006,-12.72649956,-12.56420040,-12.44439983,-12.35589981,-12.25739956, &
         -12.18089962,-12.11240005,-12.03289986,-11.93480015,-11.81389999,-11.67720032,-11.40919971,-11.07330036, &
         -10.69830036,-10.41100025,-14.34010029,-14.15550041,-14.01509953,-13.88949966,-13.65480042,-13.43579960, &
         -13.23670006,-12.97910023,-12.77519989,-12.61489964,-12.47119999,-12.32730007,-12.17119980,-12.00699997, &
         -11.70969963,-11.35060024,-10.96660042,-10.65600014,-15.41930008,-15.26990032,-15.16590023,-15.07689953, &
         -14.89920044,-14.68239975,-14.41380024,-13.96300030,-13.55609989,-13.23820019,-12.98649979,-12.77019978, &
         -12.56400013,-12.36279964,-12.01889992,-11.63490009,-11.22630024,-10.92059994,-17.23080063,-17.09090042, &
         -16.99939919,-16.92959976,-16.82430077,-16.74169922,-16.66230011,-16.49410057,-16.09469986,-15.37609959, &
         -14.65410042,-14.07880020,-13.62609959,-13.26119995,-12.75020027,-12.25749969,-11.79109955,-11.45209980, &
         -11.17109966,-11.14850044,-11.14080048,-11.14089966,-11.15730000,-11.18900013,-11.23019981,-11.29240036, &
         -11.33609962,-11.35649967,-11.34700012,-11.30399990,-11.22620010,-11.12189960,-10.90159988,-10.60239983, &
         -10.27840042,-10.12619972,-11.52420044,-11.47599983,-11.44919968,-11.43290043,-11.42099953,-11.42879963, &
         -11.44970036,-11.48649979,-11.51080036,-11.51659966,-11.49520016,-11.44289970,-11.35809994,-11.24800014, &
         -11.01669979,-10.71510029,-10.36550045,-10.18050003,-11.81890011,-11.75230026,-11.71059990,-11.68099976, &
         -11.64519978,-11.63220024,-11.63479996,-11.64830017,-11.65499973,-11.64729977,-11.61489964,-11.55410004, &
         -11.46290016,-11.34809971,-11.10849953,-10.80249977,-10.44029999,-10.22789955,-12.15250015,-12.06820011, &
         -12.01109982,-11.96700001,-11.90400028,-11.86610031,-11.84650040,-11.83150005,-11.81620026,-11.79199982, &
         -11.74639988,-11.67520046,-11.57610035,-11.45580006,-11.20800018,-10.89459991,-10.52560043,-10.28439999, &
         -12.80440044,-12.69659996,-12.61540031,-12.54609966,-12.43089962,-12.34109974,-12.27270031,-12.19359970, &
         -12.12870026,-12.06789970,-11.99429989,-11.90139961,-11.78509998,-11.65240002,-11.39009953,-11.05939960, &
         -10.68850040,-10.40359974,-13.98359966,-13.87979984,-13.79139996,-13.70489979,-13.52700043,-13.34549999, &
         -13.17129993,-12.93599987,-12.74400043,-12.59039974,-12.45110035,-12.31060028,-12.15719986,-11.99520016, &
         -11.70090008,-11.34430027,-10.96230030,-10.65270042,-15.04539967,-14.97089958,-14.91469955,-14.86209965, &
         -14.74139977,-14.57019997,-14.33740044,-13.92109966,-13.53120041,-13.22140026,-12.97410011,-12.76060009, &
         -12.55630016,-12.35649967,-12.01439953,-11.63179970,-11.22420025,-10.91909981,-16.85479927,-16.78750038, &
         -16.74150085,-16.70509911,-16.64769936,-16.59880066,-16.54520035,-16.40889931,-16.04310036,-15.35249996, &
         -14.64319992,-14.07279968,-13.62240028,-13.25860023,-12.74860001,-12.25650024,-11.79049969,-11.45160007 /)

   real(r8), parameter, dimension(nz) :: zz =     &! Tabulated heights
      (/ 0., 1., 2., 3., 5., 7., 9., 12., 15., 18., 21., 24., &
         27., 30., 35., 40., 45., 50. /)

   real(r8), parameter, dimension(nza) :: zasec = &! Tabulated sec(zen.ang)
      (/ 1., 1.3, 1.6, 2., 3., 6., 12., 50. /)

   real(r8), parameter, dimension(nalb) :: alb =  &! Tabulated albedo = 0.2 and 0.5
      (/ 0.2, 0.5 /)

   real(r8), parameter, dimension(121) :: yield = &! yield of DMS+OH rxn going to SO2 
      (/ 0.76905,  0.76953,  0.77002,  0.77050,  0.77098,  0.77146, &
         0.77195,  0.77243,  0.77291,  0.77340,  0.77388,  0.77436, &
         0.77485,  0.77533,  0.77581,  0.77629,  0.77678,  0.77726, &
         0.77774,  0.77822,  0.77870,  0.77918,  0.77967,  0.78015, &
         0.78063,  0.78111,  0.78159,  0.78207,  0.78255,  0.78303, &
         0.78351,  0.78399,  0.78447,  0.78495,  0.78543,  0.78591, &
         0.78640,  0.78688,  0.78737,  0.78786,  0.78835,  0.78884, &
         0.78934,  0.78984,  0.79035,  0.79086,  0.79137,  0.79190, &
         0.79243,  0.79297,  0.79353,  0.79409,  0.79467,  0.79526, &
         0.79587,  0.79650,  0.79715,  0.79783,  0.79853,  0.79927, &
         0.80004,  0.80084,  0.80169,  0.80258,  0.80352,  0.80452, &
         0.80558,  0.80671,  0.80790,  0.80918,  0.81055,  0.81200, &
         0.81356,  0.81522,  0.81700,  0.81890,  0.82093,  0.82309, &
         0.82540,  0.82786,  0.83047,  0.83325,  0.83618,  0.83928, &
         0.84255,  0.84598,  0.84957,  0.85332,  0.85722,  0.86126, &
         0.86543,  0.86973,  0.87412,  0.87861,  0.88316,  0.88777, &
         0.89240,  0.89705,  0.90169,  0.90631,  0.91087,  0.91537, &
         0.91978,  0.92409,  0.92829,  0.93236,  0.93630,  0.94009, &
         0.94373,  0.94721,  0.95054,  0.95370,  0.95670,  0.95954, &
         0.96223,  0.96476,  0.96714,  0.96938,  0.97148,  0.97344, &
         0.97528 /)

   real(r8), parameter, dimension(121) :: dmsrate = &! rate coeff of DMS+OH 
      (/ 0.21261E-10,0.21112E-10,0.20966E-10,0.20823E-10,0.20683E-10, &
         0.20546E-10,0.20411E-10,0.20280E-10,0.20151E-10,0.20024E-10, &
         0.19900E-10,0.19779E-10,0.19660E-10,0.19543E-10,0.19428E-10, &
         0.19316E-10,0.19206E-10,0.19098E-10,0.18991E-10,0.18887E-10, &
         0.18785E-10,0.18684E-10,0.18585E-10,0.18488E-10,0.18393E-10, &
         0.18299E-10,0.18207E-10,0.18116E-10,0.18027E-10,0.17939E-10, &
         0.17852E-10,0.17766E-10,0.17682E-10,0.17598E-10,0.17516E-10, &
         0.17434E-10,0.17354E-10,0.17274E-10,0.17194E-10,0.17115E-10, &
         0.17036E-10,0.16958E-10,0.16880E-10,0.16801E-10,0.16722E-10, &
         0.16643E-10,0.16563E-10,0.16482E-10,0.16401E-10,0.16317E-10, &
         0.16233E-10,0.16146E-10,0.16057E-10,0.15966E-10,0.15871E-10, &
         0.15774E-10,0.15673E-10,0.15567E-10,0.15458E-10,0.15343E-10, &
         0.15223E-10,0.15097E-10,0.14965E-10,0.14826E-10,0.14680E-10, &
         0.14526E-10,0.14365E-10,0.14195E-10,0.14016E-10,0.13828E-10, &
         0.13632E-10,0.13426E-10,0.13211E-10,0.12987E-10,0.12754E-10, &
         0.12514E-10,0.12265E-10,0.12009E-10,0.11747E-10,0.11479E-10, &
         0.11207E-10,0.10932E-10,0.10655E-10,0.10376E-10,0.10098E-10, &
         0.98217E-11,0.95482E-11,0.92787E-11,0.90145E-11,0.87565E-11, &
         0.85057E-11,0.82630E-11,0.80289E-11,0.78042E-11,0.75893E-11, &
         0.73845E-11,0.71899E-11,0.70058E-11,0.68321E-11,0.66687E-11, &
         0.65155E-11,0.63722E-11,0.62384E-11,0.61140E-11,0.59985E-11, &
         0.58915E-11,0.57926E-11,0.57013E-11,0.56173E-11,0.55402E-11, &
         0.54694E-11,0.54047E-11,0.53456E-11,0.52917E-11,0.52426E-11, &
         0.51981E-11,0.51578E-11,0.51214E-11,0.50886E-11,0.50591E-11, &
         0.50327E-11 /)


!##############################################################################
contains
!##############################################################################


   subroutine inisulchem()

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Initialize sulfur cycle chemistry module.
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      ! Local variables
      integer ::&
         i, l, ik, k

      !-----------------------------------------------------------------------
      !     initialize variables for jh2o2 interpolation
      !     interpolate jval data onto 1 km height grid
      do l = 1, nza*nalb
         do ik = 1, nz2-1
            k = 2
11          continue
            if( zz(k) .gt. ik-1 ) then
               jper2((l-1)*nz2+ik) = jper((l-1)*nz+k-1) &
                                     + ((ik-1)-zz(k-1))/(zz(k)-zz(k-1)) &
                                     * (jper((l-1)*nz+k)-jper((l-1)*nz+k-1))
            else
               k = k + 1
               go to 11
            endif
         end do
         jper2(l*nz2) = jper(l*nz) 
      end do
      
      do i = 1,nza1
         delang(i) = 1. / (zasec(i+1) - zasec(i))
      end do
      delalb  = 1. / (alb(2) - alb(1))
 
      !     Constants from radiation routines needed for dicor()
      dayspy  =  365.
      pie     =  4.*atan(1.)

   end subroutine inisulchem

!##############################################################################

   subroutine chemwdepdr( lchnk, nstep, lat, clat, lon, calday, dtime, &
                          pmid, pdel, zm, t, q, &
                          cldwat, cldf, concld, cldv, cmfdqr, &
                          prain, evapr, cme, rain, icwmr1, &
                          chtr, fracis, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Sulfur chemistry and wet deposition driver.  The prognostic species 
      ! are DMS, SO2, SO4, and H2O2.
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use acbnd, only: acbndget
      use history, only: outfld
      use wetdep, only: wetdepa, wetdepg
#ifdef SCYC_MASSBGT
      use massbgt, only: gamwet
#endif
      implicit none

      integer, intent(in) ::    lchnk
      integer, intent(in) ::    nstep
      integer, intent(in) ::    lat(pcols)                  ! latitude indices
      integer, intent(in) ::    lon(pcols)                  ! longitude indices
      integer, intent(in) ::    ncol

      real(r8), intent(in) ::&
         clat(pcols),              &
         calday,                   &
         dtime,                    &! time step
         pmid(pcols,pver),         &! pressure at layer midpoints
         pdel(pcols,pver),         &! delta-p between layer interfaces
         zm(pcols,pver),           &! height of layer midpoints
         t(pcols,pver),            &! air temperature (K)
         q(pcols,pver),            &! specific humidity (kg/kg)
         cldwat(pcols,pver),       &! cloud water mixing ratio
         cldf(pcols,pver),         &! total cloud fraction
         concld(pcols,pver),       &! convective cloud fraction
         cldv(pcols,pver),         &! cloudy volume which is occupied by rain or cloud water 
         cmfdqr(pcols,pver),       &! dq/dt due to moist convective rainout
         prain(pcols,pver),        &! rate of conversion of condensate to precip
         evapr(pcols,pver),        &! rate of evaporation of falling precip 
         cme(pcols,pver),          &! rate of cond-evap within the cloud
         rain(pcols,pver),         &! total precip mixing ratio
         icwmr1(pcols,pver)         ! in cloud water mixing ratio for zhang+hack scheme

!      real(r8) chtr(pcols,pver,4)   ! chemical tracer array
      real(r8), intent(inout) :: chtr(pcols,pver,4)   ! chemical tracer array
!        obuf(*)                    ! history output buffer

      real(r8), intent(out) ::&
!      real(r8) &
         fracis(pcols,pver,4)       ! fraction of transported species that are isouble 

      ! Local variables:
      integer :: i, k, m
      real(r8), dimension(pcols,pver) ::&
         so2,         &! so2 array
         so4,         &! so4 array
         dms,         &! dms
         h2o2,        &! h2o2 
         so2new,      &! so2 array
         so2tmp,      &! temporary so2 array
         so4new,      &! so4 array
         so4tmp,      &! temporary so4 array
         dmsnew,      &! dms
         dmssnk,      &! total sink of dms
         so2srcg,     &! total source of so2  
         so2srcg2,    &! dms + oh source of so2  
         so4src,      &! total source of so4  
         so4srcg,     &! gas source of so4  
         so4srca,     &! aq. source of so4  
         h2o2src,     &! total source of h2o2  
         h2o2snkg,    &! gas sink of h2o2  
         h2o2snka      ! aqueous sink of h2o2  

      real(r8), dimension(pcols,pver) ::&
         h2o2new,     &! h2o2
         h2o2tmp,     &! temporary h2o2
         h2o23d,      &! h2o2 -- from IMAGES
         oh,          &! oh -- from IMAGES
         no3,         &! no3 -- from IMAGES
         ho2im,       &! ho2 -- from IMAGES
         o3,          &! o3 -- from IMAGES
         ho2,         &! ho2 = ho2im*dicorfac
         ph,          &! ph of drops
         hion          ! Hydrogen ion concentration in drops

      real(r8) ::&
         dicorfac(pcols)  ! factors to increase HO2 diurnal averages

      integer ::&
         indcp(pcols,pver),    &! Indices of cloudy grid points
         ncldypts(pver)         ! Number of cloudy grid points

      real(r8), dimension(pcols,pver) ::&
         jh2o2,       &! photolysis rate of H2O2
         iscavt,      &
         scavt

      real(DBLKIND) :: &
         ekh2o2(pcols,pver),  &! H2O2 Henry's Law coefficient
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      real(r8) o3int(pcols), const

      ! The chemistry and wet dep schemes currently only use the in cloud water mixing
      ! ratio as a sum of icwmr1+icwmr2.  The value supplied by MATCH is that sum,
      ! and icwmr2 is set to zero.
      real(r8), dimension(pcols,pver) ::&
         icwmr2        ! in cloud water mixing ratio for hack scheme

      real(r8) :: molwt
      real(r8) :: sol_fact
!----------------------------------------------------------------------c

      do k = 1,pver
         do i = 1,ncol
            icwmr2(i,k) = 0.0
         end do
      end do

      do k = 1,pver
         do i = 1,ncol
            so2(i,k) = chtr(i,k,1)
            so4(i,k) = chtr(i,k,2)
            dms(i,k) = chtr(i,k,3)
            h2o2(i,k) = chtr(i,k,4)
         end do
      end do
      do m = 1, 4
         do k = 1, pver
            do i = 1, ncol
               fracis(i,k,m) = 1.
            end do
         end do
      end do

      call cldychmini( t, pmid, cldwat, rain, cldf, &
                       cldv, icwmr1, icwmr2, so2, so4, &
                       ph, hion, ekso2, ekh2o2, indcp, &
                       ncldypts, ncol )

      call outfld( 'PH', ph, pcols, lchnk )
      
                                                       ! Nominal concentrations
      call acbndget( lat, lon, 'H2O2', pcols, ncol, h2o23d )
      call acbndget( lat, lon, 'O3', pcols, ncol, o3 )
      call acbndget( lat, lon, 'HO2', pcols, ncol, ho2im )        ! 1.52e8
      call acbndget( lat, lon, 'NO3', pcols, ncol, no3   )        ! 1.64e7
      call acbndget( lat, lon, 'OH', pcols, ncol, oh     )        ! 1.e6

      ! alter HO2 from diurnal average to instantaneous value
      ! calculate integrate o3 (dobson units)
      call dicor( lat, calday, clat, ncol, dicorfac)
      const = 6.02e23 / (10.* 28.97 * 2.7e16*9.8) 
      o3int = 0.
      do k = 1, pver
         do i = 1, ncol
            ho2(i,k) = ho2im(i,k)*dicorfac(i)
            o3int(i) = o3int(i) + o3(i,k)*const*pdel(i,k)
         end do
      end do
      call outfld( 'O3INT', o3int, pcols, lchnk )
      call outfld( 'O3', o3, pcols, lchnk )  !"mol/mol in oxid_bnd.nc"

      call getj( clat, lat, lon, calday, zm, ncol, jh2o2 )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Toggle the chemistry and wet deposition routines !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if( mod( nstep/2, 2 ) .eq. 0 ) then

         call aqchem( so2, so4, h2o2, t, pmid, &
                      cldwat, rain, so2new, so4new, h2o2new, &
                      o3, cldf, so4srca, hion, h2o2snka, &
                      ekso2, ekh2o2, dtime, indcp, ncldypts, &
                      lat, lon, cldv, icwmr1, icwmr2, ncol )

         so2tmp = so2new
         so4tmp = so4new
         h2o2tmp = h2o2new
         call gaschem( so2tmp, so4tmp, dms, h2o2tmp, t, &
                       pmid, so2new, so4new, dmsnew, h2o2new, &
                       h2o23d, oh, no3, ho2, q, &
                       jh2o2, so4srcg, dmssnk, so2srcg, so2srcg2, &
                       h2o2src, h2o2snkg, dtime, lat, ncol )
#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'dmssnk', lat, pdel, dmssnk )
         call gamwet( 'so2srcg', lat, pdel, so2srcg )
         call gamwet( 'so4srca', lat, pdel, so4srca )
         call gamwet( 'so4srcg', lat, pdel, so4srcg )
#endif

         ! Wet deposition of SO2.
         molwt = 64._r8
         call wetdepg( lat, lon, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, prain, evapr, &
                       rain, cldwat, so2new, dtime, molwt, &
                       ekso2, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis(1,1,1), ncol )
         do k = 1,pver
            do i = 1,ncol
               so2new(i,k) = so2new(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'SO2WET', scavt, pcols, lchnk )
#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'so2wet', lat, pdel, scavt )
#endif

         ! Wet deposition of SO4.
         sol_fact = 0.3_r8
         call wetdepa( lat, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, icwmr1, prain, cme, &
                       evapr, cldwat, so4new, dtime, &
                       scavt, iscavt, cldv, fracis(1,1,2), sol_fact, ncol )
         do k = 1,pver
            do i = 1,ncol
               so4new(i,k) = so4new(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'SO4WET', scavt, pcols, lchnk )
#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'so4wet', lat, pdel, scavt )
#endif

         ! Wet deposition of H2O2.
         molwt = 34._r8
         call wetdepg( lat, lon, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, prain, evapr, &
                       rain, cldwat, h2o2new, dtime, molwt, &
                       ekh2o2, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis(1,1,4), ncol )
         do k = 1,pver
            do i = 1,ncol
               h2o2new(i,k) = h2o2new(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'H2O2WET', scavt, pcols, lchnk )

               !###############################
      else     !# toggle wetdep and chemistry #
               !###############################

         ! Wet deposition of SO2.
         molwt = 64._r8
         call wetdepg( lat, lon, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, prain, evapr, &
                       rain, cldwat, so2, dtime, molwt, &
                       ekso2, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis(1,1,1), ncol )
         do k = 1,pver
            do i = 1,ncol
               so2(i,k) = so2(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'SO2WET', scavt, pcols, lchnk )
#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'so2wet', lat, pdel, scavt )
#endif

         ! Wet deposition of SO4.
         sol_fact = 0.3_r8
         call wetdepa( lat, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, icwmr1, prain, cme, &
                       evapr, cldwat, so4, dtime, &
                       scavt, iscavt, cldv, fracis(1,1,2), sol_fact, ncol )
         do k = 1,pver
            do i = 1,ncol
               so4(i,k) = so4(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'SO4WET', scavt, pcols, lchnk )
#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'so4wet', lat, pdel, scavt )
#endif

         ! Wet deposition of H2O2.
         molwt = 34._r8
         call wetdepg( lat, lon, t, pmid, q, pdel, &
                       cldf, concld, cmfdqr, prain, evapr, &
                       rain, cldwat, h2o2, dtime, molwt, &
                       ekh2o2, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis(1,1,4), ncol )
         do k = 1,pver
            do i = 1,ncol
               h2o2(i,k) = h2o2(i,k) + scavt(i,k)*dtime
            end do
         end do
         call outfld( 'H2O2WET', scavt, pcols, lchnk )

         call aqchem( so2, so4, h2o2, t, pmid, &
                      cldwat, rain, so2new, so4new, h2o2new, &
                      o3, cldf, so4srca, hion, h2o2snka, &
                      ekso2, ekh2o2, dtime, indcp, ncldypts, &
                      lat, lon, cldv, icwmr1, icwmr2, ncol )

         so2tmp = so2new
         so4tmp = so4new
         h2o2tmp = h2o2new
         call gaschem( so2tmp, so4tmp, dms, h2o2tmp, t, &
                       pmid, so2new, so4new, dmsnew, h2o2new, &
                       h2o23d, oh, no3, ho2, q, &
                       jh2o2, so4srcg, dmssnk, so2srcg, so2srcg2, &
                       h2o2src, h2o2snkg, dtime, lat, ncol )

#ifdef SCYC_MASSBGT
         write (6,*) 'ERROR:  subroutine gamwet expects scalar for actual argument lat'
         call endrun
         call gamwet( 'dmssnk', lat, pdel, dmssnk )
         call gamwet( 'so2srcg', lat, pdel, so2srcg )
         call gamwet( 'so4srca', lat, pdel, so4srca )
         call gamwet( 'so4srcg', lat, pdel, so4srcg )
#endif

      endif     ! End of toggle

      do k = 1, pver
         do i = 1, ncol
            chtr(i,k,1) = so2new(i,k)
            chtr(i,k,2) = so4new(i,k)
            chtr(i,k,3) = dmsnew(i,k)
            chtr(i,k,4) = h2o2new(i,k)
            so4src(i,k) = so4srca(i,k) + so4srcg(i,k)
         end do
      end do
      call outfld( 'DMSSNK  ', dmssnk, pcols, lchnk )
      call outfld( 'SO2SRCG ', so2srcg, pcols, lchnk )
      call outfld( 'SO2SRCG2', so2srcg2, pcols, lchnk )
      call outfld( 'SO4SRC  ', so4src, pcols, lchnk )
      call outfld( 'SO4SRCG ', so4srcg, pcols, lchnk )
      call outfld( 'H2O2SRC ', h2o2src, pcols, lchnk )
      call outfld( 'H2O2SNKA', h2o2snka, pcols, lchnk )
      call outfld( 'H2O2SNKG', h2o2snkg, pcols, lchnk )

   end subroutine chemwdepdr

!##############################################################################

   subroutine aqchem( so2, so4, h2o2, tair, pres, &
                      qc, qr, so2new, so4new, h2o2new, &
                      o3, cldf, so4srca, hion, h2o2snka, &
                      ekso2, ekh2o2, ztodt, ind, ncpts, &
                      lat, lon, cldv, icwmr1, icwmr2, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the sources and sinks from the aqueous chemistry:
      !          S(IV) + H2O2 --> SO4 
      !          S(IV) + O3   --> SO4 
      !
      !  where concentrations are concentrations in aqueous phase.
      !
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         lat(pcols),          &! latitude indices
         lon(pcols),          &! longititude indices
         ind(pcols,pver),     &! indices for cloudy grid points
         ncpts(pver),         &! number of cloudy grid points
         ncol                  ! number of atmospheric columns used in chunk
      real(r8), intent(in) ::&
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 from aqueous chemistry 
         h2o2(pcols,pver),    &! h2o2 in kg/kg
         o3(pcols,pver),      &! o3
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         qc(pcols,pver),      &! cloud water mixing ratio
         qr(pcols,pver),      &! rain mixing ratio
         ztodt,               &! time step
         cldf(pcols,pver)     ! fraction cloud cover

      real(r8), intent(in) ::&
         cldv(pcols,pver),    &! estimate of local volume occupied by clouds
         icwmr1(pcols,pver),  &! in cloud water mixing ration for zhang scheme
         icwmr2(pcols,pver)    ! in cloud water mixing ration for hack  scheme

      real(r8), intent(inout) ::&
         hion(pcols,pver)      ! hydrogen ion concentration (mol/L)

       real(DBLKIND)           &
         theff(pcols),        &! ekso2 at cloudy grid points
         sheff,               &! temp array for theff
         ekh2o2(pcols,pver),  &! H2O2  Henry's Law coefficient
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      real(r8), intent(out) ::&
         so2new(pcols,pver),  &! so2 array
         so4new(pcols,pver),  &! so4 from aqueous chemistry 
         h2o2new(pcols,pver), &! h2o2 from aqueous chemistry 
         so4srca(pcols,pver), &! so4 production rate from achem
         h2o2snka(pcols,pver)  ! aqueous sink of H2O2

      ! Local variables:

      integer ::&
         i, k, n,    &
         i2,         &! indice for cloudy grid points
         ilw,        &! number of cloudy grid points
         ncyc         ! number of times to loop through aqchem

      real(r8), dimension(pcols,pver) ::&
         cwat,       &! cloud water
         rwat,       &! rain
         twat         ! total liquid water = cwat + rwat
      real(r8) ::&
         h2o2rate,   &! rate of S(IV) + H2O2 reaction 
         o3rate,     &! first order rate coefficient for S(IV) + O3
         molarso4,   &! total so4 in units mol SO4/L H2O
         ahso3        ! (molar)**2 concentration of HSO3-

      real(r8), dimension(pcols) ::&
         tso2,       &! so2 at cloudy grid points
         tso4,       &! so4 at cloudy grid points
         th2o2,      &! h2o2 at cloudy grid points
         to3,        &! o3 at cloudy grid points
         ttwat,      &! twat at cloudy grid points
         temp,       &! tair at cloudy grid points
         tprs,       &! pres at cloudy grid points
         tso2n,      &! so2new at cloudy grid points
         tso4n,      &! so4new at cloudy grid points
         th2o2n,     &! h2o2new at cloudy grid points
         tso4src,    &! so4srca at cloudy grid points
         th2o2snk,   &! h2o2snk at cloudy grid points
         thion,      &! hion at cloudy grid points
         tkh2o2       ! ekh2o2 at cloudy grid points

      real(r8) ::&
         molso2,     &! total so2 in units mol SO2/mol air
         molh2o2,    &! total h2o2 in units mol H2O2/mol air
         molarh2o2,  &! total h2o2 in units mol H2O2/L H2O
         molo3,      &! total o3 in units mol O3/mol air
         eko3,       &! O3   Henry's Law coefficient
         so4src,     &! source of SO4 from SO2 + H2O2
         so2snk,     &! sink of SO2 from SO2 + H2O2
         h2o2snk,    &! sink of H2O2 from SO2 + H2O2
         siv,        &! concentration of S(IV) (mol/L)
         a1, a2,     &! fraction of S(IV) that is HSO3, SO3
         totsb,      &! sum of sulfur at beg of aqchem for one point
         totse,      &! sum of sulfur at end of aqchem for one point
         weight,     &
         tc           ! temperature in deg C

      real(r8) ::&
         dtchem,     &! chemistry timestep
         patm,       &! pressure (atm)
         fa1,        &!  = Khs*K1s/[H+]
         fa2,        &!  = Khs*K1s*K2s/[H+]**2
         khs,        &! Henry's Law constant for SO2
         khsk1s,     &! Khs * K1s
         kh12s,      &! Khs * K1s * K2s
         kh2o2,      &! rate constant for H2O2 + HSO3-
         ko3hs,      &! rate constant for O3 + HSO3-
         ko3so3,     &! rate constant for O3 + SO3=
         qso2,       &! temp array for so2
         qso4,       &! temp array for so4
         qh2o2,      &! temp array for h2o2
         ahion,      &! temp array for ahion
         molrate,    &! rate of rxn in mol/mol-s
         tmprate,    &! temp rate of rxn
         sumrate,    &! h2o2rate + o3rate
         omsm

      ! Function:
!-drb       real(r8) e298, dhr, tt, ek
      real(r8) e298, dhr, tt
      real (DBLKIND) ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1./tt - 1./298.))
      !----------------------------------------------------------------------

      ncyc = 20                  !perform aqchem on smaller dt
      dtchem = ztodt/ncyc
      omsm = 1. - 1.e-10

      do k=1,pver
         do i=1,ncol

            ! Determine cloud water and rain water mixing ratios
            tc     = tair(i,k) - 273.16
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cwat(i,k) = (qc(i,k)/max(cldf(i,k), 0.01_r8) + &
            ! add suspended water from convection and do aqchem in only the liquid phase
                        icwmr1(i,k) + icwmr2(i,k)) * (1.-weight)

            if(tair(i,k) .gt. 273.16) then
               weight = 0.                    ! fraction of condensate that is ice
            else
               weight = 1.
            endif

            rwat(i,k) = qr(i,k)/max(cldv(i,k), 0.01_r8) &
                        *(1.-weight) ! aqchem in only the liquid phase

            ! Sum cwat and rwat 
            twat(i,k) = cwat(i,k) + rwat(i,k)

            !        if(twat(i,k) .gt. 0.5e-3) then
            !         tc     = tair(i,k) - 273.16
            !         weight = max(0._r8,min(-tc*0.05,1.0_r8)) 
            !         write(6,'(a,3i4,5e12.5)') 'TWAT',i,lat,k, twat(i,k), cwat(i,k),
            !     _    rwat(i,k), icwmr1(i,k)*(1.-weight), icwmr2(i,k)*(1.-weight)
            !         write (6,*) ' rwat stuff, qr, cldv, cldf ', qr(i,k), cldv(i,k),
            !     $        cldf(i,k)
            !        endif

         end do
      end do

      ! Initialize arrays
      do k=1,pver
         do i=1,ncol
            so2new(i,k) = so2(i,k)       !keep concentration
            so4new(i,k) = so4(i,k)       
            h2o2new(i,k) = h2o2(i,k)       
            so4srca(i,k) = 0.
            h2o2snka(i,k) = 0.
         end do
      end do

      ! Initialize cloudy grid point arrays
      !      write(6,*) ' Beginning of aqchem'
      do k=1,pver
         ilw = ncpts(k)

         do i2=1,ilw
            tso2(i2)   = so2(ind(i2,k),k)
            tso4(i2)   = so4(ind(i2,k),k)
            th2o2(i2)  = h2o2(ind(i2,k),k)
            to3(i2)    = o3(ind(i2,k),k)
            ttwat(i2)  = twat(ind(i2,k),k)
            temp(i2)   = tair(ind(i2,k),k)
            tprs(i2)   = pres(ind(i2,k),k)
            thion(i2)  = hion(ind(i2,k),k)
            theff(i2)  = ekso2(ind(i2,k),k)
            tkh2o2(i2) = ekh2o2(ind(i2,k),k)
            tso4src(i2) = 0.
            th2o2snk(i2) = 0.
            !        write(6,'(3i5,4e12.5,f8.3,a)') ind(i2,k), lat, k, ttwat(i2), &
            !        thion(i2), tso4(i2), tso2(i2), -log10(thion(i2)), '  begin'
         end do
                    
         ! Aqueous chemistry calculations begin here
         do i2=1,ilw
            patm = tprs(i2)/101325.
            khs    = ek(1.23_r8,3120._r8,temp(i2))
            khsk1s = khs * ek(1.3e-2_r8,2015._r8,temp(i2))
            kh12s = khsk1s * ek(6.31e-8_r8,1505._r8,temp(i2)) 
            kh2o2  = 6.2534e14*exp(-4751./temp(i2)) 
            molo3 = to3(i2) 
            eko3  = ek(1.15e-2_r8,2560._r8,temp(i2))
            ko3hs = 4.28e13*exp(-5533./temp(i2)) 
            ko3so3= 7.436e16*exp(-5280./temp(i2)) 

            !        if(lat .eq. 8 .and. k .eq. 14 .and. i2 .gt. 8) then
            !         write(6,'(i3,7e11.4)') i2, tso2(i2), th2o2(i2), thion(i2), &
            !          theff(i2), temp(i2), patm, molo3
            !         write(6,'(3x,e11.4,f8.1)') ttwat(i2), dtchem
            !        endif
            do n=1,ncyc                     !Do aqchem on a smaller dt   5/5/97

               if(n .eq. 1) then
                  qso2 = tso2(i2)
                  qso4 = tso4(i2)
                  qh2o2 = th2o2(i2)
                  ahion = thion(i2)
                  sheff = theff(i2)
               endif

               fa1 = khsk1s/max(ahion,1.e-30_r8)
               fa2 = kh12s/max((ahion*ahion),1.e-30_r8)

               ! Calculate rate of HSO3 + H2O2 reaction
               molso2  = qso2  * mwa2so2
               molh2o2 = qh2o2 * mwa2h2o2
               ahso3 = fa1 * patm*molso2
               molarh2o2= tkh2o2(i2)*patm*molh2o2

               h2o2rate = kh2o2 * ahion * ahso3 * molarh2o2/(1. + 13.*ahion)  

               ! Convert to mol/(mol-s) 
               molrate = h2o2rate*ttwat(i2)*28.97e-3
               ! check for conservation in kg/(kg-s)
               tmprate = mwa2h2o2*min(molrate/mwa2h2o2, qh2o2/dtchem)
               tmprate = mwa2so2*min(tmprate/mwa2so2, qso2/dtchem)

               h2o2rate = tmprate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3)

               ! Calculate rate of o3 oxidation
               siv = sheff * molso2 * patm         ![S(IV)]
               a1 = fa1 / sheff
               a2 = fa2 / sheff

               o3rate  = (ko3hs*a1 + ko3so3*a2) * eko3*molo3 * patm*siv

               ! Convert to mol/(mol-s) 
               molrate = o3rate*ttwat(i2)*28.97e-3
               ! check for conservation in kg/(kg-s)
               tmprate = mwa2so2*min(molrate/mwa2so2, qso2/dtchem)

               o3rate = tmprate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3)

               ! check total rate for conservation
               sumrate = h2o2rate + o3rate
               !         if(lat .eq. 8 .and. k .eq. 14 .and. i2 .gt. 8) then
               !          write(6,'(i3,5e11.4)') n, sumrate, h2o2rate, o3rate, &
               !            qso2/dtchem, qh2o2/dtchem
               !         endif
               molrate = sumrate*ttwat(i2)*28.97e-3
               molrate = mwa2so2*min(molrate/mwa2so2, qso2/dtchem)
               tmprate = molrate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3)
               h2o2rate = h2o2rate/max(sumrate, 1.e-30_r8) * tmprate
               o3rate = o3rate/max(sumrate, 1.e-30_r8) * tmprate

               ! Update so2, so4, h2o2
               so2snk = (h2o2rate + o3rate)*ttwat(i2)*1.e-3 * 64.
               h2o2snk = h2o2rate*ttwat(i2)*1.e-3 * 34.
               !         if(qh2o2 .lt. omsm*h2o2snk*dtchem) then
               !          write(6,'(a,2e12.5,3i4)') 'AQCHEM qh2o2 < h2o2snk*dt ', qh2o2, &
               !           h2o2snk*dtchem, lat, ind(i2,k), k
               !         endif
               !         if(qso2 .lt. omsm*so2snk*dtchem) then
               !          write(6,'(a,2e12.5,3i4)') 'AQCHEM qso2 < so2snk*dt ', qso2, &
               !           so2snk*dtchem, lat, ind(i2,k), k
               !         endif
               so4src = 1.5*so2snk

               qso2 = max(qso2 - so2snk*dtchem, 0._r8)
               qso4 = qso4 + so4src*dtchem
               qh2o2 = max(qh2o2 - h2o2snk*dtchem, 0._r8)

               tso4src(i2) = tso4src(i2) + so4src
               th2o2snk(i2) = th2o2snk(i2) + h2o2snk

               ! calculate new ahion, sheff
               molso2 = qso2 * mwa2so2
               molarso4 = qso4/(max(ttwat(i2), 1.e-30_r8)*96.e-3) 
               ahso3 = khsk1s * molso2 * patm
!++pjr
! rewrote this expression so we dont take the sqrt(max(a,1.e-60))
!                               but rather max(sqrt(a),1.e-30)
!               ahion = (molarso4 + sqrt(max(molarso4*molarso4 + 4.*ahso3, &
!                       1.e-60_r8)))*0.5
!--pjr
               ahion = (molarso4 + sqrt(max(molarso4*molarso4 + 4._r8*ahso3, &
                       0._r8)))*0.5
               ahion = max(ahion,1.e-30_r8)
               sheff = khs + khsk1s/ahion + kh12s/max((ahion*ahion),1.e-30_r8)
            end do    !ncyc

            tso2n(i2) = qso2
            tso4n(i2) = qso4
            th2o2n(i2)= qh2o2
            thion(i2) = ahion
            theff(i2) = sheff
            tso4src(i2) = tso4src(i2)/ncyc
            th2o2snk(i2) = th2o2snk(i2)/ncyc

         end do


         do i2=1,ilw
            ! so2new is modified for just the cloudy portion of the grid cell.  Now
            !  adjust it for the entire grid cell.
            so2new(ind(i2,k),k) = cldv(ind(i2,k),k)*tso2n(i2) + &
                                  (1.-cldv(ind(i2,k),k))*tso2(i2)
            so4new(ind(i2,k),k) = cldv(ind(i2,k),k)*tso4n(i2) + &
                                  (1.-cldv(ind(i2,k),k))*tso4(i2)
            so4srca(ind(i2,k),k) = cldv(ind(i2,k),k)*tso4src(i2) 
            h2o2snka(ind(i2,k),k)= cldv(ind(i2,k),k)*th2o2snk(i2) 
            h2o2new(ind(i2,k),k) = cldv(ind(i2,k),k)*th2o2n(i2) + &
                                   (1.-cldv(ind(i2,k),k))*th2o2(i2)
            hion(ind(i2,k),k) = thion(i2)  
            ekso2(ind(i2,k),k) = theff(i2)  
            !         write(6,'(3i5,4e12.5,f8.3,a)') ind(i2,k), lat, k, ttwat(i2), &
            !         thion(i2), tso4(i2), tso2(i2), -log10(thion(i2)), '  end'
         end do


         do i=1,ncol
            totsb = 0.5*so2(i,k) + 1./3. * so4(i,k)
            totse = 0.5*so2new(i,k) + 1./3. * so4new(i,k)
            !        if(abs(totsb-totse) .gt. .001*totsb) then
            !         write(6,*) 'Total sulfur not conserved in aqueous chemistry'
            !         write(6,*) i, k, totsb, totse
            !         stop
            !        endif

         end do

      end do

   end subroutine aqchem

!##############################################################################

   subroutine cldychmini( tair, pres, qc, qr, cldf, &
                          cldv, icwmr1, icwmr2, so2, so4, &
                          ph, hion, ekso2, ekh2o2, ind, &
                          ncpts, ncol ) 

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the pH, solubility constants of SO2 and H2O2, and the
      !  location of the cloudy grid points
      ! 
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         qc(pcols,pver),      &! cloud water mixing ratio
         qr(pcols,pver),      &! rain mixing ratio
         cldf(pcols,pver),   &! fraction cloud cover
         cldv(pcols,pver),    &! cloudy volume which is occupied by rain or cloud water
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 array
         icwmr1(pcols,pver),  &!  in cloud water mixing ration for zhang scheme
         icwmr2(pcols,pver)    !  in cloud water mixing ration for hack  scheme

      real(r8), intent(out) ::&
         ph(pcols,pver),      &! pH of drops
         hion(pcols,pver)      ! hydrogen ion concentration (mol/L)

      real(DBLKIND)           &
         theff(pcols),        &! ekso2 at cloudy grid points
         ekh2o2(pcols,pver),  &! H2O2 Henry's Law coefficient
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      integer, intent(in) ::  &
         ncol                  ! number of atmospheric columns used in chunk
      integer, intent(out) ::&
         ind(pcols,pver),     &! indices for cloudy grid points
         ncpts(pver)           ! number of cloudy grid points

      ! Local variables:

      integer :: &
         i, k,                &
         i2,                  &! index for cloudy grid points
         ilw                   ! number of cloudy grid points

      real(r8), dimension(pcols,pver) ::&
         cwat,       &! cloud water
         rwat,       &! rain
         twat         ! total liquid water = cwat + rwat

      real(r8), dimension(pcols) ::&
         tso2,       &! so2 at cloudy grid points
         tso4,       &! so4 at cloudy grid points
         ttwat,      &! twat at cloudy grid points
         temp,       &! tair at cloudy grid points
         tprs,       &! pres at cloudy grid points
         tph,        &! ph at cloudy grid points
         thion        ! hion at cloudy grid points
      real(r8)  ::   &
         molarso4,   &! total so4 in units mol SO4/L H2O
         ahso3,      &! (molar)**2 concentration of HSO3-
         molso2,     &! so2 in mol/mol units
         weight,     &
         tc           ! temperature in deg C

      ! Function:
!-drb      real(r8) e298, dhr, tt, ek
      real(r8) e298, dhr, tt
      real(DBLKIND) ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1./tt - 1./298.))
      !----------------------------------------------------------------------

      do k=1,pver
         do i=1,ncol

            ekh2o2(i,k) = ek(7.4e4_r8, 6621._r8, tair(i,k))
            ekso2(i,k)  = ek(1.23_r8,  3120._r8, tair(i,k)) 
            ! Determine cloud water and rain water mixing ratios
            tc     = tair(i,k) - 273.16
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cwat(i,k) = (qc(i,k)/max(cldf(i,k), 0.00001_r8) + &
               ! add suspended water from convection and do aqchem in only the liquid phase
                        icwmr1(i,k) + icwmr2(i,k)) * (1.-weight)

            if(tair(i,k) .gt. 273.16) then
               weight = 0.                    ! fraction of condensate that is ice
            else
               weight = 1.
            endif

            rwat(i,k) = qr(i,k)/max(cldv(i,k), 0.00001_r8) &
                        *(1.-weight) ! aqchem in only the liquid phase

            ! Sum cwat and rwat 
            twat(i,k) = cwat(i,k) + rwat(i,k)

         end do
      end do

      ! Initialize arrays
      do k=1,pver
         do i=1,ncol
            ph(i,k) = 99.
            hion(i,k) = 0.
            ind(i,k) = 0
         end do
      end do

      ! Find which grid points have liquid water
      do k=1,pver
         ncpts(k) = 0
         ilw = 0
         do i=1,ncol
            if(cldv(i,k) .ge. 0.02 .and. twat(i,k) .ge. 1.e-12) then
               ilw = ilw + 1
               ind(ilw,k) = i
            endif
         end do

         ncpts(k) = ilw       !assign number of cloudy grid points to array

         do i2=1,ilw
            tso2(i2) = so2(ind(i2,k),k)
            tso4(i2) = so4(ind(i2,k),k)
            ttwat(i2) = twat(ind(i2,k),k)
            temp(i2) = tair(ind(i2,k),k)
            tprs(i2) = pres(ind(i2,k),k)

            ! Set pH as a function of HSO3- and SO4=
            !  Say NH4+ = 1.5 SO4= and NO3- = 0.5 SO4=
            !  Thus, H+ = HSO3- + SO4=  as done in ECHAM
            molso2 = tso2(i2)  * mwa2so2

            molarso4 = tso4(i2)/(max(ttwat(i2), 1.e-30_r8)*96.e-3) 
            ahso3 = ek(1.23_r8,3120._r8,temp(i2)) * &
                    ek(1.3e-2_r8,2015._r8,temp(i2)) * &
                    molso2 * tprs(i2)/101325.
!++pjr
! rewrote this expression so we dont take the sqrt(max(a,1.e-50))
!                               but rather max(sqrt(a),1.e-25)
!            thion(i2) = (molarso4 + sqrt(max(molarso4*molarso4 + 4.*ahso3, &
!                         1.e-50_r8)))*0.5
!--pjr
            thion(i2) = (molarso4 + &
                         sqrt(max(molarso4*molarso4 + 4._r8*ahso3,0._r8)))*0.5
            thion(i2) = max(thion(i2),1.e-25_r8)
            tph(i2) = -log10(thion(i2))

            theff(i2) = ek(1.23_r8, 3120._r8, temp(i2)) * (1. + &
                        ek(1.3e-2_r8, 2015._r8, temp(i2))/thion(i2) * (1. + &
                        ek(6.31e-8_r8, 1505._r8, temp(i2))/thion(i2)))

         end do


         do i2=1,ilw
            hion(ind(i2,k),k) = thion(i2)  
            ph(ind(i2,k),k) = tph(i2)  
            ekso2(ind(i2,k),k) = theff(i2)

         end do
      end do

   end subroutine cldychmini


!##############################################################################

   subroutine coszen( calday, dodiavg, clat, lat, lon, ncol, coszrs )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Compute solar distance factor and cosine solar zenith angle usi
      ! day value where a round day (such as 213.0) refers to 0z at
      ! Greenwich longitude.
      !
      ! Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
      ! Processes in Meterology and Climatology, Elsevier Scientific
      ! Publishing Company, New York  p. 57, p. 62,63.
      ! 
      ! Author: CCM2
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         calday,              &! Calendar day, including fraction
         clat(pcols)           ! latitudes (radians)

      integer, intent(in) ::&
         lat(pcols),          &! latitude indices
         lon(pcols),          &! longitude indices
         ncol                  ! number of atmospheric columns used in chunk

      logical, intent(in) ::&
         dodiavg               ! true => do diurnal averaging

      real(r8), intent(out) ::&
         coszrs(pcols)        ! Cosine solar zenith angle

      ! Local variables
      integer :: i     ! Longitude loop index
      real(r8) ::&
         phi,     &! Greenwich calendar day + local time + long offset
         theta,   &! Earth orbit seasonal angle in radians
         delta,   &! Solar declination angle  in radians
         sinc,    &! Sine   of latitude
         cosc,    &! Cosine of latitude
         sind,    &! Sine   of declination
         cosd      ! Cosine of declination
      real(r8) ::&
         frac,    &! Daylight fraction
         arg,     &!
         tsun,    &! temporary term in diurnal averaging
         coszrsu   ! uniform cosine zenith solar angle 
      !-----------------------------------------------------------------------

      theta = 2.*pie*calday/dayspy

      ! Solar declination in radians:
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) - &
              .006758*cos(2.*theta) + .000907*sin(2.*theta) - &
              .002697*cos(3.*theta) + .001480*sin(3.*theta)

      do i=1,ncol
         ! Compute local cosine solar zenith angle,
         sinc = sin(clat(i))
         sind = sin(delta)
         cosc = cos(clat(i))
         cosd = cos(delta)
   
         ! If using diurnal averaging, then compute the average local cosine solar 
         ! zenith angle using formulas from paltridge and platt 1976  p. 57, p. 62,63.
         if (dodiavg) then
            arg = -(sinc/cosc)*(sind/cosd)
            if (arg .lt. -1.) then
               frac = 1.0
            else if (arg .gt. 1.) then
               frac = 0.0
            else
               frac = (1./pie)*acos(arg)
            endif
            tsun = pie*frac
            if (tsun .gt. 0.) then
               coszrsu =  sinc*sind + (cosc*cosd*sin(tsun))/tsun
            else
               coszrsu = 0.0
            endif
            coszrs(i) = coszrsu
         else                       ! No diurnal averaging
   
            ! Calday is the calender day for Greenwich, including fraction
            ! of day; the fraction of the day represents a local time at
            ! Greenwich; to adjust this to produce a true instantaneous time
            ! For other longitudes, we must correct for the local time change:
            ! local time based on the longitude and day of year
            ! then compute the local cosine solar zenith angle
            phi       = calday + (real(lon(i)-1)/real(plon))
            coszrs(i) = sinc*sind - cosc*cosd*cos(2.*pie*phi)
         end if
      end do

   end subroutine coszen

!##############################################################################

   subroutine dicor( lat, calday, clat, ncol, corr )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate a correction factor for the ho2 reaction rate to account for
      ! the diurnal variation of ho2
      !
      ! Author: Mary Barth
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) :: lat(pcols)             ! latitude indices

      real(r8), intent(in) :: calday                ! Calendar day, including fraction
      real(r8), intent(in) :: clat(pcols)           ! latitudes (radians)

      integer, intent(in) ::  ncol                  ! number of atmospheric columns used in chunk

      real(r8), intent(out) :: corr(pcols)          ! correction factor

      !---------------------------Local variables-----------------------------

      integer ic        ! index in chunk
      integer i         ! Longitude loop index
      real(r8) phi      ! Greenwich calendar day + local time + long offset
      real(r8) theta    ! Earth orbit seasonal angle in radians
      real(r8) delta    ! Solar declination angle  in radians
      real(r8) sinc     ! Sine   of latitude
      real(r8) cosc     ! Cosine of latitude
      real(r8) sind     ! Sine   of declination
      real(r8) cosd     ! Cosine of declination
      real(r8) coszrs       ! Cosine solar zenith angle
      real(r8) xxx      ! Work space
      real(r8) yyy      ! Work space
      integer nc        ! counter

      real(r8) ::     dicor_corr          ! h2o reaction rate correction
      integer ntimes    ! number of points to evaluate the diurnal average

      !-----------------------------------------------------------------------
      !

      do ic=1,ncol

         ! Compute solar distance factor and cosine solar zenith angle using
         ! day value where a round day (such as 213.0) refers to 0z at
         ! Greenwich longitude.
         !
         ! Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
         ! Processes in Meterology and Climatology, Elsevier Scientific
         ! Publishing Company, New York  p. 57, p. 62,63.

         theta = 2.*pie*calday/dayspy
         !
         ! Solar declination in radians:
         !
         delta = .006918 - .399912*cos(theta) + .070257*sin(theta) - &
              .006758*cos(2.*theta) + .000907*sin(2.*theta) - &
              .002697*cos(3.*theta) + .001480*sin(3.*theta)

         ! now calculate the solar zenith angle 
         ! and some useful quantities for the whole day
         !
         ! Compute local cosine solar zenith angle,
         !
         sinc = sin(clat(ic))
         sind = sin(delta)
         cosc = cos(clat(ic))
         cosd = cos(delta)
         !
         ! Calday is the calender day for Greenwich, including fraction
         ! of day; 
         ! since all we care about in this calculation is the diurnal variation
         ! of some quantities we just increment calday over one day.
         !
         xxx = 0.
         yyy = 0.
         nc  = 0
         ntimes = 48
         do i=1,ntimes
            phi       = calday + (real(i)-1)/real(ntimes)
            coszrs = sinc*sind - cosc*cosd*cos(2.*pie*phi)
            if (coszrs.gt.0) then
               nc  = nc + 1
               xxx = xxx + coszrs**2
               yyy = yyy + coszrs
            endif
         end do

         if (yyy.gt.0.) then
            dicor_corr = ntimes*xxx/(yyy**2)
         else
            dicor_corr = 1.
         endif


         !
         ! Return correction
         !
         corr(ic) = dicor_corr
      end do

    end subroutine dicor

!##############################################################################

   subroutine gaschem( so2, so4, dms, h2o2, tair, &
                       pres, so2new, so4new, dmsnew, h2o2new, &
                       h2o23d, oh, no3, ho2, h2o, &
                       jh2o2, so4srcg, dmssnk, so2srcg, so2srcg2, &
                       h2o2src, h2o2snkg, ztodt, lat, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the sources and sinks from the gas chemistry:
      !          DMS + OH     --> a*SO2 + (1-a)*MSA
      !          DMS + NO3    --> SO2 
      !          SO2 + OH + M --> SO4 + M
      !          HO2 + HO2    --> H2O2    
      !          H2O2+ OH     --> HO2 + H2O   
      !          H2O2+ hv     --> 2OH
      !
      !  where a is the yield of SO2 from DMS oxidation, MSA is methane sulfonic
      !   acid which is not carried in this model, and M is an air molecule and
      !   is accounted for in the rate coefficient expression.
      !
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) :: &
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 from gas chemistry 
         dms(pcols,pver),     &! dms
         h2o2(pcols,pver),    &! h2o2
         oh(pcols,pver),      &! oh
         no3(pcols,pver),     &! no3
         ho2(pcols,pver),     &! ho2
         h2o(pcols,pver),     &! h2o
         jh2o2(pcols,pver),   &! photolysis rate of H2O2
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         ztodt,               &! time step
         h2o23d(pcols,pver)    ! h2o2 -- 3d field from images

      integer, intent(in) ::  &
         lat(pcols),          &! latitude indices
         ncol                  ! number of atmospheric columns used in chunk

      real(r8), intent(out) ::  &
         so4new(pcols,pver),    &! so4 from gas chemistry 
         dmsnew(pcols,pver),    &! dms
         h2o2new(pcols,pver),   &! h2o2
         so4srcg(pcols,pver),   &
         so2srcg2(pcols,pver)

      real(r8), intent(inout) ::  &
         so2new(pcols,pver),      &! so2 array
         dmssnk(pcols,pver),      &
         so2srcg(pcols,pver),     &
         h2o2src(pcols,pver),     &
         h2o2snkg(pcols,pver)

      ! Local variables:

      integer ::&
         i, k, ik
      real(r8), parameter ::&         !  from NASA(1992)
         rl300=3.e-31,    &! rate constant at low pressure
         an=3.3,          &! constant:  (300/Tair)**an
         rh300=1.5e-12,   &! rate constant at high pressure
         am=0.             ! constant:  (300/Tair)**am

      real(r8) ::&
         air,             &! concentration of air (molec/cm3)
         rlt,             &! used in calculation of SO2 rate coefficient
         rht,             &! used in calculation of SO2 rate coefficient
         ax, bx,          &! used in calculation of SO2 rate coefficient
         psi,             &! used in calculation of SO2 rate coefficient
         so2rate,         &! first order rate coefficient for SO2 + OH
         so2snk,          &
         dmsrate1,        &! first order rate coefficient for DMS + OH
         dmsno3,          &! first order rate coefficient for DMS + NO3
         dmstot            ! first order rate coefficient for both rxns
      real(r8) rk2         ! rate coeff. in HO2 + HO2 rxn
      real(r8) rk3         ! rate coeff. in HO2 + HO2 rxn
      real(r8) fh2o        ! rate coeff. in HO2 + HO2 rxn
      real(r8) ho2rate     ! total rate coeff for HO2 + HO2 rxn
      real(r8) h2o2des     ! rate for H2O2*OH
      real(r8) h2o2phot    ! rate jh2o2*H2O2
      !      real(r8) tauh2o2         !relaxation time for H2O2 generation
      !      real(r8) xh2o2           !zonally avgd H2O2 in kg/kg
      !----------------------------------------------------------------------c
      
      !  Assume a 1.5 day relaxation time for H2O2 generation
      !      tauh2o2 = 129600.    !seconds

      
      do k=1,pver
         do i=1,ncol

            so2rate = 0.

            ! determine concentration of air (molec/cm3) = rho(kg/m3)*Na/MWair
            air = pres(i,k)/(287.*tair(i,k)) * 1.e-3/28.966 * 6.022e23
            !         air = 0.8*air                   !N2 concentration
 
            ! SO2 + OH rate coefficient:
            if(tair(i,k) .le. 0.) then
               write(6,*) 'TAIR bad ', i, k, lat(i), tair(i,k)
               call endrun
            endif
            rlt = rl300 * (300./tair(i,k))**an
            rht = rh300 * (300./tair(i,k))**am
            ax = rlt*air/rht 
            psi = 1./(1. + (log10(ax))**2)
            bx = rlt*air/(1. + ax)
            so2rate = bx * 0.6**psi             !rate coef. in cm3/molec-s
            ! rewrite so2rate as first order rate coefficient
            so2rate = so2rate * oh(i,k)

            so2snk = so2(i,k)*(1.-exp(-so2rate*ztodt))
            ! Update sulfur species
            so2new(i,k) = so2(i,k)  - so2snk

            so4new(i,k) = so4(i,k) + 1.5*so2snk
            so4srcg(i,k) = 1.5*so2snk/ztodt

            !   Note: 1.5 factor accounts for different molecular weights of so4 and so2
            !    mw(so4)=96;  mw(so2)=64;  mw(so4)/mw(so2)=3/2=1.5
         end do
      end do

      do k=1,pver
         do i=1,ncol

            ! DMS oxidation
            !c DMS + OH rate coefficient and SO2 yield:
            dmstot = 0.
            dmsrate1 = 0.
            dmsno3 = 0.
            ! Get dmsrate from tabulated data
            ik = int(tair(i,k)-195.)
            ik = min(max(ik,1), 121)

            if(dms(i,k) .ge. 1.e-30) then
               dmsrate1 = dmsrate(ik)*oh(i,k)

               ! DMS + NO3
               dmsno3 = 1.e-12*exp(-500.*(0.0033557 - 1./tair(i,k))) * &
                        no3(i,k)

               ! Total first order rate coeff for DMS
               dmstot = dmsrate1 + dmsno3
            endif

            ! Update sulfur species
            dmssnk(i,k) = dms(i,k)*(1.-exp(-dmstot*ztodt)) / ztodt
            dmsnew(i,k) = dms(i,k) - dmssnk(i,k)*ztodt
            so2srcg(i,k) = 1.032258 * (yield(ik) * dms(i,k)*(1.-exp(-dmsrate1*ztodt)) + &
                           dms(i,k)*(1.-exp(-dmsno3*ztodt)) ) / ztodt
            so2new(i,k) = so2new(i,k) + so2srcg(i,k)*ztodt
            so2srcg2(i,k) = 1.032258 / ztodt * &
                            yield(ik) * dms(i,k)*(1.-exp(-dmsrate1*ztodt)) 
            !   Note: 1.032258 = 64/62 accounts for different MW of so2 and dms
         end do
      end do

      do k=1,pver
         do i=1,ncol
            !-----------------------------------------------------------------------
            !c Generate H2O2 via gas phase processes -- parameterizing chemistry
            !        xh2o2 = h2o23d(i,k) * (287.*tair(i,k)/pres(i,k))*34./6.022e20
            !        h2o2new(i,k) = h2o2(i,k) + ztodt/tauh2o2 * (xh2o2 - h2o2(i,k))
            !
            ! Generate H2O2 via the HO2 + HO2 reaction
            !  HO2 concentrations are prescribed
            !  water vapor concentrations taken from q(i,k,1)
            !
            ! determine concentration of air (molec/cm3) = rho(kg/m3)*Na/MWair
            air = pres(i,k)/(287.*tair(i,k)) * 1.e-3/28.966 * 6.022e23
            rk2 = 1.7e-12*exp( 600.*(1./tair(i,k) - 0.0033557))
            rk3 = 4.9e-32*exp(1000.*(1./tair(i,k) - 0.0033557))
            fh2o = 1.4e-21*exp(2200./tair(i,k))
            ho2rate = (rk2 + rk3*air)*(1. + fh2o*h2o(i,k)) * &
                      ho2(i,k)*ho2(i,k)
            !        xh2o2 = h2o23d(i,k) * (287.*tair(i,k)/pres(i,k))*34./6.022e20

            ! Gas-phase destruction of H2O2 occurs via:
            !      H2O2 + hv --> 2OH
            !      H2O2 + OH --> HO2 + H2O
            !  These reactions seem to have the same reaction rate.  Therefore
            !   parameterize this chemistry as doubling the rate of reaction of
            !   H2O2 + OH.
      
            ! h2o2*air*28.966/34. = h2o2 in molec/cm3

            h2o2des = 1.7e-12*exp(-200.*(1./tair(i,k) - 0.0033557)) * &
                      oh(i,k) * h2o2(i,k)*air*mwa2h2o2

            h2o2phot = jh2o2(i,k) * h2o2(i,k)*air*mwa2h2o2

            h2o2src(i,k) = ho2rate/(air*mwa2h2o2)
            h2o2snkg(i,k) = (h2o2des + h2o2phot)/(air*mwa2h2o2)
!           h2o2new(i,k) = h2o2(i,k) + (ho2rate - h2o2des - h2o2phot) * &
!                           ztodt/(air*mwa2h2o2)
            h2o2new(i,k) = h2o2(i,k) + (h2o2src(i,k) - h2o2snkg(i,k)) * ztodt

            !c        h2o2new(i,k) = min(h2o2(i,k) + 
            !c     _      ho2rate * ztodt*(287.*tair(i,k)/pres(i,k))*34./6.022e20,
            !c     _      xh2o2)

         end do
      end do

   end subroutine gaschem

!##############################################################################

   subroutine getj( clat, lat, lon, calday, zm, ncol, jout )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Set h2o2 photolysis rates.
      ! 
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         clat(pcols),           &! latitudes (radians)
         calday,                &! Current calendar day = julian day + fraction
         zm(pcols,pver)          ! height of midpoint of layer

      integer, intent(in) ::&
         lat(pcols),            &! latitude indices
         lon(pcols),            &! longitude indices
         ncol                    ! number of atmospheric columns used in chunk

      real(r8), intent(out) ::&
         jout(pcols,pver)        ! photolysis rate for location

      ! Local variables
      logical dodiavg                 ! when true, diurnally avg zenith angle
      real(r8) cosza(pcols)               ! cosine of zenith angle
      real(r8) dicosza                    ! diurnally avgd cosine of zenith angle
      real(r8) disecza                    ! diurnally avgd cosine of zenith angle

      integer i, is, k
      integer ind(3)
      integer mz
      real(r8) psum
      real(r8) ratza                      !ratio of zenith angles
      real(r8) ratalb                     !ratio of albedos = (0.3-0.2)/(0.5-0.2)
      real(r8) sumrat                     !sum of the ratios
      real(r8) w1                         !1-sumrat
      !-----------------------------------------------------------------------c
      dodiavg = .true.
      call coszen( calday, dodiavg, clat, lat, lon, ncol, cosza )
! TBH:  Loops are not in optimal order for memory access...  
      do i=1,ncol
         dicosza = cosza(i)
   
         if( dicosza .lt. 0.02 ) then
            do k=1,pver
               jout(i,k) = 0.
            end do
         else
            disecza = min( 1._r8/max( dicosza, 1.e-7_r8 ), 50._r8 )
      
            ! Interpolate data onto current level and zenith angle for alb=0.3
      
            ! zenith angle:
            do is = 1, 8                            !zasec = za in table
               if( zasec(is) .gt. disecza ) go to 5 !secza = za at gridpt
            end do
            is = is - 1
5           is = max( is-1, 1 )
            ratza = (disecza - zasec(is)) * delang(is)
      
            !     albedo -- set to 0.3 (data is for 0.2 and 0.5)
            ratalb = (0.3 - alb(1)) * delalb
      
            sumrat = ratza + ratalb
            w1 = 1. - sumrat
            
            ! Set indices
            ! note albedo index is either 1 or 2; therefore only one jh2o2 will
            ! have ialb=2
            ind(1) = nz2*(is-1) 
            ind(2) = nz2*is 
            ind(3) = nz2*nza + nz2*(is-1) 
      
            do k = 1, pver
               mz = nint( min( 50._r8, zm(i,k)*0.001_r8 + 1._r8 ) )
               psum = w1*jper2(ind(1)+mz) + ratza*jper2(ind(2)+mz) + &
                      ratalb*jper2(ind(3)+mz)
               jout(i,k) = exp( psum )
            end do
         endif
      end do

   end subroutine getj

!##############################################################################

end module sulchem

