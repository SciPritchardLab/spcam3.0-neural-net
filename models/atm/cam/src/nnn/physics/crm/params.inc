c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c
c   Flags:
c

      logical CEM	! flag for Cloud Ensemble Model
      logical LES	! flag for Large-Eddy Simulation
      logical OCEAN	! flag indicating that surface is water
      logical LAND	! flag indicating that surface is land
      logical SFC_FLX_FXD ! surface laten and sensible fluxes are fixed
      logical SFC_TAU_FXD! surface drag is fixed

c   Constants:
c
      real  cp	       ! Specific heat of air at constant pressure, J/kg/K
      real  ggr        ! Gravity acceleration, m/s2
      real  lcond      ! Latent heat of condensation, J/kg
      real  lfus       ! Latent heat of fusion, J/kg
      real  lsub       ! Latent heat of sublimation, J/kg
      real  rv         ! Gas constant for water vapor, J/kg/K
      real  rgas       ! Gas constant for dry air, J/kg/K
      real  pi	       ! Just pi
      real  fac_cond, fac_fus, fac_sub
c
c  Parameters:
c
c            
      real  pres0      ! Reference surface pressure, Pa
      real  ug	       ! Velocity of the Domain's drift in x direction
      real  vg	       ! Velocity of the Domain's drift in y direction
      real  fcor       ! Coriolis parameter	
      real  fcorz      ! Vertical Coriolis parameter

      real longitude   
      real latitude    

      real uhl        ! current large-scale velocity in x near sfc
      real vhl        ! current large-scale velocity in y near sfc

c
c  Surface stuff:   	
c
        real   tabs_s	! surface temperature,K
	real   fluxt0   ! surface sensible flux, Km/s
	real   fluxq0   ! surface latent flux, m/s
	real   tau0     ! surface stress, m2/s2
	real   taux0    ! surface stress in x, m2/s2
	real   tauy0    ! surface stress in y, m2/s2
	real   z0	! roughness length
c
c  Domain top stuff:
c        
        real   gamt0	! gradient of t() at the top,K/m
        real   gamq0	! gradient of q() at the top,g/g/m
        
c
c  Misc. stuff
c
 	real  timelargescale ! time to start large-scale forcing
c
c  Microphysics stuff:
c
      real  diffelq    ! Diffusivity of water vapor in the air, m2/s
      real  therco     ! Thermal conductivity of air, J/m/s/K
      real  tbgmax     ! Maximum temperature for Bergeron proc. to occur, K
      real  tbgmin     ! Minimum temperature for Bergeron proc. to occur, K
      real  tprmin     ! Minimum temperature for snow+graupel formation, K
      real  tprmax     ! Maximum temperature for snow+graupel formation, K
      real  tgrmin     ! Minimum temperature for graupel formation, K
      real  tgrmax     ! Maximum temperature for graupel formation, K
      real  a_rain     ! Coeff. in fall relation for rain, m^(1-b_rain)/s
      real  b_rain     ! Fall speed exponent for rain
      real  a_snow     ! Coeff. in fall relation for snow, m^(1-b_snow)/s
      real  b_snow     ! Fall speed exponent for snow
      real  a_grau     ! Coeff. in fall relation for graupel,m^(1-b_grau)/s
      real  b_grau     ! Fall speed exponent for graupel
      real  qcw0       ! Threshold for autoconversion of cloud w. to rain, g/g  
      real  qci0l      ! Threshold for autoconversion of cloud ice to snow,g/g  (land)
      real  qci0o      ! Threshold for autoconversion of cloud ice to snow,g/g  (ocean)
      real  erccoef    ! Rain/Cloud water collection efficiency
      real  esccoef    ! Snow/Cloud water collection efficiency
      real  esicoef    ! Snow/cloud ice collection efficiency
      real  egccoef    ! Graupel/Cloud water collection efficiency
      real  egicoef    ! Graupel/Cloud ice collection efficiency
      real  nzeror     ! Intercept coeff. in Marshal-Palmer dist. for rain  
      real  nzeros     ! Intersept coeff. in Marshal-Palmer dist. for snow
      real  nzerog     ! Intersept coeff. in Marshal-Palmer dist. forgraupel
      real  alphaelq   ! Rate coeff. for autoconversion of cloud water
      real  betaelq    ! Rate coeff. for autoconversion of cloud ice
      real  muelq      ! Dynamic viscosity of air
      real  rhor       ! Density of water, kg/m3
      real  rhos       ! Density of snow, kg/m3
      real  rhog       ! Density of graupel, kg/m3
      real  gam3       ! Gamma function of 3
      real  gams1      ! Gamma function of (3 + b_snow)
      real  gams2      ! Gamma function of (5 + b_snow)/2
      real  gams3      ! Gamma function of (4 + b_snow)
      real  gamg1      ! Gamma function of (3 + b_grau)
      real  gamg2      ! Gamma function of (5 + b_grau)/2
      real  gamg3      ! Gamma function of (4 + b_grau)
      real  gamr1      ! Gamma function of (3 + b_rain)
      real  gamr2      ! Gamma function of (5 + b_rain)/2
      real  gamr3      ! Gamma function of (4 + b_rain)
      real  qp_threshold ! minimal rain/snow water content
!      real  vt_ice     ! terminal velocity of ice cystals	
      
      real  accrsc(nzm),accrsi(nzm),accrrc(nzm),coefice(nzm)
      real  accrgc(nzm),accrgi(nzm)
      real  evaps1(nzm),evaps2(nzm),evapr1(nzm),evapr2(nzm)
      real  evapg1(nzm),evapg2(nzm)
            
      real a_bg, a_pr, a_gr

	common /comparam/ cp,ggr,lcond,lfus,lsub,rv,rgas,pi,pres0,
     &		  ug,vg,fcor,fcorz,tabs_s,z0,uhl,vhl,
     &		  longitude, latitude,
     &		  fluxt0,fluxq0,gamt0,gamq0,tau0,taux0,tauy0,
     &            timelargescale, diffelq,therco,
     &            tbgmax, tbgmin, tprmin, tprmax, tgrmin, tgrmax,
     &            a_rain, b_rain, a_grau, b_grau, a_snow, b_snow,
     &            qcw0, qci0l, qci0o, egccoef, egicoef,
     &            erccoef, esccoef, esicoef, nzeror, nzeros, nzerog, 
     &            alphaelq, betaelq, muelq, rhor, rhos, rhog,
     &            gam3, gams1, gams2, gams3, gamr1, gamr2, gamr3,
     &            gamg1, gamg2, gamg3,
     &		  fac_cond,fac_fus,fac_sub,qp_threshold,
!     &		  fac_cond,fac_fus,fac_sub,qp_threshold,vt_ice,
     &		  accrsc,accrsi,accrgc,accrgi,accrrc,coefice,
     &		  evaps1,evaps2,evapg1,evapg2,evapr1,evapr2,
     &            a_bg, a_pr, a_gr
     
     
	common /logicflags/ CEM, LES, OCEAN, LAND, 
     &		  SFC_FLX_FXD, SFC_TAU_FXD

!tBM*   threadlocal /comparam/,/logicflags/
!$OMP   threadprivate (/comparam/, /logicflags/)
