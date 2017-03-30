	
subroutine setparm
	
!       initialize parameters:

use vars
use params
!use isccp, only : isccp_zero
!use isccpTables, only : isccp_tables_init
implicit none
	
!NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
!                dolongwave, doshortwave, dosgs, &
!                docoriolis, dosurface, dolargescale, doradforcing, &
!		nadams,fluxt0,fluxq0,tau0,tabs_s,z0,tauls,nelapse, &
!		dt, dx, dy, fcor, ug, vg, nstop, caseid, &
!		nstat, nstatfrq, nprint, nrestart, doradsimple, &
!		nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
!		donudging_uv, donudging_tq, dosmagor, doscalar,  &
!		timelargescale, longitude0, latitude0, day0, nrad, &
!		CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
!                doensemble, nensemble, doxy, dowallx, dowally, &
!                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
!                dt_cup, docup, docolumn, save2Dbin, save3Dbin, &
!                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
!	        doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
!	        dodynamicocean, ocean_type, &
!		dosolarconstant, solar_constant, zenith_angle, raddir, &
!                dotracers
	
!----------------------------
!  Set defaults:

	dodamping 	= .true.
	doupperbound   	= .false.
	docloud   	= .true.
	doprecip        = .true.
	dolongwave	= .false.
	doshortwave	= .false.
	doradsimple 	= .false.
	dosgs		= .true.
	dosmagor	= .true.
	doscalar	= .false.
	dosubsidence	= .false.
	docoriolis	= .false.
	dosurface	= .true.
	dolargescale    = .false.
	doradforcing    = .false.
	dosfcforcing    = .false.
	donudging_uv	= .false.
	donudging_tq	= .false.
	doensemble	= .false.
	doxy    	= .false.
	dowallx    	= .false.
	dowally    	= .false.
	docup		= .false.
	docolumn	= .false.
	doseasons	= .false.
	doperpetual	= .false.
	doradhomo	= .false.
	dosfchomo	= .false.
	doisccp		= .false.
	dodynamicocean 	= .false.
        dosolarconstant = .false.
        dotracers       = .false.
	CEM		= .false.
	LES		= .true.
	OCEAN		= .false.
	LAND		= .false.
	SFC_FLX_FXD	= .true.
	SFC_TAU_FXD	= .false.
		
	nadams		= 3
	dt		= crm_dt
	dt_cup		= 0.
	dx		= crm_dx
	dy		= crm_dy
	longitude0	= 0.
	latitude0	= 0.
	fcor	        = -999.
	day0		= 0.
	nrad		= 1
	ug		= 0.
	vg		= 0.
	fluxt0		= 0.
	fluxq0		= 0.
	tau0		= 0.
	z0		= 0.035
        soil_wetness    = 1.
	timelargescale  = 0.
	tauls		= 7200.
	tabs_s 		= 0.
	nstop 		= 0
	nelapse		= 999999999
	caseid		= 'les00000'
	nstat		= 1000
	nstatfrq	= 50
	nprint		= 1000
	nrestart 	= 0
	restart_sep 	= .false.
	save3Dbin	= .false.
	save2Dsep	= .false.
	save3Dsep	= .false.
	nsave3D		= 1
	nsave3Dstart	= 99999999
	nsave3Dend	= 999999999
	dogzip2D	= .false.
	dogzip3D	= .false.
	save2Dbin	= .false.
	nsave2D		= 1
	nsave2Dstart	= 99999999
	nsave2Dend	= 999999999
	nensemble	= 0
 	qnsave3D	= 0.
	ocean_type	= 0 
        raddir          = './RADDATA'

        ! Specify solar constant and zenith angle for perpetual insolation.
        ! Note that if doperpetual=.true. and dosolarconstant=.false.
        ! the insolation will be set to the daily-averaged value on day0.

        solar_constant = 685. ! Values from Tompkins & Craig, J. Climate (1998)
        zenith_angle = 51.7

!----------------------------------
!  Read namelist variables from the standard input:
!------------

!open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
!read (55,PARAMETERS)
!close(55)

!------------------------------------
!  Set parameters 


	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif

	pi = acos(-1.)
	if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)
	fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)	  
	coszrs = 0.637 ! default mean solar zenith angle
	
	if(ny.eq.1) dy=dx

	na = 1
	nb = 2
	nc = 3
	nstep = 0
        time = 0.
	dtn = dt

	a_bg = 1./(tbgmax-tbgmin)
	a_pr = 1./(tprmax-tprmin)
	a_gr = 1./(tgrmax-tgrmin)

	notopened2D = .true.
	notopened3D = .true.

!        call isccp_tables_init()   ! initialize isccp tables
!	call isccp_zero()


        rank =0
        dompi = .false.


        masterproc = rank.eq.0
end
