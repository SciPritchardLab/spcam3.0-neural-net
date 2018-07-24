c-----------------------------------------------------------------
c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c       initialize parameters:
c-----------------------------------------------------------------
	
	subroutine setparm()
	
	implicit none
	include 'common.inc'
	
	include 'params.inc'
		
c----------------------------
c  Set defaults:
c
	dodamping 	= .true.
	doupperbound   	= .false.
	docloud   	= .true.
	doprecip        = .true.
	dolongwave	= .false.
	doshortwave	= .false.
	dosgs		= .true.
	dosubsidence	= .false.
	docoriolis	= .false.
	dosurface	= .true.
	dolargescale    = .false.
	doradforcing    = .false.
	dosfcforcing    = .false.
	doradsimple	= .false.
	donudging_uv	= .false.
	donudging_tq	= .false.
	dosmagor	= .true.
	doscalar	= .false.
	doensemble	= .false.
	CEM		= .false.
	LES		= .false.
	OCEAN		= .false.
	LAND		= .false.
	SFC_FLX_FXD	= .true.
	SFC_TAU_FXD	= .false.
		
	rank = 0
	dompi = .false.
	nadams		= 3
	dt		= crm_dt
	dx		= crm_dx
	dy		= crm_dy
	dz		= 0
	longitude	= 0.
	latitude	= 0.
	fcor	        = -999.
	day0		= 0.
	nrad		= 1
	ug		= 0.
	vg		= 0.
	fluxt0		= 0.
	fluxq0		= 0.
	tau0		= 0.
	taux0		= 0.
	tauy0		= 0.
	z0		= 0.035
	timelargescale  = 0.
	tabs_s 		= 0.
	nstop 		= 0
	caseid		= 'les00000'
	nstat		= 1000
	nstatfrq	= 6
	nprint		= 1000
	nrestart 	= 0
	nsave3D		= 0
	nsave3Dstart	= 0
	nsave3Dend	= 0
	nensemble	= 0
c------------------------------------
c  Set parameters 
c

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  stop
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  stop
	endif

	pi = acos(-1.)
	if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude*pi/180.)
	fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)	  
	
	if(ny.eq.1) dy=dx

	na = 1
	nb = 2
	nc = 3

        a_bg = 1./(tbgmax-tbgmin)
        a_pr = 1./(tprmax-tprmin)
        a_gr = 1./(tgrmax-tgrmin)


	nstep = 0
	dtn = 0.

	fac_cond = lcond/cp
	fac_fus = lfus/cp
	fac_sub = lsub/cp

	notopened = .true.


	return
	end
