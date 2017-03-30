	include 'crmdims.inc'

        integer nx_gl           ! number of the global grid points in x-dir.
        integer ny_gl           ! number of the global grid points in y-dir.
        integer nz_gl           ! number of the global grid points in z-dir.
        integer npressureslabs  ! number of pressure slabs for MPI
        integer nsubdomains_x   ! number of subdomain divisions in x-direction
        integer nsubdomains_y   ! number of subdomain divisions in y-direction
c--------------------------------------------

        parameter (nx_gl = crm_nx)
        parameter (ny_gl = crm_ny)
        parameter (nz_gl = crm_nz+1)
        parameter (nsubdomains_x  = 1)
        parameter (nsubdomains_y  = 1)
        parameter (npressureslabs = 1)

c--------------------------------------------------------------
        integer nsubdomains     ! total number of subdomains
        parameter (nsubdomains = nsubdomains_x * nsubdomains_y)

        logical RUN3D           ! logical flag for 3D run
        logical RUN2D           ! logical flag for 2D run
        parameter (RUN3D = ny_gl.gt.1)
        parameter (RUN2D = ny_gl.eq.1)
!----------------------------------------------------------------------
        
        integer nx 	! number of subdomain grid points in x-dir
        integer ny	! number of subdomain grid points in y-dir
        integer nz	! number of subdomain grid points in z-dir

	parameter (nx = nx_gl/nsubdomains_x)
	parameter (ny = ny_gl/nsubdomains_y) 
	parameter (nz = nz_gl)
        
        integer nzm, nxp1, nyp1, nxp2, nyp2, nxp3, nyp3, nxp4, nyp4

	parameter (nxp1 = nx + 1)
	parameter (nyp1 = ny + 1 * YES3D)		
	parameter (nxp2 = nx + 2)
	parameter (nyp2 = ny + 2 * YES3D)		
	parameter (nxp3 = nx + 3)
	parameter (nyp3 = ny + 3 * YES3D)		
	parameter (nxp4 = nx + 4)
	parameter (nyp4 = ny + 4 * YES3D)		
	parameter (nzm = nz-1)


        integer dimx1_u,dimx2_u,dimy1_u,dimy2_u
        integer dimx1_v,dimx2_v,dimy1_v,dimy2_v
        integer dimx1_w,dimx2_w,dimy1_w,dimy2_w
        integer dimx1_s,dimx2_s,dimy1_s,dimy2_s
      parameter(dimx1_u=-1,dimx2_u=nxp3,dimy1_u=1-2*YES3D,dimy2_u=nyp2)
      parameter(dimx1_v=-1,dimx2_v=nxp2,dimy1_v=1-2*YES3D,dimy2_v=nyp3)
      parameter(dimx1_w=-1,dimx2_w=nxp2,dimy1_w=1-2*YES3D,dimy2_w=nyp2)
      parameter(dimx1_s=-2,dimx2_s=nxp3,dimy1_s=1-3*YES3D,dimy2_s=nyp3)

	
	real dx 	! grid spacing in x direction
	real dy		! grid spacing in y direction
	real dz		! grid spacing in z direction for the lowest grid layer
	real adz(nzm)   ! ratio of the grid spacing to dz for pressure levels
	real adzw(nz)	! ratio of the grid spacing to dz for w levels
        real rho(nzm)   ! air density at pressure levels,kg/m3
        real rhow(nz)   ! air density at vertical velocity levels,kg/m3
        real at, bt, ct ! coefficients for the Adams-Bashforth scheme 
        real dt		! dynamical timestep
        real dtn	! current dynamical timestep (can be smaller than dt)
        real dt3(3) 	! dynamical timesteps for three most recent time steps
        real time	! current time in sec.
	real day0	! starting day (including fraction)
	real day	! current day (including fraction)
        
        common /comgrid/ dx, dy, dz, adz, adzw, rho, rhow, 
     $		         at, bt, ct, dt, dtn, dt3, time, day, day0
     
        integer nstep	! current number of performed time steps 
        integer nstop   ! time step number to stop the integration
        integer na, nb, nc ! indeces for swapping the rhs arrays for AB scheme
        integer ncycle  ! number of subcycles over the dynamical timestep
        integer icycle  ! current subcycle 
        integer nadams	! the order of the AB scheme (should be kept at 3)        
        integer nstat	! the interval in time steps to compute statistics
        integer nstatis	! the interval between substeps to compute statistics
        integer nstatfrq! frequency of computing statistics 
        integer nprint 	! frequency of printing a listing (steps)
        integer nrestart! switch to control starting/restarting of the model
	integer nrad	! frequency of calling the radiation routines
        integer nsave3D ! frequency of writting 3D fields (steps)
	integer nsave3Dstart ! timestep to start writting 3D fields
	integer nsave3Dend   ! timestep to end writting 3D fields
        character *40 caseid! 8-symbol id-string to identify a run	
        character *40 case  ! 8-symbol id-string to identify a case-name	
        logical dostatis! flag to permit the gathering of statistics
	integer nensemble ! the number of subensemble set of perturbations
        logical notopened ! flag to see if the 2D output datafile is opened	
        common /comstep/ 
     &          nstep, na, nb, nc, nadams, ncycle, icycle, nstop,
     &          caseid, case, dostatis,
     &          nstat, nstatis, nstatfrq, nprint, nrestart, nrad, 
     &	        nsave3D, nsave3Dstart, nsave3Dend, nensemble, 
     &          notopened
c       Multitasking staff:     
          
        integer rank   ! rank of the current subdomain task (default 0) 
	integer ranknn ! rank of the "northern" subdomain task
	integer rankss ! rank of the "southern" subdomain task
	integer rankee ! rank of the "eastern"  subdomain task
	integer rankww ! rank of the "western"  subdomain task
	integer rankne ! rank of the "north-eastern" subdomain task
	integer ranknw ! rank of the "north-western" subdomain task
	integer rankse ! rank of the "south-eastern" subdomain task
	integer ranksw ! rank of the "south-western" subdomain task
	logical dompi  ! logical switch to do multitasking
	
        common /comranks/ rank, ranknn, rankss, rankee, rankww, 
     $			  rankne, ranknw, rankse, ranksw, dompi 
!tBM*	threadlocal /comstep/
!tBM*	threadlocal /comgrid/
!tBM*	threadlocal  /comranks/
!$OMP	threadprivate (/comgrid/)
!$OMP	threadprivate (/comstep/)
!$OMP	threadprivate  (/comranks/)
     			 
