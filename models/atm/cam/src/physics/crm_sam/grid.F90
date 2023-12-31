module grid

use domain

implicit none
        
integer, parameter :: nx = nx_gl/nsubdomains_x
integer, parameter :: ny = ny_gl/nsubdomains_y 
integer, parameter :: nz = nz_gl+1
integer, parameter :: nzm = nz-1
        
integer, parameter :: nsubdomains = nsubdomains_x * nsubdomains_y

logical, parameter :: RUN3D = ny_gl.gt.1
logical, parameter :: RUN2D = .not.RUN3D

integer, parameter :: nxp1 = nx + 1
integer, parameter :: nyp1 = ny + 1 * YES3D
integer, parameter :: nxp2 = nx + 2
integer, parameter :: nyp2 = ny + 2 * YES3D
integer, parameter :: nxp3 = nx + 3
integer, parameter :: nyp3 = ny + 3 * YES3D
integer, parameter :: nxp4 = nx + 4
integer, parameter :: nyp4 = ny + 4 * YES3D

integer, parameter :: dimx1_u = -1
integer, parameter :: dimx2_u = nxp3
integer, parameter :: dimy1_u = 1-2*YES3D
integer, parameter :: dimy2_u = nyp2
integer, parameter :: dimx1_v = -1
integer, parameter :: dimx2_v = nxp2
integer, parameter :: dimy1_v = 1-2*YES3D
integer, parameter :: dimy2_v = nyp3
integer, parameter :: dimx1_w = -1
integer, parameter :: dimx2_w = nxp2
integer, parameter :: dimy1_w = 1-2*YES3D
integer, parameter :: dimy2_w = nyp2
integer, parameter :: dimx1_s = -2
integer, parameter :: dimx2_s = nxp3
integer, parameter :: dimy1_s = 1-3*YES3D
integer, parameter :: dimy2_s = nyp3

integer, parameter :: ncols = nx*ny


real dx 	! grid spacing in x direction
real dy		! grid spacing in y direction
real dz		! grid spacing in z direction for the lowest grid layer

real z(nz)      ! height of the pressure levels above surface,m
real pres(nzm)  ! pressure,mb at scalar levels
real zi(nz)     ! height of the interface levels
real presi(nz)  ! pressure,mb at interface levels
real adz(nzm)   ! ratio of the grid spacing to dz for pressure levels
real adzw(nz)	! ratio of the grid spacing to dz for w levels
real grdf_x(nzm)! grid factor for eddy diffusion in x
real grdf_y(nzm)! grid factor for eddy diffusion in y
real grdf_z(nzm)! grid factor for eddy diffusion in z

real at, bt, ct ! coefficients for the Adams-Bashforth scheme 
real dt		! dynamical timestep
real dtn	! current dynamical timestep (can be smaller than dt)
real dt3(3) 	! dynamical timesteps for three most recent time steps
real dt_cup	! timestep for convection parameterization 
real time	! current time in sec.
real day0	! starting day (including fraction)
real day	! current day (including fraction)
real dtfactor   ! dtn/dt

        
integer nstep	! current number of performed time steps 
integer nstop   ! time step number to stop the integration
integer nelapse ! time step number to elapse before stoping
integer na, nb, nc ! indeces for swapping the rhs arrays for AB scheme
integer ncycle  ! number of subcycles over the dynamical timestep
integer icycle  ! current subcycle 
integer nadams	! the order of the AB scheme (should be kept at 3)        
integer nstat	! the interval in time steps to compute statistics
integer nstatis	! the interval between substeps to compute statistics
integer nstatfrq! frequency of computing statistics 
integer nprint 	! frequency of printing a listing (steps)
integer nrestart! switch to control starting/restarting of the model
logical restart_sep ! write separate restart files for sub-domains
integer nrad	! frequency of calling the radiation routines
logical save3Dbin ! save 3D data in binary format(no byte compression)
logical save3Dsep !write a separate file for snapshots for 2-D model  
logical save2Dsep !write a separate file for 2D horizontal fields for 3-D model  
integer nsave3D ! frequency of writting 3D fields (steps)
integer nsave3Dstart ! timestep to start writting 3D fields
integer nsave3Dend   ! timestep to end writting 3D fields
real    qnsave3D !threshold manimum cloud water(kg/kg) to save 3D fields
logical dogzip3D ! gzip compress a 3D output file   
logical dogzip2D ! gzip compress a 2D output file if save2Dsep=.true.   
logical save2Dbin !save 2D data in binary format, rather than compressed
integer nsave2D ! frequency of writting 2D fields (steps)
integer nsave2Dstart ! timestep to start writting 2D fields
integer nsave2Dend   ! timestep to end writting 2D fields
character *40 caseid! 8-symbol id-string to identify a run	
character *40 case  ! 8-symbol id-string to identify a case-name	
logical dostatis! flag to permit the gathering of statistics
integer nensemble ! the number of subensemble set of perturbations
logical notopened2D ! flag to see if the 2D output datafile is opened	
logical notopened3D ! flag to see if the 3D output datafile is opened	
character *256 raddir ! path to data directory containing files for CAM radiation


!   Flags:

logical CEM     ! flag for Cloud Ensemble Model
logical LES     ! flag for Large-Eddy Simulation
logical OCEAN   ! flag indicating that surface is water
logical LAND    ! flag indicating that surface is land
logical SFC_FLX_FXD  ! surface sensible flux is fixed
logical SFC_TAU_FXD  ! surface drag is fixed

!       Multitasking staff:     
          
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
logical masterproc ! .true. if rank.eq.0 

	
!   Logical switches and flags:

logical   dodamping, doupperbound, docloud, doprecip, &
          dolongwave, doshortwave, dosgs, dosubsidence, &
          docoriolis, dosurface, dolargescale, doradforcing, &
          dosfcforcing, doradsimple, donudging_uv, donudging_tq, & 
          dosmagor, doscalar, doensemble, doxy, dowallx, dowally, docup, &
          docolumn, doperpetual, doseasons, doradhomo, dosfchomo, &
          doisccp, dodynamicocean, dosolarconstant, dotracers

! For dosolarconstant simulations, allow solar constant and zenith
! angle to be specified individually
real solar_constant  ! solar constant (in W/m2)
real zenith_angle    ! zenith angle (in degrees)

common/comgrd1/dx,dy,dz,z,pres,zi,presi,adz,adzw,grdf_x,grdf_y,grdf_z, &
               at,bt,ct,dt,dtn,dt3,dt_cup,time,day0,day,dtfactor
common/comgrd2/ nstep,nstop,nelapse,na,nb,nc,ncycle,icycle, nadams, &
                CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, &
                caseid, case, dostatis, &
                nstat, nstatis, nstatfrq, nprint, nrestart, nrad
common /comgrd3/ rank, ranknn, rankss, rankee, rankww, &
                      rankne, ranknw, rankse, ranksw, dompi, masterproc
!$OMP threadprivate (/comgrd1/,/comgrd2/,/comgrd3/)


end module grid
