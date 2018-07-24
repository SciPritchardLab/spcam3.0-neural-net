#include <misc.h>
#include <preproc.h>

module initializeMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: initializeMod
! 
! !DESCRIPTION: 
! Performs land model initialization 
!
! !USES:
  use spmdMod, only : masterproc
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein 
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  public :: header    ! echo version numbers
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize
!
! !INTERFACE:
#if (defined OFFLINE) || (defined COUP_CSM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp    )
#elif (defined COUP_CAM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp    , cam_caseid , & 
                        cam_ctitle  , cam_nsrest  , cam_nstep , cam_irad  , cam_crtinic, &
                        cam_nhtfrq  , cam_mfilt   , cam_longxy, cam_latixy, cam_numlon , &
                        cam_landmask, cam_landfrac, cam_irt   )
#endif
!
! !DESCRIPTION: 
! Land model initialization.  Initialization routine for land surface
! model. Initialize land surface variables for atmosphere model if not a
! continuation run.  Model constants are set in block data
! subprograms. This subroutine:
! o Initializes run control variables via the [clmexp] namelist. 
! o Initializes active Fortran unit numbers 1 to 99, except 5 and 6, to false.
! o Reads surface data on model grid. 
! o Defines the multiple plant types and fraction areas for each surface type. 
! o Builds the appropriate subgrid <-> grid mapping indices and weights. 
! o Set up parallel processing. 
! o Initializes time constant variables.
! o Initializes history file output.
! o Initializes river basin data.
! o Initializes accumulation variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an
!   initial run.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar, only : lsmlon, lsmlat, maxpatch
    use clm_varsur, only : varsur_alloc, varsur_dealloc
    use clm_varctl, only : fsurdat, finidat, nsrest, irad, &
                           mksrf_offline_fgrid, mksrf_offline_fnavyoro
    use clmtypeInitMod, only : clmtypeInit
    use clmpoint, only : clmpoint_init
    use controlMod, only : control_init, control_print
    use mksrfdatMod, only : mksrfdat
    use surfFileMod, only : surfrd
    use pftvarcon, only : pftconrd
    use clm_mapping, only : clm_map, clm_map1d
    use histFldsMod, only : initHistFlds
    use restFileMod, only : restart
#if (defined OFFLINE)
    use atmdrvMod, only : atm_getgrid
#endif
     use accFldsMod, only : initAccFlds, initAccClmtype
#if (defined DGVM)
     use DGVMMod, only : resetTimeConstDGVM
     use EcosystemDynDGVMMod, only : EcosystemDynini
#else
     use mvegFileMod, only : monthveg_ini, interpMonthlyVeg
#endif
#if (defined BGC) 
    use DustMod, only : Dustini
#endif
#if (defined COUP_CAM)
    use time_manager, only: get_curr_date, get_nstep
#else               
    use time_manager, only: get_curr_date, get_nstep, advance_timestep, timemgr_init 
#endif
#if (defined RTM) 
    use RtmMod, only : Rtmgridini, Rtmlandini
#endif
#if (defined COUP_CSM)
    use clm_csmMod    
#endif
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(inout) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(inout) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(inout) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
#if (defined COUP_CAM)
    character(len=*),  intent(in) :: cam_caseid   !cam caseid
    character(len=*),  intent(in) :: cam_ctitle   !cam title
    character(len=*),  intent(in) :: cam_crtinic  !cam initial dataset generation frequency
    integer ,  intent(in) :: cam_irad             !cam radiation frequency
    integer ,  intent(in) :: cam_nsrest           !cam 0=initial run, > 0=continuation run
    integer ,  intent(in) :: cam_nstep            !cam current time index
    integer ,  intent(in) :: cam_nhtfrq           !cam history write freq for tape 1
    integer ,  intent(in) :: cam_mfilt            !cam number of files per tape for tape 1
    integer ,  intent(in) :: cam_irt              !cam mss retention time
    integer ,  intent(in) :: cam_numlon(:)        !cam number of longitudes 
    real(r8),  intent(in) :: cam_longxy(:,:)      !cam lon values
    real(r8),  intent(in) :: cam_latixy(:,:)      !cam lat values 
    real(r8),  intent(in) :: cam_landfrac(:,:)    !cam fractional land
    integer ,  intent(in) :: cam_landmask(:,:)    !cam land mask
#endif
!
! !CALLED FROM:
! routine program_off if cpp token OFFLINE is defined
! routine program_csm if cpp token COUP_CSM is defined
! routine atmlnd_ini in module atm_lndMod if cpp token COUP_CAM is defined
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein 
!
!EOP
!
! !LOCAL VARIABLES:
    logical  :: readini               !true if read in initial data set
    integer  :: i,j,k                 !indices
    integer  :: yr                    !current year (0, ...)
    integer  :: mon                   !current month (1 -> 12)
    integer  :: day                   !current day (1 -> 31)
    integer  :: ncsec                 !current time of day [seconds]
    integer  :: vegxy(lsmlon,lsmlat,maxpatch) !vegetation type
    real(r8) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid weights
#if (defined COUP_CSM)
    integer  :: cam_numlon(lsmlat)           !cam number of longitudes 
    real(r8) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8) :: cam_latixy(lsmlon,lsmlat)    !cam lat values 
    real(r8) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer  :: cam_landmask(lsmlon,lsmlat)  !cam land mask
#endif
    integer  :: nstep                        !current time step
!-----------------------------------------------------------------------

    ! Initialize run control variables, time manager, timestep

    call header()

    if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
#if (defined COUP_CSM) || (defined OFFLINE)
       write (6,*) 'Preset Fortran unit numbers:'
       write (6,*) '   unit  5 = standard input'
       write (6,*) '   unit  6 = standard output'
#endif
       write (6,*)
    endif

#if (defined COUP_CAM)
    call control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                       cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )
#else
    call control_init ()
#endif
    if (masterproc) call control_print()

    ! Allocate surface grid dynamic memory

    call varsur_alloc ()

#if (defined OFFLINE) || (defined COUP_CSM)
    ! Initialize time manager for initial run

    if (nsrest == 0) call timemgr_init()
#endif    

#if (defined OFFLINE)
    ! Start at nstep = 1 for an initial offline run

    if (nsrest == 0) call advance_timestep()
#endif

#if (defined RTM) 
    ! Initialize RTM river routing grid and mask

     call Rtmgridini ()
#endif

#if (defined COUP_CSM)

     ! Get grid and land mask back from flux coupler
     
     call csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask) 
#endif

     ! Read list of PFTs and their corresponding parameter values
     ! This is independent of the model resolution

     call pftconrd ()

     ! If no surface dataset name is specified then make surface dataset
     ! from original data sources. Always read surface boundary data in.
     ! This insures that bit for bit results are obtained for a run where a 
     ! surface dataset file is generated and a run where a surface dataset 
     ! is specified and read in. Set up vegetation type [veg] and weight [wt]
     ! arrays for [maxpatch] subgrid patches.
     
#if (defined OFFLINE)
     if (fsurdat == ' ') then
        call mksrfdat ()
     endif
     call surfrd (vegxy, wtxy)
#else
     if (fsurdat == ' ') then
        call mksrfdat (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask) 
     endif
     call surfrd (vegxy, wtxy, &
          cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
#endif
     
     ! Build mapping indices and weights
     
     call clm_map (vegxy, wtxy)
     
     ! Initialize clm pointer arrays
     ! Must be called after clm_map

     call clmpoint_init()

     ! Initialize values of clmtype data structures
     ! Must be called after clmpoint_init (due to pointers utilized)
	
     call clmtypeInit()

#if (defined DGVM)
     ! Initialize DGVM LPJ
     
     call EcosystemDynini ()
#else
     ! Allocate and initialize non-DGVM dynamic memory 

     call monthveg_ini  
#endif

#if (defined BGC) 
     ! Initialize dust emissions model
     
     call Dustini ()
#endif

     ! Initialize time constant variables 

     if (masterproc) write (6,*) 'Attempting to initialize time invariant variables'
     call iniTimeConst ()
     if (masterproc) write (6,*) 'Successfully initialized time invariant variables'

     ! Do 1d mapping (must be done after iniTimeConst to have values for water type) 

     call clm_map1d
     
#if (defined RTM) 
     ! Initialize river routing model(s) 

     if (masterproc) write(6,*)'Attempting to initialize RTM'
     call Rtmlandini ()
     if (masterproc) write(6,*)'Successfully initialized RTM'
#endif

#if (defined COUP_CSM)

     ! Initialize flux coupler communication

     call csm_initialize(irad,eccen, obliqr, lambm0, mvelpp)

#endif

    ! Read restart files if continuation run

    if (nsrest > 0) call restart('read')

    ! Initialize master history list. Note, routine initHistFlds will 
    ! initialize active history fields if not a restart run.

    call initHistFlds () 

#if (defined DGVM)
    ! Initialize/reset DGVM time constant variables

    if (nsrest > 0) call resetTimeConstDGVM
#endif

#if (defined OFFLINE)
    ! Read atmospheric forcing dataset one time to obtain the longitudes 
    ! and latitudes of the atmospheric dataset, as well as the edges. When
    ! coupled to atm model, these are input variables. If no
    ! atmospheric data files are provided, model uses dummy atmospheric
    ! forcing and sets atmospheric grid to land grid.

    if (masterproc) write (6,*) 'Attempting to set up atmospheric grid '
    call atm_getgrid ()
    if (masterproc) write (6,*) 'Successfully set up atmospheric grid '
#endif

    ! Get current date

    call get_curr_date(yr, mon, day, ncsec)

    ! Initialize accumulator fields to be time accumulated for various purposes.

    if (nsrest == 0) then
       call initAccFlds ()
    end if

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0 for cam and csm mode 
    ! and at nstep=1 for offline mode. This routine is also always called for a 
    ! restart run and must therefore be called after the restart file is read in

    call initAccClmtype ()

#if (!defined DGVM)
    ! Read monthly vegetation data for interpolation to daily values
    ! Note: routine mapinit needs to be called first since interpMonthlyVeg
    ! needs the array landvec%wtxy
    
    call interpMonthlyVeg (fsurdat, mon, day) 
#endif

    ! If initial run: initialize time-varying data 
    ! If continuation run: end of initialization because time varying
    ! read in from restart file
    
    if (nsrest == 0) then
       if (masterproc) write (6,*) 'Attempting to initialize time variant variables '
       if (finidat == ' ') then
          readini = .false.
       else
          readini = .true.
       end if
       call iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp)
       if (masterproc) then
          write (6,*) 'Successfully initialized time variant variables'
          write (6,*)
       endif
    endif

#if (defined COUP_CSM)
    ! Send first land model data to flux coupler. 
    
    call csm_sendalb ()
#endif
    
    ! Deallocate surface grid dynamic memory

    call varsur_dealloc ()

    ! End initialization

    if (masterproc) then
       write (6,*) 'Successfully initialized the land model'
       if (nsrest == 0) then
          write (6,*) 'begin initial run at: '
       else
          write (6,*) 'begin continuation run at:'
       end if
       write (6,*) '   nstep= ',get_nstep(), &
            ' year= ',yr,' month= ',mon,' day= ',day,' seconds= ',ncsec
       write (6,*)
       write (6,'(72a1)') ("*",i=1,60)
       write (6,*)
    endif

  end subroutine initialize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: header
!
! !INTERFACE:
  subroutine header()
!
! !DESCRIPTION: 
! Echo and save model version number
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : version
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

    version = 'CLM MODEL version 2.1'
    if ( masterproc )then
      write (6,*) trim(version)
      write (6,*)
    end if

  end subroutine header

end module initializeMod
