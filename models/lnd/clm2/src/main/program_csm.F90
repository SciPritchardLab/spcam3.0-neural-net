#include <misc.h>
#include <preproc.h>

#if (defined COUP_CSM)

! ---------------------- begin copyright notice -----------------------
!
!            NCAR Land Surface Model, version 2.0
!            copyright (c) 1996
!            University Corporation for Atmospheric Research
!            all rights reserved
!
!               ------------ ----- --- ---------- ------
!  ********** | distribution terms and conditions notice | ************
!               ------------ ----- --- ---------- ------
!
! (c) copyright 1996 University Corporation for Atmospheric Research/
! National Center for Atmospheric Research/
! Climate and Global Dynamics Division
!
! Access and use of this software shall impose the following obligations
! and understandings on the user.  The user is granted the right,
! without any fee or cost, to use, copy, modify, alter, enhance and
! distribute this software, and any derivative works thereof, and its
! supporting documentation for any purpose whatsoever, except commercial
! sales, provided that this entire notice appears in all copies of the
! software, derivative works and supporting documentation.  Further, the
! user agrees to credit UCAR/NCAR/CGD in any publications that result
! from the use of this software or in any software package that includes
! this software.  The names UCAR/NCAR/CGD, however, may not be used in
! any advertising or publicity to endorse or promote any products or
! commercial entity unless specific written permission is obtained from
! UCAR/NCAR/CGD.
!
! The CLM materials are made available with the understanding that
! UCAR/NCAR/CGD is not obligated to provide (and will not provide) the
! user with any support, consulting, training, or assistance of any kind
! with regard to the use, operation and performance of this software, nor
! to provide the user with any updates, revisions, new versions, or "bug
! fixes."
!
! This software is provided by UCAR/NCAR/CGD "as is" and any express or
! implied warranties, including but not limited to, the implied
! warranties of merchantability and fitness for a particular purpose are
! disclaimed.  In no event shall UCAR/NCAR/CGD be liable for any
! special, indirect or consequential damages or any damages whatsoever,
! including but not limited to claims associated with the loss of data
! or profits, which may result from an action in contract, negligence or
! other tortious claim that arises out of or in connection with the
! access, use or performance of this software.
!
! ---------------------- end copyright notice -------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: program_csm
!
! !INTERFACE:
PROGRAM program_csm
!
! !DESCRIPTION: 
! Driver for CLM as the land component of CCSM.
! This program is the driver for CLM to work as the land component of
! CCSM.  The flux coupler will provide all the appropriate atmospheric
! forcing for the land model to run.
! o the land surface model returns to the CCSM flux coupler surface
!   fluxes, temperatures, and albedos for only the land points on the
!   [lsmlon x lsmlat] grid 
! o the land surface model uses its own surface type data set. because
!   it processes only land points, this data set must correctly define
!   what points are land and what are not land
! o the land surface model uses its own grid dimensions (lsmlon and
!   lsmlat). currently these must equal lsmlon and lsmlat so that there
!   is a direct correspondence between the atmosphere and land grids
! o the zenith angle calculation is calculated for the
!   NEXT time step rather than the current time step. make sure
!   the calendar day is for the NEXT time step. make sure the
!   solar declination calculation is the same as in the 
!   atmospheric model but for the NEXT time step. make sure the
!   calendar day is for greenwich time (see next comment).
! o subroutine calendr: this generates a julian day (with fraction)
!   based on the time step, which is used to calculate the solar
!   zenith angle. this time must be at the greenwich meridian to
!   get the correct zenith angle. also, output from this subroutine 
!   is used to calculate the month (1, ..., 12), day (1, ..., 31), 
!   and year (00, ...) of the simulation. 
! o the land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. the net longwave radiation
!   at the surface will differ somewhat from that calculated from the
!   CCSM flux coupler because the cpl model will use the upward 
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. the net solar radiation should equal that calculated
!   from the flux coupler. if not, there is a problem.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar         !parameters
  use clm_varctl         !run control variables    
  use shr_orb_mod        !orbital parameters and routines 
  use shr_msg_mod        !csm message passing routines and variables    
#if (defined SPMD)
  use spmdMod      , only : masterproc, iam, spmd_init, mpicom
#else
  use spmdMod      , only : masterproc, iam
#endif
  use initializeMod, only : initialize
  use clm_csmMod   , only : csmstop_now, csm_setup, csm_shutdown
  use time_manager , only : advance_timestep, get_nstep 
!
! !ARGUMENTS:
    implicit none
#include "gpt.inc"
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: i,j        !loop indices 	
  integer :: nstep      !time step 
  logical :: doalb      !true if surface albedo calculation time step
!
! Earth's orbital characteristics
!
  integer  :: iyear_AD  !Year (AD) to simulate above earth's orbital parameters for
  real(r8) :: eccen     !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) :: obliq     !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) :: mvelp     !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)
!
! Orbital information after call to routine shr_orbit_params
!
  real(r8) :: obliqr    !Earth's obliquity in radians
  real(r8) :: lambm0    !Mean longitude (radians) of perihelion at the vernal equinox 
  real(r8) :: mvelpp    !Earth's moving vernal equinox longitude
  logical  :: log_print !true=> print diagnostics  

#if (!defined SPMD)
  integer :: mpicom_dummy ! needed for  call to csm_setup
#endif
!-----------------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! Initialize model
  ! -----------------------------------------------------------------
  
  ! Initialize timing library.  2nd arg 0 means disable, 1 means enable
  
  call t_setoptionf (usrsys, 0)
  call t_initializef ()

  ! Determine input/output units 

  call shr_msg_stdio ('lnd')

  ! Initialize MPI communication groups for flux coupler

#if (defined SPMD)
  call csm_setup(mpicom)
#else
  call csm_setup(mpicom_dummy)
#endif

  ! Initialize intra-MPI communication stuff or 
  ! set masterproc if not running in SPMD mode

#if (defined SPMD)
  call spmd_init()
#endif

  ! Initialize land model - initialize communication with flux coupler

  call initialize (eccen, obliqr, lambm0, mvelpp)

! -----------------------------------------------------------------
! Time stepping loop
! -----------------------------------------------------------------

  call t_startf('lnd_timeloop')

  ! begin time stepping loop

  do 

     ! doalb is true when the next time step is a radiation time step
     ! this allows for the fact that an atmospheric model may not do 
     ! the radiative calculations every time step. for example:
     !      nstep dorad doalb
     !        1     F     F
     !        2     F     T
     !        3     T     F
     ! The following expression for doalb is specific to CAM
     
     nstep = get_nstep() 
     doalb = ((irad==1 .and. nstep+1/=1) .or. (mod(nstep,irad)==0 .and. nstep+1/=1))
     
     ! Call land surface model driver 
     
     call driver (doalb, eccen, obliqr, lambm0 ,mvelpp)
     
     ! determine if time to stop
     
     if (csmstop_now) exit 
     
     ! increment time step 
     
     call advance_timestep()

  end do

  call t_stopf ('lnd_timeloop')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
! call shr_msg_finalize
  call csm_shutdown()

  stop
end program program_csm

#else

!The following is only here since empty file won't compile
subroutine program_csm_stub
  write(6,*) 'PROGRAM_CSM: this routine should not be called'
  stop 99
end subroutine program_csm_stub

#endif





