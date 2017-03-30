#include <misc.h>
#include <params.h>

subroutine inital

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define initial conditions for first run of case
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          B. Boville, April 1996
!
! $Id: inital.F90,v 1.15.2.10 2003/10/22 22:53:08 rosinski Exp $
! $Author: rosinski $
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, i1, plond, masterproc
   use inidat,       only: read_inidat
   use prognostics,  only: ps, u3, v3, t3, q3, vort, div, ptimelevels, initialize_prognostics
   use buffer,       only: initialize_buffer
   use comslt,       only: initialize_comslt
   use comspe,       only: alp, dalp
   use comsrf,       only: initialize_comsrf
   use ghg_surfvals, only: ghg_surfvals_init
   use ioFileMod,    only: getfil
   use radae,        only: initialize_radbuffer
   use phys_buffer,  only: pbuf_allocate
   use phys_grid,    only: phys_grid_init
   use time_manager, only: timemgr_init
   use filenames,    only: ncdata
#if (defined COUP_CSM)
   use ccsm_msg, only: initialize_ccsm_msg
#endif
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!---------------------------Local variables-----------------------------
!
   integer n                ! index
   character(len=256) locfn ! local filename
!
!-----------------------------------------------------------------------
!
!
! Obtain initial dataset
!
   if (masterproc) then
      call getfil (ncdata, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid_ini)
   end if
!
! Check for consistent settings on initial dataset
!
   call readinitial (ncid_ini)

! Initialize time manager.

   call timemgr_init()

! Initialize ghg surface values before default initial distributions 
! are set in inidat.
   call ghg_surfvals_init()
!
! Initialize comslt and prognostics variables 
!
   call initialize_prognostics
   call initialize_comslt
!
! Set commons
!
   call initcom
!
! Define physics data structures
!
   call phys_grid_init

#if (defined COUP_CSM)
!
! Initialize ccsm arrays (must be done after phys_grid_init where
! begchunk and endchunk are defined
!
   call initialize_ccsm_msg
#endif
!
! Initialize buffer, comsrf, and radbuffer variables 
! (which must occur after the call to phys_grid_init)
!
   call pbuf_allocate('global')
   call initialize_buffer
   call initialize_comsrf
   call initialize_radbuffer
!
! Read in initial data
!
   call read_inidat
#if (! defined PVP)
! Recover space used for ALP and DALP arrays
! (no longer needed after spectral truncations
! inside of read_inidat)
   deallocate ( alp )
   deallocate ( dalp )
#endif
!       
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!
   do n=2,ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(i1:i1+plon-1,:,:,:,n) = q3(i1:i1+plon-1,:,:,:,1)
      vort(:,:,:,n) = vort(:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
   end do
   call print_memusage ('post-inidat')

   return
end subroutine inital
