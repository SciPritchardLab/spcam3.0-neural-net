#include <misc.h>
#include <params.h>

subroutine readinitial (ncid)
!----------------------------------------------------------------------- 
! 
! Purpose: Ensure that requisite netcdf variables are on the initial dataset.
!          Set base day and date info using the "current" values from it.
! 
! Method: Issue proper netcdf wrapper calls.  Broadcast to slaves if SPMD
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use rgrid
   use time_manager, only: ic_ymd, ic_tod
#if ( defined SPMD )
   use mpishorthand
#endif
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: ncid  ! History file unit
!                                
! Local variables
!                                
   integer :: lonid            !------------------------------------------------------
   integer :: levid            ! 
   integer :: latid            ! 
   integer :: ntrkid           ! 
   integer :: ntrmid           ! 
   integer :: ntrnid           ! 
   integer :: ncdateid         ! Netcdf variable and dimension ids for variable of that
   integer :: ncsecid          ! name with "id" tacked on to the end
   integer :: hyaiid           ! 
   integer :: hybiid           ! 
   integer :: hyamid           ! 
   integer :: hybmid           ! 
   integer :: ilev             ! 
   integer :: ilevid           ! 
   integer :: rlonid           ! 
   integer :: nlonid           ! 
   integer :: wnummaxid        !------------------------------------------------------

   integer :: mlon             ! longitude dimension length from dataset
   integer :: mlev             ! level dimension length from dataset
   integer :: morec            ! latitude dimension length from dataset

#ifdef STAGGERED
   integer :: slonid           ! staggered longitude dimension length from dataset
   integer :: slatid           ! staggered latitude dimension length from dataset
#endif

!
!-----------------------------------------------------------------------
!
   if (masterproc) then
!
! Get and check dimension/date info
!
      call wrap_inq_dimid (ncid, 'lon' , lonid)
      call wrap_inq_dimid (ncid, 'lev' , levid)
      call wrap_inq_dimid (ncid, 'ilev', ilevid)
      call wrap_inq_dimid (ncid, 'lat' , latid)
!
#ifdef STAGGERED
      call wrap_inq_dimid (ncid, 'slon' , slonid)
      call wrap_inq_dimid (ncid, 'slat' , slatid)
#endif
      call wrap_inq_varid (ncid, 'ntrk'   , ntrkid)
      call wrap_inq_varid (ncid, 'ntrm'   , ntrmid)
      call wrap_inq_varid (ncid, 'ntrn'   , ntrnid)
      call wrap_inq_varid (ncid, 'date'   , ncdateid)
      call wrap_inq_varid (ncid, 'datesec', ncsecid)
      call wrap_inq_varid (ncid, 'hyai'   , hyaiid)
      call wrap_inq_varid (ncid, 'hybi'   , hybiid)
      call wrap_inq_varid (ncid, 'hyam'   , hyamid)
      call wrap_inq_varid (ncid, 'hybm'   , hybmid)
!
      call wrap_inq_dimlen (ncid, lonid , mlon)
      call wrap_inq_dimlen (ncid, levid , mlev)
      call wrap_inq_dimlen (ncid, ilevid, ilev)
      call wrap_inq_dimlen (ncid, latid , morec)
!
! Check for reduced grid info on initial dataset.  If not present, define
! variables for full grid
!
      if (nf_inq_varid (ncid, 'nlon', nlonid)  ==  nf_noerr) then
         call wrap_inq_varid (ncid, 'wnummax', wnummaxid)

         call wrap_get_var_int (ncid, nlonid, nlon)
         call wrap_get_var_int (ncid, wnummaxid, wnummax)
      else
         wnummax(:) = ptrm
         nlon(:) = plon
      end if

      call wrap_get_var_int (ncid, ncdateid, ic_ymd)
      call wrap_get_var_int (ncid, ncsecid , ic_tod)

      if (mlev /= plev.or.mlon /= plon.or.morec /= plat) then
         write(6,*)'READINITIAL: model parameters do not match initial dataset parameters'
         write(6,*)'Model Parameters:   plev = ',plev,' plon = ',plon,' plat = ',plat
         write(6,*)'Dataset Parameters: dlev = ',mlev,' dlon = ',mlon,' dlat = ',morec
         call endrun
      end if

      call wrap_get_var_realx (ncid, hyamid,hyam)
      call wrap_get_var_realx (ncid, hybmid,hybm)
      call wrap_get_var_realx (ncid, hyaiid,hyai)
      call wrap_get_var_realx (ncid, hybiid,hybi)
   end if

#if ( defined SPMD )
   call mpibcast (ic_ymd,  1,    mpiint, 0, mpicom)
   call mpibcast (ic_tod,  1,    mpiint, 0, mpicom)
   call mpibcast (nlon,    plat, mpiint, 0, mpicom)
   call mpibcast (wnummax, plat, mpiint, 0, mpicom)

   call mpibcast (hyam  ,plev ,mpir8,  0, mpicom)
   call mpibcast (hybm  ,plev ,mpir8,  0, mpicom)
   call mpibcast (hyai  ,plevp,mpir8,  0, mpicom)
   call mpibcast (hybi  ,plevp,mpir8,  0, mpicom)
#endif

   return
end subroutine readinitial
