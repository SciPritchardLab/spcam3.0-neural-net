#include <misc.h>
#include <params.h>

subroutine somint ()
!-----------------------------------------------------------------------
! Slab Ocean Model (SOM) boundary dataset time interpolation
!
! Time interpolate SOM mixed layer depth and Q flux to current time, 
! reading in new monthly data if necessary
!
!---------------------------Code history--------------------------------
!
! Modified for SOM:  B. Briegleb  April 1995
! Reviewed:          B. Briegleb  March 1996
! Reviewed:          J. Kiehl     April 1996
! Reviewed:          B. Briegleb  May   1996
! Rewritten for netcdf: J. Truesdale 5/18/97
!
!-----------------------------------------------------------------------
!
! $Id: somint.F90,v 1.2.2.3 2003/02/27 00:58:22 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------
!
! Local workspace
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ocean_data
   use phys_grid, only: get_ncols_p, begchunk, endchunk, scatter_field_to_chunk
   use pmgrid,    only: masterproc
   use time_manager, only: get_curr_date, get_curr_calday

   implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
!
! Local workspace
!
   integer cnt3(3)        ! array of counts for each dimension
   integer strt3(3)       ! array of starting indices
   integer yr             ! year number (only relevant if no cycling)
   integer mon            ! month number
   integer day            ! day number
   integer i,c            ! indices
   integer ncol           ! number of columns
   integer ntmp           ! temporary
   integer ncdate         ! current date
   integer ncsec          ! current seconds into the day

   real(r8) fact1, fact2      ! time interpolation factors
   real(r8) caldayloc         ! calendar day to locate (includes yr if no cycling)
   real(r8) calday            ! calendar day
   real(r8) deltat            ! time (days) between interpolating som data
   real(r8), parameter :: daysperyear = 365.0  ! Number of days in a year
!
!-----------------------------------------------------------------------
!
! Use year information only if a multiyear dataset
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
   if (sstcyc) then
      caldayloc = calday
   else
      yr = ncdate/10000
      caldayloc = calday + yr*daysperyear
   end if
      
   strt3(1) = 1
   strt3(2) = 1
   strt3(3) = 1
   cnt3(1)  = lonsiz
   cnt3(2)  = latsiz
   cnt3(3)  = 1
!
! If model time is past current forward som timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for sstcyc = .true. 
! interpolation between December and January (np1.eq.1).  Note that 
! np1 is never 1 when sstcyc is .false.
!
   if (caldayloc > cdaysomp .and. .not. (np1 == 1 .and. caldayloc > cdaysomm)) then
      if (sstcyc) then
         np1 = mod (np1,12) + 1
      else
         np1 = np1 + 1
      end if
      if (np1 > psomtim) then
         write(6,*)'SOMINT: Attempt to read past end of SOM dataset'
         call endrun ()
      end if
      cdaysomm = cdaysomp
      call bnddyi (date_som(np1), sec_som(np1), cdaysomp)
      if (.not.sstcyc) then
         yr = date_som(np1)/10000
         cdaysomp = cdaysomp + yr*daysperyear
      end if
      if (np1 == 1 .or. caldayloc <= cdaysomp) then
         ntmp = nm
         nm = np
         np = ntmp
         strt3(3) = np1
         if (masterproc) then
            call wrap_get_vara_realx (ncid_sst, qfluxid,  strt3, cnt3, qfluxdat)
            write(6,*)'SOMINT: Read som for date (yyyymmdd) ', date_som(np1),' sec ',sec_som(np1)
         end if
            
         call scatter_field_to_chunk (1, 1, 1, plon, qfluxdat, qfluxm(1,begchunk,np))

      else
         write(6,*)'SOMINT: Input som for date', date_som(np1), ' sec ', sec_som(np1), &
                   'does not exceed model date', ncdate, ' sec ', ncsec, ' Stopping.'
         call endrun ()
      end if
   end if
!
! Time interpolation.  Account for December-January interpolation if
! cycling som dataset.  Again note that np1 is never 1 when sstcyc is false
!
   if (np1 == 1) then                    ! Dec-Jan interpolation
      deltat = cdaysomp + daysperyear - cdaysomm
      if (caldayloc > cdaysomp) then      ! We're in December
         fact1 = (cdaysomp + daysperyear - caldayloc)/deltat
         fact2 = (caldayloc - cdaysomm)/deltat
      else                                ! We're in January
         fact1 = (cdaysomp - caldayloc)/deltat
         fact2 = (caldayloc + daysperyear - cdaysomm)/deltat
      end if
   else
      deltat = cdaysomp - cdaysomm
      fact1 = (cdaysomp - caldayloc)/deltat
      fact2 = (caldayloc - cdaysomm)/deltat
   end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
   if (abs (fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. fact2 < -1.e-6) then
      write(6,*)'SOMINT: Bad fact1 and/or fact2=',fact1,fact2
      call endrun ()
   end if

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      do i=1,ncol
         qflux(i,c)  = qfluxm(i,c,nm)  + fact2*(qfluxm(i,c,np)  - qfluxm(i,c,nm))
      end do
   end do

   return
end subroutine somint
