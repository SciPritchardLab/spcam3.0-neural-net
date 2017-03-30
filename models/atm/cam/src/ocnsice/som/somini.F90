#include <misc.h>
#include <params.h>

subroutine somini ()
!-----------------------------------------------------------------------
! Slab Ocean Model (SOM) boundary dataset initialization 
!
! Do initial read of time-varying som boundary dataset, reading two
! consecutive months on either side of the current model date.
!
!---------------------------Code history--------------------------------
!
! Original version:  L. Bath
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Hack, B. Boville, August 1992
! Modified for SOM:  B. Briegleb  April 1995
! Reviewed:          B. Briegleb  March 1996
! Reviewed:          J. Kiehl     April 1996
! Reviewed:          B. Briegleb  May   1996
! Rewritten for netcdf: J. Truesdale 5/7/97
!-----------------------------------------------------------------------
!
! $Id: somini.F90,v 1.2.2.4 2003/02/27 00:58:22 rosinski Exp $
! $Author: rosinski $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use time_manager, only: get_curr_calday, get_curr_date
   use ppgrid,       only: pcols
   use phys_grid,    only: scatter_field_to_chunk, begchunk, endchunk, get_ncols_p
   use physconst,    only: stebol
   use pmgrid,       only: masterproc
   use ocean_data
#if ( defined SPMD )
   use mpishorthand
#endif

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
   integer dateid                ! netcdf id for date variable
   integer secid                 ! netcdf id for seconds variable
   integer lonid                 ! netcdf id for longitude variable
   integer latid                 ! netcdf id for latitude variable
   integer timeid                ! netcdf id for time variable
   integer cnt3(3)               ! array of counts for each dimension
   integer strt3(3)              ! array of starting indices
   integer yr                    ! year number (only relevant if no cycling)
   integer mon                   ! month number
   integer day                   ! day number
   integer n                     ! indices
   integer ncdate                ! current date
   integer ncsec                 ! seconds in current day

   real(r8) caldayloc            ! calendar day to locate (includes yr if no cycling)
   real(r8) calday               ! calendar day
   real(r8), parameter :: daysperyear = 365.0  ! Number of days in a year

!
   call allocate_ocean ()       ! Allocate memory for boundary data
!
! Use year information only if not cycling som dataset
!
   calday = get_curr_calday()
   call get_curr_date (yr, mon, day, ncsec)

   if (sstcyc) then
      caldayloc = calday
   else
      caldayloc = calday + yr*daysperyear
   end if
!
! Read mixed layer depths and q flux info from boundary dataset.
! Get and check dimension info
!
   if (masterproc) then
      call wrap_inq_dimid (ncid_sst, 'lon', lonid)
      call wrap_inq_dimid (ncid_sst, 'time', timeid)
      call wrap_inq_dimid (ncid_sst, 'lat', latid)

      call wrap_inq_dimlen (ncid_sst, lonid, lonsiz)
      if (lonsiz /= plon) then
         write(6,*)'SOMINI: lonsiz=',lonsiz,' must = plon=',plon
         call endrun ()
      end if

      call wrap_inq_dimlen (ncid_sst, latid, latsiz)
      if (latsiz  /=  plat) then
         write(6,*)'SOMINI: latsiz=',latsiz,' must = plat=',plat
         call endrun ()
      end if

      call wrap_inq_dimlen (ncid_sst, timeid, timesiz)
!
! Check to make sure space allocated for time variables is sufficient
!
      if (timesiz > totsomsz) then
         write(6,*)'SOMINI:  Allocated space for som date variables is insufficient.'
         write(6,*)'Increase parameter totsomsz to at least', timesiz
         call endrun ()
      end if
!
      call wrap_inq_varid (ncid_sst, 'lon',     lonid)
      call wrap_inq_varid (ncid_sst, 'date',    dateid)
      call wrap_inq_varid (ncid_sst, 'datesec', secid)
      call wrap_inq_varid (ncid_sst, 'MLDANN',  mldid)
      call wrap_inq_varid (ncid_sst, 'QFLUX',   qfluxid)
      call wrap_inq_varid (ncid_sst, 'lat',     latid)
!
! Retrieve entire date and sec variables.
!
      call wrap_get_var_int (ncid_sst,dateid,date_som)
      call wrap_get_var_int (ncid_sst,secid,sec_som)
      if (sstcyc) then
         if (mod (date_som(1), 10000)/100 /= 1) then
            write(6,*)'SOMINI: When cycling som, 1st month must be 1'
            call endrun ()
         end if
         if (mod(date_som(psomtim),10000)/100 /= 12) then  
            write(6,*)'SOMINI: When cycling som, must have 12 consecutive months of data'
            call endrun ()
         end if
      end if

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = 1
   end if

#if ( defined SPMD )
   call mpibcast (date_som, totsomsz, mpiint, 0, mpicom)
   call mpibcast (sec_som,  totsomsz, mpiint, 0, mpicom)
#endif
!
! Get mixed layer depths (time-invariant)
!
   if (masterproc) then
      call wrap_get_vara_realx (ncid_sst, mldid, strt3, cnt3, mlddat)
   end if
   call scatter_field_to_chunk (1, 1, 1, plon, mlddat, mld(1,begchunk))
!
! Special code for interpolation between December and January
!
   if (sstcyc) then
      n = 12
      np1 = 1
      call bnddyi (date_som(n  ), sec_som(n  ), cdaysomm)
      call bnddyi (date_som(np1), sec_som(np1), cdaysomp)

      if (caldayloc <= cdaysomp .or. caldayloc > cdaysomm) then
         strt3(3) = n
         if (masterproc) then
            call wrap_get_vara_realx (ncid_sst, qfluxid,  strt3, cnt3, qfluxdat)
            write(6,*)'SOMINI: Read som data for date (yyyymmdd) ', date_som(n),' sec ',sec_som(n)
         end if
         
         call scatter_field_to_chunk (1, 1, 1, plon, qfluxdat, qfluxm(1,begchunk,nm))

         strt3(3) = np1
         if (masterproc) then
            call wrap_get_vara_realx (ncid_sst, qfluxid,  strt3, cnt3, qfluxdat)
            write(6,*)'SOMINI: Read som data for date (yyyymmdd) ', date_som(np1),' sec ',sec_som(np1)
         end if
         
         call scatter_field_to_chunk (1, 1, 1, plon, qfluxdat, qfluxm(1,begchunk,np))
         
         return
      end if
   end if
!
! Normal interpolation between consecutive time slices.
!
   do n=1,timesiz-1
      np1 = n + 1
      call bnddyi(date_som(n  ), sec_som(n  ), cdaysomm)
      call bnddyi(date_som(np1), sec_som(np1), cdaysomp)
      
      if (.not.sstcyc) then
         yr = date_som(n)/10000
         cdaysomm = cdaysomm + yr*daysperyear
         yr = date_som(np1)/10000
         cdaysomp = cdaysomp + yr*daysperyear
      end if

      if (caldayloc > cdaysomm .and. caldayloc <= cdaysomp) then
         strt3(3) = n
         if (masterproc) then
            call wrap_get_vara_realx (ncid_sst, qfluxid,  strt3, cnt3, qfluxdat)
         end if

         call scatter_field_to_chunk (1, 1, 1, plon, qfluxdat, qfluxm(1,begchunk,nm))

         strt3(3) = np1
         if (masterproc) then
            call wrap_get_vara_realx (ncid_sst, qfluxid,  strt3, cnt3, qfluxdat)
         end if

         call scatter_field_to_chunk (1, 1, 1, plon, qfluxdat, qfluxm(1,begchunk,np))

         if (masterproc) then
            write(6,*)'SOMINI: Read som data for dates ',date_som(n), sec_som(n),' and ', &
                      date_som(np1),sec_som(np1)
         end if
         return
      end if
   end do

   ncdate = yr*10000 + mon*100 + day

   write(6,*)'SOMINI: Failed to find dates bracketing mcdate, mcsec=', ncdate, ncsec
   call endrun ()
end subroutine somini
