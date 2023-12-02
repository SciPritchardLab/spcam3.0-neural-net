module sulbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Interpolate the GEIA sulfur emissions data (SO2, SO4, DMS) from the
   ! dataset prepared by Mary Barth.  This dataset contains seasonal average
   ! data for 1985.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On non-CRAY
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
!nf90   use netcdf

   implicit none
   save
   private
   public :: &
      sulbndini,      &! initialize sulbnd module
      sulbndint,      &! interpolate sulbnd data to requested date/time
      sulbndget        ! return latitude slice data at current date/time

   ! public module data
   
   public :: sulems_data ! full pathname for time-variant sulfer emissions dataset
   character(len=256) :: sulems_data='sulfur_emissions.nc'

   ! private module data

#include <netcdf.inc>

   integer, parameter ::&
      ndlev=2          ! number of levels in emissions data

   real(r8), allocatable, dimension(:) :: &
      time             ! time coordinate (calander days + frac)

   real(r8), dimension(plon,plat,ndlev,2) :: &
      so2in,          &! SO2 input data.  Upper and lower time bounds.
      so4in,          &! SO4 input data.  Upper and lower time bounds.
      dmsin            ! DMS input data.  Upper and lower time bounds.

   real(r8), dimension(plon,plat,ndlev) :: &
      so2,            &! SO2 interpolated data.
      so4,            &! SO4 interpolated data.
      dms              ! DMS interpolated data.

   integer ::&
      ncid,           &! ID for netCDF file
      nrec,           &! number of records (time samples)
      lotim,          &! time(lotim) .le. current time
      hitim,          &! current time .lt. time(hitim)
      loin,           &! index into input data array containing time(lotim) data
      hiin,           &! index into input data array containing time(hitim) data
      start(4),       &! start vector for netCDF hyperslabs
      count(4)         ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################

   subroutine sulbndini( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing sulfur emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that the time coordinate is increasing and represents
      ! calendar days (range = [1.,366.)).
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr
!nf90      use netcdf

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         i, istat, did, nlon, vid, recid
      character(len=*), parameter ::&
         routine_name = 'sulbndini'
      !-----------------------------------------------------------------------

      ! Open file.
      write(6,*)'Sulfer emissions dataset is: ', trim(sulems_data)
      call handle_ncerr( nf_open( trim( sulems_data ), NF_NOWRITE, ncid ), &
         routine_name//': error opening file '//trim( sulems_data ) )

      ! Check that input data is a right resolution.
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), routine_name//': ' )
      call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), routine_name//': ' )
      if ( nlon .ne. plon ) then
         write(6,*)routine_name//': number of longitudes (', nlon, ')', &
                   ' doesn''t match model resolution.'
         call endrun
      end if

      ! Get size of unlimited dimension.
      call handle_ncerr( nf_inq_unlimdim( ncid, recid ), routine_name//': ' )
      call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), routine_name//': ') 

      ! Get time coordinate.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, routine_name, 'time', nrec )
      call handle_ncerr( nf_inq_varid( ncid, 'time', vid ), &
           routine_name//': cannot find time coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, time ), &
!nf90         routine_name//': error getting time coordinate data' )
      call handle_ncerr( nf_get_var_double( ncid, vid, time ), &
         routine_name//': error getting time coordinate data' )

      ! Make sure the time coordinate looks like calander day, and is
      ! increasing.
      call chktime( time, nrec )

      ! Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      ! Read data.
      loin = 1
      hiin = 2
      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = plon
      count(2) = plat
      count(3) = ndlev
      count(4) = 1

      start(4) = lotim
      call handle_ncerr( nf_inq_varid( ncid, 'SO2', vid ), routine_name//': cannot find variable '//'SO2' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so2in(1,1,1,loin) ), &
         routine_name//': cannot read data for '//'SO2' )

      call handle_ncerr( nf_inq_varid( ncid, 'SO4', vid ), routine_name//': cannot find variable '//'SO4' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so4in(1,1,1,loin) ), &
         routine_name//': cannot read data for '//'SO4' )
      call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), routine_name//': cannot find variable '//'DMS' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmsin(1,1,1,loin) ), &
         routine_name//': cannot read data for '//'DMS' )

      start(4) = hitim
      call handle_ncerr( nf_inq_varid( ncid, 'SO2', vid ), routine_name//': cannot find variable '//'SO2' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so2in(1,1,1,hiin) ), &
         routine_name//': cannot read data for '//'SO2' )
      call handle_ncerr( nf_inq_varid( ncid, 'SO4', vid ), routine_name//': cannot find variable '//'SO4' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so4in(1,1,1,hiin) ), &
         routine_name//': cannot read data for '//'SO4' )
      call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), routine_name//': cannot find variable '//'DMS' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmsin(1,1,1,hiin) ), &
         routine_name//': cannot read data for '//'DMS' )

 

      write(6,*) routine_name//': calendar day = ',calday, ' : read data for days ',time(lotim), &
                 ' and ',time(hitim)

   end subroutine sulbndini

!#######################################################################

   subroutine sulbndint( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate sulfur data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr
!nf90      use netcdf

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current calendar day (w/ fraction), range = [1.,366.)

      ! Local variables:
      integer ::&
         oldlotim, oldhitim, &
         vid, ret
      real(r8) ::&
         dt, dt1, tint
      character(len=*), parameter ::&
         routine_name = 'sulbndint'
      !-----------------------------------------------------------------------

      ! Check to see if model time is still bounded by dataset times.
      oldlotim = lotim
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
         ! Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1
         start(4) = hitim
         call handle_ncerr( nf_inq_varid( ncid, 'SO2', vid ), routine_name//': cannot find variable '//'SO2' )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so2in(1,1,1,hiin) ), &
            routine_name//': cannot read data for '//'SO2' )
         call handle_ncerr( nf_inq_varid( ncid, 'SO4', vid ), routine_name//': cannot find variable '//'SO4' )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so4in(1,1,1,hiin) ), &
            routine_name//': cannot read data for '//'SO4' )
         call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), routine_name//': cannot find variable '//'DMS' )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmsin(1,1,1,hiin) ), &
            routine_name//': cannot read data for '//'DMS' )

         write(6,*) routine_name//': read data for day ',time(hitim)

         if ( lotim .ne. oldhitim ) then
            ! Read in new lotim data.  Replace old hitim data.
            start(4) = lotim
            call handle_ncerr( nf_inq_varid( ncid, 'SO2', vid ), routine_name//': cannot find variable '//'SO2' )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so2in(1,1,1,loin) ), &
               routine_name//': cannot read data for '//'SO2' )
            call handle_ncerr( nf_inq_varid( ncid, 'SO4', vid ), routine_name//': cannot find variable '//'SO4' )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, so4in(1,1,1,loin) ), &
               routine_name//': cannot read data for '//'SO4' )
            call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), routine_name//': cannot find variable '//'DMS' )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmsin(1,1,1,loin) ), &
               routine_name//': cannot read data for '//'DMS' )

            write(6,*) routine_name//': read data for day ',time(lotim)
         end if

      end if


      ! Linear interpolation...  Start by computing the number of days between
      !                          the lower and upper bounds, and days between
      !                          the model time and lower bound.

      if( time(hitim) .lt. time(lotim) )then
         dt = 365. - time(lotim) + time(hitim)
         if( calday .le. time(hitim) )then
            dt1 = 365. - time(lotim) + calday
         else
            dt1 = calday - time(lotim)
         end if
      else
         dt = time(hitim) - time(lotim)
         dt1 = calday - time(lotim)
      end if
      tint = dt1/dt
      call linintp( plon*plat*2, 0._r8, 1._r8, tint, so2in(1,1,1,loin), so2in(1,1,1,hiin), so2 )
      call linintp( plon*plat*2, 0._r8, 1._r8, tint, so4in(1,1,1,loin), so4in(1,1,1,hiin), so4 )
      call linintp( plon*plat*2, 0._r8, 1._r8, tint, dmsin(1,1,1,loin), dmsin(1,1,1,hiin), dms )

   end subroutine sulbndint

!#######################################################################

   subroutine sulbndget(pcols, ncol, lat, lon, semis )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return sulfur emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::     pcols
      integer, intent(in) ::     ncol
      integer, intent(in) ::     lat(pcols)             ! requested latitude indices
      integer, intent(in) ::     lon(pcols)             ! requested longitude

      real(r8), intent(out) ::&
         semis(pcols,6)  ! sulfur emissions

      ! Local variables:
      integer :: i
      !-----------------------------------------------------------------------

      do i = 1, ncol
         semis(i,1) = so2(lon(i),lat(i),1)
         semis(i,2) = so4(lon(i),lat(i),1)
         semis(i,3) = dms(lon(i),lat(i),1)
         semis(i,4) = so2(lon(i),lat(i),2)
         semis(i,5) = so4(lon(i),lat(i),2)
         semis(i,6) = dms(lon(i),lat(i),2)
      end do

   end subroutine sulbndget

!#######################################################################

end module sulbnd
