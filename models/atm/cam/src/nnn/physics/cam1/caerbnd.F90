module caerbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! This code does time interpolation for Carbon aerosol boundary data in a netCDF
   ! file.  Assumptions on the data in the netCDF file are:
   ! 1. Coordinates are ordered (lon,lat,time)
   ! 2. The time coordinate is in days, and the data is assumed to be periodic
   !    annually.
   !
   ! 
   ! Author: B. Eaton
   ! Modified 11 Jan 2001 PJR to use netcdf f90 interface so 4 byte reals can be used
   !----------------------------------------------------------------------- 

!nf90   use netcdf
   use pmgrid
#ifdef MATCH
   use precision
#else
  use shr_kind_mod,only: r8 => shr_kind_r8
#endif
   implicit none
   save
   private
   public :: &
      caerbndini,      &! initialize dmsbnd module
      caerbndint,      &! interpolate dmsbnd data to requested date/time
      caerbndget        ! return latitude slice data at current date/time

#include <netcdf.inc>
   ! private module data
   integer, parameter ::&
      nspec = 5       ! number of carbon species
   real(r8), allocatable, dimension(:) :: &
      time            ! time coordinate (calander days + frac)
   real(r8), dimension(plon,plat,nspec,2) :: &
      caerin          ! input data
   real(r8), dimension(plon,plat,nspec) :: &
      caer            ! interpolated data

   integer :: &
      ncid,          &! ID for netCDF file
      nrec,          &! number of records (time samples)
      lotim,         &! time(lotim) .le. current time
      hitim,         &! current time .lt. time(hitim)
      loin,          &! index into input data array containing time(lotim) data
      hiin,          &! index into input data array containing time(hitim) data
      start(3),      &! start vector for netCDF hyperslabs
      count(3)        ! count vector for netCDF hyperslabs

   character(len=8) ::&
      vnam(nspec)     ! Names of species in input data

!##############################################################################
contains
!##############################################################################

   subroutine caerbndini( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing carbon aerosol emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that the time coordinate is increasing and represents
      ! calendar days (range = [1.,366.)).
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr
      use ioFileMod, only: getfil
      use filenames, only: co_emis

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         did,   &
         nlon,  &
         istat, &
         recid, &! record ID
         m, mm, &
         vid
      character(len=256) :: co_emis_file

      !-----------------------------------------------------------------------

      ! Specie names in netCDF file.
      vnam(1) = 'BBOCSF'        ! biomass burning organic carbon
      vnam(2) = 'BBBCSF'        ! biomass burning black carbon
      vnam(3) = 'FFOCSF'        ! fossil fuel organic carbon
      vnam(4) = 'FFBCSF'        ! fossil fuel black carbon
      vnam(5) = 'NOCSF'         ! natural organic carbon

      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 1

      ! Get file name.  
      call getfil(co_emis, co_emis_file, 0)

      ! Open file.
!nf90      call handle_ncerr( nf90_open( co_emis_file, NF_NOWRITE, ncid ), &
!nf90         'caerbndini: error opening file '//trim(co_emis_file) )
      call handle_ncerr( nf_open( co_emis_file, NF_NOWRITE, ncid ), &
         'caerbndini: error opening file '//trim(co_emis_file) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'caerbndini: no record variables ' )
      call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
         'caerbndini: no record variables ' )

      ! Get size of unlimited dimension.
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'caerbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), 'caerbndini: ' )

      ! Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'caerbndini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'caerbndini: ' )

      call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), 'caerbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), 'caerbndini: ' )
      if ( nlon .ne. plon ) then
         write(*,*)'caerbndini: model plon = ', plon, ', dataset nlon = ', nlon
         call endrun
      end if

      ! Allocate space for time coordinate data.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, 'caerbndini', 'time', nrec )

      ! Get time coordinate.
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'time', vid ), &
!nf90         'caerbndini: cannot find time coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, time ), &
!nf90         'caerbndini: error getting time coordinate data' )

      call handle_ncerr( nf_inq_varid( ncid, 'time', vid ), &
         'caerbndini: cannot find time coordinate variable' )
      call handle_ncerr( nf_get_var_double( ncid, vid, time ), &
         'caerbndini: error getting time coordinate data' )

      ! Make sure the time coordinate looks like calander day, and is
      ! increasing.
      call chktime( time, nrec )

      ! Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      ! Read data.
      loin = 1
      hiin = 2

      start(3) = lotim
      do m = 1, nspec
!nf90         call handle_ncerr( nf90_inq_varid( ncid, vnam(m), vid ), &
!nf90            'caerbndini: cannot find variable '//vnam(m) )
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, caerin(:,:,m,loin), start, count ), &
!nf90            'caerbndini: cannot read data for '//vnam(m) )

         call handle_ncerr( nf_inq_varid( ncid, vnam(m), vid ), &
            'caerbndini: cannot find variable '//vnam(m) )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, caerin(:,:,m,loin) ), &
            'caerbndini: cannot read data for '//vnam(m) )
      end do

      start(3) = hitim
      do m = 1, nspec
!nf90         call handle_ncerr( nf90_inq_varid( ncid, vnam(m), vid ), &
!nf90            'caerbndini: cannot find variable '//vnam(m) )
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, caerin(:,:,m,hiin), start, count ), &
!nf90            'caerbndini: cannot read data for '//vnam(m) )

         call handle_ncerr( nf_inq_varid( ncid, vnam(m), vid ), &
            'caerbndini: cannot find variable '//vnam(m) )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, caerin(:,:,m,hiin) ), &
            'caerbndini: cannot read data for '//vnam(m) )
      end do

      write(*,*)'caerbndini: calendar day = ',calday, ' : read data for days ',time(lotim), &
                ' and ',time(hitim)

   end subroutine caerbndini

!#######################################################################

   subroutine caerbndint( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         oldhitim,  &
         m, vid
      real(r8) ::&
         dt, dt1, tint, &
         r0 = 0.,     &
         r1  = 1.
      !-----------------------------------------------------------------------

      ! Check to see if model time is still bounded by dataset times.
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
         ! Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1
         start(3) = hitim
         do m = 1, nspec
!nf90            call handle_ncerr( nf90_inq_varid( ncid, vnam(m), vid ), &
!nf90               'caerbndint: cannot find variable '//vnam(m) )
!nf90            call handle_ncerr( nf90_get_var( ncid, vid, caerin(:,:,m,hiin), start, count ),&
!nf90               'caerbndint: cannot read data for '//vnam(m) )

            call handle_ncerr( nf_inq_varid( ncid, vnam(m), vid ), &
               'caerbndint: cannot find variable '//vnam(m) )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, caerin(:,:,m,hiin) ),&
               'caerbndint: cannot read data for '//vnam(m) )
         end do
         write(*,*)'caerbndint: read data for day ',time(hitim)

         if ( lotim .ne. oldhitim ) then
            ! Read in new lotim data.  Replace old hitim data.
            start(3) = lotim
            do m = 1, nspec
!nf90               call handle_ncerr( nf90_inq_varid( ncid, vnam(m), vid ), &
!nf90                  'caerbndint: cannot find variable '//vnam(m) )
!nf90               call handle_ncerr( nf90_get_var( ncid, vid, caerin(:,:,m,loin), start, count ),&
!nf90                  'caerbndint: cannot read data for '//vnam(m) )

               call handle_ncerr( nf_inq_varid( ncid, vnam(m), vid ), &
                  'caerbndint: cannot find variable '//vnam(m) )
               call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, caerin(:,:,m,loin) ),&
                  'caerbndint: cannot read data for '//vnam(m) )
            end do
            write(*,*)'caerbndint: read data for day ',time(lotim)
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
      do m = 1, nspec
         call linintp( plon*plat, r0, r1, tint, caerin(1,1,m,loin), &
                       caerin(1,1,m,hiin), caer(1,1,m) )
      end do

   end subroutine caerbndint

!#######################################################################

   subroutine caerbndget( lat, lon, ncol, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return carbon aerosol emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::    ncol       ! 1st dimension of x array
      integer, intent(in) ::    lat(ncol)  ! requested latitude indeces
      integer, intent(in) ::    lon(ncol)  ! requested longitude indeces

      real(r8), intent(out) ::        x(ncol,nspec)  ! carbon emissions (kg carbon aerosol/m2/s)

      ! Local variables:
      integer ::&
         i, m
      !-----------------------------------------------------------------------

      do m = 1, nspec
         do i = 1, ncol
            x(i,m) = caer(lon(i),lat(i),m)
         end do
      end do

   end subroutine caerbndget

!#######################################################################

end module caerbnd
