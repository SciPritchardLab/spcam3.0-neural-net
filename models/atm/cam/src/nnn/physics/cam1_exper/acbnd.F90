module acbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! This code does time interpolation for 3D boundary data in a netCDF
   ! file.  Assumptions on the data in the netCDF file are:
   ! 1. Coordinates are ordered (lon,lev,lat,time)
   ! 2. The time coordinate is in days, and the data is assumed to be periodic
   !    annually.
   !
   ! Author: B. Eaton
   ! Modified 11 Jan 2000: to use F90 interface and allow 4 byte reals
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc, plon, plat, plev
!nf90   use netcdf

   implicit none
   save
   private
   public :: &
      acbndini,      &! initialize acbnd module
      acbndint,      &! interpolate acbnd data to requested date/time
      acbndget        ! return latitude slice data at current date/time

#include <netcdf.inc>

   real(r8), allocatable, dimension(:) :: &
      time            ! time coordinate (calander days + frac)
   real(r8), allocatable, dimension(:,:,:,:,:) :: &
      datin           ! input data (plon,plev,plat,2,nvar)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      dat             ! interpolated data (plon,plev,plat,nvar)

   character(len=128), allocatable, dimension(:) :: &
      varnam          ! names of variables to get from netCDF file

   integer :: &
      nvar,          &! number of variables requested
      ncid,          &! ID for netCDF file
      nrec,          &! number of records (time samples)
      lotim,         &! time(lotim) .le. current time
      hitim,         &! current time .lt. time(hitim)
      loin,          &! index into input data array containing time(lotim) data
      hiin,          &! index into input data array containing time(hitim) data
      start(4),      &! start vector for netCDF hyperslabs
      count(4)        ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################

   subroutine acbndini( ncfile, calday, xnvar, ncnam )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing annual cycle data.  Initialize
      ! arrays with the data to be interpolated to the current time.
      !
      ! It is assumed that the time coordinate is increasing and represents
      ! calendar days.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday   ! current time in calendar days + fraction.
      integer, intent(in) ::&
         xnvar    ! number of requested variables in the netCDF file
      character(len=*), intent(in) ::&
         ncfile,        &! netCDF data file
         ncnam(xnvar)    ! names of requested variables in the netCDF file

      ! Local variables:
      integer ::&
         i,     &
         istat, &
         recid, &! record ID
         vid     ! variable ID
      !-----------------------------------------------------------------------

      if ( xnvar .eq. 0 ) return

      nvar = xnvar
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      count(1) = plon
      count(2) = plev
      count(3) = plat
      count(4) = 1

      ! Allocate space for data.
      allocate( varnam(nvar), stat=istat )
      call alloc_err( istat, 'acbndini', 'varnam', nvar )
      allocate( datin(plon,plev,plat,2,nvar), stat=istat )
      call alloc_err( istat, 'acbndini', 'datin', plon*plev*plat*2*nvar )
      allocate( dat(plon,plev,plat,nvar), stat=istat )
      call alloc_err( istat, 'acbndini', 'dat', plon*plev*plat*nvar )

      do i = 1, nvar
         varnam(i) = ncnam(i)
      end do

      ! Open file.
!nf90      call handle_ncerr( nf90_open( trim( ncfile ), NF_NOWRITE, ncid ), &
!nf90         'acbndini: error opening file '//trim( ncfile ) )

      call handle_ncerr( nf_open( trim( ncfile ), NF_NOWRITE, ncid ), &
         'acbndini: error opening file '//trim( ncfile ) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'acbndini: no record variables ' )

      call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
         'acbndini: no record variables ' )

      ! Get size of unlimited dimension.
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'acbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), 'acbndini: ' )

      ! Allocate space for time coordinate data.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, 'acbndini', 'time', nrec )

      ! Get time coordinate.
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'time', vid ), &
!nf90         'acbndini: cannot find time coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, time ), &
!nf90         'acbndini: error getting time coordinate data' )

      call handle_ncerr( nf_inq_varid( ncid, 'time', vid ), &
         'acbndini: cannot find time coordinate variable' )
      call handle_ncerr( nf_get_var_double( ncid, vid, time ), &
         'acbndini: error getting time coordinate data' )

      ! Make sure the time coordinate looks like calander day, and is
      ! increasing.
      call chktime( time, nrec )

      ! Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      ! Read data.
      loin = 1
      hiin = 2

      do i = 1, nvar
!nf90         call handle_ncerr( nf90_inq_varid( ncid, varnam(i), vid ), &
!nf90            'acbndini: cannot find variable '//varnam(i) )
         call handle_ncerr( nf_inq_varid( ncid, varnam(i), vid ), &
            'acbndini: cannot find variable '//varnam(i) )

         start(4) = lotim
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, datin(:,:,:,loin,i), start, count ), &
!nf90            'acbndini: cannot read data for '//varnam(i) )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, datin(:,:,:,loin,i) ), &
            'acbndini: cannot read data for '//varnam(i) )

         start(4) = hitim
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, datin(:,:,:,hiin,i), start, count ), &
!nf90            'acbndini: cannot read data for '//varnam(i) )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, datin(:,:,:,hiin,i) ), &
            'acbndini: cannot read data for '//varnam(i) )
      end do

      write(6,*)'acbndini: calendar day = ',calday, ' : read data for days ', &
         time(lotim), ' and ',time(hitim)

   end subroutine acbndini

!#######################################################################

   subroutine acbndint( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate annual cycle surface flux data to the current time.  Read
      ! in new time samples of the input data as necessary.
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
         i,                  &
         oldhitim,           &
         vid, ret
      real(r8) ::&
         dt, dt1, tint
      !-----------------------------------------------------------------------

      if ( nvar .eq. 0 ) return

      ! Check to see if model time is still bounded by dataset times.
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
         ! Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1
         start(4) = hitim
         do i = 1, nvar
!nf90            call handle_ncerr( nf90_inq_varid( ncid, varnam(i), vid ), &
!nf90               'acbndint: cannot find variable '//varnam(i) )
!nf90            call handle_ncerr( nf90_get_var( ncid, vid, datin(:,:,:,hiin,i), start, count ), &
!nf90               'acbndint: cannot read data for '//varnam(i) )

            call handle_ncerr( nf_inq_varid( ncid, varnam(i), vid ), &
               'acbndint: cannot find variable '//varnam(i) )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, datin(:,:,:,hiin,i) ), &
               'acbndint: cannot read data for '//varnam(i) )
         end do
         write(6,*)'acbndint: read data for day ',time(hitim)

         if ( lotim .ne. oldhitim ) then
            ! Read in new lotim data.  Replace old hitim data.
            start(4) = lotim
            do i = 1, nvar
!nf90               call handle_ncerr( nf90_inq_varid( ncid, varnam(i), vid ), &
!nf90                  'acbndint: cannot find variable '//varnam(i) )
!nf90               call handle_ncerr( nf90_get_var( ncid, vid, datin(:,:,:,loin,i), start, count ), &
!nf90                  'acbndint: cannot read data for '//varnam(i) )

               call handle_ncerr( nf_inq_varid( ncid, varnam(i), vid ), &
                  'acbndint: cannot find variable '//varnam(i) )
               call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, datin(:,:,:,loin,i) ), &
                  'acbndint: cannot read data for '//varnam(i) )
            end do
            write(6,*)'acbndint: read data for day ',time(lotim)
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
      do i = 1, nvar
         call linintp( plon*plev*plat, 0._r8, 1._r8, tint, &
                       datin(1,1,1,loin,i), datin(1,1,1,hiin,i), dat(1,1,1,i) )
      end do

   end subroutine acbndint

!#######################################################################

   subroutine acbndget( lat, lon, name, ndim1, maxdim1, slice )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return annual cycle data for the requested latitudes and variable.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         ndim1,             &  ! 1st dimension of slice array
         lat(ndim1),        &  ! requested latitude indices
         lon(ndim1),        &  ! requested longititude indices
         maxdim1               ! maximum index of 1st dimension of slice array
      character(len=*), intent(in) ::&
         name       ! requested variable

      real(r8), intent(out) ::&
         slice(ndim1,plev)    ! requested data (not really a slice anymore)

      ! Local variables:
      integer ::&
         i, k, n
      !-----------------------------------------------------------------------

      if ( nvar .eq. 0 ) call endrun

      do n = 1, nvar
         if ( name .eq. trim(varnam(n)) ) then
            do k = 1, plev
               do i = 1, maxdim1
                  slice(i,k) = dat(lon(i),k,lat(i),n)
               end do
            end do
            return
         end if
      end do

      write(6,*)'getdat: ',name,' not found'
      call endrun
   end subroutine acbndget

!#######################################################################

end module acbnd
