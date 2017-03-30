module soxbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Interpolate the SOX emissions data from the GEIA-SMITH datasets.  This
   ! dataset contains seasonal average data for various years at roughly 15
   ! year intervals.  The interpolation scheme does linear interpolations to
   ! the requested calendar day in the seasonal cycle for each of
   ! the two years in the dataset that bound the requested year.  Finally a 
   ! linear interpolation is done to the requested year.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On non-CRAY
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !
   ! Modified 11 Jan 2000 PJR: work with 4 or 8 byte floats 
   !
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
!nf90   use netcdf

   implicit none
   save
   private
   public :: &
      soxbndini,      &! initialize soxbnd module
      soxbndint,      &! interpolate soxbnd data to requested date/time
      soxbndget        ! return latitude slice data at current date/time

   ! public module data
   
   public :: soxems_data ! full pathname for time-variant SOx emissions dataset
   character(len=256) :: soxems_data

   ! private module data

#include <netcdf.inc>

   integer, parameter ::&
      ndlev=2            ! number of levels in emissions data

   integer, allocatable, dimension(:) :: &
      date,             &! date coordinate (yyyymmdd format)
      yr                 ! years of annual cycle data
   real(r8), allocatable, dimension(:) :: &
      cdays              ! mid-season calendar days (fractional)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      soxyrlo,          &! SOX input data for lower bound year (plon,plat,ndlev,ntpy)
      soxyrhi            ! SOX input data for upper bound year
   real(r8), allocatable, dimension(:,:,:) :: &
      sox                ! SOX interpolated data (plon,plat,ndlev)

   integer ::&
      ntpy,             &! number of time samples per year in emissions dataset
      nyr,              &! number of years in emissions dataset
      ncid,             &! ID for netCDF file
      loyri,            &! index in yr array for lower bound year
      start(4),         &! start vector for netCDF hyperslabs
      count(4)           ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################


   subroutine soxbndini( year )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing SOX emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that each year for which data is available contains
      ! the same number of time samples.
      !
      ! It is assumed that the input data contains a "date" coordinate that is
      ! increasing and represents dates in the format yyyymmdd.
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr
      use ioFileMod, only: getfil
      use filenames, only: sox_emis

#ifdef MATCH
      use calendar
#endif

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         year    ! current year

      ! Local variables:
      integer ::&
         i, istat, did, nlon, vid, recid, nrec
      ! Externals
!      real(r8) ::&
!         caldayr
      !-----------------------------------------------------------------------
      ! Get file name.  
      call getfil(sox_emis, soxems_data, 0)

      ! Open file.
!nf90      call handle_ncerr( nf90_open( soxems_data, NF_NOWRITE, ncid ), &
!nf90         'soxbndini: error opening file '//trim(soxems_data) )
      call handle_ncerr( nf_open( soxems_data, NF_NOWRITE, ncid ), &
         'soxbndini: error opening file '//trim(soxems_data) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'soxbndini: no record variables ' )
      call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
         'soxbndini: no record variables ' )

      ! Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'soxbndini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'soxbndini: ' )
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), 'soxbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), 'soxbndini: ' )
      if ( nlon .ne. plon ) then
         write(6,*)'soxbndini: number of longitudes (', nlon, ')', &
                   ' doesn''t match model resolution.'
         call endrun
      end if

      ! Get size of unlimited dimension.
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'soxbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), 'soxbndini: ' )

      ! Get date coordinate.
      allocate( date(nrec), stat=istat )
      call alloc_err( istat, 'soxbndini', 'date', nrec )
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'date', vid ), &
!nf90         'soxbndini: cannot find date coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, date ), &
!nf90         'soxbndini: error getting date coordinate data' )
      call handle_ncerr( nf_inq_varid( ncid, 'date', vid ), &
         'soxbndini: cannot find date coordinate variable' )
      call handle_ncerr( nf_get_var_int( ncid, vid, date ), &
         'soxbndini: error getting date coordinate data' )

      ! Determine number of time samples per year.
      ntpy = 1
      do i = 2, nrec
         if ( date(i)/10000 .eq. date(1)/10000 ) ntpy = ntpy + 1
      end do

      ! Construct the years array.
      nyr = nrec/ntpy
      allocate( yr(nyr), stat=istat )
      call alloc_err( istat, 'soxbndini', 'yr', nyr )
      do i = 1, nyr
         yr(i) = date((i-1)*ntpy+1)/10000
      end do
      write(6,*)'soxbndini: years in emission dataset:',(yr(i),i=1,nyr)

      ! Construct array of calendar days for the annual cycle.
      allocate( cdays(ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'cdays', ntpy )
      do i = 1, ntpy
#ifdef MATCH         
         cdays(i) = caldayr( date(i), 0 )      ! match version
#else
         call bnddyi( date(i), 0, cdays(i) ) ! ccm version
#endif
      end do
      write(6,*)'soxbndini: calendar days in annual cycle:', (cdays(i),i=1,ntpy)

      write (6,*) 'soxbndini: searching for emissions for year ', year

      ! Check that requested year is contained in the data.
      if ( year .lt. yr(1) .or. year .gt. yr(nyr) ) then
         write(6,*)'soxbndini: requested year outside data limits'
         call endrun
      end if

      ! Find index for the data year that is the lower bound of the
      ! interval containing the input year.
      do i = 1, nyr
         if ( yr(i) .gt. year ) then
            loyri = i - 1
            exit
         end if
      end do

      allocate( soxyrlo(plon,plat,ndlev,ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'soxyrlo', plon*plat*ndlev*ntpy )
      allocate( soxyrhi(plon,plat,ndlev,ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'soxyrhi', plon*plat*ndlev*ntpy )
      allocate( sox(plon,plat,ndlev), stat=istat )
      call alloc_err( istat, 'soxbndini', 'sox', plon*plat*ndlev )

      ! Read SOx data for years surrounding initial year.
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 2
      count(4) = ntpy

!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'SOx', vid ), &
!nf90         'soxbndini: cannot find variable '//'SOx' )
      call handle_ncerr( nf_inq_varid( ncid, 'SOx', vid ), &
         'soxbndini: cannot find variable '//'SOx' )

      start(4) = (loyri-1)*ntpy + 1
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, soxyrlo, start, count ), &
!nf90         'soxbndini: cannot read data for '//'SOx' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxyrlo ), &
         'soxbndini: cannot read data for '//'SOx' )

      start(4) = start(4) + ntpy
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, soxyrhi, start, count ), &
!nf90         'soxbndini: cannot read data for '//'SOx' )
      call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxyrhi ), &
         'soxbndini: cannot read data for '//'SOx' )

      write(6,*)'soxbndini: read data for years; ', yr(loyri),' and ',yr(loyri+1)

   end subroutine soxbndini

!#######################################################################

   subroutine soxbndint( year, calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate SOX data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         year    ! current year
      real(r8), intent(in) ::&
         calday  ! calendar day (w/ fraction) in current year, range = [1.,366.)

      ! Local variables:
      integer ::&
         i, j, k, n,    &
         vid,           &
         lotim, hitim
      real(r8) ::&
         dt, dt1, tint,             &
         anncyclo(plon,plat,ndlev), &! low year data interpolated to calday
         anncychi(plon,plat,ndlev)   ! high year data interpolated to calday
      !-----------------------------------------------------------------------

      ! Check to see if model year is still bounded by dataset years.
      if ( year .gt. yr(nyr) ) then
         write(6,*)'soxbndint: requested year = ',year, ' last dataset year = ',yr(nyr)
         call endrun
      end if

      if ( year .gt. yr(loyri+1) ) then
         loyri = loyri + 1
         ! Read in new soxyrhi data.  Replace old soxyrlo data.
         soxyrlo = soxyrhi
         start(4) = start(4) + ntpy
!nf90         call handle_ncerr( nf90_inq_varid( ncid, 'SOx', vid ), &
!nf90            'soxbndint: cannot find variable '//'SOx' )
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, soxyrhi, start, count ), &
!nf90            'soxbndint: cannot read data for '//'SOx' )
         call handle_ncerr( nf_inq_varid( ncid, 'SOx', vid ), &
            'soxbndint: cannot find variable '//'SOx' )
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxyrhi ), &
            'soxbndint: cannot read data for '//'SOx' )
         write(6,*)'soxbndint: read data for year; ',yr(loyri+1)
      end if

      ! Linear interpolation...  Start by computing the number of days between
      !                          the lower and upper bounds, and days between
      !                          the model time and lower bound.

      if ( ntpy .gt. 1 ) then

         call findplb( cdays, ntpy, calday, lotim )
         hitim = mod( lotim, ntpy ) + 1

         if( cdays(hitim) .lt. cdays(lotim) )then
            dt = 365. - cdays(lotim) + cdays(hitim)
            if( calday .le. cdays(hitim) )then
               dt1 = 365. - cdays(lotim) + calday
            else
               dt1 = calday - cdays(lotim)
            end if
         else
            dt = cdays(hitim) - cdays(lotim)
            dt1 = calday - cdays(lotim)
         end if
         tint = dt1/dt

         ! Annual cycle interpolations.
         call linintp( plon*plat*2, 0._r8, 1._r8, tint, soxyrlo(1,1,1,lotim), &
                       soxyrlo(1,1,1,hitim), anncyclo )
         call linintp( plon*plat*2, 0._r8, 1._r8, tint, soxyrhi(1,1,1,lotim), &
                       soxyrhi(1,1,1,hitim), anncychi )

      else
         anncyclo(:,:,:) = soxyrlo(:,:,:,1)
         anncychi(:,:,:) = soxyrhi(:,:,:,1)
      end if

      ! Interpolate between years for which annual cycle data is present
      dt = yr(loyri+1) - yr(loyri)
      dt1 = year - yr(loyri)
      tint = dt1/dt
      call linintp( plon*plat*2, 0._r8, 1._r8, tint, anncyclo, anncychi, sox )

   end subroutine soxbndint

!#######################################################################

   subroutine soxbndget( lat, lon, ncol, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return SOX emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::         ncol          ! 1st dimension of output array
      integer, intent(in) ::         lat(ncol)    ! requested latitude index
      integer, intent(in) ::         lon(ncol)    ! requested longitude index

      real(r8), intent(out) ::&
         x(:,:)  ! SOx emissions in Tg S/m2/s
                 ! x(pcols,ndlev)

      ! Local variables:
      integer ::&
         i
      !-----------------------------------------------------------------------

      do i = 1, ncol
         x(i,1) = sox(lon(i),lat(i),1)
         x(i,2) = sox(lon(i),lat(i),2)
      end do

   end subroutine soxbndget

!#######################################################################

end module soxbnd
