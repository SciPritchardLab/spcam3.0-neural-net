c $Id: acsf.F,v 1.1.6.1 2003/07/29 00:13:45 hender Exp $
c It is assumed that the model calling this interface has been
c compiled so that 8 byte real data are being used.  On non CRAY
c machines this implies compiling with a "-r8" flag.

c This code does time interpolation for surface flux data in a netCDF
c file.  Assumptions on the data in the netCDF file are:
c 1. Flux arrays have coordinates ordered (lon,lat,time)
c 2. The time coordinate is in days, and the data is assumed to be periodic
c    annually.

      module acsf

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid,      only: plat, plev, plevp, plond, plon, masterproc
      use constituents, only: pcnst, pnats
      use error_messages, only: alloc_err, handle_ncerr

      implicit none
      save
      private
      public ::
     $  iniacsf      ! initialize acsf module
     $, acsfint      ! interpolate acsf data to requested date/time
     $, getsf        ! return latitude slice acsf data at current date/time

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c annual cycle surface flux data.

      real(r8), dimension(:), allocatable ::
     $  time     ! time coordinate (calander days + frac)

      real(r8), dimension(:,:,:,:), allocatable ::
     $  sfin     ! surface flux input data

      real(r8), dimension(:,:,:), allocatable ::
     $  sf       ! surface flux interpolated data

      character
     $  sfnam(pcnst)*32 ! names of surface flux variables

      integer
     $  nsf            ! number of surface flux variables requested
     $, sfcnst(pcnst)  ! constituent index into the sflx(plon,pcnst) array
c                      ! passed in from calling routine.
     $, ncid           ! ID for netCDF file
     $, nrec           ! number of records (time samples)
     $, lotim          ! time(lotim) .le. current time
     $, hitim          ! current time .lt. time(hitim)
     $, loin          ! index into input data array containing time(lotim) data
     $, hiin          ! index into input data array containing time(hitim) data
     $, start(3)       ! start vector for netCDF hyperslabs
     $, count(3)       ! count vector for netCDF hyperslabs


c##############################################################################
      contains
c##############################################################################

      subroutine iniacsf( ncfile, calday, sfncnam )

c Open netCDF file containing annual cycle surface flux data.  Initialize
c arrays with the data to be interpolated to the current time.

c It is assumed that the time coordinate is increasing and represents
c calendar days; range = [1.,366.).

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input argument:
      real(r8)
     $  calday               ! current time in calendar days + fraction.
      character
     $  ncfile*(*)           ! file name of netCDF file containing data
     $, sfncnam(pcnst)*32    ! sfncnam(i) is the name of the variable in the
                             ! netCDF file that contains the surface flux
                             ! data for constituent i.  If there is no data
                             ! for constituent i then set sfncnam(i) = ' '
c Local variables:
      integer
     $  i
     $, istat
     $, recid
     $, vid
c-----------------------------------------------------------------------

      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 1

c     Count requested flux variable names, and init module data
      nsf = 0
      do i = 1, pcnst
         if ( sfncnam(i)(1:1) .ne. ' ' ) then
            nsf = nsf + 1
            sfnam(nsf) = sfncnam(i)
            sfcnst(nsf) = i
         end if
      end do

      if ( nsf .eq. 0 ) return

c     Allocate space for data.
      allocate( sfin(plon,plat,2,nsf), stat=istat )
      call alloc_err( istat, 'iniacsf', 'sfin', plon*plat*2*nsf )
      allocate( sf(plon,plat,nsf), stat=istat )
      call alloc_err( istat, 'iniacsf', 'sf', plon*plat*nsf )

c     Open file.
      call handle_ncerr( nf_open( trim(ncfile), NF_NOWRITE, ncid )
     $, 'iniacsf: error opening file '//trim(ncfile) )

c     Get size of unlimited dimension.
      call handle_ncerr( nf_inq_unlimdim( ncid, recid ), 'iniacsf: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec )
     $,                  'iniacsf: ' )

c     Allocate space for time coordinate data.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, 'iniacsf', 'time', nrec )

c     Get time coordinate.
      call handle_ncerr( nf_inq_varid( ncid, 'time', vid )
     $, 'iniacsf: cannot find time coordinate variable' )
      call handle_ncerr( nf_get_var_double( ncid, vid, time )
     $, 'iniacsf: error getting time coordinate data' )

c     Make sure the time coordinate looks like calander day, and is
c     increasing.
      call chktime( time, nrec )

c     Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

c     Read data.
      loin = 1
      hiin = 2

      do i = 1, nsf

         call handle_ncerr( nf_inq_varid( ncid, sfnam(i), vid )
     $,    'iniacsf: cannot find variable '//sfnam(i) )

         start(3) = lotim
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count
     $,                                         sfin(1,1,loin,i) )
     $,    'iniacsf: cannot read data for '//sfnam(i) )


         start(3) = hitim
         call handle_ncerr( nf_get_vara_double( ncid, vid, start, count
     $,                                         sfin(1,1,hiin,i) )
     $,    'iniacsf: cannot read data for '//sfnam(i) )

         write(*,*)'iniacsf: calendar day = ',calday
         write(*,*)'          read data for days ',time(lotim)
     $,            ' and ',time(hitim)
      end do

      return
      end subroutine iniacsf

c#######################################################################

      subroutine acsfint( calday )

c Interpolate annual cycle surface flux data to the current time.  Read
c in new time samples of the input data as necessary.

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input argument:
      real(r8)
     $  calday  ! current time in calendar days + fraction.

c Local variables:
      integer
     $  i
     $, oldlotim, oldhitim
     $, vid
      real(r8)
     $  dt, dt1, tint
c-----------------------------------------------------------------------

      if ( nsf .eq. 0 ) return

c     Check to see if model time is still bounded by dataset times.
      oldlotim = lotim
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
c        Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1
         start(3) = hitim
         do i = 1, nsf
            call handle_ncerr( nf_inq_varid( ncid, sfnam(i), vid )
     $,       'acsfint: cannot find variable '//sfnam(i) )
            call handle_ncerr( nf_get_vara_double( ncid, vid, start
     $,                                     count, sfin(1,1,hiin,i) )
     $,       'acsfint: cannot read data for '//sfnam(i) )
         end do
         write(*,*)'acsfint: read data for day ',time(hitim)

         if ( lotim .ne. oldhitim ) then
c           Read in new lotim data.  Replace old hitim data.
            start(3) = lotim
            do i = 1, nsf
               call handle_ncerr( nf_inq_varid( ncid, sfnam(i), vid )
     $,          'acsfint: cannot find variable '//sfnam(i) )
               call handle_ncerr( nf_get_vara_double( ncid, vid, start
     $,                                     count, sfin(1,1,loin,i) )
     $,          'acsfint: cannot read data for '//sfnam(i) )
            end do
            write(*,*)'acsfint: read data for day ',time(lotim)
         end if

      end if


c     Linear interpolation...  Start by computing the number of days between
c                              the lower and upper bounds, and days between
c                              the model time and lower bound.

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
      do i = 1, nsf
         call linintp( plon*plat, 0._r8, 1._r8, tint, sfin(1,1,loin,i)
     $,                sfin(1,1,hiin,i), sf(1,1,i) )
      end do

      return
      end subroutine acsfint

c#######################################################################

      subroutine getsf( lat, sflx )

c Return annual cycle surface flux data for the requested latitude.

      implicit none

c Input arguments:
      integer
     $  lat        ! requested latitude index

c Output argument:
      real(r8)
     $  sflx(plond,pcnst)    ! surface fluxes

c Local variables:
      integer
     $  i, n
c-----------------------------------------------------------------------

      if ( nsf .eq. 0 ) return

      do n = 1, nsf
         do i = 1, plon
            sflx(i,sfcnst(n)) = sf(i,lat,n)
         end do
      end do

      return
      end subroutine getsf

c#######################################################################

      end module acsf
