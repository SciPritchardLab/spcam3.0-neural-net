      module volcemission

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! This code does time interpolation for volcano emissions
   !
   ! 
   ! Author: B. Eaton 
   ! Modified by P. Rasch to:
   !  1) fortran 90
   !  2) use F90 netcdf interface
   !----------------------------------------------------------------------- 

#include <params.h>

      use shr_kind_mod, only: r8 => shr_kind_r8
      use  pmgrid

!nf90      use netcdf

      implicit none

c $Id: volcemission.F,v 1.1.6.2 2003/08/06 17:48:51 rosinski Exp $

#include <netcdf.inc>

      real(r8) so2(plon,plat)

      integer  hgtidx(plon,plat)
      integer ncid    ! ID for netCDF file

      contains

c###############################################################################

      subroutine volcemisini()

c Open netCDF file containing data for SO2 emissions from non-eruptive volcanoes.
      use error_messages, only: handle_ncerr

      implicit none

c Local variables:
      integer
     $  did
     $, nlon
     $, vid
c-----------------------------------------------------------------------

c     Open file.
!nf90      call handle_ncerr( 
!nf90     $  nf90_open( 'neVolc_emissions.nc', NF90_NOWRITE, ncid )
!nf90     $, 'volcemisini: error opening file '//'neVolc_emissions.nc' )
      call handle_ncerr( 
     $  nf_open( 'neVolc_emissions.nc', NF_NOWRITE, ncid )
     $, 'volcemisini: error opening file '//'neVolc_emissions.nc' )

      ! get the record id
      !      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid),      
      !     $     'volcemisini: no record variables ' )

c     Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did )
!nf90     $,                  'volcemisini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon )
!nf90     $,                  'volcemisini: ' )

      call handle_ncerr( nf_inq_dimid( ncid, 'lon', did )
     $,                  'volcemisini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, did, nlon )
     $,                  'volcemisini: ' )
      if ( nlon .ne. plon ) then
         write(*,*)'volcemisini: change plon to ',nlon
         stop
      end if

!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'SO2', vid )
!nf90     $,    'volcemisini: cannot find variable '//'SO2' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, so2 )
!nf90     $,    'volcemisini: cannot read data for '//'SO2' )
!nf90
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'MXVOLHGTIDX', vid )
!nf90     $,    'volcemisini: cannot find variable '//'MXVOLHGTIDX' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, hgtidx )
!nf90     $,    'volcemisini: cannot read data for '//'MXVOLHGTIDX' )

      call handle_ncerr( nf_inq_varid( ncid, 'SO2', vid )
     $,    'volcemisini: cannot find variable '//'SO2' )
      call handle_ncerr( nf_get_var_double( ncid, vid, so2 )
     $,    'volcemisini: cannot read data for '//'SO2' )

      call handle_ncerr( nf_inq_varid( ncid, 'MXVOLHGTIDX', vid )
     $,    'volcemisini: cannot find variable '//'MXVOLHGTIDX' )
      call handle_ncerr( nf_get_var_int( ncid, vid, hgtidx )
     $,    'volcemisini: cannot read data for '//'MXVOLHGTIDX' )

      write(*,*)'volcemisini: read data.'

      return
      end subroutine volcemisini

c#######################################################################

      subroutine volcemist( lat, gravit, rpdel, so2tend )

c Calculate SO2 tendency due to emissions from non-eruptive volcanoes.

c N.B. Time constraints dictated that this routine be hard-wired assuming
c      the vertical grid to be the 28 levels from the NCEP reanalysis.

      implicit none

c Input arguments:
      integer
     $  lat                  ! model latitude index
      real(r8)
     $  gravit
     $, rpdel(plond,plev)    ! reciprocal of pdel

c Output argument:
      real(r8)
     $  so2tend(plond,plev)  ! tracer tendency due to volcano emissions

c Local variables:
      integer i, k
     $, hgtidx(plond)
      real(r8)
     $  masfac               ! mass scale factor
     $, so2(plond)
     $, emis(plond,plev)     ! emissions distributed in vertical
c-----------------------------------------------------------------------

      if (plev.ne.28) then
         write (6,*) ' volcemis: configured only for 28 level model '
         call abort
      endif

c     Input emissions are scaled from 1 Tg S/yr to 8 Tg S/yr.
c     The factor of 2 scales from kg S to kg SO2.
      masfac = 8. * 2.

c     Column emissions were divided as follows by Ulrike Lohmann
C   3 EMISSION LEVELS:
C    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO) 36%
C    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)         36%
C    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)        28%

c     For the current model, the volcano height index is contained in the
c     data file, post-extra eruption levels are 21-25, and 
c     explosive eruption levels are 14-16

      call getvolcemis( lat, plond, so2, hgtidx )

      do k = 1, plev
         do i = 1, plon
            emis(i,k) = 0.0
         end do
      end do

c     Test
c      do i = 1, plon
c         emis(i,plev-1) = so2(i) * masfac*.5
c         emis(i,plev) = so2(i) * masfac*.5
c      end do

c     Add emissions at volcano height index.
      do i = 1, plon
         if ( hgtidx(i) .gt. 0 ) then
            emis(i,hgtidx(i)) = so2(i) * masfac * .36
         end if
      end do

c     Add emissions at post-extra eruption levels.  The factor of 5 is for
c     splitting the 36% up among 5 levels.
      do k = 21, 25
         do i = 1, plon
            emis(i,k) = emis(i,k) + so2(i) * masfac * .36 / 5.
         end do
      end do

c     Add emissions at explosive eruption levels.  The factor of 3 is for
c     splitting the 28% up among 3 levels.
      do k = 14, 16
         do i = 1, plon
            emis(i,k) = emis(i,k) + so2(i) * masfac * .28 / 3.
         end do
      end do

c     Scale emissions from kg/m2/s to kg/kg/s
      do k = 1, plev
         do i = 1, plon
            so2tend(i,k) = emis(i,k) * gravit*rpdel(i,k)
         end do
      end do

      return
      end subroutine volcemist

c#######################################################################

      subroutine getvolcemis( lat, nd1, so2lat, hgtidxlat )

c Return SO2 emission data for the requested latitude.

      implicit none

c Input arguments:
      integer
     $  lat        ! requested latitude index
     $, nd1        ! 1st dimension of output array

c Output argument:
      real(r8)
     $  so2lat(nd1)   ! SO2 emissions in kg S/m2/s
      integer
     $  hgtidxlat(nd1)! Index for layer that contains the top of 
c                     ! the volcano.

c Local variables:
      integer
     $  i
c-----------------------------------------------------------------------

      do i = 1, plon
         so2lat(i) = so2(i,lat)
         hgtidxlat(i) = hgtidx(i,lat)
      end do

      return
      end subroutine getvolcemis

c#######################################################################

      end module volcemission
