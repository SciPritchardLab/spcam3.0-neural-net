#ifdef QVORTDAMP

#include <misc.h>
#include <params.h>

module vorticity_anncycle

 use shr_kind_mod, only: r8 => shr_kind_r8
 use comsrf
 use pmgrid
 use comspe

   implicit none


 real(r8) :: vz_anncycle(psp,plev,365)
 real(r8) :: vz_anncycle_int(psp,plev) ! linearly interpolated to model time
 real(r8) :: vz_dailymean_buffer(psp,plev,48)
 real(r8) :: vz_dailymean (psp,plev)
 integer :: prev_nstep
 logical,public :: need_to_read

contains

subroutine accumulate_daily_mean_vz

use pmgrid, only: iam,masterproc
use comspe, only: nlen, locm, numm
use time_manager, only: get_nstep

integer :: lm,m,mr,mc,n,ir,ii,mlength,k
integer :: nstep

! Note any given processor only uses a subset of
! the global array (in the scope of grcalc where we
! care about accumulating a daily mean). 

! Hnece the subset loop structures below were ripped out of grcalc

    nstep = get_nstep()
!    if (nstep .le. 48) then
     mlength = numm(iam)
    do k=1,plev
     do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
         ir = mc + 2*n - 1
         ii = ir + 1
         call storeit(vz(ir,k),ir,k,nstep)
         call storeit(vz(ii,k),ii,k,nstep)
       end do
       do n=2,nlen(m),2
         ir = mc + 2*n - 1
         ii = ir + 1
!         if (k.eq. 26 .and. lm .eq. 1 .and. n .eq. 2 .and. masterproc) then
!           write (6,*) 'MDEBUG (nstep=',nstep,'): BEFORE:'
!           write (6,*) vz_dailymean(ir,k)
!           write (6,*) vz_dailymean_buffer(ir,k,1:10)       
!         end if
         call storeit(vz(ir,k),ir,k,nstep)
         call storeit(vz(ii,k),ii,k,nstep)
!         if (k .eq. 26 .and. lm .eq. 1 .and. n .eq. 2 .and. masterproc) then
!           write (6,*) 'MDEBUG (nstep=',nstep,'): AFTER:'
!           write (6,*) vz_dailymean(ir,k)
!           write (6,*) vz_dailymean_buffer(ir,k,1:10)       
!         end if
       end do
     end do   ! lm
   end do ! k
 
end subroutine accumulate_daily_mean_vz

subroutine storeit(vzval,ir,k,nstep)
  real(r8), intent (in) :: vzval
  integer, intent (in) :: ir,k,nstep
  integer :: jj

  if (nstep .le. 48 .and. nstep .ge. 1) then ! avoid nstep = 0, invalid subscript.
     vz_dailymean_buffer(ir,k,nstep) = vzval
     vz_dailymean(ir,k) = vzval
!  elseif (nstep .ne. prev_nstep ) then  ! this may be unnecessary.
  else
    ! cycle back:
    do jj=1,47
      vz_dailymean_buffer(ir,k,jj) = vz_dailymean_buffer(ir,k,jj+1)
    end do
    vz_dailymean_buffer(ir,k,48) = vzval ! to make room for the newest value

    ! recalculate the running daily mean:
    vz_dailymean(ir,k) = 0.
    do jj=1,48
      vz_dailymean(ir,k) = vz_dailymean(ir,k) + vz_dailymean_buffer(ir,k,jj)
    end do
    vz_dailymean(ir,k) = vz_dailymean(ir,k)/48.
  end if
!  prev_nstep = nstep 
end subroutine storeit

subroutine read_vorticity_anncycle

   use ioFileMod,    only: getfil
   use filenames,    only: wind3danncycle
   use mpishorthand, only: mpir8,mpicom
 
  include 'netcdf.inc'

   integer lonsiz, latsiz, levsiz       ! Dimension sizes
    integer londimid, levdimid, latdimid ! Dimension ID's
    integer uid, vid               ! Variable ID's
  integer ndims2d                      ! number of dimensions
    integer dims2d(NF_MAX_VAR_DIMS)      ! variable shape
    integer ndims3d                      ! number of dimensions
    integer dims3d(NF_MAX_VAR_DIMS)      ! variable shape
       integer natt, ret, attlen            ! netcdf return values
    integer tmptype
    character*(NF_MAX_NAME) tmpname
!    real(r8) arr4d(plon,plev,plat,365)
!    real(r8) arr4d(plon,plat,plev,365)
    real(r8) arr3d(plon,plat,plev)
    integer :: ncid,iday
    integer strt3d(4)                    ! start lon, lev, lat, time for netcdf 3-d
    data strt3d/4*1/                     ! Only indices 2,3 will ever change
    integer cnt3d(4)                     ! lon, lat, lev, time counts for netcdf 2-d
!    data cnt3d/plon,plev,plat,365/         ! 3-d arrs: Always grab a full time slice
    data cnt3d/plon,plat,plev,1/         ! 3-d arrs: Always grab a full time slice
   character(len=256) locfn ! local filename
!

!   real(r8) u3_anncycle(plond,plev,plat,365) 
!   real(r8) v3_anncycle(plond,plev,plat,365) 
   real(r8) u3_tmp(plond,plev,plat)
   real(r8) v3_tmp(plond,plev,plat)
  real(r8)  vz_tmp(psp,plev) ! vorticity, the same as used in grcalc, I fricking hope. 

       integer :: i,j,k
   if (masterproc) then
      call getfil (wind3danncycle, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid)

! Get dimension IDs and lengths 
       call wrap_inq_dimid  (ncid, 'lat', latdimid)
       call wrap_inq_dimlen (ncid, latdimid, latsiz)
       call wrap_inq_dimid  (ncid, 'lev', levdimid)
       call wrap_inq_dimlen (ncid, levdimid, levsiz)
       call wrap_inq_dimid  (ncid, 'lon', londimid)
       call wrap_inq_dimlen (ncid, londimid, lonsiz)
       call wrap_inq_varid (ncid, 'U'   , uid)
       call wrap_inq_varid (ncid, 'V'   , vid)


       call wrap_inq_var (ncid, uid, tmpname, tmptype,ndims3d, dims3d, natt)

!       if (dims3d(1).ne.londimid .or. dims3d(2).ne.levdimid .or. &
!           dims3d(3).ne.latdimid .or. ndims3d.gt.4) then
       if (dims3d(1).ne.londimid .or. dims3d(2).ne.latdimid .or. &
           dims3d(3).ne.levdimid .or. ndims3d.gt.4) then
              write(6,*)'INIDAT: Bad number of dims or ordering on 3d fld'
              call endrun
       end if

     write (6,*) '======== READING ANNUAL CYCLE OF U,V ========'
     write (6,*) 'and calculating spectral vorticity annual cycle'

     do iday = 1,365
       strt3d(4) = iday  !  A HUGE OMISSION TO HAVE MISSED THIS!!!
       call wrap_get_vara_realx(ncid, uid, strt3d, cnt3d, arr3d)

       do i=1,plon
         do k=1,plev
           do j=1,plat
             u3_tmp(i,k,j) = arr3d(i,j,k)
           end do
         end do
       end do 

       call wrap_get_vara_realx(ncid, vid, strt3d, cnt3d, arr3d)
       do i=1,plon
         do k=1,plev
           do j=1,plat
             v3_tmp(i,k,j) = arr3d(i,j,k)
           end do
         end do
       end do 
!       write (6,*) 'HEY u3(lon=100,lev=20,lat=33,day=24) = ',u3_anncycle(100,20,33,24)
       
!write (6,*) 'HEY v3(lon=100,lev=20,lat=33,day=24) = ',v3_anncycle(100,20,33,24)

       call vorticity_spetru (u3_tmp,v3_tmp,vz_tmp)
       vz_anncycle(:psp,:plev,iday) = vz_tmp(:psp,:plev)
       if (mod (iday,20) .eq. 0) then
         write (6,*) dble(iday)/dble(365)*100., ' percent done.'
       endif
     end do ! iday
     end if ! masterproc
#ifdef SPMD
     write (6,*) '==== Broadcasting over MPI ==='
     do iday = 1,365
       vz_tmp (:psp,:plev) = vz_anncycle(:psp,:plev,iday)
       call mpibcast (vz_tmp,psp*plev,mpir8,0,mpicom)
       vz_anncycle(:psp,:plev,iday) = vz_tmp(:psp,:plev)
       if (mod (iday,20) .eq. 0) then
         write (6,*) dble(iday)/dble(365)*100., ' percent done.'
       endif
     end do 
!    call mpibcast   (vz_anncycle   ,psp*plev*365,mpir8,0 , mpicom)
!    above was too big message size
#endif

end subroutine read_vorticity_anncycle

subroutine vorticity_spetru (mmyu3, mmyv3, mmyvz)
 
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use constituents, only: pcnst, pnats
  use pspect
  use comspe
  use rgrid,        only: nlon, nmmax
  use commap,       only: w, xm, rsq, cs
  use dynconst,     only: ez, ra, rearth

#include <comctl.h>
#include <comfft.h>  

 real(r8), intent(in):: mmyu3   (plond,plev,plat) ! Fourier -> spec. coeffs. for u-wind
 real(r8), intent(in):: mmyv3   (plond,plev,plat) ! Fourier -> spec. coeffs. for v-wind
 real(r8), intent(out)::mmyvz(psp,plev) ! vorticity, the same as used in grcalc, I fricking hope. 

#if ( ! defined USEFFTLIB )
  real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
  real(r8) work((plon+1)*pcray) ! Workspace for fft
#endif

  real(r8) alpn (pspt)          ! alp*rsq*xm*ra
  real(r8) dalpn(pspt)          ! dalp*rsq*ra
  real(r8) tmp1                 ! vector temporary
  real(r8) tmp2                 ! vector temporary
  real(r8) tmpr                 ! vector temporary (real)
  real(r8) tmpi                 ! vector temporary (imaginary)
  real(r8) phialpr,phialpi      ! phi*alp (real and imaginary)
  real(r8) phdalpr,phdalpi      ! phi*dalp (real and imaginary)
  real(r8) zwalp                ! zw*alp
  real(r8) zwdalp               ! zw*dalp
  real(r8) psdalpr,psdalpi      ! alps (real and imaginary)*dalp
  real(r8) psalpr,psalpi        ! alps (real and imaginary)*alp
  real(r8) zrcsj                ! ra/(cos**2 latitude)
  real(r8) zw                   ! w**2
  real(r8) filtlim              ! filter function
  real(r8) ft                   ! filter multiplier for spectral coefficients
  real (r8) :: locv3 (plond,plev,plat)
  real (r8) :: locu3 (plond,plev,plat)
  real(r8) zsqcs

  integer ir,ii                 ! indices complex coeffs. of spec. arrs.
  integer i,k                   ! longitude, level indices
  integer irow                  ! latitude pair index
  integer latm,latp             ! symmetric latitude indices
  integer lat                   ! index
  integer m                     ! longitudinal wavenumber index (non-PVP)
  integer n                     ! latitudinal wavenumber index (non-PVP)
  integer nspec                 ! index
  integer mr,mc                 ! spectral indices



#if (defined PVP )
   write (6,*) 'HEY THIS ISNT IMPLEMENTED HERE (Pritch)'
   stop
#endif


! Zero spectral arrays
  mmyvz  (:psp,:plev) = 0. 

! Compute the quantities which are transformed to spectral space:
!   1. u = u*sqrt(1-mu*mu),   u * cos(phi)
!   2. v = v*sqrt(1-mu*mu),   v * cos(phi)
  do lat=1,plat
     irow = lat
     if (lat.gt.plat/2) irow = plat - lat + 1
     zsqcs = sqrt(cs(irow))
     do k=1,plev
        do i=1,nlon(lat)
           locu3(i,k,lat) = mmyu3(i,k,lat)*zsqcs
           locv3(i,k,lat) = mmyv3(i,k,lat)*zsqcs
        end do
     end do
!
! Transform grid -> fourier
!
     call fft991 (locu3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond   ,nlon(lat),plev        ,-1         )
     call fft991 (locv3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,plev        ,-1         )
   end do                    ! lat=1,plat

!
! Loop over latitude pairs
!
  do irow=1,plat/2
     latp = irow
     latm = plat - irow + 1
     zrcsj = ra/cs(irow)
     zw = w(irow)*2.
!
! Compute symmetric and antisymmetric components
     do k=1,plev
        do i=1,2*nmmax(irow)

           tmp1 = 0.5*(locu3(i,k,latm) - locu3(i,k,latp))
           tmp2 = 0.5*(locu3(i,k,latm) + locu3(i,k,latp))
           locu3(i,k,latm) = tmp1
           locu3(i,k,latp) = tmp2

           tmp1 = 0.5*(locv3(i,k,latm) - locv3(i,k,latp))
           tmp2 = 0.5*(locv3(i,k,latm) + locv3(i,k,latp))
           locv3(i,k,latm) = tmp1
           locv3(i,k,latp) = tmp2
        end do
     end do


!     
! Compute vzmn
!
     do k=1,plev
        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=1,nlen(m),2
              zwdalp   = zw*dalp(mr+n,irow)
              zwalp    = zw*alp (mr+n,irow)
              ir       = mc + 2*n - 1
              ii       = ir + 1
              mmyvz(ir,k) = mmyvz(ir,k) + (zwdalp*locu3(2*m-1,k,latm) - &
                         xm(m)*zwalp*locv3(2*m  ,k,latp))*zrcsj
              mmyvz(ii,k) = mmyvz(ii,k) + (zwdalp*locu3(2*m  ,k,latm) + &
                         xm(m)*zwalp*locv3(2*m-1,k,latp))*zrcsj
           end do
        end do

        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=2,nlen(m),2
             zwdalp   = zw*dalp(mr+n,irow)
              zwalp    = zw*alp (mr+n,irow)
              ir       = mc + 2*n - 1
              ii       = ir + 1
              mmyvz(ir,k) = mmyvz(ir,k) + (zwdalp*locu3(2*m-1,k,latp) - &
                         xm(m)*zwalp*locv3(2*m  ,k,latm))*zrcsj
              mmyvz(ii,k) = mmyvz(ii,k) + (zwdalp*locu3(2*m  ,k,latp) + &
                         xm(m)*zwalp*locv3(2*m-1,k,latm))*zrcsj
           end do
        end do
    end do
  end do

end subroutine vorticity_spetru

subroutine ref_vorticity_int
  use time_manager, only: get_curr_date, get_curr_calday, get_nstep
  use comspe, only: nlen, locm, numm
  use pmgrid, only: masterproc

  real(r8) :: ref_day,cur_day 

  integer :: iref_right,iref_left
  real(r8) :: time_right,time_left,cur_time
  real(r8) :: ref_times(365)
  logical :: unfound,need_to_interpolate
  integer :: mlength,k,lm,m,mr,mc,n,ir,ii,iref_day


  iref_right = 0

! We want to avoid introducing jumpiness at a daily frequency so rather than 
! just pick out the current calendar from the reference data we want to time
! interpolate between successive calendar days. Assume a 0-1 calendar day and
! that the reference data exists at 0.5. 
  cur_time = get_curr_calday()  - 1.! note this is a decimal date. Subtracting
! the one makes it go from 0 to 364 as the logic below expects.  

  do iref_day = 1,365
    ref_times(iref_day) = dble(iref_day) - 0.5
  end do

  unfound = .true.
  do iref_day = 1,365
    ! find the rightmost day
    if (unfound .and. ref_times(iref_day) .ge. cur_time) then
      iref_right = iref_day
      unfound = .false. 
    end if
  end do

  if (unfound) then
     write (6,*) 'Vorticity int. -- unexpected exception!'
  end if

  need_to_interpolate = .true.
! Special case of periodic points, when cur_time > 364.5 or cur_time < 0.5
  if (cur_time .ge. 364.5) then
    ! to handle periodicity, leave cur_time as is but 
    ! Add ghost times either greater than 365 (max value)
    time_right = 365.5
    iref_right = 1
    time_left = 364.5
    iref_left = 365
  elseif (cur_time .le. 0.5) then
    ! ... or less than 0.5 (min value)
    time_right = 0.5
    iref_right = 1
    time_left = -0.5
    iref_left = 365 
  else
    time_right = ref_times(iref_right)
    iref_left = iref_right - 1
    time_left = ref_times(iref_left)
    iref_left = iref_right - 1
    if (time_left .ge. cur_time) then
      write (6,*) 'Vorticitiy int WTF -- this should never happen.'
      write (6,*) cur_time,iref_right,time_right,iref_left,time_left
      call endrun
    end if
    if (time_right .eq. cur_time) then
      need_to_interpolate = .false.
    end if
  end if

! ==========================================================

  if (need_to_interpolate) then
!vz_anncycle(psp,plev,365)
!vz_anncycle_int(psp,plev.365)
! As for accumualted model internal mean, only have to work on subset of the
! array:
    mlength = numm(iam)
    do k=1,plev
     do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
         ir = mc + 2*n - 1
         ii = ir + 1
! Do the linear interpolation:
         vz_anncycle_int(ir,k) = vz_anncycle(ir,k,iref_left) + ( vz_anncycle(ir,k,iref_right) - vz_anncycle(ir,k,iref_left) )/( time_right - time_left)*(cur_time - time_left)
         vz_anncycle_int(ii,k) = vz_anncycle(ii,k,iref_left) + ( vz_anncycle(ii,k,iref_right) - vz_anncycle(ii,k,iref_left) )/( time_right - time_left)*(cur_time - time_left)
       end do
       do n=2,nlen(m),2
         ir = mc + 2*n - 1
         ii = ir + 1
         vz_anncycle_int(ir,k) = vz_anncycle(ir,k,iref_left) + ( vz_anncycle(ir,k,iref_right) - vz_anncycle(ir,k,iref_left) )/( time_right - time_left)*(cur_time - time_left)
         vz_anncycle_int(ii,k) = vz_anncycle(ii,k,iref_left) + ( vz_anncycle(ii,k,iref_right) - vz_anncycle(ii,k,iref_left) )/( time_right - time_left)*(cur_time - time_left)

!         if (k .eq. 10 .and. lm .eq. 1 .and. n .eq. 2 .and. masterproc) then
!           write (6,*) 'MDEBUG ON INTERPOLATION'
!           write (6,*) 'time_left, cur_time,time_right=',time_left,cur_time,time_right
!           write (6,*) 'iref_left, iref_right=',iref_left,iref_right
!           write (6,*) 'vz(left), vz_int, vz(right) =',vz_anncycle(ii,k,iref_left),vz_anncycle_int(ii,k),vz_anncycle(ii,k,iref_right)
!         end if
       end do
     end do   ! lm
   end do ! k
     

  else
!    write (6,*) 'MDEBUG NO INTERPOLATION NEEDED, curtime=',cur_time
    vz_anncycle_int(:psp,:plev) = vz_anncycle(:psp,:plev,iref_right) 
  end if ! need to interpolate?

end subroutine ref_vorticity_int

end module vorticity_anncycle

#endif
