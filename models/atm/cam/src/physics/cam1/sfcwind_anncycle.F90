#ifdef FLUXDAMP

#include <misc.h>
#include <params.h>

module sfcwind_anncycle

 use shr_kind_mod, only: r8 => shr_kind_r8
 use comsrf
 use pmgrid
 use comspe
 use ppgrid
 use phys_grid

   implicit none

  real(r8), allocatable :: ubot_anncycle_phys(:,:,:) 
  real(r8), allocatable :: ubot_anncycle_int(:,:)
  real(r8), allocatable :: ubot_dailymean_buffer(:,:,:)
  real(r8), allocatable :: ubot_dailymean(:,:) 

  real(r8), allocatable :: vbot_anncycle_phys(:,:,:) 
  real(r8), allocatable :: vbot_anncycle_int(:,:)
  real(r8), allocatable :: vbot_dailymean_buffer(:,:,:)
  real(r8), allocatable :: vbot_dailymean(:,:) 

 logical,public :: need_to_read_sfcwind

contains

subroutine allocate_sfcwindanncycle
   allocate (ubot_anncycle_phys(pcols,begchunk:endchunk,365))
   allocate (ubot_anncycle_int(pcols,begchunk:endchunk))
   allocate (ubot_dailymean_buffer(pcols,48,begchunk:endchunk))
   allocate (ubot_dailymean(pcols,begchunk:endchunk))

   allocate (vbot_anncycle_phys(pcols,begchunk:endchunk,365))
   allocate (vbot_anncycle_int(pcols,begchunk:endchunk))
   allocate (vbot_dailymean_buffer(pcols,48,begchunk:endchunk))
   allocate (vbot_dailymean(pcols,begchunk:endchunk))

end subroutine allocate_sfcwindanncycle

subroutine ref_sfcwindanncycle_int
   use time_manager, only: get_curr_date, get_curr_calday, get_nstep
   use phys_grid,      only: get_ncols_p

  real(r8) :: ref_day,cur_day
  ! Note this is copied right out of vorticity_anncycle.mod
  integer :: iref_right,iref_left
  real(r8) :: time_right,time_left,cur_time
  real(r8) :: ref_times(365)
  logical :: unfound,need_to_interpolate
  integer lchnk, ncol, k, i, iref_day

  iref_right = 0
  cur_time = get_curr_calday()  - 1.! note this is a decimal date. Subtracting

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
     write (6,*) 'SFCWIND int. -- unexpected exception!'
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
      write (6,*) 'SFCWIND int WTF -- this should never happen.'
      write (6,*) cur_time,iref_right,time_right,iref_left,time_left
      call endrun
    end if
    if (time_right .eq. cur_time) then
      need_to_interpolate = .false.
    end if
  end if

! ==========================================================


  if (need_to_interpolate) then
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do i=1,ncol
        ubot_anncycle_int(i,lchnk) = ubot_anncycle_phys(i,lchnk,iref_left) + ( ubot_anncycle_phys(i,lchnk,iref_right) - ubot_anncycle_phys(i,lchnk,iref_left) ) / ( time_right - time_left)*(cur_time - time_left)
        vbot_anncycle_int(i,lchnk) = vbot_anncycle_phys(i,lchnk,iref_left) + ( vbot_anncycle_phys(i,lchnk,iref_right) - vbot_anncycle_phys(i,lchnk,iref_left) ) / ( time_right - time_left)*(cur_time - time_left)
      end do
    end do 
  else
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      ubot_anncycle_int(::ncol,lchnk) = ubot_anncycle_phys(::ncol,lchnk,iref_right) 
      vbot_anncycle_int(::ncol,lchnk) = vbot_anncycle_phys(::ncol,lchnk,iref_right) 
    end do
  end if ! need to interpolate?

end subroutine ref_sfcwindanncycle_int

subroutine read_sfcwindanncycle

   use ioFileMod,    only: getfil
   use filenames,    only: wind3danncycle
   use mpishorthand, only: mpir8,mpicom
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physconst, only: cpair
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
    integer :: ncid,iday
    integer strt3d(4)                    ! start lon, lev, lat, time for netcdf
    data strt3d/4*1/                     ! Only indices 2,3 will ever change
    integer cnt3d(4)                     ! lon, lat, lev, time counts for netcdf
    data cnt3d/plon,plat,1,1/         ! 3-d arrs: Always grab a full time slice
   character(len=256) locfn ! local filename
   real(r8), allocatable :: ubot_tmp (:,:), vbot_tmp (:,:), arr2d (:,:)
   real(r8), allocatable :: ubot_anncycle_global(:,:,:) 
   real(r8), allocatable :: vbot_anncycle_global(:,:,:) 
   integer ncol,lchnk,icol
     integer :: i,j,k
   integer lons(pcols)
   integer lats(pcols)

   allocate (ubot_anncycle_global (plon,plat,365))
   allocate (ubot_tmp  (plon,plat) ) ! dim order matches cld_tmp in inidat (template)
   allocate (vbot_anncycle_global (plon,plat,365))
   allocate (vbot_tmp  (plon,plat) ) ! dim order matches cld_tmp in inidat (template)
   allocate (arr2d(plon,plat))
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

      if (dims3d(1).ne.londimid .or. dims3d(2).ne.latdimid .or. &
        dims3d(3).ne.levdimid .or. ndims3d.gt.4) then
        write(6,*)'WINDANNCYCLE: Bad number of dims or ordering on 3d fld'
        call endrun
      end if

      write (6,*) '======== READING ANNUAL CYCLE OF U,V ========'
      write (6,*) 'and calculating lowest model level wind annual cycle'

      do iday = 1,365
        strt3d(4) = iday
        strt3d(3) = plev ! lowest model level
        call wrap_get_vara_realx(ncid, uid, strt3d, cnt3d, ubot_tmp)
        call wrap_get_vara_realx(ncid, vid, strt3d, cnt3d, vbot_tmp)
        do i=1,plon
          do j=1,plat
            ubot_anncycle_global(i,j,iday) = ubot_tmp(i,j)
            vbot_anncycle_global(i,j,iday) = vbot_tmp(i,j)
          end do
        end do 
      end do
     end if ! masterproc
#ifdef SPMD
     write (6,*) '==== Broadcasting over MPI ==='
     do iday = 1,365
       ubot_tmp (1:plon,1:plat) = ubot_anncycle_global(1:plon,1:plat,iday)
       vbot_tmp (1:plon,1:plat) = vbot_anncycle_global(1:plon,1:plat,iday)
       call scatter_field_to_chunk(1,1,1,plon,ubot_tmp,ubot_anncycle_phys(1:pcols,begchunk:endchunk,iday))
       call scatter_field_to_chunk(1,1,1,plon,vbot_tmp,vbot_anncycle_phys(1:pcols,begchunk:endchunk,iday))
       if (mod (iday,20) .eq. 0) then
         write (6,*) dble(iday)/dble(365)*100., ' percent done.'
       endif
     end do 
#else
     write (6,*) ' HEY not implemented for non-SPMD'
     call endrun
#endif

   deallocate (ubot_anncycle_global)
   deallocate (ubot_tmp)
   deallocate (arr2d)
   deallocate (vbot_anncycle_global)
   deallocate (vbot_tmp)


! Pritch - manually inspected output below compared to netcdf loaded in matlab
! and confirmed that mpi scatter is working as intended.
! (NOTE: DID THIS INDEPENDENTLY FOR SFCWIND BOTH COMPONENTS)
!   do lchnk = begchunk,endchunk
!     call get_lat_all_p (lchnk,pcols,lats)
!     call get_lon_all_p (lchnk,pcols,lons) 
!     ncol = get_ncols_p(lchnk)
!     do icol=1,ncol
!       write (6,*) 'MDEBUG ubot (iam=',iam,'lon=',lons(icol),',lat=',lats(icol),',day=133):',ubot_anncycle_phys(icol,lchnk,133)
!       write (6,*) 'MDEBUG vbot (iam=',iam,'lon=',lons(icol),',lat=',lats(icol),',day=133):',vbot_anncycle_phys(icol,lchnk,133)
!     end do
!   end do
!   call endrun

end subroutine read_sfcwindanncycle

subroutine sfcwind_interference (lchnk,ncol,ubot,vbot,clat,fluxdampfac,fluxdamp_equatoronly, flux_dylat, flux_critlat_deg, dubot, dvbot, dwindbot, ztodt)

integer, intent(in) :: lchnk,ncol
real (r8), intent(inout) :: ubot (pcols)
real (r8), intent(inout) :: vbot (pcols)
real(r8),intent(in) :: clat(pcols)  ! latitude in radians.
real(r8), intent(in) :: fluxdampfac
logical, intent(in) :: fluxdamp_equatoronly
real(r8), intent(in) :: flux_dylat
real(r8), intent(in) :: flux_critlat_deg
real (r8), intent(out) :: dubot(pcols), dvbot(pcols), dwindbot(pcols)
real(r8), intent(in) :: ztodt

real (r8) :: eff_fluxdampfac, curlat_deg
real(r8) :: ubot_orig (pcols), vbot_orig (pcols)
integer :: i

do i=1,ncol
  ! Determine effective damping based on latitudinal restriction:
  if (fluxdamp_equatoronly) then ! equatorial restriction of inteference
    ! Determine latitude:
    curlat_deg = 180./3.14159*clat(i)
    if (curlat_deg .ge. 0. ) then
      eff_fluxdampfac = ( 1. - fluxdampfac)/(1.+exp(-2.e-2*flux_dylat*(curlat_deg - flux_critlat_deg))) + fluxdampfac
    else
      eff_fluxdampfac = ( 1. - fluxdampfac)*(1. - 1./(1. + exp(-2.e-2*flux_dylat*(curlat_deg + flux_critlat_deg)))) + fluxdampfac
    end if
  else
    eff_fluxdampfac = fluxdampfac  ! if no equatorial restriction
  end if

  ubot_orig(i) = ubot(i)
  vbot_orig(i) = vbot(i)

  ! Gamma modulates instantaneous anomalies of surface wind vector with respect
  ! to the mean annual cycle:
  ubot(i) = ubot_anncycle_int(i,lchnk) + eff_fluxdampfac*(ubot_orig(i) - ubot_anncycle_int(i,lchnk)) 
  vbot(i) = vbot_anncycle_int(i,lchnk) + eff_fluxdampfac*(vbot_orig(i) - vbot_anncycle_int(i,lchnk)) 
  dubot(i) = (ubot(i) - ubot_orig(i))/ztodt
  dvbot(i) = (vbot(i) - vbot_orig(i))/ztodt
  dwindbot (i) = (sqrt(ubot(i)**2 + vbot(i)**2) - sqrt(ubot_orig(i)**2 + vbot_orig(i)**2 ))
end do ! i

end subroutine sfcwind_interference

end module sfcwind_anncycle

#endif
