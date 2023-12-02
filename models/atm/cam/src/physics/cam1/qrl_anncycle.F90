#ifdef QRLDAMP

#include <misc.h>
#include <params.h>

module qrl_anncycle

 use shr_kind_mod, only: r8 => shr_kind_r8
 use comsrf
 use pmgrid
 use comspe
 use ppgrid
 use phys_grid

   implicit none

  real(r8), allocatable :: qrl_anncycle_phys(:,:,:,:) 
  real(r8), allocatable :: qrl_anncycle_int(:,:,:)
  real(r8), allocatable :: qrl_dailymean_buffer(:,:,:,:)
  real(r8), allocatable :: qrl_dailymean(:,:,:) 

 logical,public :: need_to_read_qrl

contains

subroutine allocate_qrlanncycle
   allocate (qrl_anncycle_phys(pcols,pver,begchunk:endchunk,365))
   allocate (qrl_anncycle_int(pcols,pver,begchunk:endchunk))
   allocate (qrl_dailymean_buffer(pcols,pver,48,begchunk:endchunk))
   allocate (qrl_dailymean(pcols,pver,begchunk:endchunk))
end subroutine allocate_qrlanncycle

subroutine accumulate_dailymean_qrl (lchnk,ncol,qrl,nstep)

integer, intent(in) :: lchnk
integer, intent(in) :: ncol
real(r8), intent(in) :: qrl(pcols,pver)
integer, intent(in) :: nstep
integer k,i,jj

      do k=1,pver
        do i=1,ncol
          if (nstep .le. 48 .and. nstep .ge. 1) then
            qrl_dailymean_buffer(i,k,nstep,lchnk) = qrl(i,k)
            qrl_dailymean(i,k,lchnk) = qrl(i,k)
          else
            do jj=1,47
               ! cycle back
              qrl_dailymean_buffer(i,k,jj,lchnk) =qrl_dailymean_buffer(i,k,jj+1,lchnk) 
            end do
            qrl_dailymean_buffer(i,k,48,lchnk) = qrl(i,k)
            ! recalculate the running daily mean:
            qrl_dailymean(i,k,lchnk) = 0.
            do jj=1,48
              qrl_dailymean(i,k,lchnk) = qrl_dailymean(i,k,lchnk) +qrl_dailymean_buffer(i,k,jj,lchnk) 
            end do
            qrl_dailymean(i,k,lchnk) = qrl_dailymean(i,k,lchnk)/48.  
          end if
        end do ! icol
    end do ! pver
     
end subroutine accumulate_dailymean_qrl

subroutine ref_qrl_anncycle_int
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
     write (6,*) 'QRL int. -- unexpected exception!'
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
      write (6,*) 'QRL int WTF -- this should never happen.'
      write (6,*) cur_time,iref_right,time_right,iref_left,time_left
      call endrun
    end if
    if (time_right .eq. cur_time) then
      need_to_interpolate = .false.
    end if
  end if

! ==========================================================


  if (need_to_interpolate) then
    do k=1,pver
      do lchnk = begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        do i=1,ncol
          qrl_anncycle_int(i,k,lchnk) = qrl_anncycle_phys(i,k,lchnk,iref_left) + ( qrl_anncycle_phys(i,k,lchnk,iref_right) - qrl_anncycle_phys(i,k,lchnk,iref_left) ) / ( time_right - time_left)*(cur_time - time_left)

!         if (k .eq. 10 .and. lchnk .eq. begchunk .and. i .eq. 1 .and. masterproc) then
!           write (6,*) 'MDEBUG ON INTERPOLATION'
!           write (6,*) 'time_left,cur_time,time_right=',time_left,cur_time,time_right
!           write (6,*) 'iref_left, iref_right=',iref_left,iref_right
!           write (6,*) 'qrl(left), qrl_int, qrl(right)=',qrl_anncycle_phys(i,k,lchnk,iref_left),qrl_anncycle_int(i,k,lchnk),qrl_anncycle_phys(i,k,lchnk,iref_right)
!
!         end if

        end do
      end do 
    end do
  else
!    write (6,*) 'MDEBUG NO INTERPOLATION NEEDED, curtime=',cur_time
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      qrl_anncycle_int(::ncol,::pver,lchnk) = qrl_anncycle_phys(::ncol,::pver,lchnk,iref_right) 
    end do
  end if ! need to interpolate?

end subroutine ref_qrl_anncycle_int

subroutine read_qrl_anncycle

   use ioFileMod,    only: getfil
   use filenames,    only: qrl3danncycle
   use mpishorthand, only: mpir8,mpicom
   use phys_grid,      only: get_ncols_p, get_lat_all_p, get_lon_all_p
   use physconst, only: cpair
 
  include 'netcdf.inc'

   integer lonsiz, latsiz, levsiz       ! Dimension sizes
    integer londimid, levdimid, latdimid ! Dimension ID's
    integer qrlid
  integer ndims2d                      ! number of dimensions
    integer dims2d(NF_MAX_VAR_DIMS)      ! variable shape
    integer ndims3d                      ! number of dimensions
    integer dims3d(NF_MAX_VAR_DIMS)      ! variable shape
       integer natt, ret, attlen            ! netcdf return values
    integer tmptype
    character*(NF_MAX_NAME) tmpname
    real(r8), allocatable ::  arr3d(:,:,:)
    integer :: ncid,iday
    integer strt3d(4)                    ! start lon, lev, lat, time for netcdf 3-d
    data strt3d/1,1,1,1/                     ! Only indices 2,3 will ever change
    integer cnt3d(4)                     ! lon, lat, lev, time counts for netcdf 2-d
    data cnt3d/plon,plat,plev,1/         ! 3-d arrs: Always grab a full time slice
   character(len=256) locfn ! local filename
!

   real(r8), allocatable :: qrl_anncycle_global(:,:,:,:) 
   real(r8),allocatable :: qrl_tmp (:,:,:)
   integer :: i,j,k
   integer lons(pcols)
   integer lats(pcols)
   integer ncol,lchnk,icol




   allocate (qrl_anncycle_global (plon,plev,plat,365))
   allocate (qrl_tmp  (plon,plev,plat) ) ! dim order matches cld_tmp in inidat (template)
   allocate (arr3d(plon,plat,plev))
   if (masterproc) then
      call getfil (qrl3danncycle, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid)
! Get dimension IDs and lengths 
       call wrap_inq_dimid  (ncid, 'lat', latdimid)
       call wrap_inq_dimlen (ncid, latdimid, latsiz)
       call wrap_inq_dimid  (ncid, 'lev', levdimid)
       call wrap_inq_dimlen (ncid, levdimid, levsiz)
       call wrap_inq_dimid  (ncid, 'lon', londimid)
       call wrap_inq_dimlen (ncid, londimid, lonsiz)
       call wrap_inq_varid (ncid, 'QRL'   , qrlid)

       call wrap_inq_var (ncid, qrlid, tmpname, tmptype,ndims3d, dims3d, natt)

       if (dims3d(1).ne.londimid .or. dims3d(2).ne.latdimid .or. &
           dims3d(3).ne.levdimid .or. ndims3d.gt.4) then
              write(6,*)'INIDAT: Bad number of dims or ordering on 3d fld'
              call endrun
       end if

       write (6,*) '======== READING ANNUAL CYCLE OF QRL ========'

      do iday = 1,365
        strt3d(4) = iday
        call wrap_get_vara_realx(ncid, qrlid, strt3d, cnt3d, arr3d)
        do i=1,plon
         do k=1,plev
           do j=1,plat
             qrl_anncycle_global(i,k,j,iday) = cpair*arr3d(i,j,k)
           end do
         end do
        end do
      end do ! iday

!       write (6,*) 'HEY qrl(lon=100,lev=20,lat=33,day=24) = ',qrl_anncycle_global(100,20,33,24)

     end if ! masterproc
#ifdef SPMD
     write (6,*) '==== Broadcasting over MPI ==='
!     write (6,*) 'plev,pver=',plev,pver
     do iday = 1,365
       qrl_tmp (1:plon,1:plev,1:plat) = qrl_anncycle_global(1:plon,1:plev,1:plat,iday)
       call scatter_field_to_chunk(1,plev,1,plon,qrl_tmp,qrl_anncycle_phys(1:pcols,1:pver,begchunk:endchunk,iday))
       if (mod (iday,20) .eq. 0) then
         write (6,*) dble(iday)/dble(365)*100., ' percent done.'
       endif
     end do 
#else
     write (6,*) ' HEY not implemented for non-SPMD'
     call endrun
#endif

   deallocate (qrl_anncycle_global)
   deallocate (qrl_tmp)
   deallocate (arr3d)

! Pritch - manually inspected output below compared to netcdf loaded in matlab
! and confirmed that mpi scatter is working as intended.

!   do lchnk = begchunk,endchunk
!     call get_lat_all_p (lchnk,pcols,lats)
!     call get_lon_all_p (lchnk,pcols,lons) 
!     ncol = get_ncols_p(lchnk)
!     do icol=1,ncol
!       write (6,*) 'MDEBUG (iam=',iam,',lev=20,lon=',lons(icol),',lat=',lats(icol),',day=133):',qrl_anncycle_phys(icol,20,lchnk,133)
!     end do
!   end do

!   write (6,*) ' MDEBUG PAUSE hey made it to end of read_qrl_anncycle!'
!   call endrun

end subroutine read_qrl_anncycle

subroutine qrl_interference (lchnk,ncol,qrl,clat,myqrl,qrldampfac,qrldamp_equatoronly, qrl_dylat, qrl_critlat_deg, qrl_dailymean_interference,qrldamp_freetroponly,qrl_pbot,qrl_ptop,qrl_dp,pmid)

integer, intent(in) :: lchnk,ncol
real (r8), intent(in) :: qrl (pcols,pver)
real(r8),intent(in) :: clat(pcols)  ! latitude in radians.
real (r8), intent(out) :: myqrl (pcols,pver)
real(r8), intent(in) :: qrldampfac
logical, intent(in) :: qrldamp_equatoronly
real(r8), intent(in) :: qrl_dylat
real(r8), intent(in) :: qrl_critlat_deg
logical, intent(in) :: qrl_dailymean_interference
logical, intent(in) :: qrldamp_freetroponly
real(r8), intent(in) :: qrl_pbot
real(r8), intent(in) :: qrl_ptop
real(r8), intent(in) :: qrl_dp
real(r8), intent(in) :: pmid(pcols,pver)

integer :: k,i
real(r8) :: pcenter, pcrit, pp

real (r8) :: eff_qrldampfac, efflat_qrldampfac ! effective modulation factor (includes vertical, meridional restriction of effect.
real (r8) :: curlat_deg

myqrl (:pcols,:pver) = 0.

do i=1,ncol
  ! Determine effective damping based on latitudinal restriction:
  if (qrldamp_equatoronly) then ! equatorial restriction of inteference
    ! Determine latitude:
    curlat_deg = 180./3.14159*clat(i)
    if (curlat_deg .ge. 0. ) then
      efflat_qrldampfac = ( 1. - qrldampfac)/ (1.+exp(-2.e-2*qrl_dylat*(curlat_deg - qrl_critlat_deg))) + qrldampfac
    else
      efflat_qrldampfac = ( 1. - qrldampfac)*(1. - 1./(1. + exp(-2.e-2*qrl_dylat*(curlat_deg + qrl_critlat_deg)))) + qrldampfac 
    end if
  else
    efflat_qrldampfac = qrldampfac  ! if no equatorial restriction
  end if

  do k=1,pver
    if (qrldamp_freetroponly) then ! vertical restriction of interference
      pcenter = 0.5*(qrl_pbot + qrl_ptop)
      pcrit = 0.5*(qrl_pbot - qrl_ptop)
      pp = 1.e-2*pmid(i,k) - pcenter 
      if (pp .ge. 0.) then
        eff_qrldampfac = (1.-efflat_qrldampfac)/(1. + exp(-2.*0.0005*qrl_dp*(pp - pcrit))) + efflat_qrldampfac
      else
        eff_qrldampfac = (1.-efflat_qrldampfac)*(1. - 1./(1. + exp(-2.*0.0005*qrl_dp*(pp + pcrit)))) + efflat_qrldampfac 
      end if
    else
      eff_qrldampfac = efflat_qrldampfac
    end if

!    if (masterproc .and. k .eq. 10 .and. i .eq. 1 .and. lchnk .eq. begchunk) then
!      write (6,*) 'qrl(i,k) = ',qrl(i,k)
!      write (6,*) 'qrl_dailymean(i,k,lchnk) = ',qrl_dailymean(i,k,lchnk)
!      write (6,*) 'qrl_anncycle_int(i,k,lchnk) = ',qrl_anncycle_int(i,k,lchnk)
!    endif
    if (qrl_dailymean_interference) then
      ! Beta modulates daily mean anomalie with respect ot the daily mean annual cycle
      myqrl(i,k) = qrl(i,k) - qrl_dailymean(i,k,lchnk) + qrl_anncycle_int(i,k,lchnk) + eff_qrldampfac*(qrl_dailymean(i,k,lchnk) - qrl_anncycle_int(i,k,lchnk))
    else
      ! Beta modulates instantaneous anomalies with respect to the daily mean
      ! annual cycle. 
      myqrl(i,k) = qrl_anncycle_int(i,k,lchnk) + eff_qrldampfac*(qrl(i,k) - qrl_anncycle_int(i,k,lchnk)) 
    end if
  end do
end do

end subroutine qrl_interference

end module qrl_anncycle

#endif
