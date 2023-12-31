!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module cloudbrain_output_limiter

!----------------------------------------------------------------------- 
! 
! Author: Sungduk Yu
! Date: Thu Oct 13 21:34:27 EDT 2022
!
! Purpose: 
! Limit neural_net outputs by preset thresholds (upper and lower bonds).
! The thresholds are 3-D vars and provided by a user 
! via "cloudbrain_output_limiter.nc" that contains:
!   - "up_PHQ": upper bounds of PHQ
!   - "lo_PHQ": lower bounds of PHQ
!   - "up_TPHYSTND": upper bounds of TPHYSTND
!   - "lo_TPHYSTND": lower bounds of TPHYSTND
! The order of the dimensions of these variables  should be (lon, lev, lat).
!
! To turn on this limiter, compile with "-DCBLIMITER".
!!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plev, plon, plat, masterproc

   implicit none

   character(len=128) :: cblimiter_filename='cloudbrain_output_limiter.nc'
                         ! note that the dimension in this file should be (lat,lev,lon)
                         ! Fortran and netcdf library has a different object ordering
   real (r8), dimension(plon, plev, plat) ::  up_PHQ, up_TPHYSTND,& ! Upper bound
                                              lo_PHQ, lo_TPHYSTND   ! Lower lowerbound

contains

subroutine init_cb_limiter

   implicit none 

   include 'netcdf.inc'

   ! Local variables:
   integer :: status_RES, ncid_RES, varid_RES
   integer :: id_lon, id_lev, id_lat
   
   status_RES = nf_open(cblimiter_filename,nf_nowrite,ncid_RES)
   status_RES = nf_inq_varid(ncid_RES,"up_PHQ",varid_RES)
   status_RES = nf_get_var(ncid_RES,varid_RES,up_PHQ)
   status_RES = nf_inq_varid(ncid_RES,"up_TPHYSTND",varid_RES)
   status_RES = nf_get_var(ncid_RES,varid_RES,up_TPHYSTND)
   status_RES = nf_inq_varid(ncid_RES,"lo_PHQ",varid_RES)
   status_RES = nf_get_var(ncid_RES,varid_RES,lo_PHQ)
   status_RES = nf_inq_varid(ncid_RES,"lo_TPHYSTND",varid_RES)
   status_RES = nf_get_var(ncid_RES,varid_RES,lo_TPHYSTND)
   status_RES = nf_inq_dimid(ncid_RES,'lon',id_lon)
   status_RES = nf_inq_dimid(ncid_RES,'lev',id_lev)
   status_RES = nf_inq_dimid(ncid_RES,'lat',id_lat)
   status_RES = nf_close(ncid_RES)
   write (6,*) '[CBLIMITER] Thresholds are loaded from a netcdf file: ', cblimiter_filename
   
   ! Debugging (will be moved under BRAINDEBUG)
   if (masterproc) then
      write (6,*) '[CBLIMITER] dimid(lon,lev,lat) = ', id_lon, id_lev, id_lat
      write (6,*) '[CBLIMITER] shape(lo_PHQ) =  ', shape(lo_PHQ)
      write (6,*) '[CBLIMITER] shape(up_PHQ) =  ', shape(up_PHQ)
      write (6,*) '[CBLIMITER] lo_PHQ(10,19,40)', lo_PHQ(10,19,40)
      write (6,*) '[CBLIMITER] lo_PHQ(50,19,40)', lo_PHQ(50,19,40)
      write (6,*) '[CBLIMITER] lo_PHQ(10,25,40)', lo_PHQ(10,25,40)
      write (6,*) '[CBLIMITER] lo_PHQ(10,19,60)', lo_PHQ(10,19,60)
      write (6,*) '[CBLIMITER] up_PHQ(10,19,40)', up_PHQ(10,19,40)
      write (6,*) '[CBLIMITER] up_PHQ(50,19,40)', up_PHQ(50,19,40)
      write (6,*) '[CBLIMITER] up_PHQ(10,25,40)', up_PHQ(10,25,40)
      write (6,*) '[CBLIMITER] up_PHQ(10,19,60)', up_PHQ(10,19,60)
      write (6,*) '[CBLIMITER] lo_TPHYSTND(10,19,40)', lo_TPHYSTND(10,19,40)
      write (6,*) '[CBLIMITER] lo_TPHYSPHQ(50,19,40)', lo_TPHYSTND(50,19,40)
      write (6,*) '[CBLIMITER] lo_TPHYSPHQ(10,25,40)', lo_TPHYSTND(10,25,40)
      write (6,*) '[CBLIMITER] lo_TPHYSPHQ(10,19,60)', lo_TPHYSTND(10,19,60)
      write (6,*) '[CBLIMITER] up_TPHYSTND(10,19,40)', up_TPHYSTND(10,19,40)
      write (6,*) '[CBLIMITER] up_TPHYSPHQ(50,19,40)', up_TPHYSTND(50,19,40)
      write (6,*) '[CBLIMITER] up_TPHYSPHQ(10,25,40)', up_TPHYSTND(10,25,40)
      write (6,*) '[CBLIMITER] up_TPHYSPHQ(10,19,60)', up_TPHYSTND(10,19,60)
   end if

end subroutine init_cb_limiter

subroutine cb_limiter (ptend,lchnk,ncol)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use physics_types,      only: physics_ptend
   use physconst,       only: cpair
   use phys_grid,    only: get_lat_p, get_lon_p
   use ppgrid, only: pver

   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
     type(physics_ptend), intent(inout)   :: ptend                 ! indivdual parameterization tendencies
     integer, intent(in) :: lchnk            ! chunk index
     integer, intent(in) :: ncol

   ! Local variables:
     integer :: idx_lon,idx_lat, i
     real(r8) :: work1(1:pver), work2(1:pver)
   !-----------------------------------------------------------------------

   do i=1,ncol ! this is the loop over independent GCM columns.
      idx_lat = get_lat_p(lchnk, i)
      idx_lon = get_lon_p(lchnk, i)
      
      work1 = lo_PHQ(idx_lon,1:pver,idx_lat)
      work2 = up_PHQ(idx_lon,1:pver,idx_lat)
      ptend%q(i,1:pver,1) = merge(ptend%q(i,1:pver,1), work1, ptend%q(i,1:pver,1).ge.work1)
      ptend%q(i,1:pver,1) = merge(ptend%q(i,1:pver,1), work2, ptend%q(i,1:pver,1).le.work2)

      work1 = lo_TPHYSTND(idx_lon,1:pver,idx_lat)*cpair
      work2 = up_TPHYSTND(idx_lon,1:pver,idx_lat)*cpair
      ptend%s(i,1:pver) = merge(ptend%s(i,1:pver), work1, ptend%s(i,1:pver).ge.work1)
      ptend%s(i,1:pver) = merge(ptend%s(i,1:pver), work2, ptend%s(i,1:pver).le.work2)
   end do

end subroutine cb_limiter

end module cloudbrain_output_limiter
