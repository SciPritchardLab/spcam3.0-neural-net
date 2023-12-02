#include <misc.h>
#include <params.h>
subroutine realloc7 (vmax2d  ,vmax2dt ,vcour   )
#ifdef SPMD
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for energy and log stats
!
! Author:  J. Rosinski
! Modified: P. Worley, October 2002
!
!----------------------------------------------------------------------------
!
! $Id: realloc7.F90,v 1.5.8.4 2003/12/18 16:22:34 pworley Exp $
! $Author: pworley $
!
!----------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plat, plev, iam, numlats, beglat, endlat
  use mpishorthand
  use spmd_dyn
  use spmd_utils, only: pair, ceil2

  implicit none

#include <comsta.h>

!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype=1000 ! message passing id
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(inout)   :: vmax2d (plev,plat) ! Max. wind at each level, latitude
  real(r8), intent(inout)   :: vmax2dt(plev,plat) ! Max. truncated wind at each lvl,lat
  real(r8), intent(inout)   :: vcour  (plev,plat) ! Maximum Courant number in slice
!
!---------------------------Local workspace-----------------------------
!
  integer procid
  integer bufpos
  integer procj
  integer j,k,jstrt
  integer beglat_p,endlat_p,numlats_p,jstrt_p
!
!-----------------------------------------------------------------------
!
  jstrt = beglat - 1
  do procj=1,ceil2(npes)-1
     procid = pair(npes,procj,iam)
     if (procid.ge.0) then
        buf1(1) = beglat
        buf1(2) = numlats
        bufpos = 2
! psurf
        do j=1,numlats
           buf1(bufpos+j) = psurf(jstrt+j)
        enddo
        bufpos = bufpos + numlats
! stq
        do j=1,numlats
           buf1(bufpos+j) = stq(jstrt+j)
        enddo
        bufpos = bufpos + numlats
! rmst
        do j=1,numlats
           buf1(bufpos+j) = rmst(jstrt+j)
        enddo
        bufpos = bufpos + numlats
! rmsd
        do j=1,numlats
           buf1(bufpos+j) = rmsd(jstrt+j)
        enddo
        bufpos = bufpos + numlats
!vmax2d
        do j=beglat,endlat
           do k=1,plev
              buf1(bufpos+k) = vmax2d(k,j)
           enddo
           bufpos = bufpos + plev
        enddo
! vmax2dt
        do j=beglat,endlat
           do k=1,plev
              buf1(bufpos+k) = vmax2dt(k,j)
           enddo
           bufpos = bufpos + plev
        enddo
! vcour
        do j=beglat,endlat
           do k=1,plev
              buf1(bufpos+k) = vcour(k,j)
           enddo
           bufpos = bufpos + plev
        enddo
!
        call mpisendrecv(buf1,bufpos,mpir8,procid,msgtype, &
           buf2,bsiz  ,mpir8,procid,msgtype,mpicom)
!
        beglat_p = buf2(1)
        numlats_p = buf2(2)
        bufpos = 2
! psurf
        jstrt_p  = beglat_p - 1
        endlat_p = beglat_p - 1 + numlats_p
        do j=1,numlats_p
           psurf(jstrt_p+j) = buf2(bufpos+j)
        enddo
        bufpos = bufpos + numlats_p
! stq
        do j=1,numlats_p
           stq(jstrt_p+j) = buf2(bufpos+j)
        enddo
        bufpos = bufpos + numlats_p
! rmst
        do j=1,numlats_p
           rmst(jstrt_p+j) = buf2(bufpos+j)
        enddo
        bufpos = bufpos + numlats_p
! rmsd
        do j=1,numlats_p
           rmsd(jstrt_p+j) = buf2(bufpos+j) 
        enddo
        bufpos = bufpos + numlats_p
! vmax2d
        do j=beglat_p,endlat_p
           do k=1,plev
             vmax2d(k,j) = buf2(bufpos+k)
           enddo
           bufpos = bufpos + plev
        enddo
! vmax2dt
        do j=beglat_p,endlat_p
           do k=1,plev
              vmax2dt(k,j) = buf2(bufpos+k)
           enddo
           bufpos = bufpos + plev
        enddo
! vcour
        do j=beglat_p,endlat_p
           do k=1,plev
              vcour(k,j) = buf2(bufpos+k)
           enddo
           bufpos = bufpos + plev
        enddo
     endif
  enddo
#endif
  return
end subroutine realloc7
