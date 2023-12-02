#include <misc.h> 
#include <params.h>

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   1) After FFT preceding Legendre analysis, reallocate fftbuf
!      to decompose over wavenumber, recombining latitudes.
!   2) Before FFT following Legendre synthesis, reallocate fftbuf
!      to recombine wavenumbers, decomposing over latitude.
! 
!-----------------------------------------------------------------------
!
! $Id: realloc4.F90,v 1.3.22.3 2003/12/18 16:22:32 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

subroutine realloc4a(fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   After FFT preceding Legendre analysis, reallocate fftbuf
!   to decompose over wavenumber, combining latitudes.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
! Modified:          P. Worley, September 2002
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use spmd_utils, only: pair, ceil2
   use mpishorthand
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype  = 1000
!---------------------------Input arguments--------------------------
!
   real(r8), intent(in)  :: fftbuf_in(plond,plev,9,beglat:endlat) 
                            ! buffer used for in-place FFTs
   real(r8), intent(out) :: fftbuf_out(2*maxm,plev,9,plat) 
                            ! buffer used for reordered Fourier coefficients
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer length_r, length_l, locrm(2*maxm)
   integer bpos
   integer procj, if, k, i
   integer lat_l, lat_r
   integer numlats_r, nl_r
!     
! First, copy local data to new location
   length_l = 2*numm(iam)
   do i=1,numm(iam)
      locrm(2*i-1) = 2*locm(i,iam)-1
      locrm(2*i)   = 2*locm(i,iam)
   enddo
   do lat_l=beglat,endlat
      do if=1,8
         do k=1,plev
            do i=1,length_l
               fftbuf_out(i,k,if,lat_l) = fftbuf_in(locrm(i),k,if,lat_l)
            enddo
         enddo
      enddo
      do i=1,length_l
         fftbuf_out(i,1,9,lat_l) = fftbuf_in(locrm(i),1,9,lat_l)
      enddo
   enddo
!
! Next, get remote data
   do procj=1,ceil2(npes)-1
      procid = pair(npes,procj,iam)
      length_r = 2*numm(procid)
      if ((length_r > 0 .or. length_l > 0) .and. procid >= 0) then
!
         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo
!
         buf1(1) = endlat-beglat+1
         bpos = 1
         do lat_l=beglat,endlat
            buf1(bpos+1) = lat_l
            bpos = bpos + 1
            do if=1,8
               do k=1,plev
                  do i=1,length_r
                     buf1(bpos+i) = fftbuf_in(locrm(i),k,if,lat_l)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo
            do i=1,length_r
               buf1(bpos+i) = fftbuf_in(locrm(i),1,9,lat_l)
            enddo
            bpos = bpos+length_r
         enddo

         call mpisendrecv (buf1,bpos,mpir8,procid,msgtype, &
                           buf2,bsiz,mpir8,procid,msgtype,mpicom)

         numlats_r = buf2(1)
         bpos = 1
         do nl_r=1,numlats_r
            lat_r = buf2(bpos+1)
            bpos = bpos + 1
            do if=1,8
               do k=1,plev
                  do i=1,length_l
                     fftbuf_out(i,k,if,lat_r) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo
            do i=1,length_l
               fftbuf_out(i,1,9,lat_r) = buf2(bpos+i)
            enddo
            bpos = bpos+length_l
         enddo
      end if
!
   end do
#endif
   return
end subroutine realloc4a
!
subroutine realloc4b(fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   Before FFT following Legendre synthesis, reallocate fftbuf
!   to combine wavenumbers, decomposing over latitude.
! 
! Author:  P. Worley, September 2002
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use spmd_utils, only: pair, ceil2
   use mpishorthand
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype  = 1000
!---------------------------Input arguments--------------------------
!
   real(r8), intent(in)  :: fftbuf_in(2*maxm,8,plevp,plat) 
                            ! buffer of Fourier coefficients to be reordered
   real(r8), intent(out) :: fftbuf_out(plond,8,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer locrm(2*maxm)
   integer bpos
   integer length_l, length_r
   integer lat_l, lat_r
   integer beglat_r, endlat_r
   integer procj, if, k, i
!     
! First, copy local data to new location
   length_l = 2*numm(iam)
   do i=1,numm(iam)
      locrm(2*i-1) = 2*locm(i,iam)-1
      locrm(2*i)   = 2*locm(i,iam)
   enddo
   do lat_l=beglat,endlat
      do k=1,plev
         do if=1,8
            do i=1,length_l
               fftbuf_out(locrm(i),if,k,lat_l) = fftbuf_in(i,if,k,lat_l)
            enddo
         enddo
      enddo

      do if=1,4
         do i=1,length_l
            fftbuf_out(locrm(i),if,plevp,lat_l) = fftbuf_in(i,if,plevp,lat_l)
         enddo
      enddo
   enddo
!
! Next, get remote data
   do procj=1,ceil2(npes)-1
      procid = pair(npes,procj,iam)
      length_r = 2*numm(procid)
      if ((length_r > 0 .or. length_l > 0) .and. procid >= 0) then
!
         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo
!
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = 0
         do lat_r=beglat_r,endlat_r
            do k=1,plev
               do if=1,8
                  do i=1,length_l
                     buf1(bpos+i) = fftbuf_in(i,if,k,lat_r)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo

            do if=1,4
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,if,plevp,lat_r)
               enddo
               bpos = bpos+length_l
            enddo

         enddo

         call mpisendrecv (buf1,bpos,mpir8,procid,msgtype, &
                           buf2,bsiz,mpir8,procid,msgtype,mpicom)

         bpos = 0
         do lat_l=beglat,endlat
            do k=1,plev
               do if=1,8
                  do i=1,length_r
                     fftbuf_out(locrm(i),if,k,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo

            do if=1,4
               do i=1,length_r
                  fftbuf_out(locrm(i),if,plevp,lat_l) = buf2(bpos+i)
               enddo
               bpos = bpos+length_r
            enddo

         enddo
      end if
!
   end do
!
#endif
   return
end subroutine realloc4b

