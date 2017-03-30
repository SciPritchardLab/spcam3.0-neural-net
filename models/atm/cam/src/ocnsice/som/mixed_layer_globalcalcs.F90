#include <misc.h>

module mixed_layer_globalcalcs
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perform mixed layer global calculations for energy conservation checks.
!
! Method: 
! Gather to a master processor who does all the work
!
! Author: Jim Rosinski/Bruce Briegleb
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, begchunk, endchunk
   use rgrid,        only: nlon
   use commap,       only: w
   use phys_grid,    only: gather_chunk_to_field
   use pmgrid,       only: plon, plat, masterproc
   use comsrf,       only: landfrac_field, tsocn
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

   CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(r8) function gmean (arr, wght)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of field "arr" over non-land
!
!-----------------------------------------------------------------------
!
! Arguments
!
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk)  ! Input array 
      real(r8), intent(out) :: wght                         ! weight 
!
! Local workspace
!
      real(r8) :: arr_field(plon,plat)   ! rectangular version of arr
      real(r8) :: wt                     ! weight
      real(r8) :: mean                   ! mean value
      real(r8) :: tmp                    ! temporary
      integer :: i, j                    ! longitude, latitude indices

      call t_startf ('gmean')
      call gather_chunk_to_field (1, 1, 1, plon, arr, arr_field)

      if (masterproc) then
         mean = 0.
         wght = 0.
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_field(i,j) < 1.) then
                  tmp = wt*(1. - landfrac_field(i,j))
                  mean = mean + arr_field(i,j)*tmp
                  wght  = wght + tmp
               end if
            end do
         end do
         mean = mean/wght
      end if
#if ( defined SPMD )
      call mpibcast (mean, 1, mpir8, 0, mpicom)
      call mpibcast (wght, 1, mpir8, 0, mpicom)
#endif
      gmean = mean
      call t_stopf ('gmean')

      return
   end function gmean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(r8) function gmean_warm (arr, wght, wghtw)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of field "arr" over ocean points whose temperature is > 0 deg C
!
!-----------------------------------------------------------------------
!
! Arguments
!
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk)
      real(r8), intent(in) :: wght
      real(r8), intent(out) :: wghtw
!
! Local workspace
!
      real(r8) :: arr_field(plon,plat)
      real(r8) :: tsocn_field(plon,plat)
      real(r8) :: wt
      real(r8) :: mean
      real(r8) :: tmp
      integer :: i, j

      call t_startf ('gmean_warm')
      call gather_chunk_to_field (1, 1, 1, plon, arr, arr_field)
      call gather_chunk_to_field (1, 1, 1, plon, tsocn, tsocn_field)

      mean  = 0.
      wghtw = 0.
      if (masterproc) then
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_field(i,j) < 1.) then
                  tmp = wt*(1. - landfrac_field(i,j))
                  mean = mean + arr_field(i,j)*tmp
                  if (tsocn_field(i,j) > 0.) then
                     wghtw  = wghtw + tmp
                  end if
               end if
            end do
         end do
         mean = mean/wght
      end if

#if ( defined SPMD )
      call mpibcast (mean, 1, mpir8, 0, mpicom)
      call mpibcast (wghtw, 1, mpir8, 0, mpicom)
#endif
      gmean_warm = mean
      call t_stopf ('gmean_warm')
      return
   end function gmean_warm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_conservation (qflux_adj, qfluxgmin, wght)
!
! Arguments
!
      real(r8), intent(in) :: qflux_adj(pcols,begchunk:endchunk)
      real(r8), intent(in) :: qfluxgmin
      real(r8), intent(in) :: wght
!
! Local workspace
!
      real(r8) :: qflux_adj_field(plon,plat)
      real(r8) :: wt
      real(r8) :: qadjgm
      integer :: i, j

      call t_startf ('check_conserve')
      call gather_chunk_to_field (1, 1, 1, plon, qflux_adj, qflux_adj_field)

      if (masterproc) then
         qadjgm = 0.
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_field(i,j) < 1. ) then
                  qadjgm = qadjgm + qflux_adj_field(i,j)*wt*(1.-landfrac_field(i,j))
               end if
            end do
         end do
         qadjgm = qadjgm/wght
      end if

#if ( defined SPMD )
      call mpibcast (qadjgm, 1, mpir8, 0, mpicom)
#endif

!BPB Compare with initial value:

      if (abs (qadjgm/qfluxgmin - 1.) > 0.0000001) then
         if (masterproc) then
            write(6,*) 'Problem in check_conservation: qadjgm,qfluxgmin=',qadjgm,qfluxgmin
            write(6,*) 'Should be equal to roundoff'
         end if
         call endrun ()
      end if
      call t_stopf ('check_conserve')
      return
   end subroutine check_conservation
end module mixed_layer_globalcalcs
