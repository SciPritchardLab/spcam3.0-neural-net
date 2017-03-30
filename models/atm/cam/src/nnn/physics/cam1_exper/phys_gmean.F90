#include <misc.h>
module phys_gmean
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perform mixed layer global calculations for energy conservation checks.
!
! Method: 
! Gather to a master processor who does all the work
!
! Author: Byron Boville from SOM code by Jim Rosinski/Bruce Briegleb
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, begchunk, endchunk
   use rgrid,        only: nlon
   use commap,       only: w
   use phys_grid,    only: gather_chunk_to_field
   use pmgrid,       only: plon, plat, masterproc
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

   CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(r8) function gmean (arr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of field "arr" in the physics chunked decomposition
!
!-----------------------------------------------------------------------
!
! Arguments
!
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk)  ! Input array, chunked
!
! Local workspace
!
      real(r8) :: arr_field(plon,plat)   ! rectangular version of arr
      real(r8) :: zmean                  ! zonal mean value
      real(r8) :: tmean                  ! temp global mean value
      integer :: i, j                    ! longitude, latitude indices

      call t_startf ('gmean')
      call gather_chunk_to_field (1, 1, 1, plon, arr, arr_field)

      if (masterproc) then
         tmean = 0.
         do j=1,plat
            zmean = 0.
            do i=1,nlon(j)
               zmean = zmean + arr_field(i,j)
            end do
            tmean = tmean + zmean * 0.5*w(j)/nlon(j)
         end do
      end if
#if ( defined SPMD )
      call mpibcast (tmean, 1, mpir8, 0, mpicom)
#endif
      gmean = tmean
      call t_stopf ('gmean')

      return
    end function gmean
end module phys_gmean
