#include <misc.h>
#include <params.h>

module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: Parameters and variables related to the dynamics grid
! 
! Author: 
! 
!-----------------------------------------------------------------------

   integer, parameter :: plon   = PLON                     ! number of longitudes
   integer, parameter :: plev   = PLEV                     ! number of vertical levels
   integer, parameter :: plat   = PLAT                     ! number of latitudes
   integer, parameter :: plevp  = plev + 1                 ! plev + 1
   integer, parameter :: nxpt   = 1                        ! no. of pts outside active domain of interpolant
   integer, parameter :: jintmx = 2                        ! number of extra latitudes in polar region
   integer, parameter :: i1     = 1 + nxpt                 ! model starting longitude index
   integer, parameter :: j1     = 1 + nxpt + jintmx        ! model starting latitude index
   integer, parameter :: plond  = plon + 1 + 2*nxpt        ! slt extended domain longitude
   integer, parameter :: plond1 = plond - i1 +1            ! slt extended domain longitude starting at i1
   integer, parameter :: platd  = plat + 2*nxpt + 2*jintmx ! slt extended domain lat.
   integer, parameter :: numbnd = nxpt + jintmx            ! no.of lats passed N and S of forecast lat
   integer, parameter :: plnlv  = plon*plev                ! Length of multilevel field slice
   integer, parameter :: plndlv = plond*plev               ! Length of multilevel 3-d field slice
!
   integer :: beglat     ! beg. index for latitudes owned by a given proc
   integer :: endlat     ! end. index for latitudes owned by a given proc
   integer :: beglatex   ! extended grid beglat
   integer :: endlatex   ! extended grid endlat
   integer :: begirow    ! beg. index for latitude pairs owned by a given proc
   integer :: endirow    ! end. index for latitude pairs owned by a given proc
   integer :: numlats    ! number of latitudes owned by a given proc
   integer :: numlatsex  ! number of latitudes owned by a given proc extended grid
   integer :: iam        ! MPI task id

   logical :: masterproc            ! Flag for (iam eq 0)
   logical :: dyngrid_set = .false. ! flag indicates dynamics grid has been set
!
#if ( ! defined SPMD )
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (beglatex = 1)
   parameter (endlatex = platd)
   parameter (begirow  = 1)
   parameter (endirow  = plat/2)
   parameter (numlats  = plat)
   parameter (numlatsex= platd)
   parameter (iam      = 0)

   parameter (masterproc = .true.)
#endif
end module pmgrid

