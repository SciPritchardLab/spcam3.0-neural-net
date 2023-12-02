#include <misc.h>
#include <params.h>

      module comozp
!----------------------------------------------------------------------- 
! 
! Purpose: Variables associated with ozone dataset
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid

      implicit none

      real(r8) cdayozm  ! dataset calendar day previous month
      real(r8) cdayozp  ! dataset calendar day next month
      real(r8) cplos    ! constant for ozone path length integral
      real(r8) cplol    ! constant for ozone path length integral

      integer nm        ! Array indices for previous month ozone data
      integer np        ! Array indices for next month ozone data
      integer oznid     ! netcdf id for ozone variable
      integer lonsiz    ! size of longitude dimension on ozone dataset
      integer levsiz    ! size of level dimension on ozone dataset
      integer latsiz    ! size of latitude dimension on ozone dataset
      integer timesiz   ! size of time dimension on ozone dataset
      integer np1       ! current forward time index of ozone dataset

      type ozmixm_pters
        real(r8), dimension(:,:,:), pointer :: val  ! (pcols,levsiz,begchunk:endchunk)
      end type ozmixm_pters
      type (ozmixm_pters) :: ozmixm(2)          ! mixing ratios for lower and upper bounds
      real(r8), allocatable :: ozmix(:,:,:)     ! mixing ratio
                                                ! (pcols,levsiz,begchunk:endchunk)
      real(r8), allocatable :: pin(:)           ! ozone pressure level (levsiz)
      real(r8), allocatable :: ozlon(:)         ! Longitudes of bdy dataset (lonsiz)
      real(r8), allocatable :: ozlat(:)         ! Latitudes of bdy dataset (latsiz)

      integer, allocatable :: date_oz(:)        ! Date on ozone dataset (YYYYMMDD)
                                                ! (timesiz)
      integer, allocatable :: sec_oz(:)         ! seconds of date (0-86399)
                                                ! (timesiz)

      end module comozp
