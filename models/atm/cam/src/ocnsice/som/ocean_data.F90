! Specified Slab Ocean Model (SOM) fields from the
! time varying boundary dataset: ocean mixed layer 
! depth (mld) and ocean mixed layer Q flux (oqf)
!
! $Id: ocean_data.F90,v 1.1.4.2 2003/02/27 00:58:20 rosinski Exp $
! $Author: rosinski $
!
module ocean_data
   use infnan
   use pmgrid, only: plon, plat

   integer, parameter :: psomtim = 12
   integer, parameter :: totsomsz = 2000

   real(r8), allocatable :: mld(:,:)      ! on-node ocean mixed layer depth
   real(r8), allocatable :: qfluxm(:,:,:) ! on-node monthly ocean mixed layer qflux
   real(r8), allocatable :: qflux(:,:)    ! on-node ocean mixed layer qflux
   real(r8), allocatable :: sstm(:,:,:)   ! on-node monthly sst
   real(r8), allocatable :: sst(:,:)      ! on-node sst

   real(r8) :: mlddat(plon,plat) = inf    ! workspace to read in field
   real(r8) :: qfluxdat(plon,plat) = inf  ! workspace to read in field
   real(r8) :: cdaysomm          = inf    ! Calendar day for prv. month SOM values read in
   real(r8) :: cdaysomp          = inf    ! Calendar day for nxt. month SOM values read in
!
   integer :: nm                 = 1      ! Array index for prv month som data
   integer :: np                 = 2      ! Array index for nxt month som data
   integer :: np1                = bigint ! current forward time index of som dataset
   integer :: qfluxid            = bigint ! netcdf id for som variable
   integer :: mldid              = bigint ! netcdf id for som variable
   integer :: sstid              = bigint ! netcdf id for som variable
   integer :: lonsiz             = bigint ! size of longitude dimension on som dataset
   integer :: levsiz             = bigint ! size of level dimension on som dataset
   integer :: latsiz             = bigint ! size of latitude dimension on som dataset
   integer :: timesiz            = bigint ! size of time dimension on som dataset
   integer :: date_som(totsomsz) = bigint ! Date on som dataset (YYYYMMDD)
   integer :: sec_som(totsomsz)  = bigint ! seconds of date on som dataset (0-86399)

CONTAINS

   subroutine allocate_ocean ()
      use ppgrid, only: pcols, begchunk, endchunk

      allocate (mld (pcols,begchunk:endchunk))
      allocate (qfluxm (pcols,begchunk:endchunk,2))
      allocate (qflux  (pcols,begchunk:endchunk))
      allocate (sstm(pcols,begchunk:endchunk,2))
      allocate (sst (pcols,begchunk:endchunk))

      mld(:,:)    = inf
      qfluxm(:,:,:)  = inf
      qflux(:,:)     = inf
      sstm(:,:,:) = inf
      sst(:,:)    = inf
   end subroutine allocate_ocean
end module ocean_data
