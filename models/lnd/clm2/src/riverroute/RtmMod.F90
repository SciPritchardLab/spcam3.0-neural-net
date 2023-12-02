#include <misc.h>
#include <preproc.h>

module RtmMod

#if (defined RTM) 

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: 
! 
! !DESCRIPTION: 
! River Routing Model (U. of Texas River Transport Model) 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : lsmlon, lsmlat, rtmlon, rtmlat 
  use shr_sys_mod, only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
  save
!
! RTM grid info
!
  integer, parameter, public :: rtmloni = 1           !rtm per-proc beginning lon index
  integer, parameter, public :: rtmlonf = rtmlon      !rtm per-proc ending lon index
  integer, parameter, public :: rtmlati = 1           !rtm per-proc beginning lat index
  integer, parameter, public :: rtmlatf = rtmlat      !rtm per-proc ending lat index

  real(r8), public, pointer :: latixy_r(:,:)          !rtm latitudes  of grid cells (degrees)       
  real(r8), public, pointer :: longxy_r(:,:)          !rtm longitudes of grid cells (degrees)       
  real(r8), public, pointer :: area_r(:,:)            !rtm gridcell area (km^2)          
  integer , public, pointer :: mask_r(:,:)            !rtm landmask (land=1,ocean=0)
!
! RTM runoff for coupled communication
!
  integer , public, pointer :: ocnrof_iindx(:)    !rtm longitude index of ocean runoff point
  integer , public, pointer :: ocnrof_jindx(:)    !rtm latitude index of ocean runoff point
  real(r8), public, pointer :: ocnrof_vec(:)      !rtm runoff vector (1/2 deg grid, kg/m^2/s)
!
! RTM time averaging calculation
!
  real(r8), public :: delt_rtm                        !rtm time step
  integer , public :: ncount_rtm                      !number of time samples to average over
!
! RTM fluxes
!
  real(r8), allocatable, public :: volr(:,:)          !water volume in cell (m^3)
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Rtmgridini   ! Initialize RTM grid and land mask
  public :: Rtmlandini   ! Initialize RTM-land interpolation weights
  public :: Rtmfluxini   ! Initialize RTM fluxout 
  public :: Rtmriverflux ! Interface with RTM river routing model
  public :: UpdateGlobal ! Update global quantities
  public :: restart_rtm  ! Read/write RTM restart data
!
! !PRIVATE MEMBER FUNCTIONS:
!  private :: RTM         ! River routing model (based on U. Texas code)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !PRIVATE TYPES:
!
! RTM grid info
!
  private
  integer  :: numlon_r(rtmlat)               !number of lon points at each lat
  real(r8), dimension(4) :: rtmedge = (/ 90., 180., -90., -180. /)  !N,E,S,W edges of rtm grid
!
! land model to RTM mapping. for each rtm grid cell:
!
  integer :: mxovr_s2r                       !max number of overlapping cells
  integer :: novr_s2r(rtmlon,rtmlat)         !number    of overlapping lsm cells
  integer , pointer :: iovr_s2r(:,:,:)   !lon index of overlapping land model cells
  integer , pointer :: jovr_s2r(:,:,:)   !lat index of overlapping land model cells
  real(r8), pointer :: wovr_s2r(:,:,:)   !weight    of overlapping land model cells
!
! RTM inputs at land model grid resolution
! rtmin_ave(1,:), rtmin_glob(1,:) and rtmin_loc(1,:) correspond to totrunin
! rtmin_ave(2,:), rtmin_glob(2,:) and rtmin_loc(2,:) correspond to prec
! rtmin_ave(3,:), rtmin_glob(3,:) and rtmin_loc(3,:) correspond to evap
!
  real(r8), pointer :: rtmin_ave(:,:)  !averaging buffer for runoff, prec, evap
  real(r8), pointer :: rtmin_glob(:,:) !global rtm input
  real(r8), pointer :: rtmin(:,:)      !local rtm input
  logical :: onproc(lsmlon,lsmlat)     !true = > grid point is on processor
!
! RTM history file variables at land model grid resolution
! qch(:,1) corresponds to qchocn2
! qch(:,2) corresponds to qchan2
!
  real(r8), pointer :: qch(:,:)     !per-proc river (channel) flow and flow into ocean (m**3 H2O /s)
!
! RTM flux variables at 1/2 degree resolution
!
  integer , allocatable :: rdirc(:,:)   !rtm river flow direction (0-8) 
  real(r8), allocatable :: fluxout(:,:) !water flux out of cell (m^3/s)
  real(r8), allocatable :: ddist(:,:)   !downstream distance (m)
  real(r8), allocatable :: rivarea(:,:) !cell area (m^2)
!
! RTM grid variables at 1/2 degree resolution
!
  real(r8), allocatable :: latsh(:)     !southern edge of cells at rtm grid     
  real(r8), allocatable :: lonwh(:,:)   !western  edge of cells at rtm grid     
!
! RTM inputs at 1/2 degree resolution
!
  real(r8), allocatable :: totrunin_r(:,:) !surface runoff (mm/s)
!
! RTM outputs returned at 1/2 degree resolution
!
  real(r8), allocatable :: flxlnd_r(:,:)   !river flux (m**3/s)
  real(r8), allocatable :: flxocn_r(:,:)   !river flux to the ocean (m**3/s)
  real(r8), allocatable :: dvolrdt_r(:,:)  !change in storage (mm/s)
  real(r8), allocatable :: volrtm(:,:)     !change in storage (m**3/s)
  real(r8), allocatable :: runrtm(:,:)     !input runoff on rtm grid (m**3/s)
!
! RTM water flux into cell
!
  real(r8) :: sfluxin(rtmloni:rtmlonf,rtmlati:rtmlatf) !water flux into cell (m3/s)
!
! RTM global averaging
!
  character(len=*),parameter :: F40="('(diag) ',a17,'    date  ', &
  &     '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
  character(len=*),parameter :: F41="('(diag) ',a17,'   nstep  ', &
  &     '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
  character(len=*),parameter :: F21="('(diag) ',a17,' ----------------------', &
  &     7('----------'))"
  character(len=*),parameter :: F22="('(diag) ',a17,i8,6(d13.4))"
!
  real(r8) :: prec_global            !total precipitation (m^3/sec) 
  real(r8) :: evap_global            !total evaporation (m^3/sec)
  real(r8) :: runlnd_global          !total input runoff on land grid (m^3/sec)
  real(r8) :: runrtm_global          !total input runoff on rtm grid (m^3/sec)
  real(r8) :: ocnrtm_global          !total ocean runoff on rtm grid (m^3/sec)
  real(r8) :: volrtm_global          !total change in storage on rtm (m^3/sec)
  integer  :: ncount_global          !global counter 
  integer  :: yrold                  !old year
!
! Parameters
!
  real(r8), parameter :: effvel = 0.35    !effective velocity (m/s)
!
! Temporaries
!
  real(r8), pointer :: temp2d(:,:)   !temporary
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmgridini
!
! !INTERFACE:
  subroutine Rtmgridini
!
! !DESCRIPTION: 
! Initialize RTM grid and land mask (U. of Texas River Transport Model)
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
#if (defined SPMD)
    use spmdMod, only : mpicom, MPI_REAL8, MPI_INTEGER, masterproc
#else
    use spmdMod, only : masterproc
#endif
    use areaMod, only : celledge, cellarea
    use clm_varctl, only : frivinp_rtm
    use clm_varcon, only : re
    use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !calc dist as in hydra
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !of grid cell down stream
    integer  :: i,j,k,n                       !loop indices
    integer  :: i2,j2                         !downstream i and j
    real(r8) :: deg2rad                       !pi/180
    real(r8) :: dx                            !lon dist. between grid cells (m)
    real(r8) :: dy                            !lat dist. between grid cells (m)
    real(r8) :: tempg(rtmlon,rtmlat)          !temporary buffer
    integer  :: tempgp(0:rtmlon+1,0:rtmlat+1) !temporary buffer  
    integer  :: ier                           !error code
!-----------------------------------------------------------------------
  
    ! Allocate rtm grid variables

    allocate(latixy_r(rtmlon,rtmlat), longxy_r(rtmlon,rtmlat), &
         latsh(rtmlat+1), lonwh(rtmlon+1,rtmlat), &
         area_r(rtmlon,rtmlat), mask_r(rtmlon,rtmlat), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'latixy_r, longxy_r, latsy, lonwh, area_r, mask_r'
       call endrun
    end if

    ! Allocate rtm flux variables

    allocate (volr(rtmloni:rtmlonf,rtmlati:rtmlat),&
              rdirc(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), &  
              fluxout(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), & 
              ddist(rtmloni:rtmlonf,rtmlati:rtmlatf), &           
              rivarea(rtmloni:rtmlonf,rtmlati:rtmlatf), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'rdirec, fluxout, ddist, rivarea'
       call endrun
    end if

    ! Allocate inputs and outputs to rtm at 1/2 degree resolution 

    allocate (totrunin_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxlnd_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxocn_r(rtmloni:rtmlonf,rtmlati:rtmlatf), & 
              dvolrdt_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              volrtm(rtmloni:rtmlonf,rtmlati:rtmlatf), &   
              runrtm(rtmloni:rtmlonf,rtmlati:rtmlatf), stat=ier)   
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'totrunin_r, flxlnd_r, flxocn_r, dvolrdt_r, volrtm, runrtm'
       call endrun
    end if

    ! Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180.
    volr(:,:) = 0.
    fluxout(:,:) = 0.
    flxocn_r(:,:) = 0.
    flxlnd_r(:,:) = 0.

    ! Open and read input data (river direction file)
    ! rtm operates from south to north and from the dateline

    if (masterproc) then
       write(6,*)'Columns in RTM = ',rtmlon
       write(6,*)'Rows in RTM    = ',rtmlat
       
       open (1,file=frivinp_rtm)
       write(6,*)'opened river direction data'
       do j = 1,rtmlat
          numlon_r(j) = 0
          do i = 1,rtmlon
             read(1,*) latixy_r(i,j),longxy_r(i,j),tempg(i,j)
             if (longxy_r(i,j) /= 1.e36) numlon_r(j) = numlon_r(j) + 1
             tempgp(i,j) = nint(tempg(i,j))
          enddo
       enddo
       close(1)

       write(6,*)'closed river direction data'
       write(6,*)
    endif

#if (defined SPMD)
    call mpi_bcast(numlon_r, size(numlon_r), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast(latixy_r, size(latixy_r), MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast(longxy_r, size(longxy_r), MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast(tempgp  , size(tempgp)  , MPI_INTEGER, 0, mpicom, ier)
#endif

    ! Determine RTM celledges, areas and interpolation masks
    
    call celledge (rtmlat    , rtmlon    , numlon_r  , longxy_r  , &
                   latixy_r  , rtmedge(1), rtmedge(2), rtmedge(3), &
                   rtmedge(4), latsh     , lonwh     )

    call cellarea (rtmlat    , rtmlon    , numlon_r  , latsh     , lonwh , &
                   rtmedge(1), rtmedge(2), rtmedge(3), rtmedge(4), area_r) 

    ! Determine rtm mask, downstream distance and area
       
    do i=1,rtmlon
       tempgp(i,0)        = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,1)
       tempgp(i,rtmlat+1) = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,rtmlat)
       if (tempgp(i,0)        /= 0) tempgp(i,0)        = mod(tempgp(i,0)       +4-1,8)+1
       if (tempgp(i,rtmlat+1) /= 0) tempgp(i,rtmlat+1) = mod(tempgp(i,rtmlat+1)+4-1,8)+1
    enddo
    do j=0,rtmlat+1
       tempgp(0,j) =tempgp(rtmlon,j)
       tempgp(rtmlon+1,j)=tempgp(1,j)
    enddo
    
    ! Determine rtm river flow direction (0-8)

    do j=rtmlati-1,rtmlatf+1
       do i=rtmloni-1,rtmlonf+1
          rdirc(i,j)=tempgp(i,j)
       enddo
    enddo
    
    ! Determine rtm ocn/land mask

    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (rdirc(i,j) == 0) then
             mask_r(i,j) = 0
          else
             mask_r(i,j) = 1
          end if
       enddo
    enddo
    
    ! Determine downstream distance - instead of reading a distance file 
    ! calculate the downstream distance 
    
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          i2 = i + ioff(tempgp(i,j))
          j2 = j + joff(tempgp(i,j))
          if (i2 == 0) i2 = 2                 !avoids i2 out of bounds in the following
          if (i2 == rtmlon+1) i2 = rtmlon-1   !avoids i2 out of bounds in the following  
          dy = deg2rad * abs(latixy_r(i,j)-latixy_r(i2,j2)) * re*1000.
          dx = deg2rad * abs(longxy_r(i,j)-longxy_r(i2,j2)) * re*1000. &
               *0.5*(cos(latixy_r(i,j)*deg2rad)+cos(latixy_r(i2,j2)*deg2rad))
          ddist(i,j) = sqrt(dx*dx + dy*dy)
          rivarea(i,j)=1.e6 * area_r(i,j)     !convert into m**2
       enddo
    enddo
    
  end subroutine Rtmgridini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmlandini
!
! !INTERFACE:
  subroutine Rtmlandini
!
! !DESCRIPTION: 
! Initialize RTM-land interpolation weights 
! (U. of Texas River Transport Model)
! and variables related to runoff time averaging
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clmpoint
    use spmdMod, only : masterproc
    use areaMod, only : areaini_point, mkmxovr
    use clm_varsur, only : numlon, area, lats, lonw, landmask
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: gi,i,j,k,n                   !loop indices
    integer  :: is,js                        !land model grid indices
    integer  :: ir,jr                        !rtm grid indices 
    real(r8) :: maskone_s(lsmlon,lsmlat)     !dummy field: see below                 
    real(r8) :: maskone_r(rtmlon,rtmlat)     !dummy field: see below                 
    integer  :: ocnrof_mask(rtmlon,rtmlat)   !rtm mask for ocean points with possible nonzero runoff
    integer  :: ocnrof_num                   !number of valid ocean points with possible nonzero runoff
    integer  :: yrnew                        !year (0, ...)
    integer  :: mon                          !month (1, ..., 12) 
    integer  :: day                          !day of month (1, ..., 31)
    integer  :: ncsec                        !seconds of current date
    integer  :: begg,endg,numg               !column 1d indices
    real(r8) :: offset                       !offset for interpolation from model->rtm grid
    real(r8) :: lonw_offset(lsmlon+1,lsmlat) !longitudinal offset for interpolation from model->rtm grid            
    integer  :: novr_i2o                     !number of overlapping land model cells in given rtm cell
    integer , pointer :: iovr_i2o(:)         !lon index of overlap input cell
    integer , pointer :: jovr_i2o(:)         !lat index of overlap input cell
    real(r8), pointer :: wovr_i2o(:)         !weight    of overlap input cell
    integer  :: ier                          !error code
!-----------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! The following section allows RTM and land model to coexist at different
    ! horizontal resolutions
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(6,*)
       write(6,*) 'Initializing area-averaging interpolation for RTM.....'
    endif
    
    ! To find fraction of each land model grid cell that is land based on rtm grid.
    ! For this purpose, want all rtm grid cells to contribute to grid cell 
    ! average on land model grid, i.e., all cells used regardless of whether land 
    ! or ocean. Do this by setting [maskone_s] = 1 
    
    ! [maskone_s]=1 means all grid cells on land model grid, regardless of whether
    ! land or ocean, will contribute to rtm grid.
    
    do j = 1,lsmlat
       do i = 1,numlon(j)
          maskone_s(i,j) = 1.
       end do
    end do
    
    ! [maskone_r] = 1 means all the rtm grid is land. Used as dummy
    ! variable so code will not abort with false, non-valid error check
    
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          maskone_r(i,j) = 1.
       end do
    end do

    ! --------------------------------------------------------------------
    ! Map weights from land model grid to rtm grid
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Initializing land model -> rtm interpolation .....'
    endif
    
    ! For each rtm grid cell: get lat [jovr_s2r] and lon [iovr_s2r] indices 
    ! and weights [wovr_s2r] of overlapping atm grid cells 
    
    call mkmxovr (lsmlon, lsmlat, numlon  , lonw , lats , &
                  rtmlon, rtmlat, numlon_r, lonwh, latsh, &
                  mxovr_s2r     , novr_s2r)

    allocate(iovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             jovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             wovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             iovr_i2o(mxovr_s2r), jovr_i2o(mxovr_s2r), wovr_i2o(mxovr_s2r), &
             stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlndini: Allocation error for ',&
            'iovr_s2r, jovr_s2r, wovr_s2r, iovr_i2o, jovr_i2o, wovr_i2o'
       call endrun
    end if

    ! Shift x-grid to locate periodic grid intersections. This
    ! assumes that all lonw(1,j) have the same value for all
    ! latitudes j and that the same holds for lonwh(1,j)
       
    if (lonw(1,1) < lonwh(1,1)) then
       offset = 360.0
    else
       offset = -360.0
    end if
    do js = 1, lsmlat
       do is = 1, numlon(js) + 1
          lonw_offset(is,js) = lonw(is,js) + offset
       end do
    end do
    
    ! Determine overlap indices and weights
    
    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf
          
          call areaini_point (ir           , jr         , lsmlon  , lsmlat  , numlon   , &
                              lonw         , lonw_offset, lats    , area    , maskone_s, &
                              rtmlon       , rtmlat     , numlon_r, lonwh   , latsh    , &
                              area_r(ir,jr), maskone_r(ir,jr), novr_i2o, iovr_i2o, jovr_i2o , &
                              wovr_i2o     , mxovr_s2r)                             
          
          if (novr_i2o /= novr_s2r(ir,jr)) then
             write(6,*)'Rtmlandini error: novr_i2o= ',novr_i2o,&
                  ' not equal to  novr_s2r ',novr_s2r(ir,jr),&
                  ' at ir,jr=',ir,jr
             call endrun
          else if (novr_i2o > mxovr_s2r) then
             write(6,*)'Rtmlandini error: novr_s2r= ',novr_s2r,&
                  ' greater than mxovr_s2r= ',mxovr_s2r
             call endrun
          endif
 
          do n = 1,novr_i2o
             iovr_s2r(ir,jr,n) = iovr_i2o(n)
             jovr_s2r(ir,jr,n) = jovr_i2o(n)
             wovr_s2r(ir,jr,n) = wovr_i2o(n)
          end do

       end do
    end do
    
    if (masterproc) then
       write(6,*) 'Successfully made land model -> rtm interpolation'
       write(6,*)
    endif

#if (defined COUP_CSM)

    ! --------------------------------------------------------------------
    ! Determine which ocean cells might have runoff values. 
    ! --------------------------------------------------------------------
    
    ! First loop over all ocean points and determine which are at the 
    ! end of rivers by examining if any neighboring points are land and 
    ! if that land neighbor points into this ocean point. Next loop over all
    ! ocean points and determine which overlap with at least one land cell.
    
    ocnrof_num = 0
    ocnrof_mask(:,:) = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (mask_r(i,j) == 0) then
             if (rdirc(i  ,j-1)==1) ocnrof_mask(i,j) = 1
             if (rdirc(i-1,j-1)==2) ocnrof_mask(i,j) = 1
             if (rdirc(i-1,j  )==3) ocnrof_mask(i,j) = 1
             if (rdirc(i-1,j+1)==4) ocnrof_mask(i,j) = 1
             if (rdirc(i  ,j+1)==5) ocnrof_mask(i,j) = 1
             if (rdirc(i+1,j+1)==6) ocnrof_mask(i,j) = 1
             if (rdirc(i+1,j  )==7) ocnrof_mask(i,j) = 1
             if (rdirc(i+1,j-1)==8) ocnrof_mask(i,j) = 1
             if (ocnrof_mask(i,j) == 0) then
                do n=1,novr_s2r(i,j)
                   is = iovr_s2r(i,j,n)
                   js = jovr_s2r(i,j,n)
                   if (landmask(is,js)==1 .and. wovr_s2r(i,j,n)>0.) then
                      ocnrof_mask(i,j) = 1
                   end if
                end do
             endif
          endif
          if (ocnrof_mask(i,j) == 1) ocnrof_num = ocnrof_num +1
       enddo
    enddo

    ! allocate ocean runoff vector and indices and determine indices
    ! need to reset ocnrof_num to 0 and do the counting again because need to first
    ! first count to allocate vector and must now count to actually determine indices
    
    allocate(ocnrof_vec(ocnrof_num), ocnrof_iindx(ocnrof_num), &
         ocnrof_jindx(ocnrof_num), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlndini: Allocation error for ',&
            'ocnrof_vec, ocnrof_iindx, ocnrof_jindx'
       call endrun
    end if
    
    ocnrof_num = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (ocnrof_mask(i,j) == 1) then
             ocnrof_num = ocnrof_num + 1
             ocnrof_iindx(ocnrof_num) = i
             ocnrof_jindx(ocnrof_num) = j
             ocnrof_vec(ocnrof_num) = 0.
          endif
       end do
    enddo
    
#endif

    ! Deallocate memory for rtm grid  - needed to be done here because
    ! rtm grid information had to be sent to coupler between calls to 
    ! Rtmgridini and Rtmlandini
    
    deallocate(latsh, lonwh, iovr_i2o, jovr_i2o, wovr_i2o)

    ! Initialize rtm time averaging variables
    ! Upon restart, the following variables will get new values 
    ! from the restart file - these values are only valid for
    ! initial runs

    ncount_rtm    = 0
    ncount_global = 0
    prec_global   = 0.  
    evap_global   = 0.
    runlnd_global = 0.
    runrtm_global = 0.
    volrtm_global = 0.
    ocnrtm_global = 0.

    call get_curr_date(yrold, mon, day, ncsec)

    begg = grid1d%beg
    endg = grid1d%end
    numg = grid1d%num

    allocate (rtmin(3,numg), rtmin_glob(3,numg), rtmin_ave(3,numg), qch(2,numg), &
         temp2d(2,numg), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlandini: Allocation error for ',&
            'rtmin, rtmin_glob, rtmin_ave, qch, temp2d'
       call endrun
    end if

    rtmin(:,:) = 0.
    rtmin_glob(:,:) = 0.
    rtmin_ave(:,:) = 0.
    qch(:,:) = 0.
    temp2d(:,:) = 0.

    onproc(:,:) = .false.
    do gi = begg,endg
       i = grid1d%ixy(gi)
       j = grid1d%jxy(gi)
       onproc(i,j) = .true.
    end do
    
  end subroutine Rtmlandini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmfluxini()
!
! !INTERFACE:
  subroutine Rtmfluxini()
!
! !DESCRIPTION: 
! Initialize RTM fluxout for case of initial run when initial data is 
! read in. For restart run, RTM fluxout is read from restart dataset.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use time_manager, only : get_step_size
    use clm_varctl, only : rtm_nsteps
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j	 !indices
    integer :: delt_rtm    !delt for rtm
!-----------------------------------------------------------------------

    delt_rtm = rtm_nsteps*get_step_size()   
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          if (mask_r(i,j)==1) then
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt_rtm)
          else
             fluxout(i,j) = 0.  
          endif
       enddo
    enddo

  end subroutine Rtmfluxini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmriverflux
!
! !INTERFACE:
  subroutine Rtmriverflux ()
!
! !DESCRIPTION: 
! Interface with RTM river routing model
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clmtype
    use clmpoint
    use clm_varpar, only : lsmlon, lsmlat
    use clm_varsur, only : numlon, area, landfrac
    use clm_varctl, only : rtm_nsteps
#if (defined SPMD)
    use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master,  &
	allgather_data, mpicom
#else
    use spmdMod, only : masterproc
#endif
    use time_manager, only : get_step_size, get_nstep
    use mapxy, only : vec2xy, xy2vec
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
!
! misc variables
!
    integer  :: gi,li,ci,pi                   !indices 
    integer  :: io,jo,ir,jr,is,js             !mapping indices
    integer  :: k,n,i,j	                      !indices
    integer  :: gindex                        !index into 1d column array 
    real(r8) :: wt                            !weight
    type(gridcell_type)     , pointer :: g    ! local pointer to derived subtype
    type(landunit_type)     , pointer :: l    ! local pointer to derived subtype
    type(column_type)       , pointer :: c    !local pointer to derived subtype
    type(column_wflux_type) , pointer :: cwf  !local pointer to derived subtype
    type(atm2lnd_flux_type) , pointer :: a2lf !local pointer to derived subtype
!
! inputs to RTM at land model resolution
!
    real(r8) :: totruninxy(lsmlon,lsmlat)     !surface runoff (mm H2O /s)
    real(r8) :: precxy(lsmlon,lsmlat)         !precipitation (mm H2O /s)
    real(r8) :: evapxy(lsmlon,lsmlat)         !evaporation (mm H2O /s)
!
! outputs returned from RTM converted to land model resolution 
!
    real(r8) :: flxout_s(lsmlon,lsmlat)       !river flow (m**3)
    real(r8) :: flxocn_s(lsmlon,lsmlat)       !flow into ocean (m**3)

    character(len=16) :: nameg                !grid 1d name 
    integer :: begg,endg,numg                 !grid 1d indices
    integer :: begc,endc                      !column 1d indices
    integer :: ier                            !error status
    integer :: nstep                          !time step index
!-----------------------------------------------------------------------

#if (defined TIMING_BARRIERS)
    call t_startf ('sync_clmrtm')
    call mpi_barrier (mpicom, ier)
    call t_stopf ('sync_clmrtm')
#endif

    begc = cols1d%beg
    endc = cols1d%end

    nameg = grid1d%name
    numg = grid1d%num
    begg = grid1d%beg
    endg = grid1d%end

    ! --------------------------------------------------------------------
    ! RTM inputs 
    ! --------------------------------------------------------------------

    ! Make gridded representation of runoff from subgrid data
    ! total surface runoff = surface runoff on soils 
    ! + runoff on glaciers, wetlands, and lakes (P-E) 

    call t_startf('rtm_input1')
    rtmin(:,:) = 0.
!$OMP PARALLEL DO PRIVATE (gi,g,a2lf,gindex,li,l,ci,c,cwf)
    do gi = begg,endg
       g => gpoint(gi)%g
       a2lf => g%a2lf
       rtmin(2,gi) = a2lf%forc_rain + a2lf%forc_snow
       do li = 1, g%gps%nlandunits
          l => g%l(li)
          do ci = 1, l%lps%ncolumns
             c => l%c(ci)
             cwf => c%cwf
             rtmin(1,gi) = rtmin(1,gi) + &
                  (cwf%qflx_surf + cwf%qflx_qrgwl + cwf%qflx_drain) * c%cps%wtxy 
             rtmin(3,gi) = rtmin(3,gi) + &
                  cwf%pwf_a%qflx_evap_tot * c%cps%wtxy  
          end do
       end do
    end do
    call t_stopf('rtm_input1')

    ! --------------------------------------------------------------------
    ! Average fluxes for RTM calculation if appropriate
    ! --------------------------------------------------------------------

    if (rtm_nsteps <= 1) then

       ! RTM averaging is not done
       ! Note: in spmd mode it is assumed that each mpi process has
       ! all the values of totruninxy, precxy and evapxy


#if (defined SPMD)
       call allgather_data(rtmin, rtmin_glob, clmlevel=nameg)
#else
       rtmin_glob => rtmin
#endif           
       totruninxy(:,:) = 0.
       precxy(:,:) = 0.
       evapxy(:,:) = 0.
       do gi = 1,numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          totruninxy(i,j) = rtmin_glob(1,gi)
          precxy(i,j) = rtmin_glob(2,gi)
          evapxy(i,j) = rtmin_glob(3,gi)
       end do
       delt_rtm = get_step_size()

    else


       ! RTM averaging is done only done by master processor - however
       ! all SPMD processe will continue below

       do gi = begg,endg
          rtmin_ave(:,gi) = rtmin_ave(:,gi) + rtmin(:,gi)
       end do
       ncount_rtm = ncount_rtm + 1     
       nstep = get_nstep()

       if ((mod(nstep,rtm_nsteps)==0) .and. (nstep>1)) then
          call t_startf('rtm_input2')
#if (defined SPMD)
          call allgather_data(rtmin_ave, rtmin_glob, clmlevel=nameg)
#else
          rtmin_glob => rtmin_ave
#endif           
          totruninxy(:,:) = 0.
          precxy(:,:) = 0.
          evapxy(:,:) = 0.
          do gi = 1,numg
             i = grid1d%ixy(gi)
             j = grid1d%jxy(gi)
             totruninxy(i,j) = rtmin_glob(1,gi)/ncount_rtm
             precxy(i,j) = rtmin_glob(2,gi)/ncount_rtm
             evapxy(i,j) = rtmin_glob(3,gi)/ncount_rtm
          end do
          delt_rtm = ncount_rtm*get_step_size()   !compute delt for rtm
          ncount_rtm = 0                          !reset counter to 0
          do gi = begg,endg
             rtmin_ave(:,gi) = 0.                    !reset averager 
          end do
          call t_stopf('rtm_input2')
       else
          do gi = begg,endg
             g => gpoint(gi)%g
             g%gwf%qchocn2 = qch(1,gi)
             g%gwf%qchan2 = qch(2,gi)
          end do
          RETURN
       endif


    endif

    ! --------------------------------------------------------------------
    ! RTM runoff - only master processor performs runoff computation
    ! --------------------------------------------------------------------

    call t_startf('rtm_input3')
    ! Map from land model grid to RTM grid (intepolate to 1/2 degree resolution)
!$OMP PARALLEL DO PRIVATE (jr,ir,n,is,js,wt)
    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf
          totrunin_r(ir,jr) = 0.
          do n = 1, novr_s2r(ir,jr)
             if (wovr_s2r(ir,jr,n) > 0.) then
                is = iovr_s2r(ir,jr,n)
                js = jovr_s2r(ir,jr,n)
                wt = wovr_s2r(ir,jr,n)
                totrunin_r(ir,jr) = totrunin_r(ir,jr) + wt*totruninxy(is,js)*landfrac(is,js)
             end if
          end do
       end do
    end do
    call t_stopf('rtm_input3')
    
    ! Determine fluxes on 1/2 degree grid
    
    call t_startf('rtm_calc')
    call Rtm 
    call t_stopf('rtm_calc')

#if (defined COUP_CSM)

    ! Determine ocean runoff vector to send to coupler 

    do n = 1,size(ocnrof_vec)
       i = ocnrof_iindx(n)
       j = ocnrof_jindx(n)
       ocnrof_vec(n) = flxocn_r(i,j)/(area_r(i,j)*1000.) ! units of kg/m^2/s
    end do
    
#endif

    ! Determine ocean runoff and total runoff on land model grid and
    ! compute global input runoff on rtm grid, global ocean runoff on
    ! rtm grid and global change in storage on rtm grid
    
    call t_startf('rtm_output1')
    do js = 1, lsmlat
       do is =1, numlon(js)
          flxout_s(is,js)  = 0.
          flxocn_s(is,js)  = 0.
       end do
    end do
    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf
          do n = 1, novr_s2r(ir,jr)
             if (wovr_s2r(ir,jr,n) > 0.) then
                is = iovr_s2r(ir,jr,n)
                js = jovr_s2r(ir,jr,n)
                if (onproc(is,js)) then
                   wt = wovr_s2r(ir,jr,n)
                   flxocn_s(is,js)  = flxocn_s(is,js) + wt*flxocn_r(ir,jr)
                   flxout_s(is,js)  = max(flxout_s(is,js), wt*flxlnd_r(ir,jr))
                end if
             end if
          end do
       end do
    end do
!$OMP PARALLEL DO PRIVATE (jr,ir)
    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf
          runrtm(ir,jr) = totrunin_r(ir,jr)*1000.*area_r(ir,jr)
          volrtm(ir,jr) = dvolrdt_r(ir,jr)*1000.*area_r(ir,jr)
       end do
    end do
    call t_stopf('rtm_output1')
    
    ! Convert gridded output to vector form 

    call t_startf('rtm_output2')
    do gi = begg,endg
       g => gpoint(gi)%g
       i = g%gps%ixy
       j = g%gps%jxy
       qch(1,gi) = flxocn_s(i,j) 
       qch(2,gi) = flxout_s(i,j)
       g%gwf%qchocn2 = qch(1,gi)
       g%gwf%qchan2 = qch(2,gi)
    end do
    call t_stopf('rtm_output2')

    ! Determine global quantities and increment global counter

    if (masterproc) then
       call t_startf('rtm_output3')
       call UpdateGlobal(totruninxy, precxy, evapxy)
       call t_stopf('rtm_output3')
    endif

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtm
!
! !INTERFACE:
  subroutine Rtm 
!
! !DESCRIPTION: 
! River routing model (based on U. Texas code)
! input is totrunin_r
! input/output is fluxout, volr
! outputs are dvolrdt_r, flxocn_r, flxlnd_r
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i, j                        !loop indices
    real(r8) :: buffern(rtmlon)             !temp buffer
    real(r8) :: buffers(rtmlon)             !temp buffer
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumdvolr(rtmlat)            !global sum (m3/s)
    real(r8) :: sumrunof(rtmlat)            !global sum (m3/s)
    real(r8) :: sumdvolr_tot                !global sum (m3/s)
    real(r8) :: sumrunof_tot                !global sum (m3/s)
!-----------------------------------------------------------------------

    ! Determine fluxout at extended points and at southern and northern outer lats

    fluxout(rtmlon+1,rtmlat+1) = fluxout(1,1)
    fluxout(rtmlon+1,0)        = fluxout(1,rtmlat)
    fluxout(0,0)               = fluxout(rtmlon,rtmlat)
    fluxout(0,rtmlat+1)        = fluxout(rtmlon,1)

    do i=1,rtmlon
       fluxout(i,0)        = fluxout(i,rtmlat)
       fluxout(i,rtmlat+1) = fluxout(i,1)
       buffern(i)          = fluxout(i,rtmlat)
       buffers(i)          = fluxout(i,1)
    enddo
    do j=1,rtmlat
       fluxout(0,j)        = fluxout(rtmlon,j)
       fluxout(rtmlon+1,j) = fluxout(1,j)
    enddo
    do i=0,rtmlon+1
       fluxout(i,0)        = buffern(mod(i+rtmlon/2-1,rtmlon)+1)
       fluxout(i,rtmlat+1) = buffers(mod(i+rtmlon/2-1,rtmlon)+1)
    enddo

    ! Determine cell-to-cell transport - calculate sfluxin

!$OMP PARALLEL DO PRIVATE (I,J)
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          sfluxin(i,j) = 0.
          if (rdirc(i  ,j-1)==1) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j-1)
          if (rdirc(i-1,j-1)==2) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j-1)
          if (rdirc(i-1,j  )==3) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j  )
          if (rdirc(i-1,j+1)==4) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j+1)
          if (rdirc(i  ,j+1)==5) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j+1)
          if (rdirc(i+1,j+1)==6) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j+1)
          if (rdirc(i+1,j  )==7) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j  )
          if (rdirc(i+1,j-1)==8) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j-1)
       enddo
    enddo

    ! Loops above and below must remain separate because fluxout is updated below

    sumdvolr(:) = 0.
    sumrunof(:) = 0.
!$OMP PARALLEL DO PRIVATE (I,J,DVOLRDT)
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf

          ! calculate change in cell storage volume change units for 
          ! totrunin from kg/m2s==mm/s -> m3/s

          dvolrdt = sfluxin(i,j) - fluxout(i,j) + 0.001*totrunin_r(i,j)*rivarea(i,j)

          ! calculate flux out of a cell:
          ! land: do not permit change in cell storage volume greater than volume present
          ! make up for the difference with an extraction term (eg from aquifers)
          ! ocean: do not permit negative change in cell storage volume,
          ! because at ocean points cell storage volume equals zero
          ! water balance check (in mm/s), convert runinxy from mm/s to m/s (* 1.e-3) 
          ! and land model area from km**2 to m**2 (* 1.e6)

          if (mask_r(i,j) == 1) then         ! land points
             volr(i,j)    = volr(i,j) + dvolrdt*delt_rtm
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt_rtm)
             flxlnd_r(i,j) = fluxout(i,j)
             flxocn_r(i,j) = 0.
          else                               ! ocean points
             flxlnd_r(i,j) = 0.
             flxocn_r(i,j) = dvolrdt
          endif
          sumdvolr(j) = sumdvolr(j) + dvolrdt
          sumrunof(j) = sumrunof(j) + totrunin_r(i,j)*1000.*area_r(i,j)
          dvolrdt_r(i,j) = 1000.*dvolrdt/rivarea(i,j)

       enddo
    enddo

    ! Global water balance calculation and error check

    sumdvolr_tot = 0.
    sumrunof_tot = 0.
    do j = 1,rtmlat
       sumdvolr_tot = sumdvolr_tot + sumdvolr(j)
       sumrunof_tot = sumrunof_tot + sumrunof(j)
    end do
    if (abs((sumdvolr_tot-sumrunof_tot)/sumrunof_tot) > 0.01) then
       write(6,*) 'RTM Error: sumdvolr= ',sumdvolr_tot,&
            ' not equal to sumrunof= ',sumrunof_tot
       call endrun
    end if

  end subroutine Rtm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UpdateGlobal
!
! !INTERFACE:
  subroutine UpdateGlobal(totruninxy, precxy, evapxy)
!
! !DESCRIPTION: 
! Update Global quantitities
! input is totrunin_r
! input/output is fluxout, volr
! outputs are dvolrdt_r, flxocn_r, flxlnd_r
!
! !USES:
    use clm_varsur, only : numlon, area, landfrac
    use time_manager, only : get_nstep, get_curr_date
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: totruninxy(lsmlon,lsmlat)  !surface runoff (mm H2O /s)
    real(r8), intent(inout) :: precxy(lsmlon,lsmlat)      !precipitation (mm H2O /s)
    real(r8), intent(inout) :: evapxy(lsmlon,lsmlat)      !evaporation (mm H2O /s)
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
!
! global balance
!
    integer  :: is, js        !indices   
    integer  :: yrnew         !year (0, ...)
    integer  :: mon           !month (1, ..., 12) 
    integer  :: day           !day of month (1, ..., 31)
    integer  :: ncsec         !seconds of current date  
    integer  :: ncdate        !current date   
    real(r8) :: prec_sum      !total precipitation (m^3/sec) 
    real(r8) :: evap_sum      !total evaporation (m^3/sec)
    real(r8) :: runlnd_sum    !total input runoff on land grid (m^3/sec)
    real(r8) :: runrtm_sum    !total input runoff on rtm grid (m^3/sec)
    real(r8) :: ocnrtm_sum    !total ocean runoff on rtm grid (m^3/sec)
    real(r8) :: volrtm_sum    !total change in storage on rtm (m^3/sec)
!-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (js,is)
    do js = 1,lsmlat
       do is = 1,numlon(js)
          precxy(is,js) = precxy(is,js)*area(is,js)*1000.*landfrac(is,js)
          evapxy(is,js) = evapxy(is,js)*area(is,js)*1000.*landfrac(is,js)
          totruninxy(is,js) = totruninxy(is,js)*area(is,js)*1000.*landfrac(is,js)
       end do
    end do
    
    prec_sum   = sum(precxy)
    evap_sum   = sum(evapxy)
    runlnd_sum = sum(totruninxy)
    runrtm_sum = sum(runrtm)
    volrtm_sum = sum(volrtm)
    ocnrtm_sum = sum(flxocn_r)
    
    prec_global   = prec_global   + prec_sum
    evap_global   = evap_global   + evap_sum
    runlnd_global = runlnd_global + runlnd_sum
    runrtm_global = runrtm_global + runrtm_sum
    volrtm_global = volrtm_global + volrtm_sum
    ocnrtm_global = ocnrtm_global + ocnrtm_sum
    
    ncount_global = ncount_global + 1
    
    ! Print out diagnostics if appropriate
    
    write(6,*)
    write(6,F41)'water inst   '
    write(6,F21)'water inst   '
    write(6,F22)'water inst   ',get_nstep(), prec_sum, evap_sum, &
         runlnd_sum, runrtm_sum, volrtm_sum, ocnrtm_sum
    write(6,*)
    
    call get_curr_date(yrnew, mon, day, ncsec)
    ncdate = yrnew*10000 + mon*100 + day
    if (yrnew /= yrold) then
       prec_global   = prec_global/ncount_global
       evap_global   = evap_global/ncount_global
       runlnd_global = runlnd_global/ncount_global
       runrtm_global = runrtm_global/ncount_global
       volrtm_global = volrtm_global/ncount_global
       ocnrtm_global = ocnrtm_global/ncount_global
       ncount_global = 0
       write(6,*)
       write(6,F40)'water tavg   '
       write(6,F21)'water tavg   '
       write(6,F22)'water tavg   ',ncdate, prec_global, evap_global,&
            runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
       write(6,*)
    endif
    yrold = yrnew

  end subroutine UpdateGlobal

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_rtm
!
! !INTERFACE:
  subroutine restart_rtm (nio, flag)
!
! !DESCRIPTION: 
! Read/write RTM restart data
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8	          
    use clmpoint
    use clmtype, only : cols1d
#if (defined SPMD)
    use spmdMod, only : masterproc, scatter_data_from_master, gather_data_to_master, &
         MPI_INTEGER, MPI_REAL8, mpicom
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS: 
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: gi                            !indices 
    integer :: begg,endg,numg                !column 1d indices
    character (len=16) :: nameg              !column 1d name 
    real(r8), pointer, dimension(:,:) :: rtmin_ave_glob  !temporary for MPI scatter
    real(r8), pointer, dimension(:,:) :: qch_glob        !temporary for MPI scatter
    integer :: ier                           !error status 
!-----------------------------------------------------------------------

    ! Set shorthand for derived types

    begg = grid1d%beg
    endg = grid1d%end
    numg = grid1d%num
    nameg = grid1d%name

    ! Allocate dynamic memory

    allocate(rtmin_ave_glob(3,numg), qch_glob(2,numg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'restart_rtm: allocation error rtmin_ave_glob and qch_glob'
       call endrun
    end if

    ! Read RTM restart 

    if (flag == 'read') then
       if (masterproc) then
          read (nio) volr
          read (nio) fluxout
          read (nio) ncount_rtm 
          read (nio) rtmin_ave_glob, qch_glob
          read (nio) ncount_global, yrold, prec_global, evap_global, &
               runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
#if (defined COUP_CSM)
          read (nio) ocnrof_vec     
#endif
       endif
#if (defined SPMD)
       call mpi_bcast(volr, size(volr), MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast(fluxout, size(fluxout), MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast(ncount_rtm, 1, MPI_INTEGER, 0, mpicom, ier)
       call scatter_data_from_master (rtmin_ave, rtmin_ave_glob, clmlevel=nameg)
       call scatter_data_from_master (qch, qch_glob, clmlevel=nameg)
#if (defined COUP_CSM)
       call mpi_bcast(ocnrof_vec, size(ocnrof_vec), MPI_INTEGER, 0, mpicom, ier)
#endif
#else
       qch(:,:) = qch_glob(:,:)
       rtmin_ave(:,:) = rtmin_ave_glob(:,:)
#endif
    endif

    ! Write RTM restart 

    if (flag == 'write') then
#if (defined SPMD)
       call gather_data_to_master (rtmin_ave, rtmin_ave_glob, clmlevel=nameg)
       call gather_data_to_master (qch, qch_glob, clmlevel=nameg)
#else
       qch_glob(:,:) = qch(:,:)
       rtmin_ave_glob(:,:) = rtmin_ave(:,:)
#endif
       if (masterproc) then
          write (nio) volr
          write (nio) fluxout
          write (nio) ncount_rtm
          write (nio) rtmin_ave_glob, qch_glob
          write (nio) ncount_global, yrold, prec_global, evap_global, &
               runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
#if (defined COUP_CSM)
          write (nio) ocnrof_vec    
#endif
       endif
    endif

    ! Deallocate dynamic memory

    deallocate(rtmin_ave_glob, qch_glob)

  end subroutine restart_rtm

#endif

end module RtmMod


