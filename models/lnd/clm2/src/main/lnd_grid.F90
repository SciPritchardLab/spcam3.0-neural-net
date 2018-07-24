#include <misc.h>
#include <preproc.h>

module lnd_grid
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_grid
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varpar, only : lsmlon, lsmlat, maxpatch, maxpatch_pft
   use clm_varsur, only : numlon, landmask
#ifdef SPMD
#ifdef COUP_CAM
   use pmgrid, only : masterproc, iam
   use spmd_dyn, only : npes
#else
   use spmdMod, only : masterproc, iam, npes
#endif
#else
   use spmdMod, only : masterproc, iam, npes
#endif
!
! !PUBLIC TYPES:
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public surface_grid_init          ! initializes land surface decomposition
   public get_nclumps                ! returns the number of clumps defined
   public get_clump_cell_id_coord    ! returns clump/cell ids based on lon/lat
   public get_clump_owner_id         ! returns clump owner based on clump id
   public get_clump_ncells_proc      ! returns number of cells for process
   public get_clump_ncells_id        ! returns number of cells in clump
   public get_clump_coord_id         ! returns lon/lat coordinates based on id
   public get_clump_info             ! returns clump info 
   public get_clump_gcell_info       ! returns 1d gridcell index
   SAVE
!
! !DESCRIPTION:
! Module provides a descomposition into a clumped data structure which can be
! mapped back to atmosphere physics chunks.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!
! !PRIVATE TYPES
   private increment_gcell_info       ! updates gridcell, landunits, columns and pfts counters

   private

   type clump
      integer, dimension(:), pointer :: lon    ! cell longitudes
      integer, dimension(:), pointer :: lat    ! cell latitudes
      integer, dimension(:), pointer :: gi     ! 1d grid cell index
      integer :: owner                         ! process id owning clump
      integer :: ncells                        ! number of gridcells in clump
      integer :: nlunits                       ! number of landunits in clump
      integer :: ncolumns                      ! number of columns in clump
      integer :: npfts                         ! number of pfts in clump
   end type clump                             
   type(clump), dimension(:), allocatable, private :: clumps
   integer :: nclumps                         ! total number of clumps

   type pmulc
      integer :: clumpid                          ! clump id for (lon, lat)
      integer :: cell                             ! matching cell id 
   end type pmulc
   type(pmulc), dimension(:,:), allocatable, private :: pmulcs

#if (defined LANDUNIT_DECOMP)
   ! Landunits per clump affects how many clumps will be defined and distributed
   ! across processes.  If lpclump is set to zero, then the land grid will be
   ! decomposed into one clump per process.

   integer, parameter, private :: lpclump = 0           ! landunits per clump
#else
   ! Pfts per clump affects how many clumps will be defined and distributed
   ! across processes.  If ppclump is set to zero, then the land grid will be
   ! decomposed into one clump per process.

   integer, parameter, private :: ppclump = 0           ! pfts per clump
#endif

!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   CONTAINS

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: surface_grid_init
!
! !INTERFACE:
   subroutine surface_grid_init(wtxy)
!
! !USES
   use clmtype, only : grid1d, land1d, cols1d, pfts1d
!
   implicit none
   real(r8), intent(in) :: wtxy(lsmlon, lsmlat, maxpatch)! subgrid patch weights
!
! !LOCAL VARIABLES:
   integer :: ppc                            ! min number of pfts per clump
   integer :: lpc                            ! min number of landunits per clump
   integer :: i, j, m, p, n                  ! loop indices
   integer :: c, cc, nc, g, gi               ! loop indices
   integer, allocatable :: cpp(:)            ! clumps per process
   integer :: icell, fcell, ncells           ! clump initial, final and total gridcells
   integer :: ilunit, flunit, nlunits        ! clump initial, final and total landunits
   integer :: icolumn, fcolumn, ncolumns     ! clump initial, final and total columns
   integer :: ipft, fpft, npfts              ! clump initial, final and total pfts
   integer :: ier                            ! error codes
   integer :: nzero                          ! first clump with zero gridcells
   logical :: clumpfull                      ! true => clump is full
!
! !CALLED FROM:
! subroutine clm_map() (clm_map.F90)
!
! !DESCRIPTION:
! This subroutine initializes the land surface decomposition into a clump data
! structure.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   ! Find total global number of grid cells, landunits, columns and pfts 

   ncells = 0
   nlunits = 0
   ncolumns = 0
   npfts = 0
   do j = 1, lsmlat
      do i = 1, numlon(j)
         if (landmask(i,j) == 1) then         
            call increment_gcell_info (i, j, wtxy, ncells, nlunits, ncolumns, npfts)
         end if
      end do
   end do

   ! Error check on total number of gridcells

   if (npes > ncells) then
      write (6,*) 'surface_grid_init(): Number of processes exceeds number ', &
           'of land grid cells'
      call endrun
   end if

   ! Diagnostic output

   if (masterproc) then
      write (6,*)' Surface Grid Characteristics'
      write (6,*)'   longitude points          = ',lsmlon
      write (6,*)'   latitude points           = ',lsmlat
      write (6,*)'   minimum number of longitude points per latitude = ',minval(numlon)
      write (6,*)'   maximum number of longitude points per latitude = ',maxval(numlon)
      write (6,*)'   total number of gridcells = ',ncells
      write (6,*)'   total number of landunits = ',nlunits
      write (6,*)'   total number of columns   = ',ncolumns
      write (6,*)'   total number of pfts      = ',npfts
      write (6,*)
   endif
    
   ! Determine number of gridcell clumps and allocate dynamic memory

#if (defined LANDUNIT_DECOMP)
   if (lpclump > 0) then 
      ! Decompose by (approximately) a fixed number of landunits per clump
      ! (all andunits in a single grid cell must reside in the same clump)
      lpc = lpclump
      nclumps = nlunits / lpc
      if (nclumps*lpclump < nlunits) nclumps = nclumps + 1
   else
      ! Normal decomposition: one clump per process.
      nclumps = npes
      ! Minimum number of landunits per clump
      lpc = nlunits / nclumps             
   endif
#else
   if (ppclump > 0) then 
      ! Decompose by (approximately) a fixed number of pfts per clump
      ! (all pfts in a single grid cell must reside in the same clump)
      ppc = ppclump
      nclumps = npfts / ppc
      if (nclumps*ppclump < npfts) nclumps = nclumps + 1
   else
      ! Normal decomposition: one clump per process.
      nclumps = npes
      ! Minimum number of patches per clump
      ppc = npfts / nclumps             
   end if
#endif
   if (nclumps < npes) then
      write (6,*) 'surface_grid_init(): Number of gridcell clumps is less than ', &
         'the number of processes'
      call endrun
   end if

   allocate(clumps(1:nclumps), stat=ier)
   if (ier /= 0) then
      write (6,*) 'surface_grid_init(): allocation error for clumps'
      call endrun
   end if
   allocate(pmulcs(1:lsmlon, 1:lsmlat), stat=ier)
   if (ier /= 0) then
      write (6,*) 'surface_grid_init(): allocation error for pmulcs'
      call endrun
   end if

   ! Determine number of grid cells in each clump

   clumps(:)%ncells = 0
   clumps(:)%nlunits = 0
   clumps(:)%ncolumns = 0
   clumps(:)%npfts = 0
   c = 1
   do j = 1,lsmlat
      do i = 1,numlon(j)
         if (landmask(i,j) == 1) then
            ! udpate number of grid cells
            call increment_gcell_info (i, j, wtxy, clumps(c)%ncells, clumps(c)%nlunits, &
                 clumps(c)%ncolumns, clumps(c)%npfts)
            ! determine if clump is full
            clumpfull = .false.
#if (defined LANDUNIT_DECOMP)
            if (clumps(c)%nlunits>=lpc .and. c<nclumps) then
               clumpfull = .true. 
            endif
#else
            if (clumps(c)%npfts>=ppc .and. c<nclumps) then
               clumpfull = .true.
            endif
#endif
            if (c > nclumps) then
               write(6,*)'clump number ',c,' is greater than nclumps ',nclumps
               call endrun
            endif
            if (clumpfull) then
               c = c + 1
            endif
         end if
      end do
   end do
      
   ! If processes end up with no land gridcells, subtract appropriate
   ! number of gridcells from the first process that has more than one gridcell 
   ! and shift them up. Note -  nzero is first clump with zero gridcells.

    nzero = 0
    do nc = 1,nclumps
       if (clumps(nc)%ncells == 0) then
          nzero = nc
          exit
       endif
    end do
    if (nzero > 0) then 	

       ! Shift process clumps

       do nc = nzero,nclumps
          do cc = 1,nclumps
             if (clumps(cc)%ncells > 1) then
                clumps(cc)%ncells = clumps(cc)%ncells -1 
                clumps(nc)%ncells = clumps(nc)%ncells + 1
                exit
             end if
          end do
       end do

       ! Recalculate landunits, columns and pfts in each clump

       ncells = 0
       clumps(:)%nlunits = 0
       clumps(:)%ncolumns = 0
       clumps(:)%npfts = 0
       c = 1
       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (landmask(i,j) == 1) then
                ! udpate number of grid cells

                call increment_gcell_info (i, j, wtxy, ncells, clumps(c)%nlunits, &
                     clumps(c)%ncolumns, clumps(c)%npfts)

                ! determine if clump is full

                if (ncells == clumps(c)%ncells .and. c<nclumps) then
                   c = c + 1
                   ncells = 0
                endif
             end if
          end do
       end do
      
    end if ! end of nzero>0 if-block   

   ! Just in case not all clumps can be used (only when nclumps > npes),
   ! reset number of clumps

   if (c < npes) then
      write(6,*)'Number of clumps ',c,' cannot be less than the number of pes ',npes
      call endrun
   end if

   nclumps = c
   if (masterproc) then
      write (6,*) 'Decomposing land gridcells into ', nclumps, ' clumps '
   end if

   do nc = 1,nclumps
      if (clumps(nc)%ncells == 0 .or. clumps(nc)%nlunits == 0 .or. &
          clumps(nc)%ncolumns == 0 .or. clumps(nc)%npfts == 0) then
         write(6,*)'Invalid clump at clump index c ',nc
         write(6,*)'clumps(c)%ncells = ',clumps(nc)%ncells
         write(6,*)'clumps(c)%nlunits = ',clumps(nc)%nlunits
         write(6,*)'clumps(c)%ncolumns = ',clumps(nc)%ncolumns
         write(6,*)'clumps(c)%npfts = ',clumps(nc)%npfts
         call endrun
      end if
   end do

    ! Allocate dynamic memory

    do c = 1,nclumps
       ncells = clumps(c)%ncells
       allocate(clumps(c)%lon(ncells), stat=ier)
       if (ier /= 0) then
          write (6,*) 'surface_grid_init(): allocation error for clumps(c)%lon'
          call endrun
       end if
       allocate(clumps(c)%lat(ncells), stat=ier)
       if (ier /= 0) then
          write (6,*) 'surface_grid_init(): allocation error for clumps(c)%lat'
          call endrun
       end if
       allocate(clumps(c)%gi(ncells), stat=ier)
       if (ier /= 0) then
          write (6,*) 'surface_grid_init(): allocation error for clumps(c)%gi'
          call endrun
       end if
    end do
    
    ! Determine clump latitudes and longitudes and 1d gridcell indices
    
    pmulcs(:,:)%clumpid = 0
    pmulcs(:,:)%cell = 0  
    
    c = 1
    g = 1
    gi = 1
    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (landmask(i,j) == 1) then
             ! set clump lon,lat and 1d gridcell index
             clumps(c)%gi(g) = gi
             clumps(c)%lon(g) = i
             clumps(c)%lat(g) = j
             
             ! set clumpid and cell in pmulcs
             pmulcs(i,j)%clumpid = c
             pmulcs(i,j)%cell = g
             
             ! determine if clump is full
             if (g == clumps(c)%ncells) then
                c = c + 1
                g = 1
             else
                g = g+1
             endif
             ! increment 1d gridcell index
             gi = gi + 1
          end if
       end do
    end do
    
    ! Assign clumps to processes.  This must be done in order because 
    ! landunits, columns and pfts must be contiguous.
    
    allocate(cpp(0:npes-1), stat=ier)
    if (ier /= 0) then
       write (6,*) 'surface_grid_init(): allocation error for cpp'
       call endrun
    end if
    cpp(:) = nclumps / npes
    do p = 0,npes-1
       if (sum(cpp) == nclumps) exit
       cpp(p) = cpp(p) + 1
    end do
    
    ! Determine gridcell, landunit, column and pft beginning and
    ! ending indices for each clump
    
    icell = 1
    ilunit = 1
    icolumn = 1
    ipft = 1
    
    fcell = 0
    flunit = 0
    fcolumn = 0
    fpft = 0
    
    ncells = 0
    nlunits = 0
    ncolumns = 0
    npfts = 0

    p = 0
    n = 0
    do c = 1,nclumps
       clumps(c)%owner = p
       n = n + 1
       ncells = ncells + clumps(c)%ncells
       nlunits = nlunits + clumps(c)%nlunits
       ncolumns = ncolumns + clumps(c)%ncolumns
       npfts = npfts + clumps(c)%npfts
       if (p < iam) then
          icell = icell + clumps(c)%ncells
          ilunit = ilunit + clumps(c)%nlunits
          icolumn = icolumn + clumps(c)%ncolumns
          ipft = ipft + clumps(c)%npfts
       end if
       if (p == iam .and. n == cpp(p)) then
          fcell = ncells
          flunit = nlunits
          fcolumn = ncolumns
          fpft = npfts
       end if
       if (n == cpp(p)) then
          p = p + 1
          n = 0
       end if
    end do

    grid1d%beg = icell
    grid1d%end = fcell
    grid1d%num = ncells

    land1d%beg = ilunit
    land1d%end = flunit
    land1d%num = nlunits

    cols1d%beg = icolumn
    cols1d%end = fcolumn
    cols1d%num = ncolumns

    pfts1d%beg = ipft
    pfts1d%end = fpft
    pfts1d%num = npfts

    if (grid1d%end < grid1d%beg .or. land1d%end < land1d%beg .or. &
        cols1d%end < cols1d%beg .or.  pfts1d%end < pfts1d%beg) then
       write(6,*)'SURFACE_GRID_INIT ERROR: iam= ',iam,' invalid 1d setting'
       write(6,*)'grid1d%beg= ',grid1d%beg,' grid1d%end= ',grid1d%end,' grid1d%num= ',grid1d%num 
       write(6,*)'land1d%beg= ',land1d%beg,' land1d%end= ',land1d%end,' land1d%num= ',land1d%num 
       write(6,*)'cols1d%beg= ',cols1d%beg,' cols1d%end= ',cols1d%end,' cols1d%num= ',cols1d%num 
       write(6,*)'pfts1d%beg= ',pfts1d%beg,' pfts1d%end= ',pfts1d%end,' pfts1d%num= ',pfts1d%num 
       call endrun()
    end if

    write(6,*)'iam= ',iam,' beg gridcell= ',icell,' end gridcell= ',fcell, &
         ' total gridcells per processor= ',fcell-icell+1  
    write(6,*)'iam= ',iam,' beg landunit= ',ilunit,' end landunit= ',flunit, &
         ' total landunitsper processor = ',flunit-ilunit+1  
    write(6,*)'iam= ',iam,' beg column  = ',icolumn,' end column= ',fcolumn, &
         ' total columns per processor  = ',fcolumn-icolumn+1 
    write(6,*)'iam= ',iam,' beg pft     = ',ipft,' end pft= ',fpft, &
         ' total pfts per processor     = ',fpft-ipft+1  

    deallocate(cpp)

   end subroutine surface_grid_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nclumps
!
! !INTERFACE:
   integer function get_nclumps()
!
! !DESCRIPTION:
! This function returns the number of clumps defining the land model grid.
!
! !ARGUMENTS
   implicit none
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   get_nclumps = nclumps

   end function get_nclumps

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_cell_id_coord
!
! !INTERFACE:
   subroutine get_clump_cell_id_coord(lon, lat, cid, cell)
!
! !DESCRIPTION:
! This subroutine returns the id of the clump and cell corresponding to the
! longitude/latitude indices provided.
!
! !ARGUMENTS
   implicit none
   integer, intent(in)  :: lon, lat              ! longitude/latitude indices
   integer, intent(out) :: cid, cell             ! clump and cell id
!
! !CALLED FROM:
! Unused.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   cid  = pmulcs(lon,lat)%clumpid
   cell = pmulcs(lon,lat)%cell

   end subroutine get_clump_cell_id_coord

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_owner_id
!
! !INTERFACE:
   integer function get_clump_owner_id(cid)
!
! !DESCRIPTION:
! This function returns the MPI process id (rank) responsible for the clump
! identified by the clump id cid.
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: cid                    ! clump id
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   get_clump_owner_id = clumps(cid)%owner

   end function get_clump_owner_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_ncells_proc
!
! !INTERFACE:
   integer function get_clump_ncells_proc(proc)
!
! !DESCRIPTION:
! This function returns the number of cells contained within clumps owned
! by process proc.
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: proc                    ! process id
!
! !LOCAL VARIABLES:
   integer :: c                                   ! loop indices
   integer :: ncells                              ! cell counter
!
! !CALLED FROM:
! subroutine clm_map() (clm_map.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   ncells = 0
   do c = 1,nclumps
      if (clumps(c)%owner == proc) ncells = ncells + clumps(c)%ncells
   end do
   get_clump_ncells_proc = ncells

   end function get_clump_ncells_proc

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_ncells_id
!
! !INTERFACE:
   integer function get_clump_ncells_id(cid)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: cid                    ! clump id
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !DESCRIPTION:
! This function returns the number of cells contained within the clump
! identified by the clump id cid.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   get_clump_ncells_id = clumps(cid)%ncells

   end function get_clump_ncells_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_coord_id
!
! !INTERFACE:
   subroutine get_clump_coord_id(cid, ncells, lons, lats)
!
! !DESCRIPTION:
! This subroutine returns the first ncells longitutde/latitude indices for the
! clump identified by the clump id cid.  In practice, this ncells is always the
! number of cells in clump cid: clumps(cid)%ncells.
!
! !ARGUMENTS:
   implicit none
   integer, intent(in)  :: cid                    ! clump id
   integer, intent(in)  :: ncells                 ! number of grid cells
   integer, intent(out) :: lons(ncells)           ! longitude indices
   integer, intent(out) :: lats(ncells)           ! latitude indices
!
! !LOCAL VARIABLES:
   integer :: i                                   ! loop index
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   do i = 1,ncells
      lons(i) = clumps(cid)%lon(i)
      lats(i) = clumps(cid)%lat(i)
   end do

   end subroutine get_clump_coord_id

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_gridinfo
!
! !INTERFACE:
   subroutine get_clump_gridinfo(cid, cell, gi)
!
! !DESCRIPTION:
!
! !ARGUMENTS:
   implicit none
   integer, intent(in)  :: cid                    ! clump id
   integer, intent(in)  :: cell                   ! clump cell id
   integer, intent(out) :: gi                 ! 1d grid index
!
! !LOCAL VARIABLES:
   integer :: i                                   ! loop index
!
! !CALLED FROM:
! subroutines alltoall_clump_to_chunk_init(), alltoall_clump_to_chunk(), and
! alltoall_chunk_to_clump() in module lp_coupling (lp_coupling.F90)
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   gi = clumps(cid)%gi(cell)

   end subroutine get_clump_gridinfo

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: increment_gcell_info
!
! !INTERFACE:
   subroutine increment_gcell_info (i, j, wtxy, ngridcells, nlunits, ncolumns, npfts)
!
! !DESCRIPTION:
! Increment number of gridcells and gridcell components
!
! !ARGUMENTS
   implicit none
   integer , intent(in) :: i             ! longitude index
   integer , intent(in) :: j             ! latitude index 
   real(r8), intent(in) :: wtxy(lsmlon, lsmlat, maxpatch) ! subgrid pft weights
   integer , intent(inout) :: ngridcells ! number of gridcells 
   integer , intent(inout) :: nlunits    ! number of landunits
   integer , intent(inout) :: ncolumns   ! number of columns 
   integer , intent(inout) :: npfts      ! number of pfts 
!
! !CALLED FROM:
! subroutines surface_grid_init int hits module
!
! !REVISION HISTORY:
! 2002.09.11  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: m             ! loop index
   integer :: nveg          ! number of vegetated pfts      
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   ! increment total number of gridcells
   ngridcells = ngridcells + 1
   
   ! increment total number of pfts
   do m = 1,maxpatch
      if (wtxy(i,j,m) > 0.0) npfts = npfts + 1
   end do

   ! increment total number of landunits
   do m = maxpatch_pft+1, maxpatch
      if (wtxy(i,j,m) > 0.) nlunits = nlunits + 1
   end do
   nveg = 0
   do m = 1, maxpatch_pft                           
      if (wtxy(i,j,m) > 0.) nveg = nveg + 1
   end do
   if (nveg > 0) nlunits = nlunits + 1

   ! increment total number of columns
#if (defined NOCOMPETE)
   ! with no competition for one column, each pft has its own column

   ncolumns = npfts
#else
   ! with all vegetated pfts competing on one column, the number of
   ! columns is the same as the number of landunits 

   ncolumns = nlunits
#endif
          
 end subroutine increment_gcell_info

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_info
!
! !INTERFACE:
   subroutine get_clump_info(proc, ncells, nlunits, ncolumns, npfts)
!
! !DESCRIPTION:
! Determine beggining grid cell index and total number of grid cells in clump
!
! !ARGUMENTS
   implicit none
   integer, intent(in)  :: proc                ! process id
   integer, intent(out) :: ncells              ! number of gridcells in clump
   integer, intent(out) :: nlunits             ! number of landunit in clump
   integer, intent(out) :: ncolumns            ! number of column in clump
   integer, intent(out) :: npfts               ! number of pft in clump
!
! !CALLED FROM:
! subroutine clm_mapp in module clm_mapping.F90
!
! !REVISION HISTORY:
! 2002.09.11  Mariana  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: c       ! clump index
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   npfts = 0
   nlunits = 0
   ncolumns = 0
   ncells = 0
   do c = 1,nclumps
     if (clumps(c)%owner == proc) then
         ncells = ncells + clumps(c)%ncells
         nlunits = nlunits + clumps(c)%nlunits
         ncolumns = ncolumns + clumps(c)%ncolumns
         npfts = npfts + clumps(c)%npfts
      end if
   end do

   end subroutine get_clump_info

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clump_gcell_info
!
! !INTERFACE:
   subroutine get_clump_gcell_info (cid, cell, gi)
!
! !ARGUMENTS
   implicit none
   integer, intent(in)  :: cid                    ! clump id
   integer, intent(in)  :: cell                   ! clump cell id
   integer, intent(out) :: gi                     ! 1d gridcell index
!
! !LOCAL VARIABLES:
   integer :: i                                   ! loop index
!
! !CALLED FROM:
! subroutines alltoall_clump_to_chunk_init(), alltoall_clump_to_chunk(), and
! alltoall_chunk_to_clump() in module lp_coupling (lp_coupling.F90)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
! 2002.11.17  Mariana Vertenstein  Creation.
!
!EOP
!------------------------------------------------------------------------------
! lnd_grid.F90,v 1.1.2.5 2003/09/10 16:36:53 mvertens Exp
! mvertens
!------------------------------------------------------------------------------

   gi = clumps(cid)%gi(cell)

   end subroutine get_clump_gcell_info

end module lnd_grid
