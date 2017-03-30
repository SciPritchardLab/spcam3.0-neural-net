#include <misc.h>
#include <preproc.h>

module clm_mapping

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clm_mapping
!
! !DESCRIPTION: 
! Initializes sub-grid mapping
!                              
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined COUP_CAM)
  use clmtype
#endif
  use clm_varpar, only : lsmlon, lsmlat, maxpatch, maxpatch_pft
  use clm_varsur, only : numlon, landmask, area, latixy, longxy
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public clm_map   ! Initialize sub-grid mapping 
  public clm_map1d ! Initialize 1d mapping arrays
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein 
!
!EOP
!----------------------------------------------------------------------- 

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map
!
! !INTERFACE:
  subroutine clm_map (vegxy, wtxy) 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
!
! !USES
#if (!defined COUP_CAM)
    use clmtype
#endif
    use lnd_grid, only : surface_grid_init, get_clump_info
#if (defined COUP_CAM)
    use lp_coupling, only : lp_coupling_init
#endif 
#if (defined SPMD)
    use spmdMod, only : masterproc, spmd_init_arrays, gather_data_to_master, &
         proc_gridtot, proc_landtot, proc_coltot, proc_pfttot, npes
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) !PFT type 
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,n,gi,li !indices
    integer  :: ier           !error status 
    integer  :: proc          !SPMD processor id
    integer  :: nlunits       !temporary index
    integer  :: ncolumns      !temporary index
    integer  :: npfts         !temporary index
    integer  :: begg,endg     !beginning and ending 1d grid indices
    integer  :: begl,endl     !beginning and ending 1d land unit indices
    integer  :: begc,endc     !beginning and ending 1d column indices
    integer  :: begp,endp     !beginning and ending 1d pft indices
    integer  :: gindex        !index into global 1d grid array
    integer  :: lindex        !index into global 1d land unit array
    integer  :: cindex        !index into global 1d land column array
    integer  :: pindex        !index into global 1d land pft array
    integer  :: nveg          !number of vegetated patches in a gridcell
    real(r8) :: wtveg         !weight (relative to gridcell) of vegetated landunit
    type(gridcell_type)       , pointer :: g   !pointer to derived subtype
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(model_pstate_type)   , pointer :: mps !pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps !pointer to derived subtype
!------------------------------------------------------------------------

#if (defined SPMD)
    ! Initialize per-processor spmd arrays 

    call spmd_init_arrays
#endif

    call surface_grid_init(wtxy)

#if (defined SPMD)
    ! Determine total subgrid components for each process

    do proc=0,npes-1
       call get_clump_info (proc, proc_gridtot(proc), proc_landtot(proc), &
            proc_coltot(proc), proc_pfttot(proc))
    end do
#endif

#if (defined COUP_CAM)
    ! Initialize mapping between the atmosphere physics chunks and
    ! and the land gridcell clumps

    call lp_coupling_init()
#endif

    ! Set up shorthand for processor beginning and ending 1d indices
    begg = grid1d%beg
    endg = grid1d%end
    begl = land1d%beg
    endl = land1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begp = pfts1d%beg
    endp = pfts1d%end

    ! --------------------------------------------------------------------
    ! Allocate dynamic memory for the CLM derived type hierarchy
    ! For now, the vegetated patches will all be gathered on a single landunit,
    ! with each vegetated type having its own column on that landunit.  The
    ! special patches (urban, lake, wetland, glacier) each get 
    ! their own landunit having a single column and 1 non-vegetated pfts
    ! --------------------------------------------------------------------

    ! Assign local pointer for simpler referencing
    mps => clm%mps
    mps%ngridcells = endg - begg + 1  

    ! Allocate for the array of gridcells to be handled by this process
    ! and for the arrays of area, weight and type for each gridcell.
    allocate(clm%g (mps%ngridcells), &
             clm%ga(mps%ngridcells), &
             clm%gw(mps%ngridcells), &
             clm%gt(mps%ngridcells), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_mapping(): gridcell allocation error'
       call endrun
    end if

    gi = 0
    do j = 1, lsmlat
       do i = 1, numlon(j)
          if (landmask(i,j)==1) then
             gi = gi + 1                           
             if (gi>=begg .and. gi<=endg) then
                clm%g(gi-begg+1)%gps%ixy = i !longitude index
                clm%g(gi-begg+1)%gps%jxy = j !latitude index       
             end if
          endif
       end do
    end do

    ! Set beggining 1d global indices
    ! (subtract 1 so that first index is beggrid for grid index, etc)
    gindex = begg - 1  
    lindex = begl - 1
    cindex = begc - 1
    pindex = begp - 1

    !Loop over all land gridcells
    !Assign local pointer for simpler referencing
    !Note - nlunits below refers to the number of landunits in a gridcell,
    !above it referred to the total number of landunits per processor

    do gi = 1,mps%ngridcells

       ! Set up pointers 
       g => clm%g(gi)					
       gps => g%gps

       ! Increment 1d global index
       gindex = gindex + 1
       gps%index1d = gindex

       ! Get 2d grid indices
       i = gps%ixy
       j = gps%jxy

       ! Set area, weight, and type information for this gridcell.
       ! For now there is only one type of gridcell, value = 1
       ! Still need to resolve the calculation of area for the gridcell
       ! Note that the area, wt, and type are stored in an array at
       ! the higher level in hierarchy, so that the higher level has information
       ! about all members of the next lower level, and the same values
       ! are stored as scalars in each member of the lower level, so
       ! that each member has information about their own area, weight, and type.
       clm%ga(gi) = area(i,j)
!      clm%gw(gi) = clm%ga(gi)/mps%area  
       clm%gt(gi) = 1
       gps%area  = clm%ga(gi)
       gps%wt    = clm%gw(gi)
       gps%itype = clm%gt(gi)

       ! Assign pointer to higher level
       gps%mps => mps

       ! Count the number of special landunit types 
       ! (urban, lake, wetland and glacier) in the gridcell.
       ! Count the number of vegetated patches in the gridcell.
       ! Add one additional landunit if there are any vegetated patches.
       ! Find the total weight (relative to gridcell) for the landunit
       ! containing the vegetated patches.  
       ! This code assumes that wtxy values sum to 1.0 for all patches 
       ! on a land gridcell.
       nlunits = 0
       do m = maxpatch_pft+1, maxpatch
          if (wtxy(i,j,m) > 0.) nlunits = nlunits + 1
       end do
       nveg = 0
       wtveg = 0.
       do m = 1, maxpatch_pft                           
          if (wtxy(i,j,m) > 0.) then
             nveg = nveg + 1
             wtveg = wtveg + wtxy(i,j,m)
          end if
       end do
       if (nveg > 0) nlunits = nlunits + 1
       gps%nlandunits = nlunits

       ! Allocate memory for the landunits on this gridcell
       ! and for the arrays of area, weight and type for each landunit.
       allocate(g%l(nlunits), g%la(nlunits), g%lw(nlunits), g%lt(nlunits), stat=ier)
       if (ier /= 0) then
          write (6,*) 'clm_mapping(): landunit allocation error'
          call endrun
       end if

       ! Determine all landunits for this gridcell.
       li = 0

       ! The first landunit contains all the vegetated patches (if any)
       if (nveg > 0) then
#if (defined NOCOMPETE)
          call landunit_veg_noncompete(nveg, wtveg, wtxy, vegxy, &
                                       g, i, j, gi, li, lindex, cindex, pindex)
#else
          call landunit_veg_compete(nveg, wtveg, wtxy, vegxy, &
                                    g, i, j, gi, li, lindex, cindex, pindex)
#endif
       end if   

       ! The other landunits are special landunits (urban, lake, wetland, glacier).
       do m = maxpatch_pft+1, maxpatch
          if (wtxy(i,j,m) > 0.) then
             call landunit_special(wtxy, g, i, j, m, gi, li, lindex, cindex, pindex)
          end if
       end do

    end do  ! end loop grid cells

  end subroutine clm_map

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_veg_compete
!
! !INTERFACE:
  subroutine landunit_veg_compete (nveg, wtveg, wtxy, vegxy, &
                                   g, i, j, gi, &
                                   li, lindex, cindex, pindex)
!
! !DESCRIPTION: 
! Initialize vegetated landunit with competition
!
! !USES
    use clm_varcon, only : istsoil
#if (!defined COUP_CAM)
    use clmtype
#endif
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nveg  !number of vegetated patches in gridcell
    real(r8), intent(in) :: wtveg !weight relative to gridcell of veg landunit
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) !PFT type 
    type(gridcell_type), pointer :: g  !pointer to derived subtype
    integer, intent(in)    :: i        !2d longitude index
    integer, intent(in)    :: j        !2d latitude index
    integer, intent(in)    :: gi       !gridcell index
    integer, intent(inout) :: li       !landunit index
    integer, intent(inout) :: lindex   !index into global 1d land unit array
    integer, intent(inout) :: cindex   !index into global 1d land column array
    integer, intent(inout) :: pindex   !index into global 1d land pft array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m,ci,pi       !indices
    integer  :: ncolumns      !temporary index
    integer  :: npfts         !temporary index
    integer  :: ier           !error status 
    real(r8) :: weight        !temporary weight
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
!------------------------------------------------------------------------

    ! Increment landunit
    li = li + 1
    
    ! Set up pointers
    l => g%l(li)
    lps => l%lps
    
    ! Set special landunit flag to false
    lps%ifspecial = .false.
    
    ! Increment 1d global land index
    lindex = lindex + 1
    lps%index1d = lindex
    
    ! Set area, weight, and type information for this landunit
    ! the weights for each landunit on the grid cell must add
    ! to one when summed over the grid cell
    weight = wtveg
    g%la(li) = clm%ga(gi) * weight
    g%lw(li) = weight
    g%lt(li) = istsoil
    lps%area  = g%la(li)
    lps%wt    = g%lw(li)
    lps%itype = g%lt(li)
    
    ! Set grid index and weight (relative to grid cell) 
    lps%ixy = i
    lps%jxy = j
    lps%wtxy = lps%area / g%gps%area
    
    ! Assign pointers to higher level
    l%g => g 
    l%a2ls => g%a2ls
    l%a2lf => g%a2lf
    l%lps%gps => g%gps
    
    ! Allocate memory for the array of columns on this landunit
    ! For competing pfts on vegetated landunit, only one column 
    ncolumns = 1
    lps%ncolumns = ncolumns
    allocate(l%c(ncolumns), l%ca(ncolumns), l%cw(ncolumns), l%ct(ncolumns), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_mapping(): vegetated column noncompete allocation error'
       call endrun
    end if
    
    ! Set only one column on landunit
    ci = 1

    ! Set up pointers
    c => l%c(ci) 
    cps => c%cps
    
    ! Increment 1d global column index
    cindex = cindex + 1
    cps%index1d = cindex
    
    ! Set column area, weight (relative to landunit), and type information. 
    ! for this column. For now all columns have the same type, value = 1.
    ! With all pfts on one column, and this as the only column
    ! on the vegetated landunit, the weight is 1.0
    weight = 1.0
    l%ca(ci) = g%la(li) * weight
    l%cw(ci) = weight
    l%ct(ci) = 1
    cps%area  = l%ca(ci)
    cps%wt    = l%cw(ci)
    cps%itype = 1
    
    ! Set grid index and weight (relative to grid cell) 
    cps%ixy = i
    cps%jxy = j
    cps%wtxy = cps%area / g%gps%area

    ! Assign pointers to higher level
    c%l => l 
    c%a2ls => g%a2ls
    c%a2lf => g%a2lf
    c%cps%lps => l%lps
    
    ! Allocate memory for the array of pfts on this column
    ! and for the arrays of area, weight and type for each pft.
    ! With competition on, all pfts are on one column.
    npfts = nveg
    cps%npfts = npfts
    allocate(c%p(npfts), c%pa(npfts), c%pw(npfts), c%pt(npfts), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_mapping(): pft allocation error'
       call endrun
    end if
    
    ! Loop through the pfts for this column.
    ! For now there is only one pft for each column
    ! but later there will be the possibility of multiple
    ! pfts on each column
    pi = 0
    do m = 1,maxpatch_pft
       if (wtxy(i,j,m) > 0.) then
          ! Set pft index with respect to column
          pi = pi+1
          
          ! Set up pointers
          p => c%p(pi)
          pps => p%pps
          
          ! Increment 1d global pft index
          pindex = pindex + 1
          pps%index1d = pindex
          
          ! Set pft area, weight (relative to column), and type information for this pft.
          ! the weight for each pft is the wtxy(i,j,m) value
          ! normalized to the total vegetated patch weight.
          weight = wtxy(i,j,m) / wtveg
          c%pa(pi) = l%ca(ci) * weight
          c%pw(pi) = weight
          c%pt(pi) = vegxy(i,j,m)
          pps%area    = c%pa(pi)
          pps%wt      = c%pw(pi)
          pps%itype   = c%pt(pi)
          
          ! Set grid index, weight (relative to grid cell) 
          ! Set m index (needed for laixy, etc. reference)
          pps%mxy = m
          pps%ixy = i
          pps%jxy = j
          pps%wtxy = pps%area / g%gps%area
          
          ! Assign pointers to higher level
          p%c => c 
          p%a2ls => g%a2ls
          p%a2lf => g%a2lf
          pps%cps => c%cps

       endif ! non-zero weight for this pft
    enddo ! loop through maxpatch_pft
    
  end subroutine landunit_veg_compete
  
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_veg_noncompete
!
! !INTERFACE:
  subroutine landunit_veg_noncompete (nveg, wtveg, wtxy, vegxy, &
                                      g, i, j, gi, &
                                      li, lindex, cindex, pindex)
!
! !DESCRIPTION: 
! Initialize vegetated landunit without competition
!
! !USES
    use clm_varcon, only : istsoil
#if (!defined COUP_CAM)
    use clmtype
#endif

!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nveg  !number of vegetated patches in gridcell
    real(r8), intent(in) :: wtveg !weight relative to gridcell of veg landunit
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) !PFT type 
    type(gridcell_type), pointer :: g  !pointer to derived subtype
    integer, intent(in)    :: i        !2d longitude index
    integer, intent(in)    :: j        !2d latitude index
    integer, intent(in)    :: gi       !gridcell index
    integer, intent(inout) :: li       !land unit index
    integer, intent(inout) :: lindex   !index into global 1d land unit array
    integer, intent(inout) :: cindex   !index into global 1d land column array
    integer, intent(inout) :: pindex   !index into global 1d land pft array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m, ci,pi      !indices
    integer  :: ncolumns      !temporary index
    integer  :: npfts         !temporary index
    integer  :: ier           !error status 
    real(r8) :: weight        !temporary weight
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
!------------------------------------------------------------------------

    ! Increment land unit 
    li = li + 1

    ! Set up pointers
    l => g%l(li)
    lps => l%lps
    
    ! Set special landunit flag to false
    lps%ifspecial = .false.
    
    ! Increment 1d global land index
    lindex = lindex + 1
    lps%index1d = lindex
    
    ! Set area, weight, and type information for this landunit
    ! the weights for each landunit on the grid cell must add
    ! to one when summed over the grid cell
    weight = wtveg
    g%la(li) = clm%ga(gi) * weight
    g%lw(li) = weight
    g%lt(li) = istsoil
    lps%area  = g%la(li)
    lps%wt    = g%lw(li)
    lps%itype = g%lt(li)
    
    ! Set grid index and weight (relative to grid cell) 
    lps%ixy = i
    lps%jxy = j
    lps%wtxy = lps%area / g%gps%area
    
    ! Assign pointers to higher level
    l%g => g
    l%a2ls => g%a2ls
    l%a2lf => g%a2lf
    l%lps%gps => g%gps
    
    ! Allocate memory for the array of columns on this landunit
    ! For non-competing pfts on vegetated landunit, ncolumns = nveg
    ncolumns = nveg
    lps%ncolumns = ncolumns
    allocate(l%c(ncolumns), l%ca(ncolumns), l%cw(ncolumns), l%ct(ncolumns), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_mapping(): vegetated column compete allocation error'
       call endrun
    end if
    
    ! Loop through regular (vegetated) patches, assign one column for each
    ! vegetated patch with non-zero weight. The weights for each column on
    ! the vegetated landunit must add to one when summed over the landunit,
    ! so the wtxy values are taken relative to the total wtveg
    ci = 0
    do m = 1, maxpatch_pft                           
       if (wtxy(i,j,m) > 0.) then
          ! Increment number of columns on landunit
          ci = ci + 1

          ! Set up pointers
          c => l%c(ci) 
          cps => c%cps
          
          ! Increment 1d global column index
          cindex = cindex + 1
          cps%index1d = cindex
          
          ! Set column area, weight (relative to landunit), and type information. 
          ! for this column. For now all columns have the same type as the
          ! associated pft.
          weight = wtxy(i,j,m) / wtveg
          l%ca(ci) = g%la(li) * weight
          l%cw(ci) = weight
          l%ct(ci) = vegxy(i,j,m)  
          cps%area  = l%ca(ci)
          cps%wt    = l%cw(ci)
          cps%itype = vegxy(i,j,m) 
          
          ! Set grid index and weight (relative to grid cell) 
          cps%ixy = i
          cps%jxy = j
          cps%wtxy = cps%area / g%gps%area
          
          ! Assign pointers to higher level
          c%l => l
          c%a2ls => g%a2ls
          c%a2lf => g%a2lf
          c%cps%lps => l%lps
          
          ! Allocate memory for the array of pfts on this column
          ! and for the arrays of area, weight and type for each pft.
          ! Without competition on, each column in the vegetated landunit is 
          ! allocated a single pft. 
          npfts = 1
          cps%npfts = npfts
          allocate (c%p(npfts), c%pa(npfts), c%pw(npfts), c%pt(npfts), stat=ier)
          if (ier /= 0) then
             write (6,*) 'clm_mapping(): pft allocation error'
             call endrun
          end if
          
          ! Each column has its own pft
          pi = 1

          ! Set up pointers
          p => c%p(pi)
          pps => p%pps
          
          ! Increment 1d global pft index
          pindex = pindex + 1
          pps%index1d = pindex
          
          ! Set area, weight (relative to column) and type information for this pft
          ! For now, a single pft per column, so weight = 1
          ! pft type comes from the m dimension of wtxy()
          weight = 1.0/npfts
          c%pa(pi) = l%ca(ci) * weight
          c%pw(pi) = weight
          c%pt(pi) = vegxy(i,j,m)
          pps%area  = c%pa(pi)
          pps%wt    = c%pw(pi)
          pps%itype = c%pt(pi)
          
          ! Set grid index, weight (relative to grid cell) 
          ! and m index (needed for laixy, etc. reference)
          pps%mxy = m
          pps%ixy = i
          pps%jxy = j
          pps%wtxy = pps%area / g%gps%area
          
          ! Assign pointers to higher level
          p%c => c
          p%a2ls => g%a2ls
          p%a2lf => g%a2lf
          pps%cps => c%cps
          
       end if   ! end if non-zero weight	
    end do   ! end loop through the possible vegetated patch indices

  end subroutine landunit_veg_noncompete

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_special
!
! !INTERFACE:
  subroutine landunit_special (wtxy, g, i, j, m, gi, &
	                       li, lindex, cindex, pindex)
!
! !DESCRIPTION: 
! Initialize special landunits (urban, lake, wetland, glacier)
!
! !USES
#if (!defined COUP_CAM)
    use clmtype
#endif
    use pftvarcon, only : noveg
    use clm_varcon, only : istice, istwet, istdlak, isturb
    use clm_varpar, only : npatch_lake, npatch_wet, npatch_urban, npatch_glacier
!
! !ARGUMENTS:
    implicit none
    type(gridcell_type), pointer :: g   !pointer to derived subtype
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
    integer, intent(in) :: i            !2-dim longitude index
    integer, intent(in) :: j            !2-dim latitude index
    integer, intent(in) :: m            !2-dim PFT patch index  
    integer, intent(inout) :: gi        !gridcell index
    integer, intent(inout) :: li        !landunit index
    integer, intent(inout) :: lindex    !index into global 1d land unit array
    integer, intent(inout) :: cindex    !index into global 1d land column array
    integer, intent(inout) :: pindex    !index into global 1d land pft array
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ci,pi         !indices
    integer  :: ncolumns      !temporary index
    integer  :: npfts         !temporary index
    integer  :: ier           !error status 
    real(r8) :: weight        !temporary weight
    type(landunit_type)       , pointer :: l   !pointer to derived subtype
    type(column_type)         , pointer :: c   !pointer to derived subtype
    type(pft_type)            , pointer :: p   !pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps !pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps !pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps !pointer to derived subtype
!------------------------------------------------------------------------

    ! Increment landunit for this gridcell
    li = li + 1

    ! Set up pointers
    l => g%l(li)  
    lps => l%lps
    
    ! Set special landunit flag to true
    lps%ifspecial = .true.

    ! Increment 1d global land index (increment total landunits)
    lindex = lindex + 1
    lps%index1d = lindex
    
    ! Set area, weight, and type information for this landunit.
    weight = wtxy(i,j,m)
    g%la(li) = clm%ga(gi) * weight
    g%lw(li) = weight

    ! Define soil type
    if (m == npatch_lake) then         !deep lake (from pctlak)
       g%lt(li) = istdlak 
    else if (m == npatch_wet) then     !wetland (from pctwet)
       g%lt(li) = istwet
    else if (m == npatch_glacier) then !glacier (from pctgla)
       g%lt(li) = istice
    else if (m == npatch_urban) then   !urban (from pcturb) 
       g%lt(li) = isturb
    else                               !error
       write(6,*)'special landunit are currently only:', & 
            ' deep lake, wetland, glacier or urban)'
       call endrun()
    endif
    lps%area = g%la(li)
    lps%wt   = g%lw(li)
    lps%itype= g%lt(li)
    
    ! Set grid index and weight (relative to grid cell) 
    lps%ixy = i
    lps%jxy = j
    lps%wtxy = lps%area / g%gps%area
    
    ! Assign pointers to higher level
    l%g => g
    l%a2ls => g%a2ls
    l%a2lf => g%a2lf
    lps%gps => g%gps
    
    ! Allocate memory for the array of columns on this landunit
    ! For the special landunits, there is only one column 
    ! For now there is only one type of column, value = 1.
    ! Later, the age classes will be implemented on different
    ! columns within the same landunit, so the column type
    ! will correspond to an age class
    ncolumns = 1
    lps%ncolumns = ncolumns
    allocate(l%c(ncolumns), l%ca(ncolumns), l%cw(ncolumns), l%ct(ncolumns), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_mapping(): column allocation error'
       call endrun
    end if
    
    ! Loop through columns for this landunit
    ! to set area, weight, and type information for this landunit
    ! We know that there is only one column for the special
    ! landunits, but the loop is included for consistency.
    do ci = 1,ncolumns
       
       ! Set up pointers
       c => l%c(ci) 
       cps => c%cps
       
       ! Increment 1d global column index
       cindex = cindex + 1
       cps%index1d = cindex
       
       ! Set area, weight, and type information for this column
       ! For now all columns have the same type, value = 1
       weight = 1.0/ncolumns
       l%ca(ci) = g%la(li) * weight
       l%cw(ci) = weight
       l%ct(ci) = 1
       cps%area  = l%ca(ci)
       cps%wt    = l%cw(ci)
       cps%itype = l%ct(ci)
       
       ! Set grid index and weight (relative to grid cell) 
       cps%ixy = i
       cps%jxy = j
       cps%wtxy = cps%area / g%gps%area
       
       ! Assign pointers to higher level
       c%l => l
       c%a2ls => g%a2ls
       c%a2lf => g%a2lf
       cps%lps => l%lps
       
       ! Set one non-vegetated pft for this column 
       npfts = 1
       cps%npfts = npfts
       allocate(c%p(npfts), c%pa(npfts), c%pw(npfts), c%pt(npfts), stat=ier)
       if (ier /= 0) then
          write (6,*) 'clm_mapping(): column allocation error'
          call endrun
       end if
       
       ! Set up pointers for non-vegetated pft
       pi = 1
       p => c%p(pi) 
       pps => p%pps
       
       ! Increment 1d global pft index
       pindex = pindex + 1
       pps%index1d = pindex
       
       ! Set area, weight (relative to column), and type information 
       ! for this non-vegetated pft
       weight = 1.0/npfts
       c%pa(pi) = l%ca(ci) * weight
       c%pw(pi) = weight
       c%pt(pi) = noveg
       pps%area    = c%pa(pi)
       pps%wt      = c%pw(pi)
       pps%itype   = c%pt(pi)
       
       ! Set grid index, weight (relative to grid cell) and  
       ! m index (needed for laixy, etc. reference)
       pps%mxy = m
       pps%ixy = i
       pps%jxy = j
       pps%wtxy = pps%area / g%gps%area
       
       ! Assign pointers to higher level
       p%c => c
       p%a2ls => g%a2ls
       p%a2lf => g%a2lf
       pps%cps => c%cps
       
    end do   ! end loop through ncolumns

  end subroutine landunit_special

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map1d
!
! !INTERFACE:
  subroutine clm_map1d()
!
! !DESCRIPTION: 
! Set up 1d array of weights and indices for xy mapping
! Note: if DGVM is defined, weights are updated in DGVM mode 
!
! !USES:
#if (!defined COUP_CAM)
    use clmtype
#endif
    use clmpoint, only : cpoint, ppoint
#if (defined SPMD)
    use spmdMod, only : masterproc, spmd_init_arrays, gather_data_to_master, &
         MPI_INTEGER, mpicom
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: gi,li,ci,pi            !indices
    integer  :: gindex                 !index into global 1d gridcell array
    integer  :: lindex                 !index into global 1d landunit array
    integer  :: cindex                 !index into global 1d column array
    integer  :: pindex                 !index into global 1d pft array
    integer  :: ixy,jxy                !gridcell indices
    real(r8) :: londeg                 !gridcell longitude (degrees)
    real(r8) :: latdeg                 !gridcell latitude (degrees)
    integer  :: begc,endc              !beginning and ending 1d column indices
    integer  :: begp,endp              !beginning and ending 1d pft indices
    integer  :: numg                   !total numger of grid cells
    integer  :: numl                   !total numger of land units
    integer  :: numc                   !total numger of columns
    integer  :: nump                   !total numger of pfts
    integer  :: ityplun                !landunit water type (soil,lake,wetland,glacier,urban)
    integer  :: ier                    !error status 
    type(gridcell_type), pointer :: g  !pointer to derived subtype
    type(landunit_type), pointer :: l  !pointer to derived subtype
    type(column_type)  , pointer :: c  !pointer to derived subtype
    type(pft_type)     , pointer :: p  !pointer to derived subtype
#if (defined SPMD)
    real(r8), pointer :: rloc(:)       !temporaries for mpi gather
    integer , pointer :: iloc(:)       !temporaries for mpi gather
    real(r8), pointer :: rglob(:)      !temporaries for mpi gather
    integer , pointer :: iglob(:)      !temporaries for mpi gather
#endif
!-------------------------------------------------------------------------

    grid1d%name = 'gridcell'
    land1d%name = 'landunit'
    cols1d%name = 'column'
    pfts1d%name = 'pft'

    numg = grid1d%num
    numl = land1d%num
    numc = cols1d%num
    nump = pfts1d%num

    allocate (grid1d%wtxy(numg), &
              grid1d%ixy(numg), &
              grid1d%jxy(numg), &
              grid1d%latdeg(numg), &
              grid1d%londeg(numg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): grid1d allocation error'
       call endrun
    end if

    allocate (land1d%wtxy(numl), &
              land1d%ixy(numl), &
              land1d%jxy(numl), &
              land1d%gindex(numl), &
              land1d%latdeg(numl), &
              land1d%londeg(numl), &
              land1d%ityplun(numl), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): land1d allocation error'
       call endrun
    end if

    allocate (cols1d%wtxy(numc), &
              cols1d%wtlnd(numc), &
              cols1d%ixy(numc), &
              cols1d%jxy(numc), &
              cols1d%gindex(numc), &
              cols1d%lindex(numc), &
              cols1d%latdeg(numc), &
              cols1d%londeg(numc), &
              cols1d%ityplun(numc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): cols1d allocation error'
       call endrun
    end if

    allocate (pfts1d%wtxy(nump), &
              pfts1d%wtlnd(nump), &
              pfts1d%wtcol(nump), &
              pfts1d%ixy(nump), &
              pfts1d%jxy(nump), &
              pfts1d%mxy(nump), &
              pfts1d%gindex(nump), &
              pfts1d%lindex(nump), &
              pfts1d%cindex(nump), &
              pfts1d%latdeg(nump), &
              pfts1d%londeg(nump), &
              pfts1d%ityplun(nump), &
              pfts1d%itypveg(nump), stat=ier)
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): pfts1d allocation error'
       call endrun
    end if

    do gi = 1, clm%mps%ngridcells
       g => clm%g(gi)					
       gindex = g%gps%index1d 
       ixy = g%gps%ixy
       jxy = g%gps%jxy
       latdeg = latixy(ixy,jxy)
       londeg = longxy(ixy,jxy)
       grid1d%wtxy(gindex) = 1._r8
       grid1d%ixy(gindex) = ixy
       grid1d%jxy(gindex) = jxy
       grid1d%latdeg(gindex) = latdeg
       grid1d%londeg(gindex) = londeg

       do li = 1, g%gps%nlandunits
          l => g%l(li)
          lindex = l%lps%index1d 
          ityplun = l%lps%itype
          land1d%ixy(lindex) = ixy
          land1d%jxy(lindex) = jxy
          land1d%latdeg(lindex) = latdeg
          land1d%londeg(lindex) = londeg
          land1d%gindex(lindex) = gindex
          land1d%ityplun(lindex) = ityplun
          land1d%wtxy(lindex) = l%lps%area / g%gps%area

          do ci = 1, l%lps%ncolumns
             c => l%c(ci) 
             cindex = c%cps%index1d 
             cols1d%ixy(cindex) = ixy
             cols1d%jxy(cindex) = jxy
             cols1d%latdeg(cindex) = latdeg
             cols1d%londeg(cindex) = londeg
             cols1d%gindex(cindex) = gindex
             cols1d%lindex(cindex) = lindex
             cols1d%ityplun(cindex) = ityplun
             cols1d%wtxy(cindex) = c%cps%area / g%gps%area
             cols1d%wtlnd(cindex) = c%cps%area / l%lps%area

             do pi = 1, c%cps%npfts
                p => c%p(pi)
                pindex = p%pps%index1d 
                pfts1d%ixy(pindex) = ixy
                pfts1d%jxy(pindex) = jxy
                pfts1d%mxy(pindex) = p%pps%mxy
                pfts1d%latdeg(pindex) = latdeg
                pfts1d%londeg(pindex) = londeg
                pfts1d%gindex(pindex) = gindex
                pfts1d%lindex(pindex) = lindex
                pfts1d%cindex(pindex) = cindex
                pfts1d%ityplun(pindex) = ityplun
                pfts1d%itypveg(pindex) = p%pps%itype
                pfts1d%wtxy(pindex) = p%pps%area / g%gps%area
                pfts1d%wtlnd(pindex) = p%pps%area / l%lps%area
                pfts1d%wtcol(pindex) = p%pps%area / c%cps%area
             end do
          end do
       end do
    end do

#if (defined SPMD)

    ! gridcell gather 

    allocate(rglob(numg), iglob(numg))
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): rglob,iglob (numg) allocation error'
       call endrun
    end if

    rloc => grid1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%wtxy(:) = rglob(:)

    iloc => grid1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=grid1d%name)
    call mpi_bcast(iglob, size(iglob), MPI_INTEGER, 0, mpicom, ier)
    grid1d%ixy(:) = iglob(:)

    iloc => grid1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=grid1d%name)
    call mpi_bcast(iglob, size(iglob), MPI_INTEGER, 0, mpicom, ier)
    grid1d%jxy(:) = iglob(:)

    rloc => grid1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%latdeg(:) = rglob(:)

    rloc => grid1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=grid1d%name)
    if (masterproc) grid1d%londeg(:) = rglob(:)

    deallocate(rglob, iglob)

    ! landunit gather

    allocate(rglob(numl),iglob(numl))
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): rglob,iglob (numl) allocation error'
       call endrun
    end if

    rloc => land1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%wtxy(:) = rglob(:)

    iloc => land1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%ixy(:) = iglob(:)

    iloc => land1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%jxy(:) = iglob(:)

    iloc => land1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%gindex(:) = iglob(:)

    rloc => land1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%latdeg(:) = rglob(:)

    rloc => land1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=land1d%name)
    if (masterproc) land1d%londeg(:) = rglob(:)

    iloc => land1d%ityplun
    call gather_data_to_master(iloc, iglob, clmlevel=land1d%name)
    if (masterproc) land1d%ityplun(:) = iglob(:)

    deallocate(rglob, iglob)

    ! column gather

    allocate(rglob(numc), iglob(numc))
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): rglob,iglob (numc) allocation error'
       call endrun
    end if

    rloc => cols1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%wtxy(:) = rglob(:)

    rloc => cols1d%wtlnd
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%wtlnd(:) = rglob(:)

    iloc => cols1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%ixy(:) = iglob(:)

    iloc => cols1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%jxy(:) = iglob(:)

    iloc => cols1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%gindex(:) = iglob(:)

    iloc => cols1d%lindex
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%lindex(:) = iglob(:)

    rloc => cols1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%latdeg(:) = rglob(:)

    rloc => cols1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%londeg(:) = rglob(:)

    iloc => cols1d%ityplun
    call gather_data_to_master(iloc, iglob, clmlevel=cols1d%name)
    if (masterproc) cols1d%ityplun(:) = iglob(:)

    deallocate(rglob, iglob)

    ! pft gather 

    allocate(rglob(nump),iglob(nump))
    if (ier /= 0) then
       write (6,*) 'clm_map1d(): rglob,iglob (nump) allocation error'
       call endrun
    end if

    rloc => pfts1d%wtxy
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtxy(:) = rglob(:)

    rloc => pfts1d%wtlnd
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtlnd(:) = rglob(:)

    rloc => pfts1d%wtcol
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%wtcol(:) = rglob(:)

    iloc => pfts1d%ixy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%ixy(:) = iglob(:)

    iloc => pfts1d%jxy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%jxy(:) = iglob(:)

    iloc => pfts1d%mxy
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%mxy(:) = iglob(:)

    iloc => pfts1d%gindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%gindex(:) = iglob(:)

    iloc => pfts1d%lindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%lindex(:) = iglob(:)

    iloc => pfts1d%cindex
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%cindex(:) = iglob(:)

    rloc => pfts1d%latdeg
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%latdeg(:) = rglob(:)

    rloc => pfts1d%londeg
    call gather_data_to_master(rloc, rglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%londeg(:) = rglob(:)

    iloc => pfts1d%ityplun
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%ityplun(:) = iglob(:)

    iloc => pfts1d%itypveg
    call gather_data_to_master(iloc, iglob, clmlevel=pfts1d%name)
    if (masterproc) pfts1d%itypveg(:) = iglob(:)

    deallocate(rglob,iglob)

#endif

#if (defined NOCOMPETE)
    ! When competition mode is not on, it is assumed that each column has 
    ! one pft and indices into the 1d column vector must equal indices 
    ! into the 1d pft vector 

    begc = cols1d%beg
    endc = cols1d%end
    begp = pfts1d%beg
    nump = pfts1d%num
    endp = pfts1d%end
    numc = cols1d%num

    if (numc /= nump) then
       write(6,*)'total number of columns must equal total number of pfts'
       write(6,*)'total number of columns = ',numc
       write(6,*)'total number of pfts = ',nump
       call endrun
    endif
    if (begc /= begp) then
       write(6,*)'beginning index of columns must equal beginning index of pfts'
       write(6,*)'beginning index of columns = ',begc
       write(6,*)'beginning index of pfts = ',begp
       call endrun
    endif
    if (endc /= endp) then
       write(6,*)'ending index of columns must equal ending index of pfts'
       write(6,*)'ending index of columns = ',endc
       write(6,*)'ending index of pfts = ',endp
       call endrun
    endif
    do pi = begp,endp
       p => ppoint(pi)%p
       c => cpoint(pi)%c
       if (p%pps%index1d /= c%cps%index1d) then
          write(6,*)'competition is not invoked: column indices must match pft indices'
          write(6,*)'indexp= ',p%pps%index1d,' indexc= ',c%cps%index1d
          call endrun
       end if
    end do
#endif

  end subroutine clm_map1d

end module clm_mapping





