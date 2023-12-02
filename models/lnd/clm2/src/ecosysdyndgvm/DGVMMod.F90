#include <misc.h>
#include <preproc.h>

module DGVMMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: DGVMMod
!
! !DESCRIPTION: 
! Module containing routines to drives the annual portion of lpj 
! (called once per year), reset variables related to lpj, 
! and initialize/Reset time invariant dgvm variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: lpj       ! drives the annual portion of lpj, called once per year
  public :: lpjreset1 ! resets variables related to lpj
  public :: lpjreset2 ! resets variables related to lpj
  public :: resetTimeConstDGVM ! Initialize/Reset time invariant dgvm variables
  public :: resetWeightsDGVM   ! Reset DGVM subgrid weights and areas
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lpj
!
! !INTERFACE:
  subroutine lpj ()
!
! !DESCRIPTION: 
! Drives the annual portion of lpj, called once per year
!
! !USES:
    use clmtype
    use clmpoint
    use globals
    use clm_varpar, only : maxpatch_pft, lsmlon , lsmlat
    use clm_varcon, only : spval
    use time_manager, only : get_curr_date, get_ref_date
    use histFileDGVMMod, only : ncid, beg3d, len3d, beg4d, len4d, &
         afirefrac_id, acfluxfire_id, bmfm_id, afmicr_id, &
         histcrt_dgvm, histwrt_dgvm
    use system_messages
    use ReproductionMod, only : Reproduction
    use TurnoverMod, only : Turnover
    use AllocationMod, only : Allocation
    use LightMod, only : Light
    use MortalityMod, only : Mortality
    use FireMod, only : Fire
    use EstablishmentMod, only : Establishment
    use KillMod, only : Kill
#if (defined SPMD)
    use spmdMod, only : masterproc, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: gi,li,ci,pi,i,j,m  !indices
    integer  :: begg, endg, numg   !1d gridcell beginning, ending, size indices
    integer  :: begp, endp, nump   !1d pft beginning, ending, size indices
    integer  :: ivt                !vegetation type index
    integer  :: yr                 !year (0 -> ...)
    integer  :: mon                !month (1 -> 12)
    integer  :: day                !day (1 -> 31)
    integer  :: sec                !seconds 
    integer  :: kyr                !used in routine climate20' below
    integer  :: ncdate             !current date
    integer  :: nbdate             !base date (reference date)
    integer  :: ier                !error status
    real(r8) :: wscal
    real(r8) :: mort_max
    real(r8) :: acfluxfire
    real(r8) :: afirefrac
    real(r8) :: bminc
    real(r8), pointer :: afirefrac_grid(:)         !for history write
    real(r8), pointer :: acfluxfire_grid(:)        !for history write
    real(r8), pointer :: bmfm_grid(:,:)            !for history write
    real(r8), pointer :: afmicr_grid(:,:)          !for history write
    real(r8), pointer :: rglob1d(:)                !temporary for history write
    real(r8), pointer :: rglob2d(:,:)              !temporary for history write
    real(r8), pointer :: rbuf2d_xy(:,:)            !temporary for history write
    real(r8), pointer :: rbuf3d_xy(:,:,:)          !temporary for history write
    type(gridcell_type)       , pointer :: g       !local pointer to derived subtype
    type(landunit_type)       , pointer :: l       !local pointer to derived subtype
    type(column_type)         , pointer :: c       !local pointer to derived subtype
    type(pft_type)            , pointer :: p       !local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps     !local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps     !local pointer to derived subtype
    type(pft_epc_type)        , pointer :: pepc    !local pointer to derived subtype
    type(pft_dgvstate_type)   , pointer :: pdgvs   !local pointer to derived subtype
    type(pft_dgvepc_type)     , pointer :: pdgvepc !local pointer to derived subtype
!-----------------------------------------------------------------------

    ! Determine per-process beginning and ending indices

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num

    begg = grid1d%beg
    endg = grid1d%end
    numg = grid1d%num

    ! Allocate and initialize dynamic memory for land output variables

    allocate (afirefrac_grid(begg:endg), &	
         acfluxfire_grid(begg:endg), &	
         bmfm_grid(maxpatch_pft,begg:endg), &	
         afmicr_grid(maxpatch_pft,begg:endg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'lpj: allocation error for afirefrac_grid ',&
            ' acfluxfire_grid, bmfm_grid, afmicr_grid '
       call endrun
    end if

    afirefrac_grid(:) = 0.0 
    acfluxfire_grid(:) = 0.0 
    bmfm_grid(:,:) = 0.0
    afmicr_grid(:,:) = 0.0

    ! kyr used in routine climate20' below
    ! NB: at end of first year, kyr = 2

    call get_curr_date (yr, mon, day, sec)  
    ncdate = yr*10000 + mon*100 + day
    call get_ref_date (yr, mon, day, sec)
    nbdate = yr*10000 + mon*100 + day
    kyr = ncdate/10000 - nbdate/10000 + 1   
    write(6,*)'LPJ: ncdate= ',ncdate,' nbdate= ',nbdate,' kyr= ',kyr

    ! Loop over pfts
    do pi = begp,endp
       
       ! Assign local pointers to pft-level derived subtypes
       p => ppoint(pi)%p
       pps => p%pps
       pepc => p%pepc
       pdgvs => p%pdgvs
       pdgvepc => p%pdgvepc

       ! *************************************************************************
       ! S. Levis version of LPJ's routine climate20 - 'Returns' tmomin20 and agdd20 
       ! for use in routine bioclim, which I have placed in routine Establishment
       ! Instead of 20-yr running mean of coldest monthly temperature,
       ! use 20-yr running mean of minimum 10-day running mean
       
       if (kyr == 2) pdgvs%tmomin20 = pdgvs%t_mo_min
       if (kyr == 2) pdgvs%agdd20 = pdgvs%agdd
       pdgvs%tmomin20 = (19.0 * pdgvs%tmomin20 + pdgvs%t_mo_min) / 20.0
       pdgvs%agdd20   = (19.0 * pdgvs%agdd20   + pdgvs%agdd    ) / 20.0
       ! *************************************************************************
       
       bminc = pdgvs%bm_inc                        ![gC/m2 patch] for output
       pdgvs%bm_inc = pdgvs%bm_inc * pdgvs%fpcgrid ![gC/m2 cell vegetated area]
       
       ! Determine grid values of npp and microbial respiration
       
       m = pps%mxy                      !m index
       if (m <= maxpatch_pft) then
          gi = pfts1d%gindex(pi)
          bmfm_grid(m,gi) = bminc                          ![gC/m2 patch] for output 
          afmicr_grid(m,gi) = pdgvs%afmicr * pdgvs%fpcgrid ![gC/m2 cell veg'd area]
       endif
       
       if (pps%itype > 0) then
          
          ! Returns updated bm_inc, litterag
          
          call Reproduction (pdgvs%bm_inc, pdgvs%litterag, pdgvs%present)
          
          ! Returns turnover_ind and updated litterag,bg, l,s,h,rm_ind
          
          call Turnover (pdgvepc%l_turn, pdgvepc%s_turn, pdgvepc%r_turn, &
                         pdgvs%litterag, pdgvs%litterbg, pdgvs%lm_ind  , &
                         pdgvs%sm_ind  , pdgvs%hm_ind  , pdgvs%rm_ind  , &
                         pdgvs%nind    , pdgvs%present , pdgvs%turnover_ind)
          
          ! Returns updated litterag, bg, and present
          
          call Kill (pdgvs%bm_inc, pdgvs%litterag, pdgvs%litterbg, &
                     pdgvs%lm_ind, pdgvs%sm_ind  , pdgvs%hm_ind  , &
                     pdgvs%rm_ind, pdgvs%nind    , pdgvs%present , &
                     pdgvepc%tree)
          
          ! Returns lai_ind, lai_inc, updates crownarea, htop, 
          ! l, s, h, rm_ind and litterag,bg
          
          if (pdgvs%annpsnpot > 0.0) then
             wscal = pdgvs%annpsn/pdgvs%annpsnpot
          else
             wscal = 1.0
          end if
          
          call Allocation (pdgvepc       , pdgvs%nind    , pdgvs%bm_inc   , &
                           pdgvepc%tree  , pps%htop      , pepc%sla       , &
                           pdgvs%lm_ind  , pdgvs%sm_ind  , pdgvs%hm_ind   , &
                           pdgvs%rm_ind  , pdgvs%present , pdgvs%crownarea, &
                           pdgvs%litterag, pdgvs%litterbg, pdgvs%lai_ind  , &
                           pdgvs%fpcgrid , pdgvs%fpcinc  , wscal          )
       end if
    end do  !end of pft loop

    ! Returns lm,rm_ind, fpcgrid, nind, litterag,bg via modules
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    call Light ()

    ! Loop over pfts
    do pi = begp,endp
       
       ! Assign local pointers to pft-level derived subtypes
       p => ppoint(pi)%p
       pps => p%pps
       pepc => p%pepc
       pdgvs => p%pdgvs
       pdgvepc => p%pdgvepc
       
       ! Obtain updated present, nind, litterag and bg
       
       ivt = pps%itype
       if (ivt > 0) then
          
          if (ivt==3 .or. ivt == 8) then
             mort_max = 0.03 !testing diff values (add to Ecosystemini?)
          else
             mort_max = 0.01 !original value for all pfts
          end if
          
          call Mortality (pdgvs%bm_inc  , pdgvs%nind   , pdgvs%turnover_ind, &
                          pdgvs%lm_ind  , pdgvs%sm_ind , pdgvs%hm_ind      , &
                          pdgvs%rm_ind  , pepc%sla     , pdgvs%litterag    , &
                          pdgvs%litterbg, pdgvs%present, pdgvepc%tree         , &
                          pdgvs%agddtw  , mort_max       )
          
          ! Returns updated litterag and nind
          
          call Fire (pdgvs%firelength, pdgvs%litterag, pdgvs%present, &
                     pdgvepc%tree    , pdgvepc%resist, pdgvs%nind   , &
                     pdgvs%lm_ind    , pdgvs%sm_ind  , pdgvs%hm_ind , &
                     pdgvs%rm_ind    , afirefrac     , pdgvs%fpcgrid, &
                     acfluxfire)
          
       else
          
          afirefrac = 0.0
          acfluxfire = 0.0
          
       end if
       
       ! Determine grid point values
       
       gi = pfts1d%gindex(pi)
       afirefrac_grid(gi) = afirefrac_grid(gi) + afirefrac*pdgvs%fpcgrid
       acfluxfire_grid(gi) = acfluxfire_grid(gi) + acfluxfire

    end do  !end of pft loop
             
    ! Returns updated present, nind, *m_ind, crownarea, fpcgrid, htop, litter*g, itypveg
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    call Establishment ()

    ! -----------------------------------------------------------------------
    ! Write grid info, time info and first set of variables to  to DGVM
    ! history file
    ! -----------------------------------------------------------------------

    ! Create DGVM history file

    if (masterproc) call histcrt_dgvm()

    ! Write grid and time information to DGVM history file

    call histwrt_dgvm()	

    ! First write of first set of output fields to DGVM history file
    ! Rest of fields are written out in lpjreset2.F90
    ! First map land or patch points to xy arrays and then write to history file

    ! Allocate dynamic memory for history writes

#if (defined SPMD)
    allocate(rglob1d(numg), rglob2d(maxpatch_pft,numg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'lpj: allocation error for rglob1d, rglob2d'
       call endrun
    end if
#endif
    if (masterproc) then
       allocate (rbuf2d_xy(lsmlon,lsmlat), rbuf3d_xy(lsmlon,lsmlat,maxpatch_pft), &
            stat=ier)
       if (ier /= 0) then
          write (6,*) 'lpj: allocation error for rbuf2d_xy, rbuf3d_xy'
          call endrun
       end if
    endif

    ! Write out afirefrac

#if (defined SPMD)
    call gather_data_to_master(afirefrac_grid, rglob1d, clmlevel=grid1d%name) 
#else
    rglob1d => afirefrac_grid
#endif
    if (masterproc) then
       rbuf2d_xy(:,:)= spval
       do gi = 1, grid1d%num
          i = grid1d%ixy(gi)    
          j = grid1d%jxy(gi)    
          rbuf2d_xy(i,j) = rglob1d(gi)
       end do
       call wrap_put_vara_realx (ncid, afirefrac_id, beg3d, len3d, rbuf2d_xy)
    endif

    ! Write out acfluxfire

#if (defined SPMD)
    call gather_data_to_master(acfluxfire_grid, rglob1d, clmlevel=grid1d%name) 
#else
    rglob1d => acfluxfire_grid
#endif
    if (masterproc) then
       rbuf2d_xy(:,:)= spval
       do gi = 1, grid1d%num
          i = grid1d%ixy(gi)    
          j = grid1d%jxy(gi)    
          rbuf2d_xy(i,j) = rglob1d(gi)
       end do
       call wrap_put_vara_realx (ncid, acfluxfire_id, beg3d, len3d, rbuf2d_xy)
    endif

    ! Write out bmfm

#if (defined SPMD)
    call gather_data_to_master(bmfm_grid, rglob2d, clmlevel=grid1d%name) 
#else
    rglob2d => bmfm_grid
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:) = spval
       do gi = 1, grid1d%num
          i = grid1d%ixy(gi)    
          j = grid1d%jxy(gi)    
          do m = 1, maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2d(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, bmfm_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out afmicr

#if (defined SPMD)
    call gather_data_to_master(afmicr_grid, rglob2d, clmlevel=grid1d%name) 
#else
    rglob2d => afmicr_grid
#endif           
    if (masterproc) then
       rbuf3d_xy(:,:,:) = spval
       do gi = 1, grid1d%num
          i = grid1d%ixy(gi)    
          j = grid1d%jxy(gi)    
          do m = 1, maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2d(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, afmicr_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Deallocate dynamic memory

    deallocate(afirefrac_grid, acfluxfire_grid, bmfm_grid, afmicr_grid)
    if (masterproc) deallocate(rbuf2d_xy, rbuf3d_xy)
#if (defined SPMD)
    deallocate(rglob1d, rglob2d)
#endif

  end subroutine lpj

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lpjreset1
!
! !INTERFACE:
  subroutine lpjreset1 (caldayp1, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION: 
! Resets variables related to lpj!
!
! !USES:
    use clmtype
    use clmpoint
    use SurfaceAlbedoMod, only : SurfaceAlbedo
    use EcosystemDynDGVMMod, only : EcosystemDyn
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: caldayp1 !calendar day at Greenwich (1.00, ..., 365.99) for nstep+1
    real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ci,pi                  !indices 
    integer :: begc,endc              !column 1d indices
    integer :: begp,endp              !pft 1d indices
    type(gridcell_type), pointer :: g !local pointer to derived subtype
    type(landunit_type), pointer :: l !local pointer to derived subtype
    type(column_type)  , pointer :: c !local pointer to derived subtype
    type(pft_type)     , pointer :: p !local pointer to derived subtype
    type(pft_pstate_type)  , pointer :: pps   !local pointer to derived subtype
    type(pft_dgvstate_type), pointer :: pdgvs !local pointer to derived subtype
!-----------------------------------------------------------------------

    ! Set 1d indices

    begc = cols1d%beg
    endc = cols1d%end

    begp = pfts1d%beg
    endp = pfts1d%end

    ! Reset a few variables here at the very end of the year

    do pi = begp,endp
       
       ! Assign local pointers to pft-level derived subtypes
       p => ppoint(pi)%p
       pps => p%pps
       pdgvs => p%pdgvs

       ! Reset variables
       pdgvs%annpsn     = 0.
       pdgvs%annpsnpot  = 0.
       pdgvs%bm_inc     = 0.
       pdgvs%afmicr     = 0.
       pdgvs%firelength = 0.
       pdgvs%agddtw     = 0.
       pdgvs%agdd       = 0.
       pdgvs%t10min     = 1.0e+36
       pdgvs%t_mo_min   = 1.0e+36

    end do  !end of pft loop

    ! call EcosystemDyn and SurfaceAlbedo because need information
    ! for first timestep of next year

    do ci = begc,endc
       c => cpoint(ci)%c

       if (.not. c%cps%lps%lakpoi) then
          call EcosystemDyn (c, .false., .true.)
       endif
       
       call SurfaceAlbedo (c, caldayp1, eccen, obliqr, lambm0, mvelpp)

    end do   

    call ResetTimeConstDGVM ()

  end subroutine lpjreset1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lpjreset2
!
! !INTERFACE:
  subroutine lpjreset2 ()
!
! !DESCRIPTION: 
! Resets variables related to lpj
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varctl, only : archive_dir, mss_wpass, mss_irt
    use clm_varpar, only : maxpatch_pft, lsmlon , lsmlat
    use clm_varcon, only : cpliq, cpice, istsoil, spval !<- for diagnostic
    use clm_varsur, only : landmask
    use histFileDGVMMod, only : ncid, dgvm_fn, beg3d, len3d, beg4d, len4d, &
         fpcgrid_id, itypveg_id, lmind_id, rmind_id, &
         smind_id, hmind_id, nind_id, begwater_id, endwater_id, &
         begenergy_id, endenergy_id
    use mapxy, only : vec2xy
    use fileutils, only : set_filename, putfil
    use time_manager, only : get_nstep
    use system_messages
#if (defined SPMD)
    use spmdMod, only : masterproc, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
    use shr_sys_mod, only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: gi,li,ci,pi,i,j,m                 !indices
    integer  :: ier                               !error status 
    real(r8) :: weight                            !temporary for weight update
    real(r8) :: cv(-nlevsno+1:nlevsoi)            !temporary for diagnostics
    real(r8), pointer :: rbuf3d_xy(:,:,:)         !history output
    integer , pointer :: ibuf3d_xy(:,:,:)         !history output
    real(r8), pointer :: water(:)                 !diagnostic column water 
    real(r8), pointer :: energy(:)                !diagnostic column energy
    real(r8) :: water_xy(lsmlon,lsmlat)           !diagnostic gridcell water
    real(r8) :: energy_xy(lsmlon,lsmlat)          !diagnostic gridcell energy
    type(gridcell_type)  , pointer :: g           !local pointer to derived subtype
    type(landunit_type)  , pointer :: l           !local pointer to derived subtype
    type(column_type)    , pointer :: c           !local pointer to derived subtype
    type(pft_type)       , pointer :: p           !local pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps    !local pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps    !local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps    !local pointer to derived subtype
    type(column_wstate_type)  , pointer :: cws    !local pointer to derived subtype
    type(column_estate_type)  , pointer :: ces    !local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps    !local pointer to derived subtype
    type(pft_dgvstate_type)   , pointer :: pdgvs  !local pointer to derived subtype
    real(r8), pointer, dimension(:,:) :: rbuf2dg  !pointer to memory to be allocated
    integer , pointer, dimension(:,:) :: ibuf2dg  !pointer to memory to be allocated
    real(r8), pointer, dimension(:,:) :: rglob2dg !pointer to memory to be allocated
    integer , pointer, dimension(:,:) :: iglob2dg !pointer to memory to be allocated
#if (defined SPMD)
    real(r8), pointer :: rbuf1dc(:)              !temporary column pointer 
#endif
    character(len=256) :: rem_dir                 !remote (archive) directory
    character(len=256) :: rem_fn                  !remote (archive) filename
    integer :: begg,endg,numg                     !gridcell 1d indices
    integer :: begl,endl,numl                     !landunit 1d indices
    integer :: begc,endc,numc                     !column 1d indices
    integer :: begp,endp,nump                     !pft 1d indices
!-----------------------------------------------------------------------

    ! Set 1d indices

    begg = grid1d%beg
    endg = grid1d%end
    numg = grid1d%num

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num

    begl = land1d%beg
    endl = land1d%end
    numl = land1d%num

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num

    !------------------------------------------------------------------
    ! Compute 1d energy and water balance
    !------------------------------------------------------------------

    allocate (water(numc), energy(numc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'lpjreset2: allocation error for water,energy'; call endrun
    end if

    do ci = begc,endc
       c => cpoint(ci)%c
       cps => c%cps
       cws => c%cws
       ces => c%ces
       if (c%l%lps%itype == istsoil) then
          energy(ci) = 0.
          water(ci) = cws%pws_a%h2ocan + cws%h2osno
          do j=1, nlevsoi
             water(ci) = water(ci) + cws%h2osoi_ice(j) + cws%h2osoi_liq(j)
             cv(j) = cps%csol(j)*(1-cps%watsat(j))*cps%dz(j) +   &
                    (cws%h2osoi_ice(j)*cpice + cws%h2osoi_liq(j)*cpliq)
             energy(ci) = energy(ci) + cv(j)*ces%t_soisno(j)
          end do
          if (cps%snl+1 == 1 .AND. cws%h2osno > 0.) then
             cv(1) = cv(1) + cpice*cws%h2osno
             energy(ci) = energy(ci) + cv(1)*ces%t_soisno(1)
          end if
          if (cps%snl+1 < 1)then
             do j = cps%snl+1, 0
                cv(j) = cpliq*cws%h2osoi_liq(j) + cpice*cws%h2osoi_ice(j)
                energy(ci) = energy(ci) + cv(j)*ces%t_soisno(j)
             enddo
          end if
       else
          water(ci) = spval
          energy(ci) = spval
       endif
    end do

#if (defined SPMD)
    allocate(rbuf1dc(begc:endc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'lpjreset2: allocation error for rbuf1dc'; call endrun
    end if
    rbuf1dc(begc:endc) = water(begc:endc)
    call gather_data_to_master(rbuf1dc, water, clmlevel=cols1d%name) 
    rbuf1dc(begc:endc) = energy(begc:endc)
    call gather_data_to_master(rbuf1dc, energy, clmlevel=cols1d%name) 
    deallocate(rbuf1dc)
#endif

    ! compute xy water and energy balance use non-updated weights 

    if (masterproc) then
       water_xy(:,:) = 0.0
       call vec2xy ('column', water, 0._r8, water_xy)
       where (landmask(:,:) /= 1)
          water_xy(:,:) = spval
       endwhere
       call wrap_put_vara_realx (ncid, begwater_id, beg3d, len3d, water_xy)
       energy_xy(:,:) = 0.0
       call vec2xy ('column', energy, 0._r8, energy_xy)
       where (landmask(:,:) /= 1)
          energy_xy(:,:) = spval
       endwhere
       call wrap_put_vara_realx (ncid, begenergy_id, beg3d, len3d, energy_xy)
    endif

    !------------------------------------------------------------------
    ! Reset subgrid weights and areas
    !------------------------------------------------------------------

    ! In CLM2 with satellite data, the number of veg patches is determined once
    ! and is less than maxpatch_pft (4) in some cells.
    ! In LSM with LPJ, the number of veg patches could be dynamic. Until we
    ! implement it as such, we will make all grid cells have 10 veg patches.

    ! Determine new subgrid weights and areas

    call resetWeightsDGVM()

    !------------------------------------------------------------------
    ! Write out more variables to DGVM history output
    ! Note that checking on ifspecial status below guarantees
    ! that the m index will always lie between 1 and maxpatch_pft
    ! (rather than 1 and maxpatch) 
    !------------------------------------------------------------------

    ! Allocate dynamic memory

    allocate(rbuf2dg(maxpatch_pft,begg:endg), ibuf2dg(maxpatch_pft,begg:endg), &
         stat=ier)
    if (ier /= 0) then
       write (6,*) 'lpjreset2: allocation error for rbuf2dg, ibuf2dg)'; call endrun
    end if
#if (defined SPMD)    
    if (masterproc) then
       allocate(rglob2dg(maxpatch_pft,numg), iglob2dg(maxpatch_pft,numg), stat=ier)
       if (ier /= 0) then
          write (6,*) 'lpjreset2: allocation error for rglob2dg, iglob2dg'; call endrun
       end if
    end if
#endif
    if (masterproc) then
       allocate (rbuf3d_xy(lsmlon,lsmlat,maxpatch_pft), ibuf3d_xy(lsmlon,lsmlat,maxpatch_pft), &
            stat=ier)
       if (ier /= 0) then
          write (6,*) 'lpjreset2: allocation error for rbuf3d_xy, ibuf3d_xy'; call endrun
       end if
    endif

    ! Write out pctxy 

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%fpcgrid*100.
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, fpcgrid_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out pftxy 

    ibuf2dg(:,:) = 0
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   ibuf2dg(m,gi) = p%pps%itype
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(ibuf2dg, iglob2dg, clmlevel=grid1d%name) 
#else
    iglob2dg => ibuf2dg
#endif
    if (masterproc) then
       ibuf3d_xy(:,:,:) = 9999
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             ibuf3d_xy(i,j,m) = iglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_int (ncid, itypveg_id, beg4d, len4d, ibuf3d_xy)
    endif

    ! Write out lm_ind

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%lm_ind
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, lmind_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out rm_ind

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%rm_ind
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, rmind_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out sm_ind

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%sm_ind
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, smind_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out hm_ind

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%hm_ind
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, hmind_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Write out nind

    rbuf2dg(:,:) = 0.
    do gi = begg,endg
       g => gpoint(gi)%g
       do li = 1,g%gps%nlandunits
          l => g%l(li)
          lps => l%lps
          if (.not. lps%ifspecial) then
             do ci = 1,l%lps%ncolumns
                c => l%c(ci)
                cps => c%cps
                do pi = 1,cps%npfts 
                   p => c%p(pi)
                   m = p%pps%mxy
                   rbuf2dg(m,gi) = p%pdgvs%nind
                end do
             end do
          end if
       end do
    end do
#if (defined SPMD)
    call gather_data_to_master(rbuf2dg, rglob2dg, clmlevel=grid1d%name) 
#else
    rglob2dg => rbuf2dg
#endif
    if (masterproc) then
       rbuf3d_xy(:,:,:)= spval
       do gi = 1, numg
          i = grid1d%ixy(gi)
          j = grid1d%jxy(gi)
          do m=1,maxpatch_pft
             rbuf3d_xy(i,j,m) = rglob2dg(m,gi)
          end do
       end do
       call wrap_put_vara_realx (ncid, nind_id, beg4d, len4d, rbuf3d_xy)
    endif

    ! Compute xy water and energy balance using updated weights 

    if (masterproc) then
       water_xy(:,:) = 0.0
       call vec2xy ('column', water, 0._r8, water_xy)
       where (landmask(:,:) /= 1)
          water_xy(:,:) = spval
       endwhere
       call wrap_put_vara_realx (ncid, endwater_id, beg3d, len3d, water_xy)

       energy_xy(:,:) = 0.0
       call vec2xy ('column', energy, 0._r8, energy_xy)
       where (landmask(:,:) /= 1)
          energy_xy(:,:) = spval
       endwhere
       call wrap_put_vara_realx (ncid, endenergy_id, beg3d, len3d, energy_xy)
    endif

    ! Deallocate dynamic memory

    deallocate (rbuf2dg, ibuf2dg)
#if (defined SPMD)
    if (masterproc) then
       deallocate (rglob2dg,iglob2dg)
    endif
#endif
    if (masterproc) then
       deallocate (rbuf3d_xy, ibuf3d_xy)
    endif
    deallocate (water, energy)

    !------------------------------------------------------------------
    ! Close and archive netcdf DGVM history file
    !------------------------------------------------------------------

    if (masterproc) then
       call wrap_close(ncid)
       write(6,*)'(LPJRESET2): Finished writing clm2 DGVM history dataset ',&
            trim(dgvm_fn), 'at nstep = ',get_nstep()
       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/hist/'
          rem_fn = set_filename(rem_dir, dgvm_fn)
          call putfil (dgvm_fn, rem_fn, mss_wpass, mss_irt, .true.)
       endif
    endif

  end subroutine lpjreset2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: resetTimeConstDGVM
!
! !INTERFACE:
  subroutine resetTimeConstDGVM()
!
! !DESCRIPTION: 
! Initialize/reset time invariant DGVM variables
!
! !USES:
    use clmtype
    use clmpoint
    use clm_varpar, only : nlevsoi
    use clm_varcon, only : spval
    use pftvarcon, only : roota_par, rootb_par, z0mr, displar, &
                          dleaf, rhol, rhos, taul, taus, xl,   &
                          qe25, vcmx25, mp, c3psn, noveg
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: m,ib,lev,pi                   !indices
    integer :: ivt                           !vegetation type index
    integer :: begp, endp                    !1d pft beginning, ending indices
    type(pft_type)   , pointer :: p	     !local pointer to derived subtype
    type(column_type), pointer :: c          !local pointer to derived subtype
    type(pft_pstate_type)   , pointer :: pps !local pointer to derived subtype
    type(column_pstate_type), pointer :: cps !local pointer to derived subtype
!-----------------------------------------------------------------------

    ! Determine per-process beginning and ending indices

    begp = pfts1d%beg
    endp = pfts1d%end

    ! Loop over mpfts
    do pi = begp,endp
       
       ! Assign local pointer for simpler referencing
       p => ppoint(pi)%p
       pps => p%pps
       c => p%c
       cps => c%cps
       
       ! Vegetation type
       ivt = pps%itype
                
       ! Reset pointer to appropriate array element of ecophysiological constants 
       p%pepc => pftcon(ivt)
       p%pdgvepc => dgv_pftcon(ivt)
                
       ! Initialize root fraction (computing from surface, d is depth in meter):
       ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
       ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with 
       ! beta & d_obs given in Zeng et al. (1998).
       if (ivt /= noveg) then
          do lev = 1, nlevsoi-1
             pps%rootfr(lev) = .5*( exp(-roota_par(ivt)*cps%zi(lev-1))  &
                                  + exp(-rootb_par(ivt)*cps%zi(lev-1))  &
                                  - exp(-roota_par(ivt)*cps%zi(lev  ))  &
                                  - exp(-rootb_par(ivt)*cps%zi(lev  )) )
          end do
          pps%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*cps%zi(nlevsoi-1))  &
                                   + exp(-rootb_par(ivt)*cps%zi(nlevsoi-1)) )
       else
          pps%rootfr(1:nlevsoi) = spval
       endif
              
    end do !end pft level initialization

  end subroutine resetTimeConstDGVM

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: resetWeightstDGVM
!
! !INTERFACE:
  subroutine resetWeightsDGVM()
!
! !DESCRIPTION: 
! Reset DGVM weights 
!
! !USES:
    use clmtype
    use clmpoint
#if (defined SPMD)
    use spmdMod, only : masterproc, gather_data_to_master
#else
    use spmdMod, only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!  subroutine lpjreset2 in this module: as part of the DGVM calculation
!  subroutine restart_dgvm in module restartDGVMMod: if the restart file is read
!  subroutine inicrd in module inicFileMod: if the initial file is read
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: gi,pi,ci,li                       !indices 
    integer :: begg,endg,numg                    !gridcell 1d indices
    integer :: begl,endl,numl                    !landunit 1d indices
    integer :: begc,endc,numc                    !column 1d indices
    integer :: begp,endp,nump                    !pft 1d indices
    integer :: ier                               !return error code	
    real(r8):: sumwt                             !consistency check
    real(r8):: pftsum                            !temporary used for column pft averaging 
    type(pft_type)     , pointer :: p	         !local pointer to derived subtype
    type(column_type)  , pointer :: c            !local pointer to derived subtype
    type(landunit_type), pointer :: l            !local pointer to derived subtype
    type(gridcell_type), pointer :: g            !local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps   !local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps   !local pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps   !local pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps   !local pointer to derived subtype
    type(pft_dgvstate_type)   , pointer :: pdgvs !local pointer to derived subtype
#if (defined SPMD)
    real(r8), pointer :: rbuf1dc(:)              !temporary column pointer 
    real(r8), pointer :: rglob1dc(:)             !temporary column pointer 
    real(r8), pointer :: rbuf1dp(:)              !temporary pft pointer 
    real(r8), pointer :: rglob1dp(:)             !temporary pft pointer 
#endif
!-----------------------------------------------------------------------

    ! Set 1d subgrid beginning and ending indices

    begg = grid1d%beg
    endg = grid1d%end
    numg = grid1d%num

    begl = land1d%beg
    endl = land1d%end
    numl = land1d%num

    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num

    ! Determine new pft properties

    do pi = begp,endp
       p => ppoint(pi)%p  ! pointer to pft
       c => p%c           ! pointer to column containing pft
       l => c%l           ! pointer to landunit containing pft
       g => l%g           ! pointer to gridcell containing pft
       pps => p%pps
       cps => c%cps
       lps => l%lps
       gps => g%gps
       pdgvs => p%pdgvs
       if (.not. lps%ifspecial) then
          ! Determine pft weight relative to column and relative to landunit
#if (defined NOCOMPETE)
          ! One column per pft - column and pft areas are equal 
          pps%wt = 1.
          pfts1d%wtcol(pi) = 1.
          pfts1d%wtlnd(pi) = pdgvs%fpcgrid
#else
          ! One column per landunit - column and landunit areas are equal 
          ! So weight relative to column and weight relative to landunit
          ! are identical
          pps%wt = pdgvs%fpcgrid  
          pfts1d%wtcol(pi) = pps%wt
          pfts1d%wtlnd(pi) = pps%wt
#endif       
          ! Determine new pft area 
          pps%area = pfts1d%wtlnd(pi) * lps%area  

          ! Determine new pft weight relative to grid cell
          pps%wtxy = pps%area / gps%area
          pfts1d%wtxy(pi) = pps%wtxy
       end if
    end do

    ! Determine new column properties

    do ci = begc,endc
       c => cpoint(ci)%c
       l => c%l
       g => l%g
       cps => c%cps
       lps => l%lps
       gps => g%gps
       if (.not. lps%ifspecial) then
          ! Note: if compete is not defined , cps%npts is equal to one
          do pi = 1,cps%npfts 
             p => c%p(pi)
             pps => p%pps
             pdgvs => p%pdgvs
#if (defined NOCOMPETE)
             ! Determine new column weight relative to landunit 
             ! When competition is not on, each column has only one pft
             cps%wt = pdgvs%fpcgrid  
             cols1d%wtlnd(ci) = cps%wt
             
             ! Determine new column area 
             cps%area = cols1d%wtlnd(ci) * lps%area 

             ! Determine new column weight relative to gridcell
             cps%wtxy = cps%area / gps%area
             cols1d%wtxy(ci) = cps%area / gps%area
#endif
             ! Determine new column pft properties
             c%pa(pi) = pps%area
             c%pw(pi) = pps%wt
             c%pt(pi) = pps%itype
          end do
       end if
    end do

#if (defined NOCOMPETE)
    ! Determine new landunit properties - column areas and
    ! weights relative to landunit 

    do li = begl,endl
       l => lpoint(li)%l
       lps => l%lps
       if (.not. lps%ifspecial) then
          do ci = 1,lps%ncolumns
             c => l%c(ci)
             cps => c%cps
             l%ca(ci) = cps%area
             l%cw(ci) = cps%wt   
          end do
       end if
    end do
#endif

    ! Consistency check - add up all the pft weights for a given gridcell
    ! and make sure they are not greater than 1.

    do gi = begg,endg
       g => gpoint(gi)%g
       gps => g%gps
       sumwt = 0._r8
       do li = 1,gps%nlandunits
          l => g%l(li)
          lps => l%lps
          do ci = 1,lps%ncolumns
             c => l%c(ci)
             cps => c%cps
             do pi = 1,cps%npfts
                p => c%p(pi)
                pps => p%pps
                sumwt = sumwt + pps%wtxy
             end do
          end do
       end do
       if (abs(sumwt - 1.0) > 1.0e-6) then
          write(6,*) 'Error in resetWeightsDGVM: sumwt of pfts for grid cell ',&
               'i,j = ',grid1d%ixy(gi),grid1d%jxy(gi),' not equal to 1'
          write(6,*) 'sum of pft weights for gridcell =', sumwt
          call endrun
       end if
    end do

    ! Determine average over all column pfts for h2ocan using new weights
    ! This is needed by begwb in DriverInitMod.F90.

    do ci = begc,endc
       c => cpoint(ci)%c
       pftsum = 0.0
       do pi = 1,c%cps%npfts
          pftsum = pftsum + c%p(pi)%pws%h2ocan * c%pw(pi)
       end do
       c%cws%pws_a%h2ocan = pftsum
    end do

#if (defined SPMD) && (defined NOCOMPETE)
    ! On master processor need to determine column weights for all processor
    ! columns. This is needed for calls to map from 1d vector to xy grid
    ! and for 1d history output. 

    if (masterproc) then
       allocate (rglob1dc(numc), stat=ier)
       if (ier /= 0) then
          write (6,*) 'lpjreset2: allocation error for rglob1dc'; call endrun
       end if
    endif
    allocate (rbuf1dc(begc:endc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'resetWeightsDGVM: allocation error for rbuf1dc'; call endrun
    end if

    rbuf1dc(begc:endc) = cols1d%wtxy(begc:endc)
    call gather_data_to_master (rbuf1dc, rglob1dc, clmlevel=cols1d%name) 
    if (masterproc) cols1d%wtxy(1:numc) = rglob1dc(1:numc)

    rbuf1dc(begc:endc) = cols1d%wtlnd(begc:endc)
    call gather_data_to_master (rbuf1dc, rglob1dc, clmlevel=cols1d%name) 
    if (masterproc) cols1d%wtlnd(1:numc) = rglob1dc(1:numc)

    if (masterproc) deallocate(rglob1dc)
    deallocate (rbuf1dc)
#endif

#if (defined SPMD)
    ! On master processor need to determine pft weights for all processor
    ! pfts. This is needed for calls to map from 1d vector to xy grid
    ! and for 1d history output. 

    if (masterproc) then
       allocate (rglob1dp(nump), stat=ier)
       if (ier /= 0) then
          write (6,*) 'resetWeightsDGVM: allocation error for rglob1dp'; call endrun
       end if
    endif
    allocate (rbuf1dp(begp:endp), stat=ier)
    if (ier /= 0) then
       write (6,*) 'resetWeigtsDGVM: allocation error for rbuf1dp'; call endrun
    end if

    rbuf1dp(begp:endp) = pfts1d%wtxy(begp:endp)
    call gather_data_to_master (rbuf1dp, rglob1dp, clmlevel=pfts1d%name) 
    if (masterproc) pfts1d%wtxy(1:nump) = rglob1dp(1:nump)

    rbuf1dp(begp:endp) = pfts1d%wtcol(begp:endp)
    call gather_data_to_master (rbuf1dp, rglob1dp, clmlevel=pfts1d%name) 
    if (masterproc) pfts1d%wtcol(1:nump) = rglob1dp(1:nump)

    rbuf1dp(begp:endp) = pfts1d%wtlnd(begp:endp)
    call gather_data_to_master (rbuf1dp, rglob1dp, clmlevel=pfts1d%name) 
    if (masterproc) pfts1d%wtlnd(1:nump) = rglob1dp(1:nump)

    deallocate (rbuf1dp)
    if (masterproc) deallocate(rglob1dp)
#endif

  end subroutine resetWeightsDGVM

#endif

end module DGVMMod




