#include <misc.h>
module dp_coupling
!BOP
!
! !MODULE: dp_coupling --- dynamics-physics coupling module
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use rgrid,         only: nlon
   use pmgrid,        only: plon, plat, plev, twod_decomp, iam,    &
                            beglat, endlat, beglev, endlev,        &
                            beglatxy, endlatxy, beglonxy, endlonxy
   use ppgrid,        only: pcols, pver
   use phys_grid
   use phys_buffer,   only: pbuf_fld, pbuf_size_max
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: ppcnst, qmin
   use physconst,     only: cpair, gravit, rair, zvir
   use geopotential,  only: geopotential_t
   use check_energy,  only: check_energy_timestep_init
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC d_p_coupling, p_d_coupling

!
! !DESCRIPTION:
!
!      This module provides 
!
!      \begin{tabular}{|l|l|} \hline \hline
!        d\_p\_coupling    &  dynamics output to physics input \\ \hline
!        p\_d\_coupling    &  physics output to dynamics input \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.10.01   Lin        Various revisions
!   01.03.26   Sawyer     Added ProTeX documentation
!   01.06.27   Mirin      Separate noncoupling coding into new routines
!   01.07.13   Mirin      Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   03.03.28   Boville    set all physics_state elements, add check_energy_timestep_init
!   03.08.13   Sawyer     Removed ghost N1 region in u3sxy
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: d_p_coupling --- convert dynamics output to physics input
!
! !INTERFACE: 
  subroutine d_p_coupling(ps,   u3s,   v3s,  pt,    coslon,  sinlon,     &
                          q3,   omga, phis,  pe,    peln,    pk,         &
                          pkz,  phys_state,  phys_tend, pbuf, full_phys, &
                          psxy,  u3sxy,       v3sxy,   ptxy,             &
                          q3xy, omgaxy,phisxy,      pexy,  pelnxy,       &
                          pkxy, pkzxy )

! !USES:
    use dynamics_vars, only: ng_d, ng_s
#if defined (SPMD)
    use mpishorthand, only : mpicom
    use spmd_dyn, only : ijk_xy_to_yz
    use parutilitiesmodule, only : sumop, parcollective
    use mod_irreg, only: mp_sendirr, mp_recvirr
#endif
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
!
    real(r8), intent(in) :: ps (plon, beglat:endlat)                      ! surface pressure
    real(r8), intent(inout) :: u3s(plon, beglat-ng_d:endlat+ng_s, beglev:endlev)    ! u-wind on d-grid
    real(r8), intent(in) :: v3s(plon, beglat-ng_s:endlat+ng_d, beglev:endlev)       ! v-wind on d-grid
    real(r8), intent(in) :: pt (plon, beglat-ng_d:endlat+ng_d, beglev:endlev)       ! Virtual pot temp
    real(r8), intent(in) :: q3 (plon, beglat-ng_d:endlat+ng_d, beglev:endlev, ppcnst) ! constituents
    real(r8), intent(in) :: omga(plon, beglev:endlev, beglat:endlat)      ! vertical velocity
    real(r8), intent(in) :: phis(plon, beglat:endlat)                     ! surface geopotential
    real(r8), intent(in) :: pe  (plon, beglev:endlev+1, beglat:endlat)    ! this fv's pint
    real(r8), intent(in) :: peln(plon, beglev:endlev+1, beglat:endlat)    ! log(pe)
    real(r8), intent(in) :: pk  (plon, beglat:endlat, beglev:endlev+1)    ! pe**cappa
    real(r8), intent(in) :: pkz (plon, beglat:endlat, beglev:endlev)      ! f-v mean of pk
    real(r8), intent(in) :: coslon(plon)                                  ! cosine of longitude
    real(r8), intent(in) :: sinlon(plon)                                  ! sin of longitudes
    logical,  intent(in) :: full_phys

! xy-decomposed instanciations below:
    real(r8), intent(in) :: psxy (beglonxy:endlonxy, beglatxy:endlatxy)                      ! surface pressure
    real(r8), intent(in) :: u3sxy(beglonxy:endlonxy, beglatxy:endlatxy, plev)       ! u-wind on d-grid
    real(r8), intent(in) :: v3sxy(beglonxy:endlonxy, beglatxy:endlatxy, plev)       ! v-wind on d-grid
    real(r8), intent(in) :: ptxy (beglonxy:endlonxy, beglatxy:endlatxy, plev)       ! Virtual pot temp
    real(r8), target, intent(in) :: q3xy (beglonxy:endlonxy, beglatxy:endlatxy, plev, ppcnst) ! constituents
    real(r8), intent(in) :: omgaxy(beglonxy:endlonxy, plev, beglatxy:endlatxy)      ! vertical velocity
    real(r8), intent(in) :: phisxy(beglonxy:endlonxy, beglatxy:endlatxy)            ! surface geopotential
    real(r8), intent(in) :: pexy  (beglonxy:endlonxy, plev+1, beglatxy:endlatxy)    ! this fv's pint
    real(r8), intent(in) :: pelnxy(beglonxy:endlonxy, plev+1, beglatxy:endlatxy)    ! log(pe)
    real(r8), intent(in) :: pkxy  (beglonxy:endlonxy, beglatxy:endlatxy, plev+1)    ! pe**cappa
    real(r8), intent(in) :: pkzxy (beglonxy:endlonxy, beglatxy:endlatxy, plev)      ! f-v mean of pk

! !OUTPUT PARAMETERS:

    type(physics_state), intent(out), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(out), dimension(begchunk:endchunk) :: phys_tend
    type(pbuf_fld),    intent(inout), dimension(pbuf_size_max)     :: pbuf

! !DESCRIPTION:
!
!   Coupler for converting dynamics output variables into physics 
!   input variables
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.07.13   AAM        Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   02.05.02   Sawyer     u3s made inout due to ghosting in d2a3dikj
!   03.08.05   Sawyer     Removed pe11k, pe11kln (for defunct Rayl fric)
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
    integer :: i,ib,j,k,m,lchnk      ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: blksiz                ! number of columns in 2D block
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8) :: qmavl                ! available q at level pver-1
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: pic(pcols)           ! ps**cappa
    real(r8), allocatable :: u3(:, :, :)       ! u-wind on a-grid
    real(r8), allocatable :: v3(:, :, :)       ! v-wind on a-grid
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer
                                     ! transpose buffers
!---------------------------End Local workspace-------------------------

    if (twod_decomp .eq. 1) then

!-----------------------------------------------------------------------
! Transform dynamics staggered winds to physics grid (D=>A)
!-----------------------------------------------------------------------

       allocate (u3(beglonxy:endlonxy, plev, beglatxy:endlatxy))
       allocate (v3(beglonxy:endlonxy, plev, beglatxy:endlatxy))

       call d2a3dikj(u3sxy,  v3sxy,   u3,    v3,   plon,  plat, plev,    &
                     beglatxy, endlatxy,  0, 0, 0, 0, 0,                 &
                     beglonxy, endlonxy,    coslon,      sinlon)

!-----------------------------------------------------------------------
! Copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------

       if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, PIC)

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)
             phys_state(lchnk)%ncol  = ncol
             phys_state(lchnk)%lchnk = lchnk

             do i=1,ncol
                phys_state(lchnk)%ps(i)   = psxy(lons(i),lats(i))
                phys_state(lchnk)%phis(i) = phisxy(lons(i),lats(i))
                pic(i) = pkxy(lons(i),lats(i),pver+1)
             enddo

             do k=1,plev
                do i=1,ncol
                   phys_state(lchnk)%u    (i,k) = u3(lons(i),k,lats(i))
                   phys_state(lchnk)%v    (i,k) = v3(lons(i),k,lats(i))
                   phys_state(lchnk)%omega(i,k) = omgaxy(lons(i),k,lats(i))
                   if (full_phys) then
                      phys_state(lchnk)%t    (i,k) = ptxy(lons(i),lats(i),k) / (1. + zvir*q3xy(lons(i),lats(i),k,1))
                      phys_state(lchnk)%exner(i,k) = pic(i) / pkzxy(lons(i),lats(i),k) 
                   else
                      phys_state(lchnk)%t    (i,k) = ptxy(lons(i),lats(i),k) * pkzxy(lons(i),lats(i),k)
                   end if
                end do
             end do

             do k=1,plev+1
                do i=1,ncol
!
! edge-level pressure arrays: copy from the arrays computed by dynpkg
!
                   phys_state(lchnk)%pint  (i,k) = pexy  (lons(i),k,lats(i))
                   phys_state(lchnk)%lnpint(i,k) = pelnxy(lons(i),k,lats(i))
                end do
             end do

!
! Copy constituents
!
             do m=1,ppcnst
                do k=1,plev
                   do i=1,ncol
                      phys_state(lchnk)%q(i,k,m) = q3xy(lons(i),lats(i),k,m)
                   end do
                end do
             end do
 
          end do   ! begchunk:endchunk loop

       else

          tsize = 7 + ppcnst

          blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
          allocate(bpter(blksiz,0:plev))
          allocate(bbuffer(tsize*block_buf_nrecs))
          allocate(cbuffer(tsize*chunk_buf_nrecs))

          call block_to_chunk_send_pters(iam+1,blksiz,plev+1,tsize,bpter)

          ib = 0
          do j=beglatxy,endlatxy
             do i=beglonxy,endlonxy
                ib = ib + 1

                bbuffer(bpter(ib,0))   = pexy(i,plev+1,j)
                bbuffer(bpter(ib,0)+1) = pelnxy(i,plev+1,j)
                bbuffer(bpter(ib,0)+2) = psxy(i,j)
                bbuffer(bpter(ib,0)+3) = phisxy(i,j)

                do k=1,plev

                   bbuffer(bpter(ib,k))   = pexy(i,k,j)
                   bbuffer(bpter(ib,k)+1) = pelnxy(i,k,j)
                   bbuffer(bpter(ib,k)+2) = u3    (i,k,j)
                   bbuffer(bpter(ib,k)+3) = v3    (i,k,j)
                   bbuffer(bpter(ib,k)+4) = omgaxy(i,k,j)
                   if (full_phys) then
                      bbuffer(bpter(ib,k)+5) = ptxy(i,j,k) / (1. + zvir*q3xy(i,j,k,1))
                      bbuffer(bpter(ib,k)+6) = pkxy(i,j,pver+1) / pkzxy(i,j,k) 
                   else
                      bbuffer(bpter(ib,k)+6) = ptxy(i,j,k) * pkzxy(i,j,k)
                   end if

                   do m=1,ppcnst
                      bbuffer(bpter(ib,k)+6+m) = q3xy(i,j,k,m)
                   end do

                end do
             end do
          end do

          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             phys_state(lchnk)%ncol  = ncol
             phys_state(lchnk)%lchnk = lchnk

             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

             do i=1,ncol

                phys_state(lchnk)%pint  (i,pver+1) = cbuffer(cpter(i,0))
                phys_state(lchnk)%lnpint(i,pver+1) = cbuffer(cpter(i,0)+1)
                phys_state(lchnk)%ps(i)            = cbuffer(cpter(i,0)+2)
                phys_state(lchnk)%phis(i)          = cbuffer(cpter(i,0)+3)

                do k=1,plev

                   phys_state(lchnk)%pint  (i,k) = cbuffer(cpter(i,k))
                   phys_state(lchnk)%lnpint(i,k) = cbuffer(cpter(i,k)+1)
                   phys_state(lchnk)%u     (i,k) = cbuffer(cpter(i,k)+2)
                   phys_state(lchnk)%v     (i,k) = cbuffer(cpter(i,k)+3)
                   phys_state(lchnk)%omega (i,k) = cbuffer(cpter(i,k)+4)
                   if (full_phys) then
                      phys_state(lchnk)%t    (i,k) = cbuffer(cpter(i,k)+5)
                      phys_state(lchnk)%exner(i,k) = cbuffer(cpter(i,k)+6)
                   else
                      phys_state(lchnk)%t    (i,k) = cbuffer(cpter(i,k)+6)
                   end if

                   do m=1,ppcnst
                      phys_state(lchnk)%q(i,k,m) = cbuffer(cpter(i,k)+6+m)
                   end do

                end do
             end do

          end do   ! begchunk:endchunk loop

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif

    else

!-----------------------------------------------------------------------
! Transform dynamics staggered winds to physics grid (D=>A)
!-----------------------------------------------------------------------

       allocate (u3(plon, plev, beglat:endlat))
       allocate (v3(plon, plev, beglat:endlat))

       call d2a3dikj(u3s,    v3s,     u3,    v3,   plon,  plat, plev,    &
                     beglat, endlat,  ng_d,  ng_d, ng_s,  ng_s, ng_d,    &
                     1, plon, coslon, sinlon)

!-----------------------------------------------------------------------
! Copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------

       if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, PIC)

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)
             phys_state(lchnk)%ncol = ncol
             phys_state(lchnk)%lchnk = lchnk

             do i=1,ncol
                phys_state(lchnk)%ps(i)   = ps(lons(i),lats(i))
                phys_state(lchnk)%phis(i) = phis(lons(i),lats(i))
                pic(i) = pk(lons(i),lats(i),pver+1)
             enddo

             do k=1,plev
                do i=1,ncol
                   phys_state(lchnk)%u    (i,k) = u3(lons(i),k,lats(i))
                   phys_state(lchnk)%v    (i,k) = v3(lons(i),k,lats(i))
                   phys_state(lchnk)%omega(i,k) = omga(lons(i),k,lats(i))
                   if (full_phys) then
                      phys_state(lchnk)%t    (i,k) = pt(lons(i),lats(i),k) / (1. + zvir*q3(lons(i),lats(i),k,1))
                      phys_state(lchnk)%exner(i,k) = pic(i) / pkz(lons(i),lats(i),k) 
                   else
                      phys_state(lchnk)%t    (i,k) = pt(lons(i),lats(i),k) * pkz(lons(i),lats(i),k)
                   end if
                end do
             end do

             do k=1,plev+1
                do i=1,ncol
!
! edge-level pressure arrays: copy from the arrays computed by dynpkg
!
                   phys_state(lchnk)%pint  (i,k) = pe  (lons(i),k,lats(i))
                   phys_state(lchnk)%lnpint(i,k) = peln(lons(i),k,lats(i))
                end do
             end do

!
! Copy constituents
!
             do m=1,ppcnst
                do k=1,plev
                   do i=1,ncol
                      phys_state(lchnk)%q(i,k,m) = q3(lons(i),lats(i),k,m)
                   end do
                end do
             end do
 
          end do   ! begchunk:endchunk loop

       else

          tsize = 7 + ppcnst

          allocate(bpter(plon,0:plev))
          allocate(bbuffer(tsize*block_buf_nrecs))
          allocate(cbuffer(tsize*chunk_buf_nrecs))

          do j=beglat,endlat

             call block_to_chunk_send_pters(j,plon,plev+1,tsize,bpter)

             do i=1,nlon(j)

                bbuffer(bpter(i,0))   = pe(i,plev+1,j)
                bbuffer(bpter(i,0)+1) = peln(i,plev+1,j)
                bbuffer(bpter(i,0)+2) = ps(i,j)
                bbuffer(bpter(i,0)+3) = phis(i,j)

                do k=1,plev

                   bbuffer(bpter(i,k))   = pe(i,k,j)
                   bbuffer(bpter(i,k)+1) = peln(i,k,j)
                   bbuffer(bpter(i,k)+2) = u3    (i,k,j)
                   bbuffer(bpter(i,k)+3) = v3    (i,k,j)
                   bbuffer(bpter(i,k)+4) = omga(i,k,j)
                   if (full_phys) then
                      bbuffer(bpter(i,k)+5) = pt(i,j,k) / (1. + zvir*q3(i,j,k,1))
                      bbuffer(bpter(i,k)+6) = pk(i,j,pver+1) / pkz(i,j,k) 
                   else
                      bbuffer(bpter(i,k)+6) = pt(i,j,k) * pkz(i,j,k)
                   end if

                   do m=1,ppcnst
                      bbuffer(bpter(i,k)+6+m) = q3(i,j,k,m)
                   end do

                end do

             end do

          end do   ! beglat:endlat loop

          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             phys_state(lchnk)%ncol  = ncol
             phys_state(lchnk)%lchnk = lchnk

             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

             do i=1,ncol

                phys_state(lchnk)%pint  (i,pver+1) = cbuffer(cpter(i,0))
                phys_state(lchnk)%lnpint(i,pver+1) = cbuffer(cpter(i,0)+1)
                phys_state(lchnk)%ps(i)            = cbuffer(cpter(i,0)+2)
                phys_state(lchnk)%phis(i)          = cbuffer(cpter(i,0)+3)

                do k=1,plev

                   phys_state(lchnk)%pint  (i,k) = cbuffer(cpter(i,k))
                   phys_state(lchnk)%lnpint(i,k) = cbuffer(cpter(i,k)+1)
                   phys_state(lchnk)%u     (i,k) = cbuffer(cpter(i,k)+2)
                   phys_state(lchnk)%v     (i,k) = cbuffer(cpter(i,k)+3)
                   phys_state(lchnk)%omega (i,k) = cbuffer(cpter(i,k)+4)
                   if (full_phys) then
                      phys_state(lchnk)%t    (i,k) = cbuffer(cpter(i,k)+5)
                      phys_state(lchnk)%exner(i,k) = cbuffer(cpter(i,k)+6)
                   else
                      phys_state(lchnk)%t    (i,k) = cbuffer(cpter(i,k)+6)
                   end if

                   do m=1,ppcnst
                      phys_state(lchnk)%q(i,k,m) = cbuffer(cpter(i,k)+6+m)
                   end do

                end do
             end do

          end do   ! begchunk:endchunk loop

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif

    endif

!
! Evaluate derived quantities
!
    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do k=1,plev
          do i=1,ncol
             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1./phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%pmid (i,k) = 0.5*(phys_state(lchnk)%pint(i,k) + phys_state(lchnk)%pint(i,k+1))
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
          end do
       end do

! Attempt to remove negative constituents in bottom layer only by moving from next level
! This is a BAB kludge to avoid masses of warning messages for cloud water and ice, since
! the vertical remapping operator currently being used for cam is not strictly monotonic 
! at the endpoints.
       do m=1,ppcnst
          do i=1,ncol
             if (phys_state(lchnk)%q(i,pver,m) < qmin(m)) then
! available q in 2nd level
                qmavl = phys_state(lchnk)%q (i,pver-1,m) - qmin(m)
! required q change in bottom level rescaled to mass fraction in 2nd level
                dqreq = (qmin(m) - phys_state(lchnk)%q(i,pver,m))                         &
                      * phys_state(lchnk)%pdel(i,pver) / phys_state(lchnk)%pdel(i,pver-1)
                qbot   = phys_state(lchnk)%q(i,pver  ,m)
                qbotm1 = phys_state(lchnk)%q(i,pver-1,m)
                if (dqreq < qmavl) then
                   phys_state(lchnk)%q(i,pver  ,m) = qmin(m)
                   phys_state(lchnk)%q(i,pver-1,m) = phys_state(lchnk)%q(i,pver-1,m) - dqreq
                   if (dqreq>1.e-14) print *, 'dpcoup dqreq', m, lchnk, i, qbot, qbotm1, dqreq
                else 
                   print *, 'dpcoup cant adjust', m, lchnk, i, qbot, qbotm1, dqreq
                end if
             end if
          end do
       end do
                   
!
! Compute initial geopotential heights
       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                            phys_state(lchnk)%t     , phys_state(lchnk)%q(1,1,1), rair,  gravit,  zvir    , &
                            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                                      + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
          end do
       end do

! Compute energy and water integrals of input state
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf)

    end do

    deallocate (u3)
    deallocate (v3)

!EOC
  end subroutine d_p_coupling
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: p_d_coupling --- convert physics output to dynamics input
!
! !INTERFACE: 
  subroutine p_d_coupling(phys_state, phys_tend, full_phys,        &
                          adiabatic,  t3, t3old,   q3old,          &
                          q3, pt,   dudt,      dvdt,    pkz,     &
                          q3xy, ptxy, dudtxy,    dvdtxy,  pkzxy,   &
                          dtime, u3s, v3s,       u3sxy,   v3sxy,   &
                          zvir,       cappa,     ptop,     pk,     &
                          peln,       ps,                          &
                          pe,         pexy,      delp,     delpxy )

! !USES:
    use dynamics_vars, only : ng_d, ng_s, yzt, q3t, yzt2
    use history, only: outfld
    use time_manager, only: is_first_step, is_first_restart_step
    use constituents, only: tendnam
#if defined ( SPMD )
    use spmd_dyn, only : uxy_to_u, vxy_to_v, ikj_xy_to_yz, qxy_to_q,     &
                         m_ttrans, q_ttrans, r_ttrans,    &
                         ptxy_to_pt, ijk_xy_to_yz, qxy3_to_q3, rxy_to_r
    use mod_irreg, only: mp_sendirr, mp_recvirr
#endif
!-----------------------------------------------------------------------
    implicit none

! Variables ending in xy are xy-decomposition instanciations.

! !INPUT PARAMETERS:
    type(physics_state), intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(in), dimension(begchunk:endchunk) :: phys_tend
    logical,  intent(in) :: full_phys
    logical,  intent(in) :: adiabatic
    real(r8), intent(in) :: dtime
    real(r8), intent(in) :: zvir
    real(r8), intent(in) :: cappa
    real(r8), intent(in) :: ptop
! xy-decomposition instanciation immediately below
    real(r8), intent(in) :: pkzxy(beglonxy:endlonxy,beglatxy:endlatxy,plev)
    real(r8), intent(in) :: u3sxy(beglonxy:endlonxy,beglatxy:endlatxy,plev)
    real(r8), intent(in) :: v3sxy(beglonxy:endlonxy,beglatxy:endlatxy,plev)

! !INPUT/OUTPUT PARAMETERS:
    real(r8), intent(inout) :: u3s(plon,beglat-ng_d:endlat+ng_s,beglev:endlev)
    real(r8), intent(inout) :: v3s(plon,beglat-ng_s:endlat+ng_d,beglev:endlev)
    real(r8), intent(inout) :: pkz(plon,beglat:endlat,beglev:endlev)
    real(r8), intent(inout) :: pe(plon,beglev:endlev+1,beglat:endlat)
! xy-decomposition instantiation immediately below
    real(r8), intent(inout) :: pexy(beglonxy:endlonxy,plev+1,beglatxy:endlatxy) ! work variable

! !OUTPUT PARAMETERS:
    real(r8), intent(out) :: pt(plon,beglat-ng_d:endlat+ng_d,beglev:endlev)
    real(r8), intent(out) :: t3(plon,beglev:endlev,beglat:endlat) !  Temperature
    real(r8), intent(inout) :: t3old(plon,beglev:endlev,beglat:endlat) !  Temperature
    real(r8), intent(out) :: q3(plon,beglat-ng_d:endlat+ng_d,beglev:endlev,ppcnst) ! constituents
    real(r8), intent(inout) :: q3old(plon,beglat-ng_d:endlat+ng_d,beglev:endlev,ppcnst) ! constituents
! xy-decomposition instantiation immediately below
    real(r8), intent(out) :: ptxy(beglonxy:endlonxy,beglatxy:endlatxy,plev)
    real(r8), intent(out) :: q3xy(beglonxy:endlonxy,beglatxy:endlatxy,plev,ppcnst)

    real(r8), intent(out) :: dudt(plon,beglev:endlev,beglat:endlat) ! U-velocity tendency
    real(r8), intent(out) :: dvdt(plon,beglev:endlev,beglat:endlat) ! V-velocity tendency
! xy-decomposition instantiation immediately below
    real(r8), intent(out) :: dudtxy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8), intent(out) :: dvdtxy(beglonxy:endlonxy,plev,beglatxy:endlatxy)

    real(r8), intent(out) :: pk(plon,beglat:endlat,beglev:endlev+1)
    real(r8), intent(out) :: peln(plon,beglev:endlev+1,beglat:endlat)
    real(r8), intent(out) :: delp(plon,beglat:endlat,beglev:endlev)
! xy-decomposition instanciation immediately below
    real(r8), intent(out) :: delpxy(beglonxy:endlonxy,beglatxy:endlatxy,plev) ! work variable

    real(r8), intent(out) :: ps(plon,beglat:endlat,beglev)

! !DESCRIPTION:
!
!   Coupler for converting physics output variables into dynamics input variables
!
! !REVISION HISTORY:
!   00.06.01   Boville    Creation
!   01.06.08   AAM        Compactified
!   01.07.13   AAM        Some support for multi-2D decompositions
!   02.03.01   Worley     Support for nontrivial physics remapping
!   02.08.06   Sawyer     T3 added -- updated to current temperature
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
    integer :: i, ib, k, m, j, lchnk  ! indices
    integer :: ncol                   ! number of columns in current chunk
    integer :: lats(pcols)            ! array of latitude indices
    integer :: lons(pcols)            ! array of longitude indices
    integer :: blksiz                 ! number of columns in 2D block
    integer :: tsize                  ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for unpacking data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
    integer :: iqa, iqb, iqc, iqd, mq     ! used for tracer transpose grouping

    real(r8) :: dt5,dtinv
    real(r8), allocatable :: ttend(:,:) ! temperature tendency
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer               ! transpose buffers
!---------------------------End Local workspace-------------------------

! -------------------------------------------------------------------------
! Copy temperature, tendencies and constituents to dynamics data structures
! For adiabatic case, compute transposes only (2-D decomposition)
! -------------------------------------------------------------------------

    if (twod_decomp .eq. 1) then

! -------------------------------------------------------------------------
! Copy onto xy decomposition, then transpose to yz decomposition
! -------------------------------------------------------------------------

       if (.not. adiabatic) then

          if (local_dp_map) then

!$omp parallel do private(lchnk, i, k, ncol, m, lons, lats)

             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                call get_lon_all_p(lchnk, ncol, lons)
                call get_lat_all_p(lchnk, ncol, lats)

                do k = 1, plev
                   do i = 1, ncol
                      dvdtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dvdt(i,k)
                      dudtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dudt(i,k)
                      ptxy  (lons(i),lats(i),k) = phys_state(lchnk)%t(i,k)
                      delpxy(lons(i),lats(i),k) = phys_state(lchnk)%pdel(i,k)
                   enddo
                enddo

                do m=1,ppcnst
                   do k=1,plev
                      do i=1,ncol
                         q3xy(lons(i),lats(i),k,m) = phys_state(lchnk)%q(i,k,m)
                      end do
                   end do
                end do

             enddo

          else

             tsize = 4 + ppcnst

             blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
             allocate(bpter(blksiz,0:plev))
             allocate(bbuffer(tsize*block_buf_nrecs))
             allocate(cbuffer(tsize*chunk_buf_nrecs))

             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)

                call chunk_to_block_send_pters(lchnk,pcols,plev+1,tsize,cpter)

                do i=1,ncol

                   do k=1,plev

                      cbuffer(cpter(i,k))   = phys_tend(lchnk)%dvdt(i,k)
                      cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt(i,k)
                      cbuffer(cpter(i,k)+2) = phys_state(lchnk)%t(i,k)
                      cbuffer(cpter(i,k)+3) = phys_state(lchnk)%pdel(i,k)

                      do m=1,ppcnst
                         cbuffer(cpter(i,k)+3+m) = phys_state(lchnk)%q(i,k,m)
                      end do

                   end do
  
                end do

             end do

             call transpose_chunk_to_block(tsize, cbuffer, bbuffer)

             call chunk_to_block_recv_pters(iam+1,blksiz,plev+1,tsize,bpter)

             ib = 0
             do j=beglatxy,endlatxy
                do i=beglonxy,endlonxy
                   ib = ib + 1

                   do k=1,plev
                      dvdtxy(i,k,j) = bbuffer(bpter(ib,k))
                      dudtxy(i,k,j) = bbuffer(bpter(ib,k)+1)
                      ptxy  (i,j,k) = bbuffer(bpter(ib,k)+2)
                      delpxy(i,j,k) = bbuffer(bpter(ib,k)+3)

                      do m=1,ppcnst
                         q3xy(i,j,k,m) = bbuffer(bpter(ib,k)+3+m)
                      end do

                   enddo
                enddo
             enddo

             deallocate(bpter)
             deallocate(bbuffer)
             deallocate(cbuffer)

          endif
       endif

#if defined (SPMD)
! Transpose from xy to yz decomposition
       call t_startf('transpose_total')
       call t_startf('transpose_bck2')

       if (.not. adiabatic) then
! Transpose dudt and dvdt
          call mp_sendirr( dudtxy, ikj_xy_to_yz%SendDesc,        &
                           ikj_xy_to_yz%RecvDesc, dudt )
          call mp_recvirr( dudt, ikj_xy_to_yz%RecvDesc )
          call mp_sendirr( dvdtxy, ikj_xy_to_yz%SendDesc,        &
                           ikj_xy_to_yz%RecvDesc, dvdt )
          call mp_recvirr( dvdt, ikj_xy_to_yz%RecvDesc )
          if (.not. full_phys) then
! Transpose pkz
! WS 02.08.06 : pkz needed to determine potential temperature,
!               pt contains temperature needed for t3
             call mp_sendirr( pkzxy, ijk_xy_to_yz%SendDesc,      &
                              ijk_xy_to_yz%RecvDesc, pkz )
             call mp_recvirr( pkz, ijk_xy_to_yz%RecvDesc )
          endif
          call mp_sendirr( delpxy, ijk_xy_to_yz%SendDesc,      &
                           ijk_xy_to_yz%RecvDesc, delp )
          call mp_recvirr( delp, ijk_xy_to_yz%RecvDesc )
       endif

! Transpose pt
       call mp_sendirr( ptxy, ijk_xy_to_yz%SendDesc,             &
                        ijk_xy_to_yz%RecvDesc, yzt )

       call mp_recvirr( yzt, ijk_xy_to_yz%RecvDesc )
! Transpose u3s
       call mp_sendirr( u3sxy, ijk_xy_to_yz%SendDesc,              &
                        ijk_xy_to_yz%RecvDesc, yzt2 )

!$omp parallel do private(i,j,k)
       do k=beglev,endlev
          do j = beglat,endlat
             do i=1,plon
                pt(i,j,k) = yzt(i,j,k)
             enddo
          enddo
       enddo

       call mp_recvirr( yzt2, ijk_xy_to_yz%RecvDesc )
! Transpose v3s
       call mp_sendirr( v3sxy, ijk_xy_to_yz%SendDesc,            &
                        ijk_xy_to_yz%RecvDesc, yzt )

!$omp parallel do private(i,j,k)
       do k=beglev,endlev
          do j = beglat,endlat
             do i=1,plon
                u3s(i,j,k) = yzt2(i,j,k)
             enddo
          enddo
       enddo

       call mp_recvirr( yzt, ijk_xy_to_yz%RecvDesc )
!$omp parallel do private(i,j,k)
       do k=beglev,endlev
          do j = beglat,endlat
             do i=1,plon
                v3s(i,j,k) = yzt(i,j,k)
             enddo
          enddo
       enddo

       call t_stopf('transpose_bck2')

       call t_startf('transpose_qbck')

! Transpose q3
       iqa = 1
       do mq = 1, q_ttrans
          iqb = iqa + m_ttrans - 1
          call mp_sendirr( q3xy(:,:,:,iqa:iqb), qxy_to_q%SendDesc,    &
                           qxy_to_q%RecvDesc, q3t(:,:,:,iqa:iqb) )
          if (mq .gt. 1) then
            iqc = iqa - m_ttrans
            iqd = iqc + m_ttrans - 1
!$omp parallel do private(i,j,k,m)
            do m=iqc,iqd
               do k=beglev,endlev
                  do j = beglat,endlat
                     do i=1,plon
                        q3(i,j,k,m) = q3t(i,j,k,m)
                     enddo
                  enddo
               enddo
            enddo
          endif
          call mp_recvirr( q3t(:,:,:,iqa:iqb), qxy_to_q%RecvDesc )
          iqa = iqa + m_ttrans
       enddo
       if (r_ttrans .ne. 0) then
          iqb = iqa + r_ttrans - 1
          call mp_sendirr( q3xy(:,:,:,iqa:iqb), rxy_to_r%SendDesc,     &
                           rxy_to_r%RecvDesc, q3t(:,:,:,iqa:iqb) )
       endif
       if (q_ttrans .ne. 0) then
            iqc = iqa - m_ttrans
            iqd = iqc + m_ttrans - 1
!$omp parallel do private(i,j,k,m)
          do m=iqc,iqd
             do k=beglev,endlev
                do j = beglat,endlat
                   do i=1,plon
                      q3(i,j,k,m) = q3t(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       endif
       if (r_ttrans .ne. 0) then
          call mp_recvirr  ( q3t(:,:,:,iqa:iqb), rxy_to_r%RecvDesc )
!$omp parallel do private(i,j,k,m)
          do m=iqa,iqb
             do k=beglev,endlev
                do j = beglat,endlat
                   do i=1,plon
                      q3(i,j,k,m) = q3t(i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       endif

       call t_stopf('transpose_qbck')
       call t_stopf('transpose_total')
#endif

    else

! -------------------------------------------------------------------------
! Copy onto yz decomposition
! -------------------------------------------------------------------------

       if (.not. adiabatic) then

          if (local_dp_map) then

!$omp parallel do private(lchnk, i, k, ncol, m, lons, lats)

             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                call get_lon_all_p(lchnk, ncol, lons)
                call get_lat_all_p(lchnk, ncol, lats)

                do k = 1, plev
                   do i = 1, ncol
                      dvdt(lons(i),k,lats(i)) = phys_tend(lchnk)%dvdt(i,k)
                      dudt(lons(i),k,lats(i)) = phys_tend(lchnk)%dudt(i,k)
                      pt  (lons(i),lats(i),k) = phys_state(lchnk)%t(i,k)
                      delp(lons(i),lats(i),k) = phys_state(lchnk)%pdel(i,k)
                   enddo
                enddo

                do m=1,ppcnst
                   do k=1,plev
                      do i=1,ncol
                         q3(lons(i),lats(i),k,m) = phys_state(lchnk)%q(i,k,m)
                      end do
                   end do
                end do

             enddo

          else

             tsize = 4 + ppcnst
             allocate(bpter(plon,0:plev))
             allocate(bbuffer(tsize*block_buf_nrecs))
             allocate(cbuffer(tsize*chunk_buf_nrecs))

             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)

                call chunk_to_block_send_pters(lchnk,pcols,plev+1,tsize,cpter)

                do i=1,ncol

                   do k=1,plev

                      cbuffer(cpter(i,k))   = phys_tend(lchnk)%dvdt(i,k)
                      cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt(i,k)
                      cbuffer(cpter(i,k)+2) = phys_state(lchnk)%t(i,k)
                      cbuffer(cpter(i,k)+3) = phys_state(lchnk)%pdel(i,k)

                      do m=1,ppcnst
                         cbuffer(cpter(i,k)+3+m) = phys_state(lchnk)%q(i,k,m)
                      end do

                   end do

                end do

             end do

             call transpose_chunk_to_block(tsize, cbuffer, bbuffer)

             do j=beglat,endlat
 
                call chunk_to_block_recv_pters(j,plon,plev+1,tsize,bpter)

                do i=1,nlon(j)

                   do k=1,plev

                      dvdt(i,k,j) = bbuffer(bpter(i,k))
                      dudt(i,k,j) = bbuffer(bpter(i,k)+1)
                      pt  (i,j,k) = bbuffer(bpter(i,k)+2)
                      delp(i,j,k) = bbuffer(bpter(i,k)+3)

                      do m=1,ppcnst
                         q3(i,j,k,m) = bbuffer(bpter(i,k)+3+m)
                      end do

                   end do

                end do

             end do

             deallocate(bpter)
             deallocate(bbuffer)
             deallocate(cbuffer)

          endif

       endif

    endif

    if (.not. adiabatic) then
! WS: 02.08.06: Update t3 to temperature; is this necessary for adiabatic??
!$omp parallel do private(i,j,k)
       do k=beglev,endlev
          do j = beglat,endlat
             do i=1,plon
                t3(i,k,j) = pt(i,j,k)
             enddo
          enddo
       enddo
       if(.not.is_first_step())then
          allocate (ttend  (plon, beglev:endlev))
          dtinv=1./dtime
          do j=beglat,endlat
             do k=beglev,endlev
                do i=1,plon
                   ttend(i,k)=(t3(i,k,j)-t3old(i,k,j))*dtinv
                enddo
             enddo
             call outfld('TTEND',ttend,plon,j)
	     do m=1,ppcnst
                do k=beglev,endlev
                do i=1,plon
                   ttend(i,k)=(q3(i,j,k,m)-q3old(i,j,k,m))*dtinv
                enddo
                enddo
                call outfld(tendnam(m),ttend,plon,j)
	     enddo
          enddo
	  deallocate(ttend)
       endif
       do j=beglat,endlat
          do k=beglev,endlev
             do i=1,plon
                t3old(i,k,j)=t3(i,k,j)
             enddo
          enddo
	  do m=1,ppcnst
             do k=beglev,endlev
             do i=1,plon
                q3old(i,j,k,m)=q3(i,j,k,m)
             enddo
             enddo
	  enddo
       enddo

       if (.not. full_phys) then
! WS 02.08.06 : this section of code is not not specific to 1D decomposition
!$omp parallel do private(i, j, k)
          do k=1,plev
             do j=beglat,endlat
                do i=1,plon
                   pt(i,j,k) = pt(i,j,k) / pkz(i,j,k)
                enddo
             enddo
          enddo
       endif

! -------------------------------------------------------------------------
! Update u3s and v3s from tendencies dudt and dvdt.
! -------------------------------------------------------------------------
       dt5 = .5*dtime
       call uv3s_update(dudt, u3s, dvdt, v3s, dt5,                     &
                        plon, plat, plev, beglat, endlat,              &
                        ng_d, ng_s, ng_s, ng_d, beglev, endlev)

    endif

! -------------------------------------------------------------------------
! Compute pt, q3, pe, delp, ps, peln, pkz and pk.
! For 2-D decomposition, delp is transposed to delpxy, pexy is computed
!  from delpxy (and ptop), and pexy is transposed back to pe.
! Note that pt, q3, delp and pe are input parameters as well.
! For ideal or adiabatic physics, fewer quantities are updated.
! -------------------------------------------------------------------------
    call p_d_adjust(pe,   pt,   q3,   delp,   ps,              &
                    peln, pk,   pkz,  zvir,    cappa,  delpxy, &
                    pexy, ptop, full_phys  )
 
!EOC
  end subroutine p_d_coupling
!-----------------------------------------------------------------------
end module dp_coupling
