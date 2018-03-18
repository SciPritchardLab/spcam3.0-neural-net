!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,        only: pcols, pver
   use rgrid,         only: nlon
   use pmgrid
   use phys_buffer,   only: pbuf_fld, pbuf_size_max
   use phys_grid
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: pcnst, ppcnst
   use physconst,     only: cpair, gravit, rair, zvir
   use geopotential,  only: geopotential_t
   use check_energy,  only: check_energy_timestep_init
   use pmgrid, only: masterproc
#ifdef QVORTDAMP
   use comslt, only: qfcst,dqfcst
#endif
   use time_manager, only: is_first_step, is_first_restart_step
                           

   implicit none

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(ps, t3, u3, v3, q3, &
#ifdef QVORTDAMP
                          u3aux,v3aux, &
#endif
                          omga, phis, phys_state, phys_tend, pbuf)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
!------------------------------------------------------------------------------
    use physconst,     only: cappa
    use history, only: outfld
    use time_manager, only: get_step_size,get_nstep

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (plond, beglat:endlat)            ! surface pressure
    real(r8), intent(in) :: t3  (plond, plev, beglatex:beglatex+numlats-1)  ! temperature
    real(r8), intent(in) :: u3  (plond, plev, beglatex:beglatex+numlats-1)  ! u-wind component
    real(r8), intent(in) :: v3  (plond, plev, beglatex:beglatex+numlats-1)  ! v-wind component
    real(r8), intent(in) :: q3  (plond, plev, ppcnst, beglatex:beglatex+numlats-1) ! constituents
#ifdef QVORTDAMP
    real(r8), intent(in) :: u3aux  (plond, plev, beglatex:beglatex+numlats-1)  ! u-wind component
    real(r8), intent(in) :: v3aux  (plond, plev, beglatex:beglatex+numlats-1)  ! v-wind component
#endif
    real(r8), intent(in) :: omga(plond, plev, beglat:endlat)      ! vertical velocity
    real(r8), intent(in) :: phis(plond, beglat:endlat)            ! Surface geopotential

    type(physics_state), intent(out), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(out), dimension(begchunk:endchunk) :: phys_tend
    type(pbuf_fld),    intent(inout), dimension(pbuf_size_max)     :: pbuf
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,k,j,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

        real(r8) :: ztodt
    real(r8) :: ttend(pcols,pver),qtend (pcols,pver),qitend(pcols,pver),qctend(pcols,pver)
    real (r8) :: stend(pcols,pver),utend (pcols,pver),vtend(pcols,pver),omgatend(pcols,pver)
    real (r8) :: t_debug(pcols,pver) !SR: debug output to check what t3 is
    integer :: nstep
#ifdef PVBUDGET
  real(r8):: pv1(pcols,pver,begchunk:endchunk)
  real(r8):: pv2(pcols,pver,begchunk:endchunk)
  real(r8):: pv3(pcols,pver,begchunk:endchunk)
  real(r8):: pvtend(pcols,pver,begchunk:endchunk)
#endif
  real(r8) :: aux
  real(r8):: u3loc(pcols,pver)
  real(r8):: u3auxloc(pcols,pver)
#ifdef QVORTDAMP
!  real(r8):: auxqfcst(pcols,pver)
  real(r8):: auxdqfcst(pcols,pver)
#endif
    nstep = get_nstep()
    ztodt = get_step_size()
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

#ifdef PVBUDGET
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, TTEND,UTEND,VTEND,OMGATEND,QTEND,QCTEND,QITEND,PV1,PV2,PV3,PVTEND,AUX,STEND,U3LOC,U3AUXLOC)
#else
#ifdef QVORTDAMP
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, TTEND,UTEND,VTEND,OMGATEND,QTEND,QCTEND,QITEND,AUX,STEND,U3LOC,U3AUXLOC,AUXDQFCST)
#else
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, TTEND,UTEND,VTEND,OMGATEND,QTEND,QCTEND,QITEND,AUX,STEND)
#endif
#endif
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          phys_state(lchnk)%ncol  = ncol
          phys_state(lchnk)%lchnk = lchnk

          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = ps  (lons(i),lats(i))
             phys_state(lchnk)%phis (i)     = phis(lons(i),lats(i))
          end do

          do k=1,plev
             do i=1,ncol
               ! MSP: store the influence of the SLD dycore relative to 
               ! the remembered state from the 
               ! physics package state structure:
#ifdef QVORTDAMP
               u3loc(i,k) = u3 (lons(i),k,lats(i))
               u3auxloc(i,k) = u3aux (lons(i),k,lats(i))
               auxdqfcst(i,k) = dqfcst(lons(i),k,lats(i))/ztodt
#endif
               t_debug(i,k) = t3 (lons(i),k,lats(i))
               ttend(i,k) = (t3 (lons(i),k,lats(i)) - phys_state(lchnk)%t(i,k))/ztodt
               utend(i,k) = (u3  (lons(i),k,lats(i)) - phys_state(lchnk)%u    (i,k))/ztodt
               vtend(i,k) = (v3  (lons(i),k,lats(i)) - phys_state(lchnk)%v    (i,k) )/ztodt
               omgatend(i,k) = (omga(lons(i),k,lats(i)) - phys_state(lchnk)%omega (i,k))/ztodt
               ! For the vapor species we need to wait until after the conversion
               ! to moist mixing ratios so just store the old state for now:
                qtend (i,k) = phys_state(lchnk)%q(i,k,1)
                qctend (i,k) = phys_state(lchnk)%q(i,k,2)
                qitend (i,k) = phys_state(lchnk)%q(i,k,3)

                phys_state(lchnk)%u    (i,k)   = u3  (lons(i),k,lats(i))
                phys_state(lchnk)%t    (i,k)   = t3  (lons(i),k,lats(i))
                phys_state(lchnk)%v    (i,k)   = v3  (lons(i),k,lats(i))
                phys_state(lchnk)%omega(i,k)   = omga(lons(i),k,lats(i))
                phys_state(lchnk)%q(i,k,1)     = q3  (lons(i),k,1,lats(i))
             end do
          end do

          ! convert constituents (except specific humidity) from dry to moist basis
          do m=2,ppcnst
             do k=1,plev
                do i=1,ncol
                   phys_state(lchnk)%q(i,k,m) = q3(lons(i),k,m,lats(i))*(1. - q3(lons(i),k,1,lats(i)))
                end do
             end do
          end do
         if (.not. is_first_step() .and. .not. is_first_restart_step()) then
           do k=1,plev
             do i=1,ncol
               qtend(i,k) = -(qtend(i,k) - phys_state(lchnk)%q(i,k,1))/ztodt
               qctend(i,k) = -(qctend(i,k) - phys_state(lchnk)%q(i,k,2))/ztodt
               qitend(i,k) = -(qitend(i,k) - phys_state(lchnk)%q(i,k,3))/ztodt
             end do
           end do
           call outfld('DBGT4',t_debug,pcols,lchnk)
           call outfld('TTEND',ttend,pcols,lchnk)
           call outfld('UTEND',utend,pcols,lchnk)
           call outfld('VTEND',vtend,pcols,lchnk)
           call outfld('OMGATEND',omgatend,pcols,lchnk)
           call outfld('QTEND',qtend,pcols,lchnk)
           call outfld('QCTEND',qctend,pcols,lchnk)
           call outfld('QITEND',qitend,pcols,lchnk)
#ifdef QVORTDAMP
           call outfld ('XXXU3',u3loc,pcols,lchnk)
           call outfld ('XXXU3AUX',u3auxloc,pcols,lchnk)
           call outfld ('DQFCST',auxdqfcst,pcols,lchnk)
!           call outfld ('QFCST',auxqfcst,pcols,lchnk)
#endif
         end if

       end do

    else

       tsize = 4 + ppcnst

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

       do j=beglat,endlat

          call block_to_chunk_send_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

             bbuffer(bpter(i,0))   = ps  (i,j)
             bbuffer(bpter(i,0)+1) = phis(i,j)

             do k=1,plev

                bbuffer(bpter(i,k))   = t3  (i,k,j)
                bbuffer(bpter(i,k)+1) = u3  (i,k,j)
                bbuffer(bpter(i,k)+2) = v3  (i,k,j)
                bbuffer(bpter(i,k)+3) = omga(i,k,j)
                bbuffer(bpter(i,k)+4) = q3  (i,k,1,j)

                ! convert constituents (except specific humidity) from dry to moist basis
                do m=2,ppcnst
                   bbuffer(bpter(i,k)+3+m) = q3(i,k,m,j)*(1. - q3(i,k,1,j))
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

             phys_state(lchnk)%ps   (i)     = cbuffer(cpter(i,0))
             phys_state(lchnk)%phis (i)     = cbuffer(cpter(i,0)+1)

             do k=1,plev

                phys_state(lchnk)%t    (i,k)   = cbuffer(cpter(i,k))
                phys_state(lchnk)%u    (i,k)   = cbuffer(cpter(i,k)+1)
                phys_state(lchnk)%v    (i,k)   = cbuffer(cpter(i,k)+2)
                phys_state(lchnk)%omega (i,k)   = cbuffer(cpter(i,k)+3)

                do m=1,ppcnst
                   phys_state(lchnk)%q (i,k,m) = cbuffer(cpter(i,k)+3+m)
                end do

             end do

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)

    endif

!-----------------------------------------------------------------------
! Fill auxilliary arrays in physics data structure
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS, AUX, STEND, ZTODT)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol

! pressure arrays
       call plevs0(ncol, pcols, pver, &
                   phys_state(lchnk)%ps,   phys_state(lchnk)%pint,    &
                   phys_state(lchnk)%pmid, phys_state(lchnk)%pdel)

! log(pressure) arrays and Exner function
       do k=1,pver+1
          do i=1,ncol
             phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
          end do
       end do
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%rpdel(i,k)  = 1./phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

! Compute initial geopotential heights
       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                            phys_state(lchnk)%t     , phys_state(lchnk)%q(1,1,1), rair,  gravit,  zvir    , &
                            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
             aux = cpair*phys_state(lchnk)%t(i,k) + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
             stend(i,k) = (aux - phys_state(lchnk)%s(i,k))/ztodt
             phys_state(lchnk)%s(i,k) = aux
#ifdef PVBUDGET
             pvtend(i,k,lchnk) = phys_state(lchnk)%pv(i,k) ! last updated after tphysac.
#endif
          end do
       end do

! Compute energy and water integrals of input state
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf)

       if (nstep .gt. 2 ) then
         call outfld('STEND',stend,pcols,lchnk)
       end if

    end do

#ifdef PVBUDGET
  ! Mike Pritchard, calculate newest PV
    call calculate_physics_PV (phys_state,pv1,pv2,pv3) ! updates state%pv

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL )

    do lchnk=begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      pvtend(:ncol,:pver,lchnk) = (phys_state(lchnk)%pv(:ncol,:pver) - pvtend(:ncol,:pver,lchnk))/ztodt
       call outfld ('PVTEND',pvtend(:,:,lchnk),pcols,lchnk) 
    end do
#endif


    return
  end subroutine d_p_coupling

!===============================================================================
  subroutine p_d_coupling(phys_state, phys_tend, t2, fu, fv, flx_net, q3)
!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend), intent(in), dimension(begchunk:endchunk) :: phys_tend

    real(r8), intent(out) :: t2(plond, plev, beglat:endlat)        ! temp tendency
    real(r8), intent(out) :: fu(plond, plev, beglat:endlat)        ! u wind tendency
    real(r8), intent(out) :: fv(plond, plev, beglat:endlat)        ! v wind tendency
    real(r8), intent(out) :: flx_net(plond,beglat:endlat)          ! net flux
    real(r8), intent(out) :: q3(plond, plev, ppcnst, beglatex:beglatex+numlats-1) ! constituents
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,j,k,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
!-----------------------------------------------------------------------

    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         call get_lon_all_p(lchnk, ncol, lons)
         call get_lat_all_p(lchnk, ncol, lats)

         do k=1,plev
            do i=1,ncol
               t2(lons(i),k,lats(i))   = phys_tend(lchnk)%dTdt (i,k)
               fu(lons(i),k,lats(i))   = phys_tend(lchnk)%dudt (i,k)
               fv(lons(i),k,lats(i))   = phys_tend(lchnk)%dvdt (i,k)
               q3(lons(i),k,1,lats(i)) = phys_state(lchnk)%q(i,k,1)
            end do
         end do

         do i=1,ncol
            flx_net(lons(i),lats(i))   = phys_tend(lchnk)%flx_net(i)
         end do

         ! convert constituents (except specific humidity) from moist to dry basis
         do m=2,ppcnst
            do k=1,plev
               do i=1,ncol
                  q3(lons(i),k,m,lats(i)) = phys_state(lchnk)%q(i,k,m) /     &
                                            (1. - phys_state(lchnk)%q(i,k,1))
               end do
            end do
         end do
      end do

    else

       tsize = 3 + ppcnst

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,plev+1,tsize,cpter)

          do i=1,ncol

             cbuffer(cpter(i,0)) = phys_tend(lchnk)%flx_net(i)

             do k=1,plev

                cbuffer(cpter(i,k))   = phys_tend(lchnk)%dTdt (i,k)
                cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt (i,k)
                cbuffer(cpter(i,k)+2) = phys_tend(lchnk)%dvdt (i,k)
                cbuffer(cpter(i,k)+3) = phys_state(lchnk)%q(i,k,1)

                ! convert constituents (except specific humidity) from moist to dry basis
                do m=2,ppcnst
                   cbuffer(cpter(i,k)+2+m) = phys_state(lchnk)%q(i,k,m) /     &
                                             (1. - phys_state(lchnk)%q(i,k,1))
                end do

             end do

          end do

       end do

       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)

       do j=beglat,endlat

          call chunk_to_block_recv_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

             flx_net(i,j) = bbuffer(bpter(i,0))

             do k=1,plev

                t2(i,k,j) = bbuffer(bpter(i,k))
                fu(i,k,j) = bbuffer(bpter(i,k)+1)
                fv(i,k,j) = bbuffer(bpter(i,k)+2)

                do m=1,ppcnst
                   q3(i,k,m,j) = bbuffer(bpter(i,k)+2+m)
                end do

             end do

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)

    endif

    return
  end subroutine p_d_coupling
end module dp_coupling
