#include <misc.h>
#include <params.h>

module inidat
!BOP
!
! !MODULE: inidat --- dynamics-physics coupling module
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use comsrf

! !PUBLIC MEMBER FUNCTIONS:
   public read_inidat, copy_inidat

! !PUBLIC DATA MEMBERS:
   real(r8), allocatable :: ps_tmp(:,:)
   real(r8), allocatable :: phis_tmp(:,:)
   real(r8), allocatable :: landfrac_tmp(:,:)                  
   real(r8), allocatable :: tsocn_tmp(:,:)                  
   real(r8), allocatable :: icefrac_tmp(:,:)                  
   real(r8), allocatable :: landm_tmp(:,:)
   real(r8), allocatable :: pblht_tmp(:,:)
   real(r8), allocatable :: tpert_tmp(:,:)
   real(r8), allocatable :: qpert_tmp(:,:)
   real(r8), allocatable :: sgh_tmp(:,:)
   real(r8), allocatable :: tsice_tmp(:,:)
   real(r8), allocatable :: tsice_rad_tmp(:,:)                   
   real(r8), allocatable :: tbot_tmp(:,:)                  
   real(r8), allocatable :: tssub_tmp(:,:,:)
   real(r8), allocatable :: sicthk_tmp(:,:)
   real(r8), allocatable :: snowhice_tmp(:,:)                

   real(r8) zgsint_tmp

   logical read_tsicerad
   logical read_tbot
   logical read_pblh
   logical read_tpert
   logical read_qpert
   logical read_cloud
   logical read_qcwat
   logical read_tcwat
   logical read_lcwat

!
! !DESCRIPTION:
!
!      This module provides 
!
!      \begin{tabular}{|l|l|} \hline \hline
!        read\_inidat    &   \\ \hline
!        copy\_inidat    &   \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   YY.MM.DD   ?????      Creation
!   00.06.01   Grant      First attempt at modifying for LRDC
!   01.10.01   Lin        Various revisions
!   01.01.15   Sawyer     Bug fixes for SPMD mode
!   01.03.26   Sawyer     Added ProTeX documentation
!   02.04.04   Sawyer     Removed comspe
!   03.08.31   Mirin      Eliminated many 3D temporary global arrays
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: read_inidat --- read initial dataset
!
! !INTERFACE: 
   subroutine read_inidat

! !USES:
      use pmgrid
      use pspect
      use rgrid
      use commap
      use prognostics
      use physconst,    only: gravit
      use history,      only: fillvalue
      use constituents, only: pcnst, pnats, cnst_name, qmin, cnst_read_iv
      use chemistry,    only: chem_implements_cnst, chem_init_cnst
      ! TBH:  combine modules aerosol_intr and aerosols?  
      use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
      use test_tracers, only: test_tracers_implements_cnst, test_tracers_init_cnst
      use cldcond,      only: cldcond_implements_cnst, cldcond_init_cnst
      use phys_buffer,  only: pbuf, pbuf_times, pbuf_get_fld_idx
      use phys_grid
#if ( defined SPMD )
      use mpishorthand
      use spmd_dyn, only : comm_y, comm_z
      use parutilitiesmodule, only: parcollective2d, BCSTOP
#endif

      implicit none

      include 'netcdf.inc'

!------------------------------Commons----------------------------------

#include <comctl.h>
#include <comqfl.h>
#include <comlun.h>
#include <perturb.h>

! !DESCRIPTION:
!
!   Read initial dataset and spectrally truncate as appropriate.
!
! !REVISION HISTORY:
!
!   00.06.01   Grant      First attempt at modifying for LRDC
!   00.10.01   Lin        Various revisions
!   01.01.15   Sawyer     Bug fixes for SPMD mode
!   01.03.09   Eaton      Modifications
!   01.03.26   Sawyer     Added ProTeX documentation
!   03.08.31   Mirin      Eliminated many 3D temporary global arrays
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer i,j,k,m,lat,lon   ! grid and constituent indices
      integer ihem              ! hemisphere index
      real(r8) pdelb(plond,plev)! pressure diff between interfaces
      real(r8) pertval          ! perturbation value
      real(r8) zgssum           ! partial sums of phis
      integer ii, ic, n, mq, ml1, ml2
!
! Netcdf related variables
!
      integer lonsiz, latsiz, levsiz ! Dimension sizes
      integer londimid, levdimid, latdimid ! Dimension ID's
      integer tid                 ! Variable ID's
      integer tracid(pcnst+pnats) ! Variable ID's
      integer phisid, sghid, psid ! Variable ID's
      integer landmid
      integer pblhtid
      integer tpertid
      integer qpertid
      integer cldid
      integer qcwatid
      integer tcwatid
      integer lcwatid
#if ( ! defined COUP_CSM )
      integer ts1id, ts2id, ts3id, ts4id,tsiceid,tsice_rad_id ! Variable ID's
      integer tbotid             ! Variable ID's
#endif
#if ( defined COUP_SOM )
      integer sicid
      integer icefracid
      integer tsocnid
#endif
      integer snowhiceid           ! Variable ID's
      integer landfracid        ! Variable ID's
      integer usid, vsid
      integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
      integer strt3d(4)         ! start lon, lev, lat, time for netcdf 3-d
      data strt2d/3*1/          ! Only index 2 will ever change
      data strt3d/4*1/          ! Only indices 2,3 will ever change

      integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
      integer cnt3d(4)          ! lon, lat, lev, time counts for netcdf 2-d
      data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
      data cnt3d/plon,plev,plat,1/ ! 3-d arrs: Always grab a full time slice

      integer ndims2d           ! number of dimensions
      integer dims2d(NF_MAX_VAR_DIMS) ! variable shape
      integer ndims3d           ! number of dimensions
      integer dims3d(NF_MAX_VAR_DIMS) ! variable shape
      integer tmptype
      integer natt, ret, attlen ! netcdf return values
      logical phis_hires        ! true => PHIS came from hi res topo
      real(r8), allocatable :: arrxyz(:,:,:)
      real(r8), allocatable :: arrxzy(:,:,:)
      real(r8), allocatable :: uv_local(:,:,:)
      real(r8), pointer, dimension(:,:,:,:) :: cld
      real(r8), pointer, dimension(:,:,:,:) :: tcwat
      real(r8), pointer, dimension(:,:,:,:) :: qcwat
      real(r8), pointer, dimension(:,:,:,:) :: lcwat

      character*(NF_MAX_NAME) tmpname
      character*256 text
      character*80 trunits      ! tracer untis

      integer istat
      integer slatid, slatdimid, slatsiz
      integer slonid, slondimid, slonsiz
      integer cnt3dus(4)        ! index counts for netcdf U staggered grid
      integer cnt3dvs(4)        ! index counts for netcdf V staggered grid

!      data cnt3dus/plon,plev,splat,1/ ! 3-d arrs: Always grab a full time slice
! SJL
      integer platm1
      parameter (platm1=plat-1)
      data cnt3dus/plon,plev,platm1,1/ ! temporary patch
      data cnt3dvs/plon,plev,plat,1/ ! 3-d arrs: Always grab a full time slice

!
!-----------------------------------------------------------------------
!     August 2003 revision described below (Mirin)
!-----------------------------------------------------------------------
! This routine has the master process reading in global arrays from disk
!   and storing them in temporary arrays. The data is then scattered to
!   the other processes in copy_inidat. Originally, all the data was
!   read in, and subsequently all of it was scattered.
! Because of the large volume of temporary global data, particularly at
!   high resolution, the procedure was modified to use fewer temporary
!   global arrays. This required interleaving reading and scattering of
!   data. Scattering of the 3D arrays (and a few others) is now part of
!   read_inidat; the remaining quantities are scattered in copy_inidat.
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
! Allocate memory for temporary arrays
!-----------------------------------------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter 
!
      allocate ( ps_tmp(plond,plat), stat=istat )
      if ( istat /= 0 ) call endrun
      allocate ( arrxyz(plond,plat,plev), stat=istat )
      if ( istat /= 0 ) call endrun
      allocate ( arrxzy(plond,plev,plat), stat=istat )             
      if ( istat /= 0 ) call endrun
      allocate( uv_local(plond,beglat:endlat,beglev:endlev), stat=istat )
      if ( istat /= 0 ) call endrun
      allocate ( phis_tmp(plond,plat), stat=istat )        
      if ( istat /= 0 ) call endrun
      allocate ( landm_tmp(plond,plat), stat=istat )                  
      if ( istat /= 0 ) call endrun
      allocate ( sgh_tmp(plond,plat), stat=istat )                 
      if ( istat /= 0 ) call endrun
      allocate ( tsice_tmp(plond,plat), stat=istat )                   
      if ( istat /= 0 ) call endrun
      allocate ( tsice_rad_tmp(plond,plat), stat=istat )                   
      if ( istat /= 0 ) call endrun
      allocate ( tbot_tmp    (plond,plat), stat=istat )                
      if ( istat /= 0 ) call endrun
      allocate ( tssub_tmp(plond,plevmx,plat), stat=istat )         
      if ( istat /= 0 ) call endrun
      allocate ( sicthk_tmp(plond,plat), stat=istat )               
      if ( istat /= 0 ) call endrun
      allocate ( snowhice_tmp(plond,plat), stat=istat )                
      if ( istat /= 0 ) call endrun
      allocate ( landfrac_tmp(plond,plat), stat=istat )                
      if ( istat /= 0 ) call endrun
      allocate ( pblht_tmp(plond,plat), stat=istat )
      if ( istat /= 0 ) call endrun
      allocate ( tpert_tmp(plond,plat), stat=istat )
      if ( istat /= 0 ) call endrun
      allocate ( qpert_tmp(plond,plat), stat=istat )
      if ( istat /= 0 ) call endrun
#if ( defined COUP_SOM )
      allocate ( icefrac_tmp(plond,plat) )                
      allocate ( tsocn_tmp(plond,plat) )                
#endif
!
!-----------------------------------------------------------------------
! Read in input variables
!-----------------------------------------------------------------------

!
! logical flags to track which "extra" fields are indeed in IC file
!
      read_tsicerad = .false.
      read_tbot     = .false.
      read_pblh     = .false.
      read_tpert    = .false.
      read_qpert    = .false.
      read_cloud    = .false.
      read_qcwat    = .false.
      read_tcwat    = .false.
      read_lcwat    = .false.
      qpert_tmp(:plon,:) = 0.

      if (masterproc) then
!
! Get dimension IDs and lengths 
!
         call wrap_inq_dimid  (ncid_ini, 'lat', latdimid)
         call wrap_inq_dimlen (ncid_ini, latdimid, latsiz)
         call wrap_inq_dimid  (ncid_ini, 'lev', levdimid)
         call wrap_inq_dimlen (ncid_ini, levdimid, levsiz)
         call wrap_inq_dimid  (ncid_ini, 'lon', londimid)
         call wrap_inq_dimlen (ncid_ini, londimid, lonsiz)
         call wrap_inq_dimid  (ncid_ini, 'slat', slatdimid)
         call wrap_inq_dimlen (ncid_ini, slatdimid, slatsiz)
         call wrap_inq_dimid  (ncid_ini, 'slon', slondimid)
         call wrap_inq_dimlen (ncid_ini, slondimid, slonsiz)

!
! Get variable id's 
! Check that all tracer units are in mass mixing ratios
!
!         call wrap_inq_varid (ncid_ini, 'U'   , uid)
!         call wrap_inq_varid (ncid_ini, 'V'   , vid)

         call wrap_inq_varid (ncid_ini, 'slat', slatid)
         call wrap_inq_varid (ncid_ini, 'slon', slonid)
         call wrap_inq_varid (ncid_ini, 'US'  , usid)
         call wrap_inq_varid (ncid_ini, 'VS'  , vsid)

         call wrap_inq_varid (ncid_ini, 'T'   , tid)
         call wrap_inq_varid (ncid_ini, 'PS'  , psid)
         call wrap_inq_varid (ncid_ini, 'PHIS', phisid)
         call wrap_inq_varid (ncid_ini, 'SGH' , sghid)
         if (nf_inq_varid (ncid_ini, 'LANDM_COSLAT', landmid) /= nf_noerr) then
            write(6,*)'INIDAT: LANDM_COSLAT not found on initial dataset.'
            write(6,*)'        Need to run definesurf to create it.'
            write(6,*)'        This field became a requirement as of cam2_0_2_dev43'
            call endrun ()
         end if

         if ( nf_inq_varid(ncid_ini, 'PBLH', pblhtid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'PBLH' , pblhtid)
            read_pblh  = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'TPERT', tpertid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'TPERT', tpertid)
            read_tpert = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'QPERT', qpertid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'QPERT', qpertid)
            read_qpert = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'CLOUD', cldid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'CLOUD', cldid  )
            read_cloud = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'QCWAT', qcwatid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'QCWAT', qcwatid)
            read_qcwat = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'TCWAT', tcwatid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'TCWAT', tcwatid)
            read_tcwat = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'LCWAT', lcwatid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'LCWAT', lcwatid)
            read_lcwat = .true.
         end if

#if ( ! defined COUP_CSM )
!
! For land-fraction check if the variable name LANDFRAC is on the dataset if not assume FLAND
!
         if ( nf_inq_varid(ncid_ini, 'LANDFRAC', landfracid ) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'LANDFRAC', landfracid)
         else
            call wrap_inq_varid (ncid_ini, 'FLAND', landfracid)
         end if
         if ( nf_inq_varid(ncid_ini, 'TBOT', tbotid) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'TBOT'    , tbotid)
            read_tbot     = .true.
         end if
         if ( nf_inq_varid(ncid_ini, 'TSICERAD', tsice_rad_id) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'TSICERAD', tsice_rad_id)
            read_tsicerad = .true.
         end if
         call wrap_inq_varid (ncid_ini, 'TSICE', tsiceid)
         call wrap_inq_varid (ncid_ini, 'TS1', ts1id)
         call wrap_inq_varid (ncid_ini, 'TS2', ts2id)
         call wrap_inq_varid (ncid_ini, 'TS3', ts3id)
         call wrap_inq_varid (ncid_ini, 'TS4', ts4id)
         call wrap_inq_varid (ncid_ini, 'SNOWHICE', snowhiceid)
#if ( defined COUP_SOM )
         call wrap_inq_varid (ncid_ini, 'SICTHK', sicid)
         call wrap_inq_varid (ncid_ini, 'ICEFRAC', icefracid)
         call wrap_inq_varid (ncid_ini, 'TSOCN', tsocnid)
#endif
#endif
!
! Guard:  Check that "Q" is on IC file
!
         if (cnst_read_iv(1) .and. &
              nf_inq_varid(ncid_ini, cnst_name(1), tracid(1) ) /= NF_NOERR ) then
            write(6,*) 'Error:  ',cnst_name(1), ' not found on IC file'
            call endrun
         end if

         do m=1,pcnst+pnats
            if (cnst_read_iv(m) .and. &
               nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
               call wrap_inq_varid (NCID_INI,cnst_name(m), tracid(m))
               call wrap_get_att_text (NCID_INI,tracid(m),'units', trunits)
               if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
                  write(6,*)'INIDAT: tracer units for tracer = ', &
                            cnst_name(m),' must be in KG/KG'
                  call endrun
               endif
            end if
         end do
!
! Check dimension ordering for one 2-d and one 3-d field.
! Assume other arrays of like rank will have dimensions ordered the same.
!
         call wrap_inq_var (ncid_ini, psid, tmpname, tmptype, &
                            ndims2d, dims2d, natt)
         if (dims2d(1).ne.londimid .or. dims2d(2).ne.latdimid .or. &
             ndims2d.gt.3) then
            write(6,*)'INIDAT: Bad number of dims or ordering on 2d fld'
            call endrun
         end if
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
         ret = nf_inq_attlen (ncid_ini, phisid, 'from_hires', attlen)
         if (ret.eq.NF_NOERR .and. attlen.gt.256) then
            write(6,*)'INIDAT: from_hires attribute length is too long'
            call endrun
         end if
         ret = nf_get_att_text (ncid_ini, phisid, 'from_hires', text)
         if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
            phis_hires = .true.
!            write(6,*)'INIDAT: Will filter input PHIS: attribute ', &
!                      'from_hires is true'
         else
            phis_hires = .false.
!            write(6,*)'INIDAT: Will not filter input PHIS: attribute ', &
!                      'from_hires is either false or not present'
         end if
!
! Read in 2d fields.  
! For stand alone run: get surface temp and 4 (sub)surface temp fields
! For stand alone run with slab-ocean: get sea-ice thickness and snow cover
!
         do j=1,plat
            strt2d(2) = j
            if (ideal_phys .or. aqua_planet) then
               do i=1,nlon(j)
                  phis_tmp(i,j) = 0.
                  sgh_tmp (i,j) = 0.
               end do
            else
               call wrap_get_vara_realx (ncid_ini, phisid, strt2d, cnt2d, &
                                         phis_tmp(1,j))
               call wrap_get_vara_realx (ncid_ini, sghid , strt2d, cnt2d, &
                                         sgh_tmp(1,j))
            endif
            call wrap_get_vara_realx (ncid_ini, landmid, strt2d, cnt2d, &
                                      landm_tmp(1,j))
            call wrap_get_vara_realx (ncid_ini, psid, strt2d, cnt2d, &
                                      ps_tmp(1,j))
#if ( ! defined COUP_CSM )
            if (aqua_planet) then
               do i=1,nlon(j)
                  landfrac_tmp(i,j) = 0.
                  tbot_tmp(i,j) = 0.
               end do
            else
               call wrap_get_vara_realx (ncid_ini, landfracid, strt2d, cnt2d, &
                                         landfrac_tmp(1,j))
               if (read_tbot) then
                  call wrap_get_vara_realx (ncid_ini, tbotid    , strt2d, cnt2d, tbot_tmp    (1,j))
               end if
            endif
            call wrap_get_vara_realx (ncid_ini, tsiceid, strt2d, cnt2d, &
                                           tsice_tmp(1,j))
            call wrap_get_vara_realx (ncid_ini, ts1id, strt2d, cnt2d, &
                                           tssub_tmp(1,1,j))
            call wrap_get_vara_realx (ncid_ini, ts2id, strt2d, cnt2d, &
                                           tssub_tmp(1,2,j))
            call wrap_get_vara_realx (ncid_ini, ts3id, strt2d, cnt2d, &
                                           tssub_tmp(1,3,j))
            call wrap_get_vara_realx (ncid_ini, ts4id, strt2d, cnt2d, &
                                           tssub_tmp(1,4,j))
            if (read_tsicerad) then
               call wrap_get_vara_realx (ncid_ini, tsice_rad_id, strt2d, cnt2d, tsice_rad_tmp(1,j))
            end if
            if (read_pblh) then
               call wrap_get_vara_realx (ncid_ini, pblhtid, strt2d, cnt2d, pblht_tmp(1,j  ))
            end if
            if (read_tpert) then
               call wrap_get_vara_realx (ncid_ini, tpertid, strt2d, cnt2d, tpert_tmp(1,j  ))
            end if
            if (read_qpert) then
               call wrap_get_vara_realx (ncid_ini, qpertid, strt2d, cnt2d, qpert_tmp(1,j))
            end if
!
! Set sea-ice thickness and snow cover:
!
#if ( defined COUP_SOM )
            call wrap_get_vara_realx(ncid_ini, sicid, strt2d, cnt2d, sicthk_tmp(1,j))
            call wrap_get_vara_realx(ncid_ini, icefracid, strt2d, cnt2d, icefrac_tmp(1,j))
            call wrap_get_vara_realx(ncid_ini, tsocnid, strt2d, cnt2d, tsocn_tmp(1,j))
#endif
            call wrap_get_vara_realx(ncid_ini, snowhiceid, strt2d, cnt2d, snowhice_tmp(1,j))
#endif
         end do
!
! Read in 3d fields.  

! Staggered grid variables and transpose

! Read in u3

         call wrap_get_vara_realx(ncid_ini,usid, strt3d, cnt3dus, arrxzy)

         do k = 1, plev
!
! SJL: initialize j=1 because later on arrxyz will be copied to u3s using f90 array syntax
               do i = 1, plon
                  arrxyz(i,1,k) = fillvalue
               enddo
!
            do j = 1, plat-1
               do i = 1, plon
                  arrxyz(i,j+1,k) = arrxzy(i,k,j)
               enddo
            enddo
         enddo

      endif                     ! end of if-masterproc

! Scatter u3

#if ( defined SPMD )
      call scatter( arrxyz, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               u3s(i,j,k) = uv_local(i,j,k)
            enddo
         enddo
      enddo
#else
!$omp parallel do private(i, j, k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               u3s(i,j,k) = arrxyz(i,j,k)
            enddo
         enddo
      enddo
#endif

      if (masterproc) then

! Read in v3

         call wrap_get_vara_realx(ncid_ini,vsid, strt3d, cnt3dvs, arrxzy)

         do k = 1, plev
            do j = 1, plat
               do i = 1, plon
                  arrxyz(i,j,k) = arrxzy(i,k,j)
               enddo
            enddo
         enddo

      endif                     ! end of if-masterproc

! Scatter v3

#if ( defined SPMD )
      call scatter( arrxyz, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               v3s(i,j,k) = uv_local(i,j,k)
            enddo
         enddo
      enddo
#else
!$omp parallel do private(i, j, k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               v3s(i,j,k) = arrxyz(i,j,k)
            enddo
         enddo
      enddo
#endif

      if (masterproc) then

! Read in t3

         call wrap_get_vara_realx(ncid_ini, tid, strt3d, cnt3d, arrxzy)

!         
! Add random perturbation to temperature if required
!
         if (pertlim.ne.0.0) then
            write(6,*)'INIDAT: Adding random perturbation bounded by +/-', &
                      pertlim,' to initial temperature field'
            do lat=1,plat
               do k=1,plev
                  do i=1,nlon(lat)
                     call random_number (pertval)
                     pertval = 2.*pertlim*(0.5 - pertval)
                     arrxzy(i,k,lat) = arrxzy(i,k,lat)*(1. + pertval)
                  end do
               end do
            end do
         endif
!$omp parallel do private(k)
         do k = 1, plev
            call xpavg(arrxzy(1,k,   1), plon)
            call xpavg(arrxzy(1,k,plat), plon)
         enddo

      endif                     ! end of if-masterproc

! Scatter t3

#if ( defined SPMD )
      call scatter( arrxzy, strip3dxzy, t3, mpicom )
#else
!$omp parallel do private(i, j, k, ic)
      do j=beglat,endlat
         do k=beglev,endlev
            do i=1,plon
               t3(i,k,j) = arrxzy(i,k,j)
            enddo
         enddo
      enddo
#endif

! IMPORTANT - the following block of coding must be kept adjacent to the t3
!    coding, as it assumes that ARRXZY contains t3 (for .not. read_tbot)

      if (masterproc) then
         if(.not. read_tbot) then
            tbot_tmp   (:plon,:) = arrxzy(:plon,plev,:)
         endif
      endif

! Read tcwat

! IMPORTANT - the following block of coding must be kept adjacent to the t3
!    coding, as it assumes that ARRXZY contains t3 (for .not. read_tcwat)

      if (masterproc) then

         if (read_tcwat) then
            call wrap_get_vara_realx(ncid_ini, tcwatid, strt3d, cnt3d, arrxzy)
         else
            write(6,*) 'Warning:  TCWAT not found on IC file; initialized with T'
         end if

      endif

! Scatter tcwat

      m = pbuf_get_fld_idx('TCWAT')
      tcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,arrxzy,tcwat(:,:,:,1))

! Read cld

      if (masterproc) then

         if (read_cloud) then
            call wrap_get_vara_realx(ncid_ini, cldid  , strt3d, cnt3d, arrxzy)
         else
            arrxzy    (:plon,:,:) = 0.
            write(6,*) 'Warning:  CLOUD not found on IC file; initialized to 0.'
         end if

      endif

! Scatter cld

      m = pbuf_get_fld_idx('CLD')
      cld => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,arrxzy,cld(:,:,:,1))

! Read qcwat

      if (masterproc) then

         if (read_qcwat) then
            call wrap_get_vara_realx(ncid_ini, qcwatid, strt3d, cnt3d, arrxzy)
         else
! Search among constituents
            mq = 0
            do m = 1, pcnst+pnats
              if (cnst_name(m) == 'Q') mq = m
            enddo
            if (mq .eq. 0) then
               write (6,*) 'Q not found - exiting'
               call endrun()
            endif
            m = mq
            write(6,*) 'Warning:  QCWAT not found on IC file; initialized with ',cnst_name(m)
            if (cnst_read_iv(m) .and. &
               nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arrxzy)
            else
! Initialize instead of input
               write(6,*) 'Warning:  Not reading ',cnst_name(m), ' from IC file.'
               arrxzy = 0.
               if (cldcond_implements_cnst(cnst_name(m))) then
                  call cldcond_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "cldcond_init_cnst"'
               else if (chem_implements_cnst(cnst_name(m))) then
                  call chem_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "chem_init_cnst"'
               else if (test_tracers_implements_cnst(cnst_name(m))) then
                  call test_tracers_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "test_tracers_init_cnst"'
               else
                  write(6,*) '          ', cnst_name(m), ' set to 0.'
               end if
            endif
            do lat=1,plat
               call qneg3('INIDAT  ', lat, nlon(lat), plon, plev, 1, &
                          qmin(m), arrxzy(:plon,:plev,lat))
            end do
            do k = 1, plev
               call xpavg(arrxzy(:,k,   1), plon)
               call xpavg(arrxzy(:,k,plat), plon)
            enddo
         end if

      endif

! Scatter qcwat

      m = pbuf_get_fld_idx('QCWAT')
      qcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,arrxzy,qcwat(:,:,:,1))

! Read lcwat

      if (masterproc) then
         if (read_lcwat) then
            call wrap_get_vara_realx(ncid_ini, lcwatid, strt3d, cnt3d, arrxzy)
         else
! Search among constituents
            ml1 = 0
            ml2 = 0
            do m = 1, pcnst+pnats
              if (cnst_name(m) == 'CLDLIQ') ml1 = m
              if (cnst_name(m) == 'CLDICE') ml2 = m
            enddo
            if (ml1 .eq. 0) then
               write (6,*) 'CLDLIQ not found - exiting'
               call endrun()
            endif
            if (ml2 .eq. 0) then
               write (6,*) 'CLDICE not found - exiting'
               call endrun()
            endif
            m = ml1
            write(6,*) 'Warning:  LCWAT not found on IC file; ',cnst_name(m),' added to LCWAT instead'
            if (cnst_read_iv(m) .and. &
               nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arrxzy)
            else
! Initialize instead of input
               write(6,*) 'Warning:  Not reading ',cnst_name(m), ' from IC file.'
               arrxzy = 0.
               if (cldcond_implements_cnst(cnst_name(m))) then
                  call cldcond_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "cldcond_init_cnst"'
               else if (chem_implements_cnst(cnst_name(m))) then
                  call chem_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "chem_init_cnst"'
               else if (test_tracers_implements_cnst(cnst_name(m))) then
                  call test_tracers_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "test_tracers_init_cnst"'
               else
                  write(6,*) '          ', cnst_name(m), ' set to 0.'
               end if
            endif
            do lat=1,plat
               call qneg3('INIDAT  ', lat, nlon(lat), plon, plev, 1, &
                          qmin(m), arrxzy(:plon,:plev,lat))
            end do
            do k = 1, plev
               call xpavg(arrxzy(:,k,   1), plon)
               call xpavg(arrxzy(:,k,plat), plon)
            enddo

! Store CLDLIQ in arrxyz
            do k = 1, plev
              do lat = 1, plat
                do lon = 1, plon
                  arrxyz(lon,lat,k) = arrxzy(lon,k,lat)
                enddo
              enddo
            enddo

            m = ml2
            write(6,*) 'Warning:  LCWAT not found on IC file; ',cnst_name(m),' added to LCWAT instead'
            if (cnst_read_iv(m) .and. &
               nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arrxzy)
            else
! Initialize instead of input
               write(6,*) 'Warning:  Not reading ',cnst_name(m), ' from IC file.'
               arrxzy = 0.
               if (cldcond_implements_cnst(cnst_name(m))) then
                  call cldcond_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "cldcond_init_cnst"'
               else if (chem_implements_cnst(cnst_name(m))) then
                  call chem_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "chem_init_cnst"'
               else if (test_tracers_implements_cnst(cnst_name(m))) then
                  call test_tracers_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "test_tracers_init_cnst"'
               else
                  write(6,*) '          ', cnst_name(m), ' set to 0.'
               end if
            endif
            do lat=1,plat
               call qneg3('INIDAT  ', lat, nlon(lat), plon, plev, 1, &
                          qmin(m), arrxzy(:plon,:plev,lat))
            end do
            do k = 1, plev
               call xpavg(arrxzy(:,k,   1), plon)
               call xpavg(arrxzy(:,k,plat), plon)
            enddo

! Add CLDLIQ and CLDICE
            do k = 1, plev
              do lat = 1, plat
                do lon = 1, plon
                  arrxzy(lon,k,lat) = arrxyz(lon,lat,k) + arrxzy(lon,k,lat)
                enddo
              enddo
            enddo
         end if

      endif

! Scatter lcwat

      m = pbuf_get_fld_idx('LCWAT')
      lcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,arrxzy,lcwat(:,:,:,1))

      if (pbuf_times > 1) then
         do n = 2, pbuf_times
            cld  (:,:,:,n) = cld  (:,:,:,1)
            tcwat(:,:,:,n) = tcwat(:,:,:,1)
            qcwat(:,:,:,n) = qcwat(:,:,:,1)
            lcwat(:,:,:,n) = lcwat(:,:,:,1)
         end do
      end if

      if (masterproc) then
!
! Compute integrals of mass, moisture, and geopotential height
!
!gg  Integrals of mass and moisture should be unnecessary in Lin-Rood dynamics
!gg  because they are conserved. What's left is the global geopotential...
!gg  Dunno if that's necessary or not, so I left it in.

         zgsint_tmp = 0.
         do lat = 1, plat
!              
! Accumulate average mass of atmosphere
!
            zgssum = 0.
            do i=1,nlon(lat)
               zgssum = zgssum + phis_tmp(i,lat)
            end do
            zgsint_tmp = zgsint_tmp + w(lat)*zgssum/nlon(lat)
         end do                  ! end of latitude loop
!
! Normalize average height
!
         zgsint_tmp = zgsint_tmp*.5/gravit
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
!
! SJL:
         tmass0 = 98222./gravit
!
! WS: Average pole information moved here: treat the global arrays
!
!-----------------------------------------------------------
! Average T, PS, PHIS and Q at the poles.       The initial
! conditions *should* already have these variables averaged,
! but do it here for safety -- no harm if it's already done.
!-----------------------------------------------------------

         call xpavg(phis_tmp(1, 1), plon)
         call xpavg(phis_tmp(1,plat), plon)
         call xpavg(ps_tmp(1,   1), plon)
         call xpavg(ps_tmp(1,plat), plon)

         if (ideal_phys) tmass0 = 100000./gravit
!        write(6,800) tmassf_tmp,tmass0,qmassf_tmp
!        write(6,810) zgsint_tmp
800      format('INIDAT: MASS OF INITIAL DATA BEFORE CORRECTION = ' &
                ,1p,e20.10,/,' DRY MASS WILL BE HELD = ',e20.10,/, &
                ' MASS OF MOISTURE AFTER REMOVAL OF NEGATIVES = ',e20.10) 
810      format(/69('*')/'INIDAT: Globally averaged geopotential ', &
                'height = ',f16.10,' meters'/69('*')/)

      endif                     ! end of if-masterproc

! Scatter ps and phis

#if ( defined SPMD )
      if (myid_z .eq. 0) then
         call scatter( ps_tmp, strip2d, ps, comm_y )
         call scatter( phis_tmp, strip2d, phis, comm_y )
      endif
      if (twod_decomp .eq. 1) then
         call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, ps ) 
         call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, phis ) 
      endif
#else
      ps(:,:) = ps_tmp(:,:)
      phis(:,:) = phis_tmp(:,:)
#endif

! Initialize constituents

      do m = 1, pcnst+pnats

         if (masterproc) then

            if (cnst_read_iv(m) .and. &
               nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arrxzy)
            else
               write(6,*) 'Warning:  Not reading ',cnst_name(m), ' from IC file.'
               arrxzy = 0.
               if (cldcond_implements_cnst(cnst_name(m))) then
                  call cldcond_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "cldcond_init_cnst"'
               else if (chem_implements_cnst(cnst_name(m))) then
                  call chem_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "chem_init_cnst"'
               else if (test_tracers_implements_cnst(cnst_name(m))) then
                  call test_tracers_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "test_tracers_init_cnst"'
               else if (aerosol_implements_cnst(cnst_name(m))) then
                  call aerosol_init_cnst(cnst_name(m), arrxzy)
                  write(6,*) '          ', cnst_name(m), ' initialized by "aerosol_init_cnst"'
               else
                  write(6,*) '          ', cnst_name(m), ' set to 0.'
               end if
            endif

!$omp parallel do private(lat)
            do lat=1,plat
               call qneg3('INIDAT  ', lat, nlon(lat), plon, plev, 1, &
                          qmin(m), arrxzy(:plon,:plev,lat))
            end do

            do k = 1, plev
               call xpavg(arrxzy(:,k,   1), plon)
               call xpavg(arrxzy(:,k,plat), plon)
            enddo

            do k = 1, plev
               do j = 1, plat
                  arrxyz(:,j,k) = arrxzy(:,k,j)
               enddo
            enddo

         end if   ! masterproc

         call scatter_q_field_to_block(arrxyz, m)

      end do
!
!-----------------------------------------------------------------------
! Copy temporary arrays to model arrays
!-----------------------------------------------------------------------
!
      call copy_inidat
!
!-----------------------------------------------------------------------
! Deallocate memory for temporary arrays
!-----------------------------------------------------------------------
!
      deallocate ( ps_tmp )
      deallocate ( arrxyz )
      deallocate ( arrxzy )
      deallocate ( uv_local )
      deallocate ( phis_tmp )        
      deallocate ( landm_tmp )
      deallocate ( sgh_tmp )
      deallocate ( tsice_tmp )
      deallocate ( tsice_rad_tmp )
      deallocate ( tbot_tmp )
      deallocate ( tssub_tmp )
      deallocate ( sicthk_tmp )
      deallocate ( snowhice_tmp )
      deallocate ( pblht_tmp )
      deallocate ( tpert_tmp )
      deallocate ( qpert_tmp )
      deallocate ( landfrac_tmp )
#if ( defined COUP_SOM )
      deallocate ( icefrac_tmp )
      deallocate ( tsocn_tmp )
#endif
!
      return
!EOC
   end subroutine read_inidat
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: copy_inidat --- Copy temporary arrays to model arrays 
!
! !INTERFACE: 
   subroutine copy_inidat

! !USES:
      use prognostics
      use buffer
      use phys_grid
      use phys_buffer, only: pbuf, pbuf_times, pbuf_get_fld_idx
#if ( defined SPMD )
      use mpishorthand
#endif

      implicit none
!------------------------------Commons----------------------------------
#include <comqfl.h>

! !DESCRIPTION:
!
! Copy temporary arrays to model arrays 
! note that the use statements below contain the definitions
! of the model arrays
!
! !REVISION HISTORY:
!
!   00.06.01   Grant      First attempt at modifying for LRDC
!   00.10.01   Lin        Various revisions
!   00.12.02   Sawyer     Use PILGRIM to scatter data sets
!   01.03.26   Sawyer     Added ProTeX documentation
!   03.08.31   Mirin      Eliminated many 3D temporary global arrays
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8), allocatable :: tmpchunk3d(:,:,:)             
      real(r8), allocatable :: tmpchunk(:,:)             
      integer i, j, ic, k, m
      integer n
      integer ncol
      real(r8) :: pmx, pmn

!-----------------------------------------------------------------------

      allocate ( tmpchunk(pcols,begchunk:endchunk) )
      allocate ( tmpchunk3d(pcols,plevmx,begchunk:endchunk) )
      tmpchunk(:,:) = 0.
      tmpchunk3d(:,:,:) = 0.
! physics variables
      call scatter_field_to_chunk(1,1,1,plond,landfrac_tmp,landfrac(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,landm_tmp,landm(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,sgh_tmp,sgh(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,tsice_tmp,tsice(1,begchunk))
#if ( defined COUP_SOM )
      call scatter_field_to_chunk(1,1,1,plond,sicthk_tmp,sicthk(1,begchunk))
#endif
      call scatter_field_to_chunk(1,1,1,plond,snowhice_tmp,snowhice(1,begchunk))

      call scatter_field_to_chunk(1,plevmx,1,plond,tssub_tmp,tmpchunk3d)
      do i =begchunk,endchunk
         ncol = get_ncols_p(i)
         surface_state2d(i)%tssub(:ncol,:) = tmpchunk3d(:ncol,:,i)
      end do

#if ( defined COUP_SOM )
      call scatter_field_to_chunk(1,1,1,plond,sicthk_tmp,sicthk(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,icefrac_tmp,icefrac(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,tsocn_tmp,tsocn(1,begchunk))

! define an initial ocean fraction and non-land ice fraction
! The 1st "where" stmt used to be done in update_srf_fractions (dev45)

      do i =begchunk,endchunk
         ncol = get_ncols_p(i)
         where (icefrac(:ncol,i) + landfrac(:ncol,i) > 1.0)
            icefrac(:ncol,i) = 1. - landfrac(:ncol,i)
         end where

         where (landfrac(:ncol,i) < 1.)
            aice(:ncol,i) = icefrac(:ncol,i)/(1. - landfrac(:ncol,i))
         elsewhere
            aice(:ncol,i) = 0.
         end where
         ocnfrac(:ncol,i) = 1. - landfrac(:ncol,i) - icefrac(:ncol,i)
      enddo
      write(6,*)'INIDAT: ocnfrac=',ocnfrac(1,begchunk)
!
! Master needs global landfrac
!
      call gather_chunk_to_field(1,1,1,plon,landfrac,landfrac_field)
!      write(6,*)'INIDAT iam=',iam,' landfrac=',landfrac
!      write(6,*)'INIDAT iam=',iam,' landfrac_field=',landfrac_field
!
!JR Could read in Focn from initial dataset if available
      Focn(:,:) = 0.
#else
      Focn(:,:) = inf
      frzmlt(:,:) = 0.  ! needs to be 0, otherwise test in tstm always true
      tsocn(:,:) = inf
#endif

! cloud and cloud water initialization should be done in their own packages.  Do it
! here for now since moving it will change answers.

      if (masterproc) then
!
! Arbitrarily initialize all "extra" fields that couldn't be found on the IC file
!
         if(.not. read_pblh) then
            pblht_tmp  (:plon,:) = 0.
            write(6,*) 'Warning:  PBLH not found on IC file; initialized to 0.'
         end if
         if(.not. read_tpert) then
            tpert_tmp  (:plon,:) = 0.
            write(6,*) 'Warning:  TPERT not found on IC file; initialized to 0.'
         end if
         if(.not. read_qpert) then
            qpert_tmp  (:plon,:) = 0.
            write(6,*) 'Warning:  QPERT not found on IC file; initialized to 0.'
         end if
         if(.not. read_tbot) then
!  tbot_tmp now initialized in read_inidat
            write(6,*) 'Warning:  TBOT not found on IC file; initialized with lowest level of T'
         end if
         if(.not. read_tsicerad) then
            tsice_rad_tmp(:plon,:) = tsice_tmp(:plon,:)
            write(6,*) 'Warning:  TSICERAD not found on IC file; initialized with TSICE'
         end if
      endif

      call scatter_field_to_chunk(1,1,1,plond,tbot_tmp,tmpchunk)
      do i =begchunk,endchunk
         ncol = get_ncols_p(i)
         surface_state2d(i)%tbot(:ncol) = tmpchunk(:ncol,i)
      end do
      call scatter_field_to_chunk(1,          1,1,plond,tsice_rad_tmp,tsice_rad(1,begchunk))
      call scatter_field_to_chunk(1,          1,1,plond,pblht_tmp,pblht(1  ,begchunk  ))
      call scatter_field_to_chunk(1,          1,1,plond,tpert_tmp,tpert(1  ,begchunk  ))

! Qpert for only the first constituent is initialized and scattered; for all other
!   constituents, qpert is set to zero.
      call scatter_field_to_chunk(1,          1,1,plond,qpert_tmp,tmpchunk(1,begchunk ))
      qpert(:,:,:) = 0.
      do i =begchunk,endchunk
        qpert(:,1,i) = tmpchunk(:,i)
      enddo
!
! Global integerals
!
      if (masterproc) then
         zgsint = zgsint_tmp
      endif
!
#if ( defined SPMD )
      call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif
      deallocate ( tmpchunk )
      deallocate ( tmpchunk3d)
!EOC
   end subroutine copy_inidat

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: scatter_q_field_to_block --- scatter a 3D constituent array to prognostic array q3
!
! !INTERFACE: 
   subroutine scatter_q_field_to_block(xyz, cnst_idx)

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid,       only: plon, plat, plev, beglat, endlat, &
                              beglev, endlev, strip3dxyz
      use prognostics,  only: q3
#if ( defined SPMD )
      use mpishorthand, only: mpicom
#endif

      implicit none

! !INPUT PARAMETERS:
      real(r8), dimension(plon,plat,plev), intent(in) :: &
         xyz        ! 3D constituent field
      integer, intent(in) :: &
         cnst_idx   ! constituent index in prognostic array q3

! !DESCRIPTION:
!
! Scatter a 3D constituent array from the master processor to the 
! q3 array in the prognostics module.
!
! !REVISION HISTORY:
!
!   02.07.31   Eaton      Initial version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer j, k
      real(r8), allocatable :: xyz_local(:,:,:)

!-----------------------------------------------------------------------

#if ( defined SPMD )

      allocate( xyz_local(plon,beglat:endlat,beglev:endlev) )
      call scatter( xyz, strip3dxyz, xyz_local, mpicom )
      do k=beglev,endlev
         do j=beglat,endlat
            q3(:,j,k,cnst_idx) = xyz_local(:,j,k)
         enddo
      enddo
      deallocate( xyz_local )

#else

      do j=beglat,endlat
         do k=beglev,endlev
            q3(:,j,k,cnst_idx) = xyz(:,j,k)
         enddo
      enddo

#endif

!EOC
   end subroutine scatter_q_field_to_block

!-----------------------------------------------------------------------
end module inidat

