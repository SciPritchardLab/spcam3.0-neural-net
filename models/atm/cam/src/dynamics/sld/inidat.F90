#include <misc.h>
#include <params.h>

module inidat

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use comsrf
   use pmgrid
   use comspe


   real(r8), allocatable :: ps_tmp(:,:)                   
   real(r8), allocatable :: u3_tmp(:,:,:)             
   real(r8), allocatable :: v3_tmp(:,:,:)             
   real(r8), allocatable :: t3_tmp(:,:,:)             
   real(r8), allocatable :: q3_tmp(:,:,:,:) 
   real(r8), allocatable :: tl_tmp(:,:,:) 
   real(r8), allocatable :: tm_tmp(:,:,:) 
   real(r8), allocatable :: ql_tmp(:,:,:) 
   real(r8), allocatable :: qm_tmp(:,:,:) 
   real(r8), allocatable :: phis_tmp(:,:)        
   real(r8), allocatable :: phisl_tmp(:,:)        
   real(r8), allocatable :: phism_tmp(:,:)        
   real(r8), allocatable :: landfrac_tmp(:,:)                  
   real(r8), allocatable :: tsocn_tmp(:,:)                  
   real(r8), allocatable :: icefrac_tmp(:,:)                  
   real(r8), allocatable :: landm_tmp(:,:)                  
   real(r8), allocatable :: pblht_tmp(:,:)                  
   real(r8), allocatable :: tpert_tmp(:,:)                  
   real(r8), allocatable :: qpert_tmp(:,:,:)                  
   real(r8), allocatable :: cld_tmp(:,:,:)             
   real(r8), allocatable :: qcwat_tmp(:,:,:)             
   real(r8), allocatable :: tcwat_tmp(:,:,:)             
   real(r8), allocatable :: lcwat_tmp(:,:,:) 
   real(r8), allocatable :: sgh_tmp(:,:)                 
   real(r8), allocatable :: tsice_tmp(:,:)                   
   real(r8), allocatable :: tsice_rad_tmp(:,:)                   
   real(r8), allocatable :: tbot_tmp(:,:)                  
   real(r8), allocatable :: tssub_tmp(:,:,:)         
   real(r8), allocatable :: dpsl_tmp(:,:)                
   real(r8), allocatable :: dpsm_tmp(:,:)                 
   real(r8), allocatable :: div_tmp(:,:,:)             
   real(r8), allocatable :: sicthk_tmp(:,:)               
   real(r8), allocatable :: snowhice_tmp(:,:)                
   real(r8) tmassf_tmp
   real(r8) qmass1_tmp 
   real(r8) qmass2_tmp 
   real(r8) zgsint_tmp 
   real(r8) qmassf_tmp 
   
   logical read_tsicerad
   logical read_tbot
   logical read_pblh
   logical read_tpert
   logical read_qpert
   logical read_cloud
   logical read_qcwat
   logical read_tcwat
   logical read_lcwat

contains

   subroutine read_inidat

!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id: inidat.F90,v 1.22.4.15 2003/12/12 22:54:58 rosinski Exp $
! $Author: rosinski $
!
!-----------------------------------------------------------------------

    use pspect
    use rgrid
    use commap
    use physconst,    only: rair, gravit
    use constituents, only: pcnst, pnats, cnst_name, qmin, cnst_read_iv
    use chemistry,    only: chem_implements_cnst, chem_init_cnst
    ! TBH:  combine modules aerosol_intr and aerosols?  
    use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
    use test_tracers, only: test_tracers_implements_cnst, test_tracers_init_cnst
    use cldcond,      only: cldcond_implements_cnst, cldcond_init_cnst

    implicit none

    include 'netcdf.inc'

#include <comctl.h>
#include <comhyb.h>
#include <comqfl.h>
#include <comlun.h>
#include <perturb.h>

!---------------------------Local workspace-----------------------------
!
    integer i,j,k,m,lat,irow             ! grid and constituent indices
    integer ihem                         ! hemisphere index
    real(r8) phi(2,psp/2)                ! used in spectral truncation of phis
    real(r8) pdelb(plond,plev)           ! pressure diff between interfaces
!                                        ! using "B" part of hybrid grid only
    real(r8) hyad (plev)                 ! del (A)
    real(r8) pssum                       ! surface pressure sum
    real(r8) dotproda                    ! dot product
    real(r8) dotprodb                    ! dot product
    real(r8) pertval                     ! perturbation value
    real(r8) zgssum                      ! partial sums of phis
    real(r8) tmp1                        ! tmp space
    integer ii                           ! index
!
! Netcdf related variables
!
    integer lonsiz, latsiz, levsiz       ! Dimension sizes
    integer londimid, levdimid, latdimid ! Dimension ID's
    integer uid, vid, tid                ! Variable ID's
    integer tracid(pcnst+pnats)          ! Variable ID's
    integer phisid, sghid, psid          ! Variable ID's
    integer landmid
    integer pblhtid
    integer tpertid
    integer qpertid
    integer cldid
    integer qcwatid
    integer tcwatid
    integer lcwatid
#if ( ! defined COUP_CSM )
    integer ts1id, ts2id, ts3id, ts4id,tsiceid,tsice_rad_id   ! Variable ID's
    integer tbotid             ! Variable ID's
#endif
#if ( defined COUP_SOM )
      integer sicid
      integer icefracid
      integer tsocnid
#endif
    integer snowhiceid           ! Variable ID's#endif
    integer landfracid        ! Variable ID's

    integer strt2d(3)                    ! start lon, lat, time indices for netcdf 2-d
    integer strt3d(4)                    ! start lon, lev, lat, time for netcdf 3-d
    data strt2d/3*1/                     ! Only index 2 will ever change
    data strt3d/4*1/                     ! Only indices 2,3 will ever change

    integer cnt2d(3)                     ! lon, lat, time counts for netcdf 2-d
    integer cnt3d(4)                     ! lon, lat, lev, time counts for netcdf 2-d
    data cnt2d/plon,1,1/                 ! 2-d arrs: Always grab only a "plon" slice
    data cnt3d/plon,plev,plat,1/         ! 3-d arrs: Always grab a full time slice

    integer ndims2d                      ! number of dimensions
    integer dims2d(NF_MAX_VAR_DIMS)      ! variable shape
    integer ndims3d                      ! number of dimensions
    integer dims3d(NF_MAX_VAR_DIMS)      ! variable shape
    integer tmptype
    integer natt, ret, attlen            ! netcdf return values
    logical phis_hires                   ! true => PHIS came from hi res topo
    real(r8) arr3d(plon,plev,plat)
    character*(NF_MAX_NAME) tmpname
    character*256 text
    character*80 trunits                 ! tracer untis

!
!-----------------------------------------------------------------------
! Allocate memory for temporary arrays
!-----------------------------------------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter 
!
    allocate ( ps_tmp(plond,plat) )                   
    allocate ( u3_tmp(plond,plev,plat) )             
    allocate ( v3_tmp(plond,plev,plat) )             
    allocate ( t3_tmp(plond,plev,plat) )             
    allocate ( q3_tmp(plond,plev,pcnst+pnats,plat) ) 
    allocate ( tl_tmp(plond,plev,plat) )
    allocate ( tm_tmp(plond,plev,plat) )
    allocate ( ql_tmp(plond,plev,plat) )
    allocate ( qm_tmp(plond,plev,plat) )
    allocate ( phis_tmp(plond,plat) )        
    allocate ( phisl_tmp(plond,plat) )        
    allocate ( phism_tmp(plond,plat) )        
    allocate ( landm_tmp(plond,plat) )                  
    allocate ( sgh_tmp(plond,plat) )                 
    allocate ( tsice_tmp(plond,plat) )                   
    allocate ( tsice_rad_tmp(plond,plat) )                   
    allocate ( tbot_tmp    (plond,plat) )                
    allocate ( tssub_tmp(plond,plevmx,plat) )         
    allocate ( dpsl_tmp(plond,plat) )                
    allocate ( dpsm_tmp(plond,plat) )                 
    allocate ( div_tmp(plond,plev,plat) )             
    allocate ( sicthk_tmp(plond,plat) )               
    allocate ( snowhice_tmp(plond,plat) )                
    allocate ( landfrac_tmp(plond,plat) )                
    allocate ( pblht_tmp(plond,plat) )                  
    allocate ( tpert_tmp(plond,plat) )                  
    allocate ( qpert_tmp(plond,pcnst+pnats,plat) )                  
    allocate ( cld_tmp  (plond,plev,plat) )             
    allocate ( qcwat_tmp(plond,plev,plat) )             
    allocate ( tcwat_tmp(plond,plev,plat) )             
    allocate ( lcwat_tmp(plond,plev,plat) ) 
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
    qpert_tmp(:plon,:,:) = 0.

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
!
! Get variable id's 
! Check that all tracer units are in mass mixing ratios
!
       call wrap_inq_varid (ncid_ini, 'U'   , uid)
       call wrap_inq_varid (ncid_ini, 'V'   , vid)
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
             call wrap_get_att_text (NCID_INI,tracid(m),'units',trunits)
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
       call wrap_inq_var (ncid_ini, uid, tmpname, tmptype,ndims3d, dims3d, natt)
       if (dims3d(1).ne.londimid .or. dims3d(2).ne.levdimid .or. &
           dims3d(3).ne.latdimid .or. ndims3d.gt.4) then
              write(6,*)'INIDAT: Bad number of dims or ordering on 3d fld'
              call endrun
       end if
       call wrap_inq_var (ncid_ini, psid, tmpname, tmptype,ndims2d,dims2d ,natt)
       if (dims2d(1).ne.londimid .or. dims2d(2).ne.latdimid .or. ndims2d.gt.3) then
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
          write(6,*)'INIDAT: Will filter input PHIS: attribute from_hires is true'
       else
          phis_hires = .false.
          write(6,*)'INIDAT: Will not filter input PHIS: attribute ', &
                    'from_hires is either false or not present'
       end if
!
! Read in 2d fields.  
! For stand alone run: get surface temp and 4 (sub)surface temp fields
!
       do j=1,plat
          strt2d(2) = j
          if (ideal_phys .or. aqua_planet) then
             do i=1,nlon(j)
                phis_tmp(i,j) = 0.
                sgh_tmp (i,j) = 0.
             end do
          else
             call wrap_get_vara_realx (ncid_ini, phisid, strt2d, cnt2d, phis_tmp(1,j))
             call wrap_get_vara_realx (ncid_ini, sghid , strt2d, cnt2d, sgh_tmp (1,j))
          end if
          call wrap_get_vara_realx(ncid_ini, landmid, strt2d, cnt2d, landm_tmp(1,j  ))
          call wrap_get_vara_realx(ncid_ini, psid   , strt2d, cnt2d, ps_tmp   (1,j  ))
#if ( ! defined COUP_CSM )
          if (aqua_planet) then
             do i=1,nlon(j)
                landfrac_tmp(i,j) = 0.
             end do
          else
             call wrap_get_vara_realx (ncid_ini, landfracid, strt2d, cnt2d, landfrac_tmp(1,j))
             if (read_tbot) then
                call wrap_get_vara_realx (ncid_ini, tbotid    , strt2d, cnt2d, tbot_tmp    (1,j))
             end if
          endif
          call wrap_get_vara_realx (ncid_ini, tsiceid, strt2d, cnt2d, tsice_tmp(1,j))
          call wrap_get_vara_realx (ncid_ini, ts1id, strt2d, cnt2d, tssub_tmp (1,1,j))
          call wrap_get_vara_realx (ncid_ini, ts2id, strt2d, cnt2d, tssub_tmp (1,2,j))
          call wrap_get_vara_realx (ncid_ini, ts3id, strt2d, cnt2d, tssub_tmp (1,3,j))
          call wrap_get_vara_realx (ncid_ini, ts4id, strt2d, cnt2d, tssub_tmp (1,4,j))

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
             call wrap_get_vara_realx (ncid_ini, qpertid, strt2d, cnt2d, qpert_tmp(1,1,j))
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
! Copies are done instead of reading directly into 
! prognostic arrays for netcdf performance.
! Constituents not read from initial file are initialized by the package
! that implements them.
!
       call wrap_get_vara_realx(ncid_ini, uid, strt3d, cnt3d, arr3d)
       u3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

       call wrap_get_vara_realx(ncid_ini, vid, strt3d, cnt3d, arr3d)
       v3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

       call wrap_get_vara_realx(ncid_ini, tid, strt3d, cnt3d, arr3d)
       t3_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)

       do m=1,pcnst+pnats
          if (cnst_read_iv(m) .and. &
             nf_inq_varid(ncid_ini, cnst_name(m), tracid(m) ) == NF_NOERR ) then
             call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arr3d)
             q3_tmp(:plon,:plev,m,:plat) = arr3d
          else
             write(6,*) 'Warning:  Not reading ',cnst_name(m), ' from IC file.'
             arr3d = 0.
             if (cldcond_implements_cnst(cnst_name(m))) then
                call cldcond_init_cnst(cnst_name(m), arr3d)
                write(6,*) '          ', cnst_name(m), ' initialized by "cldcond_init_cnst"'
             else if (chem_implements_cnst(cnst_name(m))) then
                call chem_init_cnst(cnst_name(m), arr3d)
                write(6,*) '          ', cnst_name(m), ' initialized by "chem_init_cnst"'
             else if (test_tracers_implements_cnst(cnst_name(m))) then
                call test_tracers_init_cnst(cnst_name(m), arr3d)
                write(6,*) '          ', cnst_name(m), ' initialized by "test_tracers_init_cnst"'
             else if (aerosol_implements_cnst(cnst_name(m))) then
                call aerosol_init_cnst(cnst_name(m), arr3d)
                write(6,*) '          ', cnst_name(m), ' initialized by "aerosol_init_cnst"'
             else
                write(6,*) '          ', cnst_name(m), ' set to 0.'
             end if
             q3_tmp(:plon,:plev,m,:plat) = arr3d
          endif
       end do
       do lat=1,plat
          call qneg3('INIDAT  ',lat   ,nlon(lat),plond   ,plev    , &
                     pcnst+pnats,qmin ,q3_tmp(1,1,1,lat))
       end do

!
       if (read_cloud) then
          call wrap_get_vara_realx(ncid_ini, cldid  , strt3d, cnt3d, arr3d)
          cld_tmp  (:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)
       end if

       if (read_qcwat) then
          call wrap_get_vara_realx(ncid_ini, qcwatid, strt3d, cnt3d, arr3d)
          qcwat_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)
       end if

       if (read_tcwat) then
          call wrap_get_vara_realx(ncid_ini, tcwatid, strt3d, cnt3d, arr3d)
          tcwat_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)
       end if

       if (read_lcwat) then
          call wrap_get_vara_realx(ncid_ini, lcwatid, strt3d, cnt3d, arr3d)
          lcwat_tmp(:plon,:plev,:plat) = arr3d(:plon,:plev,:plat)
       end if
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
                   t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                end do
             end do
          end do
       endif
!
!-----------------------------------------------------------------------
! Spectrally truncate ps and its derivatives (dpsl and dpsm), phis, 
! u, v, t, divergence (div).
!-----------------------------------------------------------------------
!
       arr3d = q3_tmp(:plon,:plev,1,:plat) ! save q
       call spetru (ps_tmp   ,phis_tmp  ,u3_tmp  ,v3_tmp  ,t3_tmp   , &
                    q3_tmp   ,div_tmp   ,dpsl_tmp,dpsm_tmp,tl_tmp   , &
                    tm_tmp   ,ql_tmp    ,qm_tmp  ,phi     ,phisl_tmp, &
                    phism_tmp,phis_hires)
!
! For sld do not use the spectrally truncated q3
       q3_tmp(:plon,:plev,1,:plat) = arr3d
!
! Compute ln(Ps*) (Ritchie & Tanguay, 1995) in spectral space
!
       tmp1 = 1./(rair*t0(plev))
       do ii = 1,psp/2
          i = 2*ii - 1
          lnpstar(i  ) = -phi(1,ii)*tmp1
          lnpstar(i+1) = -phi(2,ii)*tmp1
       end do
!
!-----------------------------------------------------------------------
! Integrals of mass, moisture and geopotential height
!-----------------------------------------------------------------------
!
! Compute pdel from "A" portion of hybrid vertical grid
!
       do k=1,plev
          hyad(k) = hyai(k+1) - hyai(k)
       end do
       do k=1,plev
          do i=1,plon
             pdela(i,k) = hyad(k)*ps0
          end do
       end do
!        
! Initialize mass and moisture integrals for summation
! in a third calculation loop (assures bit-for-bit compare
! with non-random history tape).
!
       tmassf_tmp = 0.
       qmass1_tmp = 0.
       qmass2_tmp = 0.
       zgsint_tmp = 0.
!
! Compute integrals of mass, moisture, and geopotential height
!
       do irow = 1,plat/2
          do ihem=1,2
             if (ihem.eq.1) then
                lat = irow
             else
                lat = plat - irow + 1
             end if
!              
! Accumulate average mass of atmosphere
!
             call pdelb0 (ps_tmp(1,lat),pdelb   ,nlon(lat))
             pssum  = 0.
             zgssum = 0.
             do i=1,nlon(lat)
                pssum  = pssum  + ps_tmp  (i,lat)
                zgssum = zgssum + phis_tmp(i,lat)
             end do
             tmassf_tmp = tmassf_tmp + w(irow)*pssum/nlon(lat)
             zgsint_tmp = zgsint_tmp + w(irow)*zgssum/nlon(lat)
!
! Calculate global integrals needed for water vapor adjustment
!
             do k=1,plev
                dotproda = 0.
                dotprodb = 0.
                do i=1,nlon(lat)
                   dotproda = dotproda + q3_tmp(i,k,1,lat)*pdela(i,k)
                   dotprodb = dotprodb + q3_tmp(i,k,1,lat)*pdelb(i,k)
                end do
                qmass1_tmp = qmass1_tmp + w(irow)*dotproda/nlon(lat)
                qmass2_tmp = qmass2_tmp + w(irow)*dotprodb/nlon(lat)
             end do
          end do
       end do                  ! end of latitude loop
!
! Normalize average mass, height
!
       tmassf_tmp = tmassf_tmp*.5/gravit
       qmass1_tmp = qmass1_tmp*.5/gravit
       qmass2_tmp = qmass2_tmp*.5/gravit
       zgsint_tmp = zgsint_tmp*.5/gravit
       qmassf_tmp = qmass1_tmp + qmass2_tmp
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
!
       tmass0 = 98222./gravit
       if (ideal_phys) tmass0 = 100000./gravit
       write(6,800) tmassf_tmp,tmass0,qmassf_tmp
       write(6,810) zgsint_tmp
800    format('INIDAT: MASS OF INITIAL DATA BEFORE CORRECTION = ' &
              ,1p,e20.10,/,' DRY MASS WILL BE HELD = ',e20.10,/, &
              ' MASS OF MOISTURE AFTER REMOVAL OF NEGATIVES = ',e20.10) 
810    format(/69('*')/'INIDAT: Globally averaged geopotential ', &
              'height = ',f16.10,' meters'/69('*')/)
!
! Compute and apply an initial mass fix factor which preserves horizontal
! gradients of ln(ps).
!
     if (adiabatic .or. ideal_phys) then
        fixmas = tmass0/tmassf_tmp
     else
        fixmas = (tmass0 + qmass1_tmp)/(tmassf_tmp - qmass2_tmp)
     end if
     do lat=1,plat
        do i=1,nlon(lat)
           ps_tmp(i,lat) = ps_tmp(i,lat)*fixmas
        end do
     end do
    endif                     ! end of if-masterproc
!
!-----------------------------------------------------------------------
! Copy temporary arrays to model arrays
!-----------------------------------------------------------------------
!
    call copy_inidat


!-----------------------------------------------------------------------
! Deallocate memory for temporary arrays
!-----------------------------------------------------------------------
!
    deallocate ( ps_tmp )
    deallocate ( u3_tmp )
    deallocate ( v3_tmp )
    deallocate ( t3_tmp )
    deallocate ( q3_tmp )
    deallocate ( tl_tmp )
    deallocate ( tm_tmp )
    deallocate ( ql_tmp )
    deallocate ( qm_tmp )
    deallocate ( phis_tmp )        
    deallocate ( phisl_tmp )        
    deallocate ( phism_tmp )        
    deallocate ( landm_tmp )
    deallocate ( sgh_tmp )
    deallocate ( tsice_tmp )
    deallocate ( tsice_rad_tmp )
    deallocate ( tbot_tmp )
    deallocate ( tssub_tmp )
    deallocate ( dpsl_tmp )
    deallocate ( dpsm_tmp )
    deallocate ( div_tmp )
    deallocate ( sicthk_tmp )
    deallocate ( snowhice_tmp )
    deallocate ( pblht_tmp )
    deallocate ( tpert_tmp )
    deallocate ( qpert_tmp )
    deallocate ( cld_tmp )
    deallocate ( qcwat_tmp )
    deallocate ( tcwat_tmp )
    deallocate ( lcwat_tmp )
    deallocate ( landfrac_tmp )
#if ( defined COUP_SOM )
      deallocate ( icefrac_tmp )
      deallocate ( tsocn_tmp )
#endif
!
    return
  end subroutine read_inidat

!*********************************************************************C

  subroutine copy_inidat
!-----------------------------------------------------------------------
!
! Purpose:
! Copy temporary arrays to model arrays 
! note that the use statements below contain the definitions
! of the model arrays
!
!-----------------------------------------------------------------------

    use prognostics
    use buffer
    use phys_grid
    use constituents, only: cnst_get_ind
    use phys_buffer, only: pbuf, pbuf_times, pbuf_get_fld_idx
#if ( defined SPMD )
    use mpishorthand
    use spmd_dyn, only: npes, compute_gsfactors
#endif

    implicit none

#include <comqfl.h>
!
!---------------------------Local workspace-----------------------------
!
    real(r8), allocatable :: tmpchunk3d(:,:,:)             
    real(r8), allocatable :: tmpchunk(:,:)             
    integer, parameter :: iend = i1+plon-1 ! last "real model" i-index in extended grid
    integer, parameter :: jend = j1+plat-1 ! last "real model" j-index in extended grid
    integer begj, endj
    integer n,i,j, m, c
    integer :: ixcldice, ixcldliq
    real(r8), pointer, dimension(:,:,:,:) :: qcwat, lcwat, tcwat, cld

#if ( defined SPMD )
    integer :: numperlat         ! number of values per latitude band
    integer :: numsend(0:npes-1) ! number of items to be sent
    integer :: numrecv           ! number of items to be received
    integer :: displs(0:npes-1)  ! displacement array
#endif
!
!-----------------------------------------------------------------------
!
#ifdef HADVTEST
!
!JR Overwrite fields for flat-earth solid-body rotation
!
    call hadvtest_init
#endif

    begj = beglatex + numbnd
    endj = begj + numlats - 1

!PW Dynamics fields
#if ( defined SPMD )
    numperlat = plond
    call compute_gsfactors (numperlat, numrecv, numsend, displs)

    call mpiscatterv (ps_tmp    ,numsend, displs, mpir8,ps    (1,beglat,1) ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (phis_tmp  ,numsend, displs, mpir8,phis  (1,beglat)   ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (phisl_tmp ,numsend, displs, mpir8,phisl (1,beglat)   ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (phism_tmp ,numsend, displs, mpir8,phism (1,beglat)   ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (dpsl_tmp  ,numsend, displs, mpir8,dpsl  (1,beglat)   ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (dpsm_tmp  ,numsend, displs, mpir8,dpsm  (1,beglat)   ,numrecv, mpir8,0,mpicom)

    numperlat = plndlv
    call compute_gsfactors (numperlat, numrecv, numsend, displs)

    call mpiscatterv (u3_tmp    ,numsend, displs, mpir8,u3    (i1,1,begj,1),numrecv, mpir8,0,mpicom)
    call mpiscatterv (v3_tmp    ,numsend, displs, mpir8,v3    (i1,1,begj,1),numrecv, mpir8,0,mpicom)
    call mpiscatterv (t3_tmp    ,numsend, displs, mpir8,t3    (i1,1,begj,1),numrecv, mpir8,0,mpicom)

    call mpiscatterv (div_tmp   ,numsend, displs, mpir8,div   (1,1,beglat,1) ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (tl_tmp    ,numsend, displs, mpir8,tl    (1,1,beglat) ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (tm_tmp    ,numsend, displs, mpir8,tm    (1,1,beglat) ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (ql_tmp    ,numsend, displs, mpir8,ql    (1,1,beglat) ,numrecv, mpir8,0,mpicom)
    call mpiscatterv (qm_tmp    ,numsend, displs, mpir8,qm    (1,1,beglat) ,numrecv, mpir8,0,mpicom)

    numperlat = plndlv*(pcnst+pnats)
    call compute_gsfactors (numperlat, numrecv, numsend, displs)

    call mpiscatterv (q3_tmp, numsend, displs, mpir8, q3(i1,1,1,begj,1), numrecv, mpir8, 0, mpicom)

    call mpibcast   (lnpstar   ,psp,mpir8,0 , mpicom)

#else

    ps    (:,:,1) = ps_tmp    (:,:)
    phis  (:,:)   = phis_tmp  (:,:)
    phisl (:,:)   = phisl_tmp (:,:)
    phism (:,:)   = phism_tmp (:,:)
    dpsl  (:,:)   = dpsl_tmp  (:,:)
    dpsm  (:,:)   = dpsm_tmp  (:,:)

    u3    (i1:iend,:,j1:jend,1) = u3_tmp(:plon,:plev,:plat)
    v3    (i1:iend,:,j1:jend,1) = v3_tmp(:plon,:plev,:plat)
    t3    (i1:iend,:,j1:jend,1) = t3_tmp(:plon,:plev,:plat)

    div   (:,:,:,1) = div_tmp   (:,:,:)
    tl    (:,:,:) = tl_tmp    (:,:,:)
    tm    (:,:,:) = tm_tmp    (:,:,:)
    ql    (:,:,:) = ql_tmp    (:,:,:)
    qm    (:,:,:) = qm_tmp    (:,:,:)

    q3    (i1:iend,:plev,:pcnst+pnats,j1:jend,1) = q3_tmp(:plon,:plev,:pcnst+pnats,:plat)
#endif
    dpsmm1(:,:)   = dpsm (:,:)
    dpsmp1(:,:)   = dpsm (:,:)
    dpslm1(:,:)   = dpsl (:,:)
    dpslp1(:,:)   = dpsl (:,:)
    tlm1  (:,:,:) = tl   (:,:,:)
    tmm1  (:,:,:) = tm   (:,:,:)
    ed1   (:,:,:) = 0.

!PW Physics fields
    allocate ( tmpchunk(pcols,begchunk:endchunk) )
    allocate ( tmpchunk3d(pcols,plevmx,begchunk:endchunk) )

    call scatter_field_to_chunk(1,1,1,plond,landfrac_tmp,landfrac(1,begchunk))
    call scatter_field_to_chunk(1,1,1,plond,landm_tmp,landm(1,begchunk))
    call scatter_field_to_chunk(1,1,1,plond,sgh_tmp,sgh(1,begchunk))
    call scatter_field_to_chunk(1,1,1,plond,tsice_tmp,tsice(1,begchunk))
    call scatter_field_to_chunk(1,1,1,plond,snowhice_tmp,snowhice(1,begchunk))

    call scatter_field_to_chunk(1,plevmx,1,plond,tssub_tmp,tmpchunk3d)
    do i =begchunk,endchunk
       surface_state2d(i)%tssub(:,:) = tmpchunk3d(:,:,i)
    end do

#if ( defined COUP_SOM )
      call scatter_field_to_chunk(1,1,1,plond,sicthk_tmp,sicthk(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,icefrac_tmp,icefrac(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,tsocn_tmp,tsocn(1,begchunk))

! define an initial ocean fraction and non-land ice fraction
! The 1st "where" stmt used to be done in update_srf_fractions (dev45)

      do c=begchunk,endchunk
         where (icefrac(:pcols,c) + landfrac(:pcols,c) > 1.0)
            icefrac(:pcols,c) = 1. - landfrac(:pcols,c)
         end where

         where (landfrac(:pcols,c) < 1.)
            aice(:pcols,c) = icefrac(:pcols,c)/(1. - landfrac(:pcols,c))
         elsewhere
            aice(:pcols,c) = 0.
         end where
         ocnfrac(:pcols,c) = 1. - landfrac(:pcols,c) - icefrac(:pcols,c)
      end do

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
!
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
            qpert_tmp  (:plon,1,:) = 0.
            write(6,*) 'Warning:  QPERT not found on IC file; initialized to 0.'
         end if
         if(.not. read_cloud) then
            cld_tmp    (:plon,:,:) = 0.
            write(6,*) 'Warning:  CLOUD not found on IC file; initialized to 0.'
         end if
         if(.not. read_qcwat) then
            qcwat_tmp  (:plon,:,:) = q3_tmp(:plon,:,1,:)
            write(6,*) 'Warning:  QCWAT not found on IC file; initialized with Q'
         end if
         if(.not. read_tcwat) then
            tcwat_tmp  (:plon,:,:) = t3_tmp(:plon,:,:)
            write(6,*) 'Warning:  TCWAT not found on IC file; initialized with T'
         end if
         if(.not. read_lcwat) then
            call cnst_get_ind('CLDICE', ixcldice)
            call cnst_get_ind('CLDLIQ', ixcldliq)
            lcwat_tmp  (:plon,:,:) = q3_tmp(:plon,:,ixcldice,:) + q3_tmp(:plon,:,ixcldliq,:)
            write(6,*) 'Warning:  LCWAT not found on IC file; initialized with CLDICE + CLDLIQ'
         end if
         if(.not. read_tbot) then
            tbot_tmp   (:plon,:) = t3_tmp(:plon,plev,:)
            write(6,*) 'Warning:  TBOT not found on IC file; initialized with lowest level of T'
         end if
         if(.not. read_tsicerad) then
            tsice_rad_tmp(:plon,:) = tsice_tmp(:plon,:)
            write(6,*) 'Warning:  TSICERAD not found on IC file; initialized with TSICE'
         end if
      endif

      call scatter_field_to_chunk(1,1,1,plond,tbot_tmp,tmpchunk)
      do i =begchunk,endchunk
         surface_state2d(i)%tbot(:) = tmpchunk(:,i)
      end do
      call scatter_field_to_chunk(1,          1,1,plond,tsice_rad_tmp,tsice_rad(1,begchunk))
      call scatter_field_to_chunk(1,          1,1,plond,pblht_tmp,pblht(1  ,begchunk  ))
      call scatter_field_to_chunk(1,          1,1,plond,tpert_tmp,tpert(1  ,begchunk  ))
      call scatter_field_to_chunk(1,pcnst+pnats,1,plond,qpert_tmp,qpert(1,1,begchunk  ))
      m = pbuf_get_fld_idx('QCWAT')
      qcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,qcwat_tmp,qcwat(:,:,:,1))
      m = pbuf_get_fld_idx('LCWAT')
      lcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,lcwat_tmp,lcwat(:,:,:,1))
      m = pbuf_get_fld_idx('TCWAT')
      tcwat => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,tcwat_tmp,tcwat(:,:,:,1))
      m = pbuf_get_fld_idx('CLD')
      cld => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
      call scatter_field_to_chunk(1,plev,1,plond,cld_tmp,cld(:,:,:,1))
!
    if (pbuf_times > 1) then
       do n = 2, pbuf_times
          cld  (:,:,:,n) = cld  (:,:,:,1)
          qcwat(:,:,:,n) = qcwat(:,:,:,1)
          lcwat(:,:,:,n) = lcwat(:,:,:,1)
          tcwat(:,:,:,n) = tcwat(:,:,:,1)
       end do
    end if
!
! Global integerals
!
    if (masterproc) then
       tmassf = tmassf_tmp
       qmass1 = qmass1_tmp
       qmass2 = qmass2_tmp
       qmassf = qmassf_tmp
       zgsint = zgsint_tmp
    endif

#if ( defined SPMD )
    call mpibcast (tmass0,1,mpir8,0,mpicom)
    call mpibcast (tmassf,1,mpir8,0,mpicom)
    call mpibcast (qmass1,1,mpir8,0,mpicom)
    call mpibcast (qmass2,1,mpir8,0,mpicom)
    call mpibcast (qmassf,1,mpir8,0,mpicom)
    call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif
    deallocate ( tmpchunk )
    deallocate ( tmpchunk3d)

  end subroutine copy_inidat

#ifdef HADVTEST
  subroutine hadvtest_init
    use pmgrid
    use rgrid
    use physconst, only:
    use commap

    implicit none

#include <comhyb.h>
#include <hadvtest.h>
!
!---------------------------Local workspace-----------------------------
!
    integer i   !
    integer k   ! - indices
    integer lat !

    real(r8) h0, u0, small_r, big_r, theta, theta_c, lambda, lambda_c
    real(r8) alfa, dlam, pie
    real(r8) pie
!
!-----------------------------------------------------------------------
!
! First: zero sgh and phis fields
!
    sgh_tmp(:,:) = 0.
    phis_tmp(:,:) = 0.
!
!JR Analytic IC and wind
!
    pie = acos(-1.)
!
!jr Define wind and constituent fields
!
    h0 = 1000.
    u0 = 2.*pie*rearth/(12.*86400.)
    big_r = rearth/3.
    theta_c = +60.*pie/180.   ! 60 deg north
    theta_c = -60.*pie/180.   ! 60 deg south
    theta_c = 0.              ! equator
    lambda_c = 0.             ! Greenwich
    lambda_c = 3.*pie/2.

    do lat=1,plat
       theta = clat(lat)
       do k=1,plev
          alfa = 0.
          if (k.eq.1) then
             alfa = 0.
          else if (k.eq.2) then
             alfa = 0.05
          else if (k.eq.plev-1) then
             alfa = 0.5*pie - 0.05
          else if (k.eq.plev) then
             alfa = 0.5*pie      ! blows north
          else
             alfa = (k-2)*pie/(2.*(plev-3))
          end if

          do i=1,nlon(lat)
             lambda = 2.*pie*(i-1)/nlon(lat)
!
!jr Use these settings in conjunction with theta_c to start the blob at 
!jr Greenwich
!
             usave(i,k,lat)    = u0*(cos(theta)*cos(alfa) + &
                                 sin(theta)*cos(lambda-0.5*pie)*sin(alfa))
             vsave(i,k,lat)    = -u0*sin(lambda-0.5*pie)*sin(alfa)
!
!jr Use these settings in conjunction with theta_c to start the blob at 270.
!
             usave(i,k,lat)    = u0*(cos(theta)*cos(alfa) + &
                                 sin(theta)*cos(lambda)*sin(alfa))
             vsave(i,k,lat)    = -u0*sin(lambda)*sin(alfa)
             u3_tmp(i,k,lat)   = usave(i,k,lat)
             v3_tmp(i,k,lat)   = vsave(i,k,lat)
             dlam              = lambda - lambda_c
             small_r           = rearth*acos(sin(theta_c)*sin(theta) + &
                                 cos(theta_c)*cos(theta)*cos(dlam))
             q3_tmp(i,k,1,lat) = 0.
             if (small_r .lt. big_r) then
                q3_tmp(i,k,1,lat) = h0/2.*(1. + cos(pie*small_r/big_r))
             end if
!
!jr Stick Q into T to test spectral advection (of what's in T)
!jr Or put 300 in T.
!
             t3_tmp(i,k,lat) = 300.
             t3_tmp(i,k,lat) = q3_tmp(i,k,1,lat)
          end do
       end do
!
!jr Save surface pressure for future timesteps.  Set to 1.e5 everywhere
!
       do i=1,nlon(lat)
          ps_tmp(i,lat) = ps0
          pssave(i,lat) = ps_tmp(i,lat)
       end do
    end do

    return
  end subroutine hadvtest_init
#endif

end module inidat
