#include <misc.h>
#include <params.h>

module dust_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use shr_kind_mod,only: r8 => shr_kind_r8
  use pmgrid,      only: plon, plat, masterproc
  use ppgrid,      only: pcols, pver,pverp
  use physconst,   only: mwdry, mwh2o,gravit,rair
  use constituents,only: ppcnst, cnst_add, cnst_name, advected, cnst_get_ind
  use dust,        only: dust_number,ndst,sz_nbr,dst_src_nbr
    

  implicit none

  private          ! Make default type private to the module

  save

  integer ncyear
  integer ix_dust, nx_dust

  integer, parameter :: ncnst=ndst                   ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'DSTX01', 'DSTX02', 'DSTX03', 'DSTX04'/)

!
! Public interfaces
!
  public dust_register_cnst                        ! register consituents
  public dust_implements_cnst                      ! returns true if consituent is implemented by this package
  public dust_init_cnst                            ! initialize mixing ratios if not read from initial file
  public dust_initialize                           ! initialize (history) variables
  public dust_wet_intr                             ! interface to wet deposition
  public dust_emis_intr                            ! interface to emission
  public dust_drydep_intr                          ! interface to tendency computation
  public dust_time_interp                          ! interpolate oxidants and fluxes to current time
  public dust_idx1                                 ! allow other parts of code to know where dust is

  real(r8) stk_crc(ndst) ![frc] Correction to Stokes settling velocity
  real(r8) dns_aer       ![kg m-3] Aerosol density
  real(r8) tmp1          !Factor in saltation computation (named as in Charlie's code)
  real(r8) ovr_src_snk_mss(dst_src_nbr,ndst)  
  real(r8) dmt_vwr(ndst) ![m] Mass-weighted mean diameter resolved
  real(r8)::  soil_erodibility(plon,plat)            ! soil erodibility factor

contains

!===============================================================================
  subroutine dust_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents for all aerosols
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    real(r8), parameter :: one  = 1._r8
    real(r8), parameter :: zero  = 0._r8
!-----------------------------------------------------------------------

! Set names of variables undergoing evolution
! returns m as current index for tracer numbers
    call cnst_add(cnst_names(1), advected, one, one, zero, m) 

! and store the start index of dust species used elsewhere in model retrieved by dust_idx1
    call set_dust_idx(m)

    call cnst_add(cnst_names(2), advected, one, one, zero, m)
    call cnst_add(cnst_names(3), advected, one, one, zero, m)
    call cnst_add(cnst_names(4), advected, one, one, zero, m)

    return
  end subroutine dust_register_cnst



!=======================================================================
  function dust_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: dust_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     dust_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           dust_implements_cnst = .true.
           return
        end if
     end do
  end function dust_implements_cnst


!=======================================================================
  subroutine dust_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial mass mixing ratios.
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------
    
    character(len=*), intent(in) :: name         ! constituent name
    
    real(r8), intent(out) :: q(:,:,:)            !  mass mixing ratio
!-----------------------------------------------------------------------
    
    if ( name == cnst_names(1) ) then
       q = 0._r8
    else if ( name == cnst_names(2) ) then
       q = 0._r8
    else if ( name == cnst_names(3) ) then
       q = 0._r8
    else if ( name == cnst_names(4) ) then
       q = 0._r8
    end if

  end subroutine dust_init_cnst



  function dust_idx1()
    implicit none
    integer dust_idx1
    dust_idx1 = ix_dust
  end function dust_idx1

  subroutine set_dust_idx(m)
    implicit none
    integer m
    
    ix_dust = m
    write (6,*) 'set_dust_idx: start index for dust is ', ix_dust
  end subroutine set_dust_idx

!===============================================================================
  subroutine dust_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of dust chemistry
!          (declare history variables)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    use history,    only: addfld, add_default, phys_decomp
!    use dust_atm,       only: dust_atm_ini
!    use dustbnd,    only: dustbndini
!    use dstszdst
!    use dsttvbds

    use shr_const_mod,    only: SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH

    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual

    use commap, only: clat, clon, w
    use ppgrid, only: pcols, begchunk, endchunk
!added by nmm
    use error_messages, only: alloc_err,handle_ncerr
!nf90     use netcdf
     use pmgrid
    use ioFileMod, only: getfil
    use filenames, only: soil_erod

    implicit none
#include <netcdf.inc>
!---------------------------Local workspace-----------------------------
    integer :: m
    integer :: ix             ! start index for dust species
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    character lpthncsf*200       ! local file for netCDF surface flux data
    character dummy*20
    character(len=256) :: soil_erod_file

! changes to read in soil erodibility factor
     real(r8) :: dum     !dummy variable for erf test
!     Parameters
     integer :: i,istat,did,nlon,vid,recid,nrec
    integer ::&
       ntpy,             &! number of time samples per year in emissions dataset
       nyr,              &! number of years in emissions dataset
       ncid,             &! ID for netCDF file
       loyri,            &! index in yr array for lower bound year
       start(4),         &! start vector for netCDF hyperslabs
       count(4)           ! count vector for netCDF hyperslabs

!    real(r8) :: gw(plat)
!-----------------------------------------------------------------------


!    call inisflx( SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH,  &
!         clon, clat, gw)

    calday = get_curr_calday()
    write (6,*) ' dust_initialize: positioning at calday ', calday
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    ncyear = yr
    if ( ncyear .lt. 1000 ) ncyear = ncyear + 1900

! use Sam's dust initialize subroutine:  call equivalent here:

    call Dustini()

! for soil erodibility in mobilization, apply inside CAM instead of lsm.
! read in soil erodibility factors, similar to Zender's boundary conditions

      ! Get file name.  
      call getfil(soil_erod, soil_erod_file, 0)

!nf90      call handle_ncerr( nf90_open( soil_erod_file, NF_NOWRITE, ncid ), &
!nf90         'dustbndini: error opening file '//trim(soil_erod_file) )
      call handle_ncerr( nf_open( soil_erod_file, NF_NOWRITE, ncid ), &
         'dustbndini: error opening file '//trim(soil_erod_file) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'dustbndini: no record variables ' )
      call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
         'dustbndini: no record variables ' )

      ! Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'dustbndini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'dustbndini: ' )
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), 'dustbndini: ' )
      call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), 'dustbndini: ' )
      if ( nlon .ne. plon ) then
         write(6,*)'dustbndini: number of longitudes (', nlon, ')', &
                   ' doesn''t match model resolution.'
         call endrun
      end if
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'mbl_bsn_fct_geo', vid ), &
!nf90         'dustbndini: cannot find variable '//'mbl_bsn_fct_geo' )
      call handle_ncerr( nf_inq_varid( ncid, 'mbl_bsn_fct_geo', vid ), &
         'dustbndini: cannot find variable '//'mbl_bsn_fct_geo' )
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 1
      count(4) = 1

!nf90     call handle_ncerr( nf90_get_var( ncid, vid,soil_erodibility, start, count ), &
!nf90         'dustbndini: cannot read data for '//'mbl_bsn_fct_geo' )
     call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soil_erodibility ), &
         'dustbndini: cannot read data for '//'mbl_bsn_fct_geo' )


!     Initialize time-invariant physical constants
!      call dst_cst_cmn_ini(rair,gravit,rearth)
    write (6,*) ' call to set time invariant fields here'
!     Initialize size grid
!      call dst_szdst_ini()
    write (6,*) ' call to set size grid here'
!     Initialize miscellaneous common blocks
!      call dst_msc_cmn_ini()
    write (6,*) ' call to set common blocks here'
!     Initialize time-invariant boundary data from netCDF file
!     call dst_tibds_ini('dst_bnd.nc')
!     NB: dst_tibds_ini() opens and closes the file, while
!     dst_tvbds_ini() opens the file and leaves it open
!     Initialize time-varying seasonal cycle data from netCDF file
!     call dst_tvbds_ini('dst_bnd.nc',calday)
       ix=dust_idx1()

! Set names of variable tendencies and declare them as history variables
       dummy = 'LND_MBL'
       call addfld (dummy,'frac ',1, 'A','Soil erodibility factor',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = 'RAM1'
       call addfld (dummy,'frac ',1, 'A','RAM1',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = 'airFV'
       call addfld (dummy,'frac ',1, 'A','FV',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = 'ORO'
       call addfld (dummy,'frac ',1, 'A','ORO',phys_decomp)
       call add_default (dummy, 1, ' ')
       do m = 1, 4
       dummy = trim(cnst_name(ix-1+m)) // 'SF'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' dust surface emission',phys_decomp)
       call add_default (dummy, 1, ' ')
       dummy = trim(cnst_name(ix-1+m)) // 'OD'
       call addfld (dummy,'Tau ',1, 'A',trim(cnst_name(ix-1+m))//' optical depth',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'TB'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' turbulent dry deposition flux',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'GV'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' gravitational dry deposition flux',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'DD'
       call addfld (dummy,'kg/m2/s ',1, 'A',trim(cnst_name(ix-1+m))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'DT'
       call addfld (dummy,'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' dry deposition',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'PP'
       call addfld (dummy,'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' wet deposition',phys_decomp)
       call add_default (dummy, 1, ' ')
        dummy = trim(cnst_name(ix-1+m)) // 'DV'
       call addfld (dummy,'m/s ',pver, 'A',trim(cnst_name(ix-1+m))//' deposition velocity',phys_decomp)
       call add_default (dummy, 1, ' ')
   end do


  end subroutine dust_initialize


!===============================================================================
  subroutine dust_wet_intr (state, ptend, cflx, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to wet processing of aerosols (source and sinks).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use wetdep, only: wetdepa

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: nstep
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                    ! latitude 
    real(r8), intent(in) :: cme(pcols,pver)            ! local condensation of cloud water
    real(r8), intent(in) :: prain(pcols,pver)            ! production of rain
    real(r8), intent(in) :: evapr(pcols,pver)            ! evaporation of rain
    real(r8), intent(in) :: cldn(pcols,pver)            ! cloud fraction
    real(r8), intent(in) :: cldc(pcols,pver)            ! convective cloud fraction
    real(r8), intent(in) :: cldv(pcols,pver)            ! cloudy volume undergoing scavenging

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    real(r8), intent(inout)  :: cflx(pcols,ppcnst)       ! Surface constituent flux (kg/m^2/s)

    real(r8), intent(inout) :: fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble

!
! Local variables
!
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tstate(pcols, pver, 4)            ! temporary state vector
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio
    real(r8) :: obuf(1)
    real(r8) :: calday        ! current calendar day
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: scavt(pcols, pver)
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: mm,i,k
    integer :: m                                  ! tracer index
    integer :: ixcldliq
    integer :: ixcldice
    real(r8) totcond(pcols, pver) ! total condensate
    real(r8) :: sol_fact

!-----------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    totcond = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

    do m = 1, ndst

       mm = ix_dust + m - 1
       scavt=0.
!       write(*,*) 'wet dep removed for debugging'
       ptend%lq(mm) = .TRUE.
       sol_fact = 0.3_r8
       call wetdepa( lat, state%t, state%pmid, state%q, state%pdel,  &
            cldn, cldc, cmfdqr, conicw, prain, cme,                     &
            evapr, totcond, state%q(:,:,mm), dt,            &
            scavt, iscavt, cldv, fracis(:,:,mm), sol_fact, ncol )
       ptend%q(:,:,mm)=scavt

       call outfld( trim(cnst_name(mm))//'PP', ptend%q(:,:,mm), pcols, lchnk)
!      write (6,*) ' range of ptend ix ', minval(ptend%q(:ncol,:,mm)),maxval(ptend%q(:ncol,:,mm))

    end do

!   write (6,*) ' dust_wet_intr: pcols, pcols ', pcols, pcols

    return

  end subroutine dust_wet_intr

  subroutine dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac,fv,ram1)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to dry deposition and sedimentation of dust
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: Natalie Mahowald and Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lat_all_p
    use constituents,  only: cnst_name
    use drydep_mod, only: drydepdr, setdvel, ddflux
    use dust_sediment_mod, only: dust_sediment_tend, dust_sediment_vel
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                 ! latitude 
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel--used over oceans and sea ice.
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: hflx(pcols)                  ! sensible heat flux
    real(r8), intent(in) :: prect(pcols)                     ! prect
    real(r8), intent(in) :: snowh(pcols)                     ! snow depth
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    integer, intent(in)  :: month
    real(r8), intent(in) :: wvflx(pcols)       ! water vapor flux
    real(r8), intent(in) :: fv(pcols)        ! for dry dep velocities from land model for dust
    real(r8), intent(in) :: ram1(pcols)       ! for dry dep velocities from land model for dust

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
!
! Local variables
!
    integer :: m                                  ! tracer index
    integer :: mm                                  ! tracer index
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tvs(pcols,pver)
    real(r8) :: dvel(pcols)            ! deposition velocity
    real(r8) :: sflx(pcols)            ! deposition flux
    real(r8) :: vlc_dry(pcols,pver,ndst)            ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,ndst)            ! dep velocity
    real(r8)::  vlc_trb(pcols,ndst)            ! dep velocity
    real(r8)::  dep_trb(pcols)       !kg/m2/s
    real(r8)::  dep_dry(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8)::  dep_dry_tend(pcols,pver)       !kg/kg/s (total of grav and trb)
     real(r8) :: obuf(1)
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real(r8)  pvdust(pcols,pverp)    ! sedimentation velocity in Pa
    integer :: i,k
     real(r8) :: oro(pcols)
!
!-----------------------------------------------------------------------
    ix = dust_idx1()
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol
       tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
       rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
! calculate oro--need to run match

       do i=1,ncol
          oro(i)=1.
          if(icefrac(i)>0.5) oro(i)=2.
          if(ocnfrac(i)>0.5) oro(i)=0.
       enddo
       call outfld( 'ORO', oro, pcols, lchnk )
       
!   write (6,*) ' dust drydep invoked '

!   Dry deposition of Dust Aerosols
!   #################################
!    call setdvel( ncol, landfrac, icefrac, ocnfrac, .001_r8, .001_r8, .001_r8, dvel )
!  we get the ram1,fv from the land model, but need to calculate it over oceans and ice.  
!  better if we got thse from the ocean and ice model
!  for friction velocity, we use ustar (from taux and tauy), except over land, where we use fv from the land model.
    call calcram(ncol,landfrac,icefrac,ocnfrac,obklen,ustar,ram1,state%t(:,pver),&
         state%pmid(:,pver),state%pdel(:,pver),fv)
       call outfld( 'airFV', fv(:), pcols, lchnk )
       call outfld( 'RAM1', ram1(:), pcols, lchnk )
!       call outfld( 'icefrc', icefrac(:), pcols, lchnk )

    call DustDryDep(ncol,state%t(:,:),state%pmid(:,:),ram1,fv,vlc_dry,vlc_trb,vlc_grv)

    do m=1,ndst

       mm = dust_idx1() + m - 1
! use pvdust instead (means making the top level 0)
       pvdust(:ncol,1)=0.
       pvdust(:ncol,2:pverp) = vlc_dry(:ncol,:,m)
       

       call outfld( trim(cnst_name(mm))//'DV', pvdust(:,2:pverp), pcols, lchnk )
       if(.true.) then ! use phil's method
!      convert from meters/sec to pascals/sec
!      pvdust(:,1) is assumed zero, use density from layer above in conversion
       pvdust(:ncol,2:pverp) = pvdust(:ncol,2:pverp) * rho(:ncol,:)*gravit        

!      calculate the tendencies and sfc fluxes from the above velocities
       call dust_sediment_tend( &
            ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
            state%q(:,:,mm) , pvdust  , ptend%q(:,:,mm), sflx  )
       else   !use charlie's method
         call d3ddflux(ncol, vlc_dry(:,:,m), state%q(:,:,mm),state%pmid,state%pdel, tvs,sflx,ptend%q(:,:,mm),dt)
        endif
! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i,m)/vlc_dry(i,pver,m)
          dep_grv(i)=sflx(i)*vlc_grv(i,pver,m)/vlc_dry(i,pver,m)
       enddo
       call outfld( trim(cnst_name(mm))//'DD',sflx, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'TB', dep_trb, pcols, lchnk )
       call outfld( trim(cnst_name(mm))//'GV', dep_grv, pcols, lchnk )

       call outfld( trim(cnst_name(mm))//'DT',ptend%q(:,:,mm), pcols, lchnk)
!       write (6,*) ' range of tends for dust ', mm, minval(ptend%q(:ncol,:,mm)), maxval(ptend%q(:ncol,:,mm))
!      call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)
      

    end do

! set flags for tendencies (water and 4 ghg's)
    ptend%name  = ptend%name//'+dust_drydep'
    ptend%lq(ioff+1:ioff+4) = .TRUE.

    return
  end subroutine dust_drydep_intr

  subroutine dust_time_interp

    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
!   use dustbnd, only: dustbndint

    implicit none

    real(r8) calday
    integer :: yr, mon, day, ncsec

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    write (6,*) ' dust_time_interp: interpolating dust emissions ', calday
!    call dustbndint(calday)      ! interpolate oxidants


  end subroutine dust_time_interp

  subroutine dust_emis_intr (state, ptend, dt,cflx)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all dusts.
! Notice that the mobilization is calculated in the land model (need #define BGC) and
! the soil erodibility factor is applied here.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: Phil Rasch and Natalie Mahowald
!
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual

!   use dust, only: dustsf
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    real(r8),            intent(inout) :: cflx(pcols,ppcnst)

    integer  lat(pcols)                  ! latitude index 
    integer  lon(pcols)                  ! longitude index
    integer lchnk
    integer ncol
    integer i
    integer m
    real(r8) :: soil_erod_tmp(pcols)
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    lchnk = state%lchnk
    ncol = state%ncol
    
    call get_lat_all_p(lchnk, ncol, lat)
    call get_lon_all_p(lchnk, ncol, lon)

    do i = 1, ncol
       soil_erod_tmp(i) = soil_erodibility(lon(i),lat(i))
    end do

         call outfld('LND_MBL',soil_erod_tmp,pcols,lchnk)
!
!         write(40,*) cflx(1,1)
!         write(41,*) soil_erodibility(1,lat(1))
    do m=dust_idx1(),dust_idx1()+3
! multiply by soil erodibility factor

       cflx(:ncol,m)=cflx(:ncol,m)*soil_erod_tmp(:ncol)
         
         call outfld(trim(cnst_name(m)) // 'SF',cflx(:,m),pcols,lchnk)
! this is being done inside of the vertical diffusion automatically
!         ptend%lq(m) = .true. ! tendencies for all dust on
!         ptend%q(:ncol,pver,m) = cflx(:ncol,m)*gravit/state%pdel(:ncol,pver)
      enddo

!      write(42,*) cflx(1,1)
    return
  end subroutine dust_emis_intr 

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Dustini()
!
! !INTERFACE:
!
  subroutine Dustini()
!
! !DESCRIPTION: 
!
! Compute source efficiency factor from topography
! Initialize other variables used in subroutine Dust:
! ovr_src_snk_mss(m,n) and tmp1.
! Define particle diameter and density needed by atm model
! as well as by dry dep model
! Source: Paul Ginoux (for source efficiency factor)
! Modifications by C. Zender and later by S. Levis
! Rest of subroutine from C. Zender's dust model
!
! !USES
!
    use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_RDAIR
!
! !ARGUMENTS:
!
    implicit none
!
! !REVISION HISTORY
! Created by Samual Levis
! Revised for CAM by Natalie Mahowald
!EOP
!------------------------------------------------------------------------

!------------------------------------------------------------------------
    !Local Variables
    integer  :: ci,m,n                  !indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ![frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ![frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ![frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ![frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 !Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 !Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ![m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ![m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ![m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ![m] Width of size bin
    real(r8) :: slp_crc(ndst)           ![frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ![m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ![m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ![m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ![frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ![frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     !temporary 
    real(r8) :: ln_gsd                  ![frc] ln(gsd)
    real(r8) :: gsd_anl                 ![frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ![m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ![m] Number median particle diameter
    real(r8) :: lgn_dst                 !Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ![frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ![frc] Current relative accuracy
    real(r8) :: itr_idx                 ![idx] Counting index
    real(r8) :: dns_mdp                 ![kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ![m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ![kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ![m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            !Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      !Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ![m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ![m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ![m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ![m] Size Bin widths
    
    ! constants
    real(r8) :: dmt_vma_src(dst_src_nbr) =    &     ![m] Mass median diameter
         (/ 0.832e-6 , 4.82e-6 , 19.38e-6 /)        !BSM96 p. 73 Table 2
    real(r8) :: gsd_anl_src(dst_src_nbr) =    &     ![frc] Geometric std deviation
         (/ 2.10     ,  1.90   , 1.60     /)        !BSM96 p. 73 Table 2
    real(r8) :: mss_frc_src(dst_src_nbr) =    &     ![frc] Mass fraction 
         (/ 0.036, 0.957, 0.007 /)                  !BSM96 p. 73 Table 2
    real(r8) :: dmt_grd(5) =                  &     ![m] Particle diameter grid
         (/ 0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6    ![m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0         ![kg m-3] Density of optimal saltation particles


    ! declare erf intrinsic function
    real(r8) :: dum     !dummy variable for erf test
#if (defined AIX) 
#define ERF erf
#else
#define ERF derf
    real(r8) derf
#endif
!------------------------------------------------------------------------

    ! Sanity check on erf: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

    dum = 1.0
    if (abs(0.8427-ERF(dum))/0.8427>0.001) then
       write (6,*) 'erf(1.0) = ',ERF(dum)
       write (6,*) 'Dustini: Error function error'
       call endrun
    endif
    dum = 0.0
    if (ERF(dum) /= 0.0) then
       write (6,*) 'erf(0.0) = ',ERF(dum)
       write (6,*) 'Dustini: Error function error'
       call endrun
    endif

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
       sqrt2lngsdi = sqrt(2.0) * log(gsd_anl_src(m))
       do n = 1, ndst
          lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
          lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
          ovr_src_snk_frc = 0.5 * (ERF(lndmaxjovrdmdni/sqrt2lngsdi) - &
                                   ERF(lndminjovrdmdni/sqrt2lngsdi))
          ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
       enddo
    enddo

!delete portions that have to do only with the sink

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)
    
    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if dust_number ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here dust_number = dst_nbr
    
    ! Override automatic grid with preset grid if available
    if (dust_number() == 4) then
       do n = 1, dust_number()
          dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
          dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
          dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
          dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
       end do
    else
       write (6,*) 'Dustini error: dust_number must equal to 4 with current code'
       call endrun                                      !see more comments above
    endif                                               !endif dust_number == 4

    ! Bin physical properties
    gsd_anl = 2.0      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e+3   ! [kg m-3] Aerosol density
    ! Set a fundamental statistic for each bin
    dmt_vma = 2.524e-6 ! [m] Mass median diameter analytic She84 p.75 Table1
    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)
    dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]
    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl
    do n = 1, dust_number()
       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do
       lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*SHR_CONST_PI))
       dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved
       do m = 1, sz_nbr
          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)
          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
       end do
       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do
    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)
    eps_max = 1.0e-4_r8
    dns_mdp = 100000. / (295.0*SHR_CONST_RDAIR) ![kg m-3] const prs_mdp & tpt_vrt
    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_r8 * ((295.0/273.0_r8)**1.5_r8) * 393.0_r8 / &
         (295.0+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000.*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*295.0)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, dust_number()
       slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
            (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
            gravit * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, dust_number()

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]
       do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_r8) then
             cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                  (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                  9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_r8) then
             cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                  (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_r8) then
             cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
          else
             write (6,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
             write(6,*)'Dustini error: Reynolds number too large in stk_crc_get()'
             call endrun 
          endif

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)
          vlc_grv(m) = sqrt(4.0_r8 * gravit * dmt_vwr(m) * slp_crc(m) * dns_aer / &
               (3.0_r8*cff_drg_grv(m)*dns_mdp))   
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          endif
          if (itr_idx > 20) then
             write (6,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
             goto 100                                        !to next iteration
          endif
          itr_idx = itr_idx + 1

       end do                                                !end while
100    continue   !Label to jump to when iteration does not converge
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities
    do m = 1, dust_number()
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

    return
  end subroutine Dustini

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine DustDryDep(c)
!
! !INTERFACE:
!
  subroutine DustDryDep(ncol,t,pmid,ram1,fv,vlc_dry,vlc_trb,vlc_grv)
!
! !DESCRIPTION: 
!
    ! Determine Turbulent dry deposition for dust. Calculate the turbulent 
    ! component of dust dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the dust dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CCSM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, dustini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
!
! !USES
!
    use shr_const_mod, ONLY: SHR_CONST_PI, SHR_CONST_RDAIR, SHR_CONST_BOLTZ
    use physconst, only: rair
!
! !ARGUMENTS:
!
    implicit none
!
    real(r8) :: t(pcols,pver)       !atm temperature (K)
    real(r8) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8) :: fv(pcols)           !friction velocity (m/s)
    real(r8) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: vlc_trb(pcols,ndst)  !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,ndst)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,ndst)  !dry deposn velocity (m/s)
    integer, intent(in) :: ncol
!
! !REVISION HISTORY
! Created by Sam Levis
! Modified for CAM by Natalie Mahowald
!EOP
!------------------------------------------------------------------------

!------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k          !indices
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn   ![frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(ndst) ![frc] Slip correction factor
    real(r8) :: rss_lmn(ndst) ![s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp           !temporary 

    ! constants
    real(r8),parameter::shm_nbr_xpn_lnd=-2./3. ![frc] shm_nbr_xpn over land
! needs fv and ram1 passed in from lnd model

!------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol
          rho=pmid(i,k)/rair/t(i,k)
          ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
          ! when code asks to use midlayer density, pressure, temperature,
          ! I use the data coming in from the atmosphere, ie t(i,k), pmid(i,k)
          
          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(SHR_CONST_PI*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air
          
          do m = 1, dust_number()
             slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
                  dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
             vlc_grv(i,k,m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                  gravit * slp_crc(m) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
             vlc_grv(i,k,m) = vlc_grv(i,k,m) * stk_crc(m)         ![m s-1] Correction to Stokes settling velocity
             vlc_dry(i,k,m)=vlc_grv(i,k,m)
          end do
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       do m = 1, dust_number()
          stk_nbr = vlc_grv(i,k,m) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))    ![frc] SeP97 p.965
          dff_aer = SHR_CONST_BOLTZ * t(i,k) * slp_crc(m) / &    ![m2 s-1]
               (3.0_r8*SHR_CONST_PI*vsc_dyn_atm(i,k)*dmt_vwr(m)) !SeP97 p.474
          shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972
          shm_nbr_xpn = shm_nbr_xpn_lnd                          ![frc]
          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 
          tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
          rss_lmn(m) = 1.0_r8 / (tmp*fv(i)) ![s m-1] SeP97 p.972,965
       end do

       ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)
       do m = 1, dust_number()
          rss_trb = ram1(i) + rss_lmn(m) + ram1(i)*rss_lmn(m)*vlc_grv(i,k,m) ![s m-1]
          vlc_trb(i,m) = 1.0_r8 / rss_trb                            ![m s-1]
          vlc_dry(i,k,m) = vlc_trb(i,m)  +vlc_grv(i,k,m)
       end do

    end do !end ncols loop

    return
  end subroutine DustDryDep
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Calcram
!
! !INTERFACE:
!

  subroutine  calcram(ncol,landfrac,icefrac,ocnfrac,obklen,ustar,ram1,t,pmid,pdel,fv)
!
! !DESCRIPTION: 
!  
! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
! from Seinfeld and Pandis, p.963.
!  
! Author: Natalie Mahowald
!
    implicit none
    integer, intent(in) :: ncol
    real(r8) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8) :: fv(pcols)                 ! sfc frc vel from land
    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: t(pcols)       !atm temperature (K)
    real(r8), intent(in) :: pmid(pcols)    !atm pressure (Pa)
    real(r8), intent(in) :: pdel(pcols)    !atm pressure (Pa)
      real(r8), parameter :: zzocen = 0.0001   ! Ocean aerodynamic roughness length
      real(r8), parameter :: zzsice = 0.0400   ! Sea ice aerodynamic roughness length
      real(r8), parameter :: xkar   = 0.4      ! Von Karman constant
 
! local variables
    real(r8) :: z,psi,psi0,nu,nu0,temp,ram
    integer :: i
!    write(*,*) rair,zzsice,zzocen,gravit,xkar
    do i=1,ncol
          z=pdel(i)*rair*t(i)/pmid(i)/gravit/2.0   !use half the layer height like Ganzefeld and Lelieveld, 1995
	  if(obklen(i).eq.0) then
	         psi=0.
	        psi0=0.
          else
          psi=min(max(z/obklen(i),-1.0_r8),1.0_r8)
          psi0=min(max(zzocen/obklen(i),-1.0_r8),1.0_r8)
          endif
	  temp=z/zzocen
          if(icefrac(i) > 0.5) then 
             psi0=min(max(zzsice/obklen(i),-1.0_r8),1.0_r8)
             temp=z/zzsice
	   endif
          if(psi> 0.) then
             ram=1/xkar/ustar(i)*(log(temp)+4.7*(psi-psi0))
          else
             nu=(1.00-15.000*psi)**(.25)
             nu0=(1.000-15.000*psi0)**(.25)
	    if(ustar(i).ne.0.) then
             ram=1/xkar/ustar(i)*(log(temp) &
               +log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)) &
                +2.0*(atan(nu)-atan(nu0)))
             else
	         ram=0.
	     endif
          endif
          if(landfrac(i) < 0.000000001) then
             fv(i)=ustar(i)
             ram1(i)=ram
          endif
!          write(*,*) i,pdel(i),t(i),pmid(i),gravit,obklen(i),psi,psi0,icefrac(i),nu,nu0,ram,ustar(i),&
!             log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)),2.0*(atan(nu)-atan(nu0))
       enddo
     return
      end subroutine calcram
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine d3ddflux
!
! !INTERFACE:
!
   subroutine  d3ddflux ( ncol, vlc_dry, q,pmid,pdel, tv, dep_dry,dep_dry_tend,dt)
! Description:
!Do 3d- settling deposition calculations following Zender's dust codes, Dec 02.
!
! Author: Natalie Mahowald
!
      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vlc_dry(pcols,pver)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols,pver)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    pmid(pcols,pver)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    pdel(pcols,pver)   ! delta pressure across level (Pa)
      real(r8), intent(in) ::    tv(pcols,pver)  ! midpoint virtual temperature in surface layer (K)
    real(r8),            intent(in)  :: dt             ! time step

! Output arguments:

      real(r8), intent(out) ::    dep_dry(pcols) ! flux due to dry deposition in kg /m^s/sec
      real(r8), intent(out) ::    dep_dry_tend(pcols,pver) ! flux due to dry deposition in kg /m^s/sec

! Local variables:

      real(r8) :: flux(pcols,0:pver)  ! downward flux at each level:  kg/m2/s 
      integer i,k
      do i=1,ncol
         flux(i,0)=0.
      enddo
      do k=1,pver
         do i = 1, ncol
            flux(i,k) = -min(vlc_dry(i,k) * q(i,k) * pmid(i,k) /(tv(i,k) * rair), &
                      q(i,k)*pdel(i,k)/gravit/dt)
            dep_dry_tend(i,k)=(flux(i,k)-flux(i,k-1))/pdel(i,k)*gravit  !kg/kg/s

         end do
      enddo
! surface flux:
      do i=1,ncol
         dep_dry(i)=flux(i,pver)
      enddo
      return
      end subroutine d3ddflux


end module dust_intr
