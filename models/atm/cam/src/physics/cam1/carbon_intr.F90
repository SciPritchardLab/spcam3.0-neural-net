#include <misc.h>
#include <params.h>

module carbon_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use shr_kind_mod,only: r8 => shr_kind_r8
  use pmgrid,      only: plat, plev, plevp, masterproc
  use ppgrid,      only: pcols, pver
  use physconst,   only: mwdry, mwh2o
  use constituents,only: ppcnst, cnst_add, cnst_name, advected, cnst_get_ind
  use caer,        only: caer_idx_set, caer_idx1, caer_ncomp
    

  implicit none

  private          ! Make default type private to the module

  save

  logical volcEmis
  character*4 ncnam(5)                           ! names of oxidants
  character(len=80) :: srcnam (4)                ! names of source/sink tendencies

  integer, parameter :: ncnst=4                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'OCPHO', 'BCPHO', 'OCPHI', 'BCPHI'/)

  integer ncyear
!
! Public interfaces
!
  public carbon_register_cnst                        ! register consituents
  public carbon_implements_cnst                      ! returns true if consituent is implemented by this package
  public carbon_init_cnst                            ! initialize mixing ratios if not read from initial file
  public carbon_initialize                           ! initialize (history) variables
  public carbon_wet_intr                             ! interface to wet deposition
  public carbon_emis_intr                            ! interface to emission
  public carbon_drydep_intr                          ! interface to tendency computation
  public carbon_time_interp                          ! interpolate oxidants and fluxes to current time

  real(r8), parameter:: cpso2  = 666.
  real(r8), parameter:: cpso4  = 666.
  real(r8), parameter:: cph2o2 = 666.
  real(r8), parameter:: cpdms  = 666.

  real(r8), parameter:: mwso2  = 64.
  real(r8), parameter:: mwso4  = 96.
  real(r8), parameter:: mwh2o2 = 34.
  real(r8), parameter:: mwdms  = 62.
  real(r8)  gravit

contains

!===============================================================================
  subroutine carbon_register_cnst
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

! and store the start index of carbon species used elsewhere in model retrieved by caer_idx1
    call caer_idx_set(m) 

    call cnst_add(cnst_names(2), advected, one, one, zero, m)
    call cnst_add(cnst_names(3), advected, one, one, zero, m)
    call cnst_add(cnst_names(4), advected, one, one, zero, m)

    return
  end subroutine carbon_register_cnst



!=======================================================================
  function carbon_implements_cnst(name)
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
     logical :: carbon_implements_cnst      ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     carbon_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           carbon_implements_cnst = .true.
           return
        end if
     end do
  end function carbon_implements_cnst


!=======================================================================
  subroutine carbon_init_cnst(name, q)
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

  end subroutine carbon_init_cnst


!===============================================================================
  subroutine carbon_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of carbon chemistry
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
    use caer,       only: caerini
    use caerbnd,    only: caerbndini
    use surface,    only: inisflx
    use drydep_mod, only: inidrydep
    use shr_const_mod,    only: SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use commap, only: clat, clon, w
    use ppgrid, only: pcols, begchunk, endchunk


    implicit none

!---------------------------Local workspace-----------------------------
    integer :: m, ix                                ! tracer index
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    character lpthncsf*200       ! local file for netCDF surface flux data
!    real(r8) :: gw(plat)
!-----------------------------------------------------------------------

    ix = caer_idx1() ! get the first index of advected carbon constituents
    gravit = SHR_CONST_G

!    do m = 1,plat/2
!       gw(m) = w(m)
!       gw(plat-m+1) = w(m)
!    end do
!    write (6,*) ' carbon_initialize:, oro dims ', pcols, begchunk, endchunk
!    call inisflx( SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH,  &
!         clon, clat, gw)

    calday = get_curr_calday()
    write (6,*) ' carbon_initialize: positioning at calday ', calday
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    ncyear = yr
    if ( ncyear .lt. 1000 ) ncyear = ncyear + 1900
    call caerini()
    call caerbndini( calday )

! Set names of variable tendencies and declare them as history variables
!    do m = 1, 4
!       srcnam(m) = trim(cnst_name(ix-1+m)) // 'SRC'
!       call addfld (srcnam(m),'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' source/sink',phys_decomp)
!       call add_default (srcnam(m), 1, ' ')
!    end do

#ifdef FLD_DEF1
    call addfld ('FFBCSF','kg/m2/s ',1, 'A','Fossil Fuel Black Carbon emissions',phys_decomp)
    call add_default ('FFBCSF', 1, ' ')
    call addfld ('FFOCSF','kg/m2/s ',1, 'A','Fossil Fuel Organic Carbon emissions',phys_decomp)
    call add_default ('FFOCSF', 1, ' ')
    call addfld ('BBBCSF','kg/m2/s ',1, 'A','Fossil Fuel Black Carbon emissions',phys_decomp)
    call add_default ('BBBCSF', 1, ' ')
    call addfld ('BBOCSF','kg/m2/s ',1, 'A','Fossil Fuel Organic Carbon emissions',phys_decomp)
    call add_default ('BBOCSF', 1, ' ')
#else
    call addfld ('OCPHOSF','kg/m2/s ',1, 'A','Organic Carbon Hydrophobic emissions',phys_decomp)
    call add_default ('OCPHOSF', 1, ' ')
    call addfld ('OCPHISF','kg/m2/s ',1, 'A','Organic Carbon Hydorphilic emissions',phys_decomp)
    call add_default ('OCPHISF', 1, ' ')
    call addfld ('BCPHOSF','kg/m2/s ',1, 'A','Black Carbon Hydrophobic emissions',phys_decomp)
    call add_default ('BCPHOSF', 1, ' ')
    call addfld ('BCPHISF','kg/m2/s ',1, 'A','Black Carbon Hydrophilic emissions',phys_decomp)
    call add_default ('BCPHISF', 1, ' ')
#endif
    call addfld ('BCPHIOD','kg/m2/s',1, 'A','hydrophilic BC Aerosol Optical Depth',phys_decomp)
    call add_default ('BCPHIOD', 1, ' ')
    call addfld ('BCPHOOD','kg/m2/s',1, 'A','hydrophobic BC Aerosol Optical Depth',phys_decomp)
    call add_default ('BCPHOOD', 1, ' ')

    call addfld ('OCPHIOD','kg/m2/s',1, 'A','hydrophilic OC Aerosol Optical Depth',phys_decomp)
    call add_default ('OCPHIOD', 1, ' ')
    call addfld ('OCPHOOD','kg/m2/s',1, 'A','hydrophobic OC Aerosol Optical Depth',phys_decomp)
    call add_default ('OCPHOOD', 1, ' ')

    call addfld ('BCPHIDRY','kg/m2/s ',1, 'A','hydrophilic BC dry deposition',phys_decomp)
    call add_default ('BCPHIDRY', 1, ' ')
    call addfld ('BCPHODRY','kg/m2/s ',1, 'A','hydrophobic BC dry deposition',phys_decomp)
    call add_default ('BCPHODRY', 1, ' ')
    call addfld ('BCPHIWET','kg/m2/s ',pver, 'A','hydrophilic BC wet deposition',phys_decomp)
    call add_default ('BCPHIWET', 1, ' ')

    call addfld ('OCPHIDRY','kg/m2/s ',1, 'A','hydrophilic OC dry deposition',phys_decomp)
    call add_default ('OCPHIDRY', 1, ' ')
    call addfld ('OCPHODRY','kg/m2/s ',1, 'A','hydrophobic OC dry deposition',phys_decomp)
    call add_default ('OCPHODRY', 1, ' ')
    call addfld ('OCPHIWET','kg/m2/s ',pver, 'A','hydrophilic OC wet deposition',phys_decomp)
    call add_default ('OCPHIWET', 1, ' ')

  end subroutine carbon_initialize


!===============================================================================
  subroutine carbon_wet_intr (state, ptend, cflx, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to we processing of aerosols (source and sinks).
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
    integer :: m                                  ! tracer index
    integer :: ixcldice, ixcldliq
    integer :: ioff                               ! offset for ghg indices
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
    real(r8) totcond(pcols, pver) ! sum of cldice and cldliq
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    integer :: mm
    integer :: nphob
    real(r8) :: sol_fact

!-----------------------------------------------------------------------
    ix = caer_idx1()
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol
    tstate(:ncol,:,:) = state%q(:ncol,:,ix:ix+3)

    totcond = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

!   Wet deposition of hydrophilic carbon aerosol species.
    nphob = caer_ncomp()/2
    ptend%name  = ptend%name//'+carbon_wetdep'
    do m = 1, nphob

       mm = caer_idx1() + nphob + m - 1

       ptend%lq(mm) = .TRUE.
       sol_fact = 0.3_r8
       call wetdepa( lat, state%t, state%pmid, state%q, state%pdel,  &
            cldn, cldc, cmfdqr, conicw, prain, cme,                     &
            evapr, totcond, state%q(:,:,mm), dt,            &
            ptend%q(:,:,mm), iscavt, cldv, fracis(:,:,mm), sol_fact, ncol )
       
       call outfld( trim(cnst_name(mm))//'WET', ptend%q(:,:,mm), pcols, lchnk)

    end do


    return

  end subroutine carbon_wet_intr

  subroutine carbon_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lat_all_p
    use constituents,  only: cnst_name
    use drydep_mod, only: drydepdr, setdvel, ddflux
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
    real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
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
    real(r8) :: tvs(pcols)
    real(r8) :: obuf(1)
    real(r8)  dvel(pcols)            ! deposition velocity
    real(r8)  sflx(pcols)            ! deposition flux
    real(r8), parameter :: mil  = .001_r8
!
!-----------------------------------------------------------------------
    ix = caer_idx1()
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol

    tvs(:ncol) = state%t(:ncol,pver)!*(1+state%q(:ncol,pver)

!   write (6,*) ' carbon drydep invoked '

!   Dry deposition of Carbon Aerosols
!   #################################
    call setdvel( ncol, landfrac, icefrac, ocnfrac, mil, mil, mil, dvel )
    do m = 1, caer_ncomp()
       mm = caer_idx1() + m - 1
       call ddflux( ncol, dvel, state%q(:,pver,mm), state%pmid(:,pver), tvs, sflx )
       ptend%q(:ncol,pver,mm) = sflx(:ncol)*gravit*state%rpdel(:ncol,pver)
#ifdef MATCH
       call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lat, obuf )
#else
       call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)
#endif
    end do

! set flags for tracer tendencies 
    ptend%name  = ptend%name//'+carbon_drydep'
    ptend%lq(ioff+1:ioff+4) = .TRUE.
!
! record tendencies on history files
!    do m = 1, 4
!       call outfld ('DRYDEP'//srcnam(m),ptend%q(:,:,ioff+m),pcols,lchnk)
!    end do

    return
  end subroutine carbon_drydep_intr

  subroutine carbon_time_interp

    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use caerbnd, only: caerbndint

    implicit none

    real(r8) calday
    integer :: yr, mon, day, ncsec

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

!   write (6,*) ' carbon_time_interp: interpolating carbon emissions ', calday
    call caerbndint(calday)      ! interpolate oxidants


  end subroutine carbon_time_interp

  subroutine carbon_emis_intr (state, ptend, dt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all carbons
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use caer, only: caersf
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    integer  lat(pcols)                  ! latitude index 
    integer  lon(pcols)                  ! longitude index
    integer lchnk
    integer ncol
    integer i
    integer m
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8) :: cflx(pcols,2)

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

    call caersf (ncol, 2, lat, lon, cflx) ! 2 means partition into organic and black

    call outfld('OCPHOSF', cflx(:, 1), pcols, lchnk)
    call outfld('BCPHOSF', cflx(:, 2), pcols, lchnk)

    m = caer_idx1()
    ptend%name  = ptend%name//'+carbon_emis'
    ptend%lq(m:m+3) = .true. ! tendencies for all carbon on
    ptend%q(:ncol,pver,m:m+3) = 0. ! zero all carbon tends
    ptend%q(:ncol,pver,m)   = cflx(:ncol,1)*gravit/state%pdel(:ncol,plev)  ! OCPHO source
    ptend%q(:ncol,pver,m+1) = cflx(:ncol,2)*gravit/state%pdel(:ncol,plev)  ! BCPHO source
  

    return
  end subroutine carbon_emis_intr 

end module carbon_intr
