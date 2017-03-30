#include <misc.h>
#include <params.h>

module sulfur_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use shr_kind_mod,only: r8 => shr_kind_r8
  use pmgrid,      only: plat, plev, plevp, masterproc
  use ppgrid,      only: pcols, pver
  use shr_const_mod, only: SHR_CONST_G
  use physconst,   only: mwdry, mwh2o
  use constituents,only: ppcnst, cnst_add, cnst_name, advected, cnst_get_ind
  use drydep_mod, only: drydepdr, setdvel, ddflux

  implicit none

  private          ! Make default type private to the module

  save

  logical volcEmis
  character*4 ncnam(5)                           ! names of oxidants
  character(len=80) :: srcnam (4)                ! names of source/sink tendencies
  integer ncyear

  integer :: ncnst                                   ! number of constituents
  character(len=8), allocatable :: cnst_names(:)     ! constituent names

!
! Public interfaces
!
  public sulfur_register_cnst                        ! register consituents
  public sulfur_implements_cnst                      ! returns true if consituent is implemented by this package
  public sulfur_init_cnst                            ! initialize mixing ratios if not read from initial file
  public sulfur_initialize                           ! initialize (history) variables
  public sulfur_emis_intr                            ! interface to emission computation
  public sulfur_wet_intr                             ! interface to wet chemistry and deposition
  public sulfur_drydep_intr                          ! interface to dry deposition
  public sulfur_time_interp                          ! interpolate oxidants and fluxes to current time

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
  subroutine sulfur_register_cnst
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
    use scyc, only: scyc_idx_set, scyc_ncomp

    implicit none
!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    integer :: istat
    real(r8), parameter :: zero  = 0._r8
!-----------------------------------------------------------------------

! Define constituent names
    ncnst = scyc_ncomp()
!TBH:  Should cnst_names live in scyc?  Probably...  
    allocate( cnst_names(ncnst), stat=istat )
    if ( istat /= 0 ) then
      write(6,*) 'sulfur_register_cnst: ERROR: allocate of cnst_names failed'
      call endrun
    end if
    cnst_names(1) = 'SO2'
    cnst_names(2) = 'SO4'
    cnst_names(3) = 'DMS'
    cnst_names(4) = 'H2O2'

! Set names of variables undergoing evolution
! returns m as current index for tracer numbers
    call cnst_add(cnst_names(1), advected, mwso2, cpso2, zero, m) 

! and store the start index of sulfur species used elsewhere in model retrieved by scyc_idx1
    call scyc_idx_set(m) 

    call cnst_add(cnst_names(2), advected, mwso4, cpso4, zero, m)
    call cnst_add(cnst_names(3), advected, mwdms, cpdms, zero, m)
    call cnst_add(cnst_names(4), advected, mwh2o2, cph2o2, zero, m)

    return
  end subroutine sulfur_register_cnst


!=======================================================================
  function sulfur_implements_cnst(name)
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
     logical :: sulfur_implements_cnst      ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     sulfur_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           sulfur_implements_cnst = .true.
           return
        end if
     end do
  end function sulfur_implements_cnst


!=======================================================================
  subroutine sulfur_init_cnst(name, q)
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

  end subroutine sulfur_init_cnst


!===============================================================================
  subroutine sulfur_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of sulfur chemistry
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
    use scyc,       only: scycini, scyc_idx1
    use surface,    only: inisflx
    use soxbnd,     only: soxbndini
    use dmsbnd,     only: dmsbndini
    use volcemission,     only: volcemisini
    use sulchem, only: inisulchem
    use drydep_mod, only: inidrydep
    use acsf,       only: iniacsf
    use acbnd,      only: acbndini
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual

    use shr_const_mod,    only: SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH
    use commap, only: clat, clon, w
    use ppgrid, only: pcols, begchunk, endchunk
    use ioFileMod, only: getfil
    use filenames, only: oxid

    implicit none

!---------------------------Local workspace-----------------------------
    integer :: m, ix                                ! tracer index
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    character lpthncsf*200       ! local file for netCDF surface flux data
    character(len=256) :: oxid_file  ! netcdf local filename for oxidants
!    real(r8) :: gw(plat)
!-----------------------------------------------------------------------

    gravit = SHR_CONST_G

    ix = scyc_idx1() ! get the first index of advected sulfur constituents

!    do m = 1,plat/2
!       gw(m) = w(m)
!       gw(plat-m+1) = w(m)
!    end do
!    write (6,*) ' sulfur_initialize:, oro dims ', pcols, begchunk, endchunk
!    call inisflx( SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH,  &
!         clon, clat, gw)

    calday = get_curr_calday()
    write (6,*) ' sulfur_initialize: positioning at calday ', calday
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    call scycini( .true., 'SMITH', volcEmis )  ! set some logicals in the sulfur cycle

    ncyear = yr
    if ( ncyear .lt. 1000 ) ncyear = ncyear + 1900
    call soxbndini( ncyear )
    call dmsbndini( calday )
    if ( volcEmis ) call volcemisini()
    call inisulchem()
    ncnam(1) = 'O3'
    ncnam(2) = 'OH'
    ncnam(3) = 'HO2'
    ncnam(4) = 'NO3'
    ncnam(5) = 'H2O2'
    ! Get file name.  
    call getfil(oxid, oxid_file, 0)
    call acbndini( oxid_file, calday, 5, ncnam )
    write (6,*) ' sulfur_initialize: opened oxidant file ',trim(oxid_file),' for sulfur species '

! Set names of variable tendencies and declare them as history variables
!    do m = 1, 4
!       srcnam(m) = trim(cnst_name(ix-1+m)) // 'SRC'
!       call addfld (srcnam(m),'kg/kg/s ',pver, 'A',trim(cnst_name(ix-1+m))//' source/sink',phys_decomp)
!       call add_default (srcnam(m), 1, ' ')
!    end do

    call addfld ('SO2SF','kg/m2/s ',1, 'A','SO2 emissions',phys_decomp)
    call add_default ('SO2SF', 1, ' ')
    call addfld ('TAUSUL','SO4 AOD ',1, 'A','SO4 Aerosol Optical Depth',phys_decomp)
    call add_default ('TAUSUL', 1, ' ')
    call addfld ('SO4SF','kg/m2/s ',1, 'A','SO4 emissions',phys_decomp)
    call add_default ('SO4SF', 1, ' ')
    call addfld ('DMSSF','kg/m2/s ',1, 'A','DMS emissions',phys_decomp)
    call add_default ('DMSSF', 1, ' ')
    call addfld ('SO2DRY','kg/m2/s ',1, 'A','SO2 dry deposition',phys_decomp)
    call add_default ('SO2DRY', 1, ' ')
    call addfld ('SO4DRY','kg/m2/s ',1, 'A','SO4 dry deposition',phys_decomp)
    call add_default ('SO4DRY', 1, ' ')
    call addfld ('SO2WET','kg/m2/s ',pver, 'A','SO2 wet deposition',phys_decomp)
    call add_default ('SO2WET', 1, ' ')
    call addfld ('SO4WET','kg/m2/s ',pver, 'A','SO4 wet deposition',phys_decomp)
    call add_default ('SO4WET', 1, ' ')
    call addfld ('H2O2WET','kg/m2/s ',pver, 'A','H2O2 wet deposition',phys_decomp)
    call add_default ('H2O2WET', 1, ' ')
    call addfld ('PH','pH',pver, 'A','Cloud Water pH',phys_decomp)
    call add_default ('PH', 1, ' ')
    call addfld ('O3','mol/mol',pver, 'A','Ozone',phys_decomp)
    call add_default ('O3', 1, ' ')
    call addfld ('O3INT','Dobson Units',1, 'A','Ozone',phys_decomp)
    call add_default ('O3INT', 1, ' ')
    call addfld('DMSSNK  ','kg/kg/s' ,pver, 'A', 'DMSsink ',phys_decomp)
    call add_default ('DMSSNK', 1, ' ')
    call addfld('SO2SRCG ','kg/kg/s' ,pver, 'A', 'SO2 Src G   ',phys_decomp)
    call add_default ('SO2SRCG', 1, ' ')
    call addfld('SO2SRCG2','kg/kg/s' ,pver, 'A', 'SO2 Src G2  ',phys_decomp)
    call add_default ('SO2SRCG2', 1, ' ')
    call addfld('SO4SRC  ','kg/kg/s' ,pver, 'A', 'SO4 Src Tot ',phys_decomp)
    call add_default ('SO4SRC', 1, ' ')
    call addfld('SO4SRCG ','kg/kg/s' ,pver, 'A', 'SO4 Src G   ',phys_decomp)
    call add_default ('SO4SRCG', 1, ' ')
    call addfld('H2O2SRC ','kg/kg/s' ,pver, 'A', 'H2O2Src ',phys_decomp)
    call add_default ('H2O2SRC', 1, ' ')
    call addfld('H2O2SNKA','kg/kg/s' ,pver, 'A', 'H2O2SnkA',phys_decomp)
    call add_default ('H2O2SNKA', 1, ' ')
    call addfld('H2O2SNKG','kg/kg/s' ,pver, 'A', 'H2O2SnkG',phys_decomp)
    call add_default ('H2O2SNKG', 1, ' ')

!   you could use this to read a generic surface flux from a netcdf file.
!   call iniacsf( lpthncsf, calday, cnst_name(ix) )! open netcdf files and position for todays fluxes

    call inidrydep (SHR_CONST_RDAIR, SHR_CONST_G, clat)

  end subroutine sulfur_initialize


  subroutine sulfur_emis_intr (state, ptend, dt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all sulfur species
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
!   use sulfur_intr,   only: sulfur_emis_intr
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
      use scyc, only: doscyc, scyc_idx1, add_volc_emis
      use sulemis, only: addsulemis
      use volcemission, only: volcemist
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
    real(r8) clat(pcols)                 ! latitude 
    integer lchnk
    integer ncol
    integer i
    integer m
    real(r8) astmp(pcols,pver,3)
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8) :: so2sf(pcols), so4sf(pcols), dmssf(pcols)
    integer :: k
    real(r8) :: gravit2       ! replace with gravit?  

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
    call get_rlat_all_p(lchnk, ncol, clat)

    m = scyc_idx1()
    astmp(:ncol,:pver,1:3) = state%q(:ncol,:pver,m:m+2)
    gravit2 = 9.8_r8
    call addsulemis( ncol, lat, lon, dt, gravit2, state%rpdel, state%zi, astmp)

!   add sources from dms, so2, and so4. Don't use or rely on sfc fluxes
!   because sources can be at 100m
    ptend%name  = ptend%name//'+sulfur_emiss'
    ptend%lq = .true.
    m = scyc_idx1()
    ptend%q(:ncol,:pver,m:m+2) = (astmp(:ncol,:pver,1:3)-state%q(:ncol,:pver,m:m+2))/dt

    so2sf = 0
    so4sf = 0
    dmssf = 0
    do k = 1,pver
       so2sf(:ncol) = so2sf(:ncol) + ptend%q(:ncol,k,m  )*state%pdel(:ncol,k)/gravit
       so4sf(:ncol) = so4sf(:ncol) + ptend%q(:ncol,k,m+1)*state%pdel(:ncol,k)/gravit
       dmssf(:ncol) = dmssf(:ncol) + ptend%q(:ncol,k,m+2)*state%pdel(:ncol,k)/gravit
    end do
    call outfld('SO2SF', so2sf, pcols, lchnk)
    call outfld('SO4SF', so4sf, pcols, lchnk)
    call outfld('DMSSF', dmssf, pcols, lchnk)

!        Volcano emissions.
!         if ( add_volc_emis() ) then
!            m = scyc_idx1()
!            call volcemist( lat, gravit, rpdel, so2tend )
!            call outfld( 'SO2VOLC', so2tend, pcols, lat)
!            do k = 1, pver
!               do i = 1, ncol
!                  as(i,k,m) = as(i,k,m) + dtime*so2tend(i,k)
!               end do
!            end do
!         end if
!      end if


    return
  end subroutine sulfur_emis_intr 

!===============================================================================
  subroutine sulfur_wet_intr (state, ptend, cflx, nstep, dt, lat, clat, lon, cme, prain, &
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
    use scyc,       only: scyc_idx1
    use sulchem,    only: chemwdepdr
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: nstep
    integer, intent(in) :: lat(pcols)                  ! latitude indices
    real(r8), intent(in) :: clat(pcols)                ! latitudes
    integer, intent(in) :: lon(pcols)                  ! longtitude indices
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
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    integer :: ixcldliq
    integer :: ixcldice
    real(r8) :: tstate(pcols, pver, 4)            ! temporary state vector
    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio
    real(r8) :: obuf(1)
    real(r8) :: calday        ! current calendar day
    real(r8) totcond(pcols, pver) ! total condensate
    integer :: yr, mon, day, ncsec
    integer :: ncdate

!-----------------------------------------------------------------------
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    ix = scyc_idx1()
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol
    tstate(:ncol,:,:) = state%q(:ncol,:,ix:ix+3)
    totcond = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

! compute tendencies and surface fluxes
    call chemwdepdr( lchnk, nstep, lat, clat, lon, calday, &
         dt, state%pmid, state%pdel, state%zm, state%t, state%q, &
         totcond, cldn, cldc, cldv, cmfdqr, &
         prain, evapr, cme, rainmr, conicw, &
         tstate, fracis(:,:,ix:ix+3), ncol)

    ptend%name  = ptend%name//'+sulfur_wetdep'
    ptend%q(:ncol,:,ix:ix+3) = (tstate(:ncol,:,:)-state%q(:ncol,:,ix:ix+3))/dt
    ptend%lq(ix:ix+3) = .true.

!   write (6,*) ' sulfur_wet_intr: pcols, ncol, lchnk ', pcols, ncol, lchnk 
!   write (6,*) ' range of ptend ix ', minval(ptend%q(:ncol,:,ix)),maxval(ptend%q(:ncol,:,ix))
!
! disabled here because we cant seperate the wet, dry, and chemical processes
! at this level
! record SO2,SO4,H2O2 tendencies on history files
!    do m = 1, 2
!       call outfld (trim(cnst_name(ioff+m))//'WET',ptend%q(:,:,ioff+m),pcols,lchnk)
!    end do
!    call outfld (trim(cnst_name(ioff+4))//'WET',ptend%q(:,:,ioff+m),pcols,lchnk)

    return

  end subroutine sulfur_wet_intr

  subroutine sulfur_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
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
    use scyc,       only: scyc_idx1
    use drydep_mod, only: drydepdr
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
    integer :: mm                                 ! tracer index
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: ix
    real(r8) :: tvs(pcols)
    real(r8) :: obuf(1)
    real(r8)  dvel(pcols)            ! deposition velocity
    real(r8)  sflx(pcols)            ! deposition flux
    real(r8), parameter :: mil   = .001_r8
    real(r8), parameter :: mil2  = .002_r8
    real(r8), parameter :: mil6  = .006_r8
    real(r8), parameter :: mil8  = .008_r8
!
!-----------------------------------------------------------------------
    ix = scyc_idx1()
    ioff  = ix - 1
    lchnk = state%lchnk
    ncol  = state%ncol

    tvs(:ncol) = state%t(:ncol,pver)!*(1+state%q(:ncol,pver)

!   Dry Deposition Velocities:
!   Here are the values used by Hans Feichter for SO2:
!   .6cm/s over land
!   .8cm/s over ocean
!   .1cm/s over ice or snow
!   For SO4 he uses .2cm/s everywhere

!   Dry deposition of SO2
!   #####################
    mm = scyc_idx1()
    call setdvel( ncol, landfrac, icefrac, ocnfrac, mil6, mil8, mil, dvel )
    call ddflux( ncol, dvel, state%q(:,pver,mm), state%pmid(:,pver), tvs, sflx )
#ifdef MATCH
    call outfld( trim(tracnam(mm))//'DRY', sflx, pcols, &
         lat, obuf )
#else
    call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk )
#endif
#ifdef SCYC_MASSBGT
    call aal( 'so2dry', lat, flx(1,mm) )
#endif
    ptend%name  = ptend%name//'+sulfur_drydep'
    ptend%q(:ncol,pver,mm) = sflx(:ncol)*gravit*state%rpdel(:ncol,pver)

!   Dry deposition of SO4
!   #####################
    mm = scyc_idx1()+1
    call setdvel( ncol, landfrac, icefrac, ocnfrac, mil2, mil, mil2, dvel )
    call ddflux( ncol, dvel, state%q(:,pver,mm), state%pmid(:,pver), tvs, sflx )
#ifdef MATCH
    call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols,  lat, obuf )
#else
    call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)
#endif
#ifdef SCYC_MASSBGT
    call aal( 'so4dry', lat, sflx )
#endif
    ptend%q(:ncol,pver,mm) = sflx(:ncol)*gravit*state%rpdel(:ncol,pver)

!   Dry deposition of H2O2
!   #####################
    mm = scyc_idx1()+3
    call setdvel( ncol, landfrac, icefrac, ocnfrac, mil6, mil8, mil, dvel )
    call ddflux( ncol, dvel, state%q(:,pver,mm), state%pmid(:,pver), tvs, sflx )
#ifdef MATCH
    call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lat, obuf )
#else
    call outfld( trim(cnst_name(mm))//'DRY', sflx, pcols, lchnk)
#endif
    ptend%q(:ncol,pver,mm) = sflx(:ncol)*gravit*state%rpdel(:ncol,pver)

!   write (6,*) ' sulfur drydep invoked '
!    call drydepdr( lchnk, ncol, month, dt, landfrac, icefrac, ocnfrac, &
!         snowh, prect, &
!         ts, wvflx, hflx, pblh, ustar, obklen, &
!         state%pmid(1,plev), state%rpdel, tvs, state%zi, state%u, state%v, &
!         state%t, fsds, state%q, ptend%q, obuf )


! set flags for tracer tendencies 
    ptend%name  = ptend%name//'+sulfur_drydep'
    ptend%lq(ioff+1:ioff+4) = .TRUE.
!
! record tendencies on history files
!    do m = 1, 4
!       call outfld ('DRYDEP'//srcnam(m),ptend%q(:,:,ioff+m),pcols,lchnk)
!    end do

    return
  end subroutine sulfur_drydep_intr

  subroutine sulfur_time_interp

    use acbnd,      only: acbndint
    use soxbnd,      only: soxbndint
    use dmsbnd,      only: dmsbndint
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual

    implicit none

    real(r8) calday
    integer :: yr, mon, day, ncsec

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

!   write (6,*) ' sulfur_time_interp: interpolating oxidants, sox and dms to calday ', calday
    call acbndint(calday)      ! interpolate oxidants
    call soxbndint(yr, calday) ! interpolate sox
    call dmsbndint(calday)     ! interpolate dms

  end subroutine sulfur_time_interp

end module sulfur_intr
