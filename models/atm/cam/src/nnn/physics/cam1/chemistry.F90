#include <misc.h>
#include <params.h>

module chemistry

!---------------------------------------------------------------------------------
! Module to parameterized greenhouse gas chemical loss frequencies from 
! Portmann and Solomon
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat, plev, plevp, plon, masterproc
  use ppgrid,       only: pcols, pver
  use physconst,    only: mwdry, mwch4, mwn2o, mwf11, mwf12, mwh2o, gravit
  use constituents, only: ppcnst, cnst_add, cnst_name, advected
  use ghg_surfvals, only: ch4vmr, n2ovmr, f11vmr, f12vmr

  implicit none

  private
  save
!
! Public interfaces
!
  public chem_register                             ! register consituents
  public chem_implements_cnst                      ! returns true if consituent is implemented by this package
  public chem_init_cnst                            ! initialize mixing ratios if not read from initial file
  public chem_init                                 ! initialize (history) variables
  public chem_timestep_init                        ! time interpolate chemical loss frequencies
  public chem_timestep_tend                        ! interface to tendency computation

! Namelist variable
  logical, public :: trace_gas=.false.             ! set true to activate this module

  integer, parameter :: ptrlon=01                  ! number of longitudes in input dataset
  integer, parameter :: ptrlat=36                  ! number of latitudes in input dataset
  integer, parameter :: ptrlev=56                  ! number of levels in input dataset
  integer, parameter :: ptrtim=12                  ! number of times(months) in input dataset

! Ratios of molecular weights
  real(r8), parameter :: rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
  real(r8), parameter :: rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
  real(r8), parameter :: rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
  real(r8), parameter :: rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
  real(r8), parameter :: rh2och4= mwh2o/mwch4      ! ratio of molecular weight h2o   to ch4

! Time/space interpolation of loss frequencies
  real(r8) :: tch4i  (plat,plev,ptrtim) ! input data  ch4   loss rate interp. in lat and lev
  real(r8) :: tn2oi  (plat,plev,ptrtim) ! input data  n2o   loss rate interp. in lat and lev
  real(r8) :: tcfc11i(plat,plev,ptrtim) ! input data  cfc11 loss rate interp. in lat and lev
  real(r8) :: tcfc12i(plat,plev,ptrtim) ! input data  cfc12 loss rate interp. in lat and lev
  real(r8) :: tch4m  (plat,plev,2)      ! input data  ch4   loss rate interp. in lat and lev
  real(r8) :: tn2om  (plat,plev,2)      ! input data  n2o   loss rate interp. in lat and lev
  real(r8) :: tcfc11m(plat,plev,2)      ! input data  cfc11 loss rate interp. in lat and lev
  real(r8) :: tcfc12m(plat,plev,2)      ! input data  cfc12 loss rate interp. in lat and lev
  real(r8) :: tch4   (plat,plev)        ! instantaneous ch4   loss rate 
  real(r8) :: tn2o   (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: tcfc11 (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: tcfc12 (plat,plev)        ! instantaneous ch4   loss rate
  real(r8) :: cdaytrm                   ! calendar day for previous month data
  real(r8) :: cdaytrp                   ! calendar day for next month data

  integer :: np                               ! array index for previous month tracer data
  integer :: nm                               ! array index for next month tracer data
  integer :: np1                              ! current forward time index of tracer dataset
  integer :: date_tr(ptrtim)                  ! date on tracer dataset (YYYYMMDD)
  integer :: sec_tr(ptrtim)                   ! seconds of date on tracer dataset (0-86399)

! dummy values for specific heats at constant pressure
  real(r8), parameter:: cpch4 = 666.
  real(r8), parameter:: cpn2o = 666.
  real(r8), parameter:: cpf11 = 666.
  real(r8), parameter:: cpf12 = 666.

  integer, parameter :: ncnst=4                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)
  character(len=8) :: srcnam(ncnst)                  ! names of source/sink tendencies
  integer :: ixghg                                   ! index of 1st constituent (N2O)

contains

!===============================================================================
  subroutine chem_register
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents for parameterized greenhouse gas chemistry
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------

! Set names of diffused variable tendencies and declare them as history variables

    call cnst_add(cnst_names(1), advected, mwn2o, cpn2o, 0._r8, ixghg, longname='Nitrous Oxide')
    call cnst_add(cnst_names(2), advected, mwch4, cpch4, 0._r8, m, longname='Methane')
    call cnst_add(cnst_names(3), advected, mwf11, cpf11, 0._r8, m)
    call cnst_add(cnst_names(4), advected, mwf12, cpf12, 0._r8, m)

    return
  end subroutine chem_register

!===============================================================================

  function chem_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: chem_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     chem_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           chem_implements_cnst = .true.
           return
        end if
     end do
  end function chem_implements_cnst

!===============================================================================
  subroutine chem_init
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterized greenhouse gas chemistry
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

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------

! Set names of diffused variable tendencies and declare them as history variables
    do m = 1, 4
       srcnam(m) = trim(cnst_name(ixghg-1+m)) // 'SRC'
       call addfld (srcnam(m),'kg/kg/s ',pver, 'A',trim(cnst_name(ixghg-1+m))//' source/sink',phys_decomp)
       call add_default (srcnam(m), 1, ' ')
    end do

    call chem_init_loss

    return
  end subroutine chem_init

!===============================================================================
  subroutine chem_timestep_tend (state, ptend, cflx, dt, fh2o)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
    use history,       only: outfld
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_lat_all_p
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

    type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
    real(r8), intent(inout) :: cflx(pcols,ppcnst)      ! Surface constituent flux (kg/m^2/s)
    real(r8), intent(out) :: fh2o(pcols)          ! H2O flux required to balance the H2O source from methane chemistry (kg/m^2/s)
!
! Local variables
!
    integer :: i, k, m                            ! indices
    integer :: ioff                               ! offset for ghg indices
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: lat(pcols)                         ! latitude index for S->N storage
!
!-----------------------------------------------------------------------
! Initialize output tendency structure
    call physics_ptend_init(ptend)

    ioff  = ixghg - 1
    lchnk = state%lchnk
    ncol  = state%ncol

! get latitude indices
    call get_lat_all_p(lchnk, ncol, lat)

! compute tendencies and surface fluxes
    call ghg_chem ( lchnk, ncol, lat,                                               &
         state%q(:,:,1), state%q(:,:,ixghg:ixghg+3),                              &
         ptend%q(:,:,1), ptend%q(:,:,ixghg:ixghg+3), cflx(:,ixghg:ixghg+3), dt)

! set flags for tracer tendencies (water and 4 ghg's)
    ptend%lq(1)             = .TRUE.
    ptend%lq(ioff+1:ioff+4) = .TRUE.
!
! record tendencies on history files
    do m = 1, 4
       call outfld (srcnam(m),ptend%q(:,:,ioff+m),pcols,lchnk)
    end do

! Compute water vapor flux required to make the conservation checker happy
    fh2o = 0
    do k = 1, pver
       do i = 1, ncol
          fh2o(i) = fh2o(i) + ptend%q(i,k,1)*state%pdel(i,k)/gravit
       end do
    end do

    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine ghg_chem (lchnk, ncol, lat, qh2o, qghg, dqh2o, dqghg, fghg, dt)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Apply the interpolated chemical loss rates from the input data to
! N2O, CH4, CFC11 and CFC12. Set the surface values to a constant.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: lat(pcols)             ! latitude index for S->N storage

    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: qh2o(pcols,pver)      ! mass mixing ratios of water vapor
    real(r8), intent(in) :: qghg(pcols,pver,4)    ! mass mixing ratios of greenhouse gases

    real(r8), intent(out) :: dqh2o(pcols,pver)    ! tendency of mass mixing ratios (water)
    real(r8), intent(out) :: dqghg(pcols,pver,4)  ! tendency of mass mixing ratios (ghg's)
    real(r8), intent(out) :: fghg(pcols,4)        ! Surface constituent flux (kg/m^2/s)
!
! Local variables
!
    integer i,k                                   ! loop indexes
    real(r8) xch4                                 ! new methane mass mixing ratio
    real(r8) xn2o                                 ! new nitrous oxide mass mixing ratio
    real(r8) xcfc11                               ! new cfc11 mass mixing ratio
    real(r8) xcfc12                               ! new cfc12 mass mixing ratio
!
!-----------------------------------------------------------------------
!
! Apply chemical rate coefficient using time split implicit method. The
! turn the new value back into a tendency. NOTE that water tendency is
! twice methane tendency. Water is specific humidity which is in mass
! mixing ratio units. Note that
!  o 1   => indx of n2o
!  o 2 => indx of ch4
!  o 3 => indx of cfc11
!  o 4 => indx of cfc12
!
    do k=1,pver-2
       do i=1,ncol
          xn2o         = qghg(i,k,1) / (1. + tn2o  (lat(i),k) * dt)
          xch4         = qghg(i,k,2) / (1. + tch4  (lat(i),k) * dt)
          xcfc11       = qghg(i,k,3) / (1. + tcfc11(lat(i),k) * dt)
          xcfc12       = qghg(i,k,4) / (1. + tcfc12(lat(i),k) * dt)

          dqghg(i,k,1) =(xn2o   - qghg(i,k,1)) / dt
          dqghg(i,k,2) =(xch4   - qghg(i,k,2)) / dt
          dqghg(i,k,3) =(xcfc11 - qghg(i,k,3)) / dt
          dqghg(i,k,4) =(xcfc12 - qghg(i,k,4)) / dt

          dqh2o(i,k)   = -2. * rh2och4 * dqghg(i,k,2)
       end do
    end do
!
! Set the "surface" tendencies (bottom 2 levels) to maintain specified
! tropospheric concentrations.
!
    do k = pver-1, pver
       do i=1,ncol
          dqghg(i,k,1) =((rmwn2o*n2ovmr) - qghg(i,k,1)) / dt
          dqghg(i,k,2) =((rmwch4*ch4vmr) - qghg(i,k,2)) / dt
          dqghg(i,k,3) =((rmwf11*f11vmr) - qghg(i,k,3)) / dt
          dqghg(i,k,4) =((rmwf12*f12vmr) - qghg(i,k,4)) / dt

          dqh2o(i,k)   = 0.
       end do
    end do
!
! For now set all tracer fluxes to 0
!
    do i=1,ncol
       fghg(i,1) = 0.
       fghg(i,2) = 0.
       fghg(i,3) = 0.
       fghg(i,4) = 0.
    end do

    return
  end subroutine ghg_chem

!===============================================================================
  subroutine chem_init_loss
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Do initial read of time-variant chemical loss rate frequency dataset, containing
! loss rates as a function of latitude and pressure.  Determine the two
! consecutive months between which the current date lies.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
!-----------------------------------------------------------------------
    use ioFileMod
    use commap
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use filenames, only: bndtvg

    implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>
    include 'netcdf.inc'

!
! Local variables
!
    integer dateid            ! netcdf id for date variable
    integer secid             ! netcdf id for seconds variable
    integer lonid             ! netcdf id for longitude variable
    integer latid             ! netcdf id for latitude variable
    integer levid             ! netcdf id for level variable
    integer timid             ! netcdf id for time variable
    integer tch4id            ! netcdf id for ch4   loss rate
    integer tn2oid            ! netcdf id for n2o   loss rate
    integer tcfc11id          ! netcdf id for cfc11 loss rate
    integer tcfc12id          ! netcdf id for cfc12 loss rate
    integer lonsiz            ! size of longitude dimension on tracer dataset
    integer levsiz            ! size of level dimension on tracer dataset
    integer latsiz            ! size of latitude dimension on tracer dataset
    integer timsiz            ! size of time dimension on tracer dataset
    integer j,n,k,nt          ! indices
    integer ki,ko,ji,jo       ! indices
    integer :: yr, mon, day   ! components of a date
    integer :: ncdate         ! current date in integer format [yyyymmdd]
    integer :: ncsec          ! current time of day [seconds]

    real(r8) :: calday        ! current calendar day
    real(r8) lato(plat)       ! cam model latitudes (degrees)
    real(r8) zo(plev)         ! cam model heights (m)
    real(r8) lati(ptrlat)     ! input data latitudes (degrees)
    real(r8) pin(ptrlev)      ! input data pressure values (mbars)
    real(r8) zi(ptrlev)       ! input data heights (m)

    real(r8) xch4i  (ptrlon,ptrlev,ptrlat,ptrtim) ! input ch4   loss rate coeff
    real(r8) xn2oi  (ptrlon,ptrlev,ptrlat,ptrtim) ! input n2o   loss rate coeff
    real(r8) xcfc11i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc11 loss rate coeff
    real(r8) xcfc12i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc12 loss rate coeff

    real(r8) xch4(ptrlat, ptrlev)   ! input ch4   loss rate coeff indices changed
    real(r8) xn2o(ptrlat, ptrlev)   ! input n2o   loss rate coeff indices changed
    real(r8) xcfc11(ptrlat, ptrlev) ! input cfc11 loss rate coeff indices changed
    real(r8) xcfc12(ptrlat, ptrlev) ! input cfc12 loss rate coeff indices changed

    real(r8) xch4lv(ptrlat, plev)   ! input ch4   loss rate coeff interp to cam levels
    real(r8) xn2olv(ptrlat, plev)   ! input n2o   loss rate coeff interp to cam levels
    real(r8) xcfc11lv(ptrlat, plev) ! input cfc11 loss rate coeff interp to cam levels
    real(r8) xcfc12lv(ptrlat, plev) ! input cfc12 loss rate coeff interp to cam levels

    character(len=256) :: locfn    ! netcdf local filename to open 
!
!-----------------------------------------------------------------------
!
! Initialize
!
    nm = 1
    np = 2
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
    if (masterproc) then
       write(6,*)'CHEM_INIT_LOSS: greenhouse loss rates from: ', trim(bndtvg)
       call getfil(bndtvg, locfn)
       call wrap_open(locfn, 0, ncid_trc)
       write(6,*)'CHEM_INIT_LOSS: NCOPN returns id ',ncid_trc,' for file ',trim(locfn)
!
!------------------------------------------------------------------------
! Read tracer data
!------------------------------------------------------------------------
!
! Get dimension info
!
       call wrap_inq_dimid(ncid_trc, 'lat' , latid)
       call wrap_inq_dimid(ncid_trc, 'lev' , levid)
       call wrap_inq_dimid(ncid_trc, 'lon' , lonid)
       call wrap_inq_dimid(ncid_trc, 'time', timid)

       call wrap_inq_dimlen(ncid_trc, lonid, lonsiz)
       call wrap_inq_dimlen(ncid_trc, levid, levsiz)
       call wrap_inq_dimlen(ncid_trc, latid, latsiz)
       call wrap_inq_dimlen(ncid_trc, timid, timsiz)
!
! Check dimension info
!
       if (ptrlon/=1) then
          write(6,*)'CHEM_INIT_LOSS: longitude dependence not implemented'
          call endrun
       endif
       if (lonsiz /= ptrlon) then
          write(6,*)'CHEM_INIT_LOSS: lonsiz=',lonsiz,' must = ptrlon=',ptrlon
          call endrun
       end if
       if (levsiz /= ptrlev) then
          write(6,*)'CHEM_INIT_LOSS: levsiz=',levsiz,' must = ptrlev=',ptrlev
          call endrun
       end if
       if (latsiz /= ptrlat) then
          write(6,*)'CHEM_INIT_LOSS: latsiz=',latsiz,' must = ptrlat=',ptrlat
          call endrun
       end if
       if (timsiz /= ptrtim) then
          write(6,*)'CHEM_INIT_LOSS: timsiz=',timsiz,' must = ptrtim=',ptrtim
          call endrun
       end if
!
! Determine necessary dimension and variable id's
!
       call wrap_inq_varid(ncid_trc, 'lat'    , latid)
       call wrap_inq_varid(ncid_trc, 'lev'    , levid)
       call wrap_inq_varid(ncid_trc, 'date'   , dateid)
       call wrap_inq_varid(ncid_trc, 'datesec', secid)
       call wrap_inq_varid(ncid_trc, 'TCH4'   , tch4id)
       call wrap_inq_varid(ncid_trc, 'TN2O'   , tn2oid)
       call wrap_inq_varid(ncid_trc, 'TCFC11' , tcfc11id)
       call wrap_inq_varid(ncid_trc, 'TCFC12' , tcfc12id)
!
! Obtain entire date and sec variables. Assume that will always
! cycle over 12 month data.
!
       call wrap_get_var_int(ncid_trc, dateid, date_tr)
       call wrap_get_var_int(ncid_trc, secid , sec_tr)

       if (mod(date_tr(1),10000)/100 /= 1) then
          write(6,*)'(CHEM_INIT_LOSS): error when cycling data: 1st month must be 1'
          call endrun
       end if
       if (mod(date_tr(ptrtim),10000)/100 /= 12) then
          write(6,*)'(CHEM_INIT_LOSS): error when cycling data: last month must be 12'
          call endrun
       end if
!
! Obtain input data latitude and level arrays.
!
       call wrap_get_var_realx(ncid_trc, latid, lati)
       call wrap_get_var_realx(ncid_trc, levid, pin )
!
! Convert input pressure levels to height (m).
! First convert from millibars to pascals.
!
       do k=1,ptrlev
          pin(k) = pin(k)*100.
          zi(k) = 7.0e3 * log (1.0e5 / pin(k))
       end do
!
! Convert approximate cam pressure levels to height (m).
!
       do k=1,plev
          zo (k) = 7.0e3 * log (1.0e5 / hypm(k))
       end do
!
! Convert cam model latitudes to degrees.
! Input model latitudes already in degrees.
!
       do j=1,plat
          lato(j) = clat(j)*45./atan(1.)
       end do
!
! Obtain all time samples of tracer data.
!
       call wrap_get_var_realx(ncid_trc, tch4id  , xch4i  )
       call wrap_get_var_realx(ncid_trc, tn2oid  , xn2oi  )
       call wrap_get_var_realx(ncid_trc, tcfc11id, xcfc11i)
       call wrap_get_var_realx(ncid_trc, tcfc12id, xcfc12i)
!
! Close netcdf file
!
       call wrap_close(ncid_trc)
!
!------------------------------------------------------------------------
! Interpolate tracer data to model grid
!------------------------------------------------------------------------
!
! Loop over all input times.
!
       do nt = 1, ptrtim
!
! Remove longitude and time index and switch level and latitude indices
! for the loss coefficients.
!
          do j=1,ptrlat
             do k=1,ptrlev
                xch4  (j,k) = xch4i  (1,k,j,nt)
                xn2o  (j,k) = xn2oi  (1,k,j,nt)
                xcfc11(j,k) = xcfc11i(1,k,j,nt)
                xcfc12(j,k) = xcfc12i(1,k,j,nt)
             end do
          end do
!
! Interpolate input data to model levels.
! If the CAM level is outside the range of the input data (this
! can happen only in troposphere) put zero for every latitude.
! Otherwise determine the input data levels bounding the current
! CAM level and interpolate.
!
          do ko=1,plev
             if (zo(ko) < zi(ptrlev)) then
                do j=1,ptrlat
                   xch4lv  (j,ko) = 0.0
                   xn2olv  (j,ko) = 0.0
                   xcfc11lv(j,ko) = 0.0
                   xcfc12lv(j,ko) = 0.0
                end do
                goto 50
             end if
             do ki=1,ptrlev-1
                if (zo(ko) < zi(ki) .and. zo(ko) >= zi(ki+1)) then
                   do j=1,ptrlat
                      xch4lv(j,ko) = xch4(j,ki) + (xch4(j,ki+1) - xch4(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xn2olv(j,ko) = xn2o(j,ki) + (xn2o(j,ki+1) - xn2o(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xcfc11lv(j,ko) = xcfc11(j,ki) + (xcfc11(j,ki+1) - xcfc11(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                      xcfc12lv(j,ko) = xcfc12(j,ki) + (xcfc12(j,ki+1) - xcfc12(j,ki)) &
                           / (zi(ki+1) - zi(ki)) * (zo(ko) - zi(ki))
                   end do
                   goto 50
                endif
             end do
             write (6,*) '(CHEM_INIT_LOSS): Error in vertical interpolation'
             call endrun
50           continue
          end do
!
! Interpolate input data to model latitudes.
! Determine the input data latitudes bounding the current CAM latitude and
! interpolate. Use last value from input data if the cam latitude is
! outside the range of the input data latitudes.
!
          do jo=1,plat
             if (lato(jo) <= lati(1)) then
                do k = 1, plev
                   tch4i(jo,k,nt)   = xch4lv(1,k)
                   tn2oi(jo,k,nt)   = xn2olv(1,k)
                   tcfc11i(jo,k,nt) = xcfc11lv(1,k)
                   tcfc12i(jo,k,nt) = xcfc12lv(1,k)
                end do
             else if (lato(jo) >= lati(ptrlat)) then
                do k = 1, plev
                   tch4i(jo,k,nt)   = xch4lv(ptrlat,k)
                   tn2oi(jo,k,nt)   = xn2olv(ptrlat,k)
                   tcfc11i(jo,k,nt) = xcfc11lv(ptrlat,k)
                   tcfc12i(jo,k,nt) = xcfc12lv(ptrlat,k)
                end do
             else
                do ji=1,ptrlat-1
                   if ( (lato(jo) > lati(ji)) .and. (lato(jo) <= lati(ji+1))) then
                      do k=1,plev
                         tch4i(jo,k,nt) = xch4lv(ji,k) + (xch4lv(ji+1,k) - xch4lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tn2oi(jo,k,nt) = xn2olv(ji,k) + (xn2olv(ji+1,k) - xn2olv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tcfc11i(jo,k,nt) = xcfc11lv(ji,k) + (xcfc11lv(ji+1,k) - xcfc11lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                         tcfc12i(jo,k,nt) = xcfc12lv(ji,k) + (xcfc12lv(ji+1,k) - xcfc12lv(ji,k)) &
                              / (lati(ji+1)   -  lati(ji)) * (lato(jo) - lati(ji))
                      end do
                      goto 90
                   endif
                end do
             end if
             write (6,*)'(CHEM_INIT_LOSS): Error in horizontal interpolation'
90           continue
          end do

       end do                 ! end loop over time samples
!
! Initial time interpolation between December and January
!
       calday = get_curr_calday()
       if ( is_perpetual() ) then
          call get_perp_date(yr, mon, day, ncsec)
       else
          call get_curr_date(yr, mon, day, ncsec)
       end if
       ncdate = yr*10000 + mon*100 + day

       n = 12
       np1 = 1
       call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
       call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
       if (calday <= cdaytrp .or. calday > cdaytrm) then
          do j=1,plat
             do k=1,plev
                tch4m  (j,k,nm) = tch4i  (j,k,n)
                tn2om  (j,k,nm) = tn2oi  (j,k,n)
                tcfc11m(j,k,nm) = tcfc11i(j,k,n)
                tcfc12m(j,k,nm) = tcfc12i(j,k,n)
                tch4m  (j,k,np) = tch4i  (j,k,np1)
                tn2om  (j,k,np) = tn2oi  (j,k,np1)
                tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                tcfc12m(j,k,np) = tcfc12i(j,k,np1)
             end do
          end do
          goto 10
       end if
!
! Initial normal interpolation between consecutive time slices.
!
       do n=1,timsiz-1
          np1 = n + 1
          call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
          call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
          if (calday > cdaytrm .and. calday <= cdaytrp) then
             do j=1,plat
                do k=1,plev
                   tch4m  (j,k,nm) = tch4i  (j,k,n)
                   tn2om  (j,k,nm) = tn2oi  (j,k,n)
                   tcfc11m(j,k,nm) = tcfc11i(j,k,n)
                   tcfc12m(j,k,nm) = tcfc12i(j,k,n)
                   tch4m  (j,k,np) = tch4i  (j,k,np1)
                   tn2om  (j,k,np) = tn2oi  (j,k,np1)
                   tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                   tcfc12m(j,k,np) = tcfc12i(j,k,np1)
                end do
             end do
             goto 10
          end if
       end do
       write(6,*)'CHEM_INIT_LOSS: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
       call endrun
!
! Data positioned correctly
!
10     continue

    endif                     ! end of masterproc

    return
  end subroutine chem_init_loss

!===============================================================================
  subroutine chem_timestep_init
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Time interpolate chemical loss rates to current time, reading
! in new monthly data if necessary
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------

    use commap
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
#if ( defined SPMD )
    use mpishorthand
#endif

    implicit none

#include <comctl.h>
!-----------------------------------------------------------------------
!
! Local workspace
!
    integer ntmp               ! temporary
    integer j,k                ! indices
    integer :: yr, mon, day    ! components of a date
    integer :: ncdate          ! current date in integer format [yyyymmdd]
    integer :: ncsec           ! current time of day [seconds]
    real(r8) :: calday         ! current calendar day
    real(r8) fact1, fact2      ! time interpolation factors
    real(r8) deltat            ! time (days) between interpolating ozone data

!-----------------------------------------------------------------------
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
    if ( masterproc) then
!
! If model time is past current forward timeslice, obtain the next
! timeslice for time interpolation.  Messy logic is for
! interpolation between December and January (np1.eq.1).
!
       calday = get_curr_calday()
       if ( is_perpetual() ) then
          call get_perp_date(yr, mon, day, ncsec)
       else
          call get_curr_date(yr, mon, day, ncsec)
       end if
       ncdate = yr*10000 + mon*100 + day

       if (calday > cdaytrp .and.  .not. (np1 == 1 .and. calday > cdaytrm)) then
          np1 = mod(np1,12) + 1
          if (np1 > ptrtim) then
             write(6,*)'CHEMINT: Attempt to access bad month'
             call endrun
          end if
          cdaytrm = cdaytrp
          call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
          if (np1 == 1 .or. calday <= cdaytrp) then
             ntmp = nm
             nm   = np
             np   = ntmp
             do j=1,plat
                do k=1,plev
                   tch4m  (j,k,np) = tch4i  (j,k,np1)
                   tn2om  (j,k,np) = tn2oi  (j,k,np1)
                   tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                   tcfc12m(j,k,np) = tcfc12i(j,k,np1)
                end do
             end do
          else
             write(6,*)'CHEMINT: Input data for date',date_tr(np1), &
                  ' sec ',sec_tr(np1), 'does not exceed model date', &
                  ncdate,' sec ',ncsec,' Stopping.'
             call endrun
          end if
       end if
!
! Determine factors for time interpolation.
!
       if (np1 == 1) then     ! Dec-Jan interpolation
          deltat = cdaytrp + 365. - cdaytrm
          if (calday > cdaytrp) then ! We're in December
             fact1 = (cdaytrp + 365. - calday)/deltat
             fact2 = (calday - cdaytrm)/deltat
          else                ! We're in January
             fact1 = (cdaytrp - calday)/deltat
             fact2 = (calday + 365. - cdaytrm)/deltat
          end if
       else                   ! Non Dec-Jan interpolation
          deltat = cdaytrp - cdaytrm
          fact1 = (cdaytrp - calday)/deltat
          fact2 = (calday - cdaytrm)/deltat
       end if
!
! Check sanity of time interpolation factors to within 32-bit roundoff.
!
       if (abs(fact1+fact2-1.) > 1.e-6 .or. &
            fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
            fact2 > 1.000001 .or. fact2 < -1.e-6) then
          write(6,*)'CHEMINT: Bad fact1 and/or fact2=',fact1,fact2
          call endrun
       end if
!
! Do time interpolation
!
       do j=1,plat
          do k=1,plev
             tch4(j,k)   = tch4m  (j,k,nm)*fact1 + tch4m  (j,k,np)*fact2
             tn2o(j,k)   = tn2om  (j,k,nm)*fact1 + tn2om  (j,k,np)*fact2
             tcfc11(j,k) = tcfc11m(j,k,nm)*fact1 + tcfc11m(j,k,np)*fact2
             tcfc12(j,k)  =tcfc12m(j,k,nm)*fact1 + tcfc12m(j,k,np)*fact2
          end do
       end do
    end if			! end of if-masterproc
#if ( defined SPMD )
    call mpibcast (tch4  , plev*plat,mpir8,0,mpicom)
    call mpibcast (tn2o  , plev*plat,mpir8,0,mpicom)
    call mpibcast (tcfc11, plev*plat,mpir8,0,mpicom)
    call mpibcast (tcfc12, plev*plat,mpir8,0,mpicom)
#endif

    return
  end subroutine chem_timestep_init

!===============================================================================
  subroutine chem_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial mass mixing ratios of CH4, N2O, CFC11 and CFC12.
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name

    real(r8), intent(out) :: q(plon,plev,plat)   !  mass mixing ratio
!-----------------------------------------------------------------------

    if ( name == 'N2O' ) then
       q = rmwn2o * n2ovmr
       return
    else if ( name == 'CH4' ) then
       q = rmwch4 * ch4vmr
       return
    else if ( name == 'CFC11' ) then
       q = rmwf11 * f11vmr
       return
    else if ( name == 'CFC12' ) then
       q = rmwf12 * f12vmr
       return
    end if

  end subroutine chem_init_cnst

end module chemistry
