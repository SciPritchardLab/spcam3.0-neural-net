    module test_sst_mod

    use shr_kind_mod, only: r8 => shr_kind_r8
    use comsrf, only: plevmx, oro, ts, tssub, snowh, sicthk, icefrac
    use rgrid, only: nlon
    use physconst, only: tmelt

    implicit none

    private    ! By default make data private to this module

    public initial, geticinfo, test


    contains

    subroutine geticinfo
!
! Get ORO and nlon fields from IC file
!
    use pmgrid, only: plon, plat
#include <comlun.h>
#include <comctl.h>
#include <netcdf.inc>
    integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
    data strt2d/3*1/          ! Only index 2 will ever change
    integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
    data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
    integer oroid             ! Variable ID's
    integer nlonid            ! Variable ID's
    integer j
    integer ret

    write(6,*) 'Read in ORO and nlon from IC dataset'
    print *, 'Snag ORO var-od'
    ret = nf_inq_varid (ncid_ini, 'ORO', oroid)
    if (ret/=NF_NOERR)then
      print *, 'Error: on inq on ORO'
      call endrun
    end if

    if ( .not. fullgrid )then
      print *, 'Snag nlon var-od'
      ret = nf_inq_varid (ncid_ini, 'nlon', nlonid)
      if (ret/=NF_NOERR)then
        print *, 'Error: on inq on nlon'
        call endrun
      end if
      print *, 'Get nlon'
      ret = nf_get_var_int (ncid_ini, nlonid, nlon)
      if (ret/=NF_NOERR)then
        print *, 'Error: on get of nlon'
        call endrun
      end if
    else
      nlon(:) = plon
    end if
    print *, 'Get ORO'
    do j=1,plat
       strt2d(2) = j
       ret = nf_get_vara_double(ncid_ini, oroid, strt2d, cnt2d, oro(1,j))
       if (ret/=NF_NOERR)then
         write(6,*)nf_strerror(ret)
         print *, 'Error: on get of ORO'
         call endrun
       end if
    end do

    end subroutine geticinfo

    subroutine test( j, icef )
    integer :: j
    real(r8) :: icef
!
! Test for consistencies between variables
!
    integer i

    do i = 1, nlon(j)
      icef = icef + icefrac(i,j)
      if( oro(i,j) == 2 ) then
        if ( icefrac(i,j) /= 1.0 ) then
          write(6,*) 'ERROR:: icefrac and oro inconstitent'
          call endrun
        end if
     else
        if ( icefrac(i,j) /= 0.0 ) then
          write(6,*) 'ERROR:: icefrac and oro inconstitent'
          call endrun
        end if
     end if
   end do
   end subroutine test

    subroutine initial
!
! Initialize surface data
!
    use pmgrid, only: plon, plat, beglat, endlat
    integer j, k

    do j=beglat,endlat
      snowh(1:nlon(j),j)  = 0.0
      sicthk(1:nlon(j),j) = 0.0
      ts(1:nlon(j),j)     = tmelt
      do k = 1, plevmx
        tssub(1:nlon(j),k,j)= ts(1:nlon(j),j)
      end do
    end do

    end subroutine initial

  end module test_sst_mod

  program test_sst
!
!	test_sst.F90
!
!	Unit-tester for the SST interpolation routines.
!
!	$Id: test_sst.F90,v 1.5.2.1 2002/06/15 13:49:57 erik Exp $
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ioFileMod, only: getfil
  use pmgrid, only: plond, beglat, endlat
  use ppgrid, only: begchunk, endchunk
  use rgrid, only: nlon
  use comsrf, only: initialize_comsrf, oro, ts, tssub, snowh, sicthk, icefrac
  use test_sst_mod, only: initial, geticinfo, test
  use sst_data, only: sstan, sstini, sstint
  use physconst, only: tmelt
  use phys_grid, only: get_ncols_p, phys_grid_init
  implicit none
#include <commss.h>
#include <comctl.h>
#include <comtim.h>
#include <comlun.h>
  character(len=256) :: locfn    ! netcdf local filename to open 
  integer, parameter :: ntypes = 3
  integer, parameter :: dates(ntypes) = (/ 901, 901, 19900901  /)
  logical, parameter :: cycle(ntypes) = (/ .true., .true., .false. /)
  logical, parameter :: grid(ntypes) = (/ .true., .false., .true. /)
  character(len=256), parameter :: dataset(ntypes) = (/ &
                          "T42M5079.nc                            ",  &
                          "T42M5079.r1up.nc                       ",  &
                          "pcmdi_sst_ccm_bcT42_1977_1998.012000.nc" /)
  character(len=256), parameter :: icdataset(ntypes) = (/ &
                          "SEP1.T42L26.112000.nc                ",  &
                          "SEP1.T42L26.112000.r1up.definesurf.nc",  &
                          "SEP1.T42L26.112000.nc                " /)
  real(r8) :: icef !  Ice fraction
  integer :: itype, lchnk, i, lat, ncol
  integer, parameter :: nyears = 17   ! Number of years to run over
  logical :: orohires = .true.
!
! Initialize comsrf data
!
  call phys_grid_init( plond,  -1 )
  call initialize_comsrf
!
! Loop over the types of tests to do
!
  do itype = 1, ntypes
!
! Set needed variables to use
!
    itsst    = 1
    ncid_sst = 1
    dtime    = 3600.0*12   ! 12-hour time-steps for testing
    ndbase   = 0
    nsbase   = 0
    nbsec    = 0
    sstcyc   = cycle(itype)
    nbdate   = dates(itype)
    fullgrid = grid(itype)
    bndtvs   = "/fs/cgd/csm/inputdata/atm/ccm3/" // dataset(itype)
    ncdata   = "/fs/cgd/csm/inputdata/atm/ccm3/" // icdataset(itype)
    nstep    = 0
!
! Initialization, set time info, get file and open it
!
    call calendr(nstep   ,dtime   ,ndbase  ,nsbase  ,nbdate  , &
                 nbsec   ,ndcur   ,nscur   ,ncdate  ,ncsec   , &
                 calday  )
    call getfil(bndtvs, locfn)
    call wrap_open(locfn, 0, ncid_sst)
    write(6,*)'NCOPN returns id ',ncid_sst,' for file ',trim(locfn)
    call getfil(ncdata, locfn)
    call wrap_open(locfn, 0, ncid_ini)
    write(6,*)'NCOPN returns id ',ncid_ini,' for file ',trim(locfn)
    call geticinfo
    call wrap_close(ncid_ini)
  
    call sstini( orohires )
    call sstint
    call initial
    where( oro(:,:) == 2 )
      icefrac = 1.0
    elsewhere
      icefrac = 0.0
    end where
!
! Initialize surface and sub-surface temperatures, set new
! new sea ice concentrations and compute longwave up over non-land
!
    do lat=beglat,endlat
       do i=1,nlon(lat)
          if (nint(oro(i,lat)) /= 1) then !non-land
             if (tssub(i,1,lat) <= tmelt ) then
                oro(i,lat)     = 2.
                icefrac(i,lat) = 1
             end if
          end if
        end do
    end do

!
! Loop over the number of time-steps in the desired years to run
!
    do nstep = 1, 365*nyears*3600*24/int(dtime)
       call calendr(nstep   ,dtime   ,ndbase  ,nsbase  ,nbdate  , &
                    nbsec   ,ndcur   ,nscur   ,ncdate  ,ncsec   , &
                    calday  )
       call sstint
       icef = 0.0
       do lchnk=begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call sstan (lchnk, ncol, oro(1,lchnk), ts(1,lchnk), &
                        tssub(1,1,lchnk), snowh(1,lchnk), &
                        sicthk(1,lchnk))
          call test( lchnk, icef )
       end do
       icef = icef / sum(nlon(:))
       if ( mod(ncdate,30) == 0 ) write(6,*) 'Date', ncdate, 'Icefraction: ', icef
    end do
    call wrap_close(ncid_sst)
  end do

  end program test_sst
