program fmain

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  include 'netcdf.inc'
!
! Local workspace
!
  real(r8), parameter :: fillvalue = 1.d36
  real(r8), parameter :: filter_coefficient = 0.25D0

  character*128 :: topofile
  character(len=128) :: landmfile = ' '
  character*80 :: filename
  character*80 :: arg
  character*256 :: cmdline             ! input command line
  character, allocatable :: history(:) ! history attribute

  logical verbose                      ! Add print statements
  logical make_ross                    ! Make Ross ice shelf south of -79
  logical terrain_filter               ! Execute SJ Lin's terrain filter

  integer cmdlen, hislen, totlen       ! character array lengths
  integer fileid
  integer lonid, rlonid, latid
  integer sghid, phisid, flandid, nlonid, landmid
  integer londimid, latdimid, timdimid
  integer start(4), count(4)
  integer plon, nlat, ntim
  integer i, j
  integer ret
  integer nargs                        ! input arg
  integer n                            ! index loops thru input args

  integer , allocatable :: nlon(:)
  real(r8), allocatable :: mlatcnts(:)    ! model cell center latitudes
  real(r8), allocatable :: mloncnts(:,:)  ! model cell center longitudes
  real(r8), allocatable :: sgh(:,:)
  real(r8), allocatable :: phis(:,:)
  real(r8), allocatable :: fland(:,:)
  real(r8), allocatable :: landm(:,:)

  integer iargc
  external iargc
!
! Default settings before parsing argument list
!
  filename = ' '
  topofile = ' '
  verbose = .false.
  make_ross = .true.
  terrain_filter = .false.

! parse input arguments

  nargs = iargc()
  n = 1
  cmdline = char(10) // 'definesurf '
  do while (n .le. nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1

    select case (arg)
    case ('-t')
      call getarg (n, arg)
      n = n + 1
      topofile = arg
      cmdline = trim(cmdline) // ' -t ' // trim(topofile)
    case ('-v')
      verbose = .true.
      cmdline = trim(cmdline) // ' -v'
    case ('-l')
      call getarg (n, arg)
      n = n + 1
      landmfile = arg
      cmdline = trim(cmdline) // ' -l'
    case ('-r')
      make_ross = .false.
      cmdline = trim(cmdline) // ' -r'
    case ('-z')
      terrain_filter = .true.
      cmdline = trim(cmdline) // ' -z'
    case default
      if (filename .eq. ' ') then
        filename = arg
      else
        write (6,*) 'Argument ', arg,' is not known'
        call usage_exit (' ')
      end if
      cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if (filename == ' ') then
    call usage_exit ('Must enter an input file')
  end if

  if (topofile == ' ') then
    call usage_exit ('Must enter topofile name via -t arg')
  end if

  if (landmfile == ' ') then
    call usage_exit ('Must enter landmfile name via -l arg')
  end if

  ret = nf_open (trim(filename), nf_write, fileid) 
  if (ret /= nf_noerr) then
    write(6,*)nf_strerror(ret)
    write(6,*)'Unable to open input file ', trim(filename), ' for writing'
    stop 999
  end if

  call wrap_inq_dimid (fileid, 'lon', londimid)
  call wrap_inq_dimlen (fileid, londimid, plon)
  
  call wrap_inq_dimid (fileid, 'lat', latdimid)
  call wrap_inq_dimlen (fileid, latdimid, nlat)
  
  if (nf_inq_dimid (fileid, 'time', timdimid) == NF_NOERR) then
    ret = nf_inq_dimlen (fileid, timdimid, ntim)
    if (ntim/=1) then
      write(6,*)'Size of time dimension, if present, must be one'
      call endrun
    end if
  end if
!
! Get longitude and latitude arrays for model grid.  
! If reduced grid, 2-d variable containing lon values for each lat is called "rlon".
! First allocate space for dynamic arrays now that sizes are known
!
  allocate (nlon(nlat))
  allocate (mlatcnts(nlat))
  allocate (mloncnts(plon,nlat))

  if (nf_inq_varid (fileid, 'nlon', nlonid) == nf_noerr) then
    if (nf_get_var_int (fileid, nlonid, nlon) /= nf_noerr) then
      write(6,*)'nf_get_var_int() failed for nlon'
      call endrun
    end if
  else
    nlon(:) = plon
  end if
    
  do j=1,nlat 
    if (nlon(j)<1 .or. nlon(j)>plon) then
      write(6,*)'nlon(',j,')=',nlon(j),' is invalid.'
      write(6,*)'Must be between 1 and ',plon
      call endrun
    end if
  end do

  call wrap_inq_varid (fileid, 'lat', latid)
  call wrap_get_var8 (fileid, latid, mlatcnts)

  if (nf_inq_varid (fileid, 'lon', lonid) == nf_noerr) then
     call wrap_get_var8 (fileid, lonid, mloncnts(1,1))
     do j=2,nlat
        mloncnts(:,j) = mloncnts(:,1)
     end do
  else
     call wrap_inq_varid (fileid, 'rlon', rlonid)
     call wrap_get_var8 (fileid, rlonid, mloncnts)
  end if

  if (nf_inq_varid (fileid, 'FLAND', flandid) == nf_noerr) then
     ret = nf_redef (fileid)
     if (ret/=NF_NOERR) call handle_error (ret)
     ret = nf_rename_var (fileid, flandid, 'LANDFRAC');
     if (ret/=NF_NOERR) call handle_error (ret)
     ret = nf_enddef (fileid)
     if (ret/=NF_NOERR) call handle_error (ret)
  endif
!
! Check dimensions and allocate space for variables
!
  call chkdims (fileid, 'SGH'  , sghid  , londimid, latdimid, timdimid, verbose)
  call chkdims (fileid, 'PHIS' , phisid , londimid, latdimid, timdimid, verbose)
  call chkdims (fileid, 'LANDFRAC', flandid  , londimid, latdimid, timdimid, verbose)
  call chkdims (fileid, 'LANDM_COSLAT', landmid, londimid, latdimid, timdimid, verbose)

  allocate (sgh(plon,nlat))
  allocate (phis(plon,nlat))
  allocate (fland(plon,nlat))
  allocate (landm(plon,nlat))
!
! Create output attributes
!
  ret = nf_redef (fileid)
  if (ret/=NF_NOERR) call handle_error (ret)

  if (make_ross) then
    write (6,*) 'Extending Ross ice shelf south of -79 degrees'
    call wrap_put_att_text (fileid, nf_global, 'make_ross', 4, 'true')
  else
    write (6,*) 'Not doing anything special for Ross ice shelf'
    call wrap_put_att_text (fileid, nf_global, 'make_ross', 5, 'false')
  end if

  call wrap_put_att_text   (fileid, sghid, 'long_name', 28, 'orography standard deviation')
  call wrap_put_att_text   (fileid, sghid, 'units', 1, 'M')
  call wrap_put_att_double (fileid, sghid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_text   (fileid, sghid, 'from_hires', 4, 'true')

  call wrap_put_att_text   (fileid, phisid, 'long_name', 20, 'surface geopotential')
  call wrap_put_att_text   (fileid, phisid, 'units', 5, 'M2/S2')
  call wrap_put_att_double (fileid, phisid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_text   (fileid, phisid, 'from_hires', 4, 'true')

  call wrap_put_att_text   (fileid, flandid, 'long_name', 21, 'gridbox land fraction')
  call wrap_put_att_text   (fileid, flandid, 'units', 4, 'FRACTION')
  call wrap_put_att_double (fileid, flandid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_text   (fileid, flandid, 'from_hires', 4, 'true')

  call wrap_put_att_text   (fileid, landmid, 'long_name', 71, &
                          'land ocean transition mask: ocean (0), continent (1), transition (0-1)')
  call wrap_put_att_text   (fileid, landmid, 'units', 4, 'none')
  call wrap_put_att_double (fileid, landmid, '_FillValue', nf_double, 1, fillvalue)
  call wrap_put_att_text   (fileid, landmid, 'from_hires', 4, 'true')
!
! Add to or define history attribute.
!
  cmdlen = len_trim (cmdline)

  if (nf_inq_attlen (fileid, nf_global, 'history', hislen)  == nf_noerr) then
    totlen = cmdlen + hislen
    allocate (history(totlen))
    if (nf_get_att_text (fileid, nf_global, 'history', history) /= nf_noerr) then
      call endrun
    end if
  else
    hislen = 0
    totlen = cmdlen
    allocate (history(totlen))
  end if

  do i=1,cmdlen
    history(hislen+i) = cmdline(i:i)
  end do

  call wrap_put_att_text (fileid, nf_global, 'history', totlen, history)

  ret = nf_enddef (fileid)
  if (ret/=NF_NOERR) call handle_error (ret)
!
! Determine model topographic height and variance
!
  call sghphis (plon, nlat, nlon, mlatcnts, mloncnts, topofile, &
                verbose, sgh, phis, fland)
!
! Calculate LANDM field required by cloud water.  
!
!JR Replace original resolution-dependent calculation with interpolation.
!JR
!JR  call inimland (plon, nlat, nlon, mlatcnts, mloncnts, topofile, &
!JR                 verbose, make_ross, landm)
!
  call interplandm (plon, nlat, nlon, mlatcnts, mloncnts, &
                    landmfile, landm)
!
! Define land mask.  Set all non-land points to ocean (i.e. not sea ice).
! Also: zero SGH over ocean.  This is to make things consistent with the way
! CAM and CCM3 has been run in the past, but will cause differences as compared with
! CCM3.5.5 when interpolating.
!
  do j=1,nlat
    do i=1,nlon(j)
!
! Overwrite FLAND flag as land for Ross ice shelf
!
      if (make_ross .and. mlatcnts(j) < -79.) then
        fland(i,j) = 1.
      end if
    end do
!
! Fill region outside reduced grid with flag values
!
    do i=nlon(j)+1,plon
      sgh(i,j)   = fillvalue
      phis(i,j)  = fillvalue
      fland(i,j)   = fillvalue
      landm(i,j) = fillvalue
    end do
  end do

! Do the terrain filter.
! Note: not valid if a reduced grid is used.

  if (terrain_filter) then

     write(6,*) 'Terrain filtering'

     call sm2(plon, nlat, phis, plon/12, filter_coefficient)
     call sm2(plon, nlat, sgh,  plon/12, filter_coefficient)

  endif

  write(6,*) 'Writing surface quantities'

  start(:) = 1
  count(1) = plon
  count(2) = nlat
  count(3:) = 1

  call wrap_put_vara8 (fileid, sghid,   start, count, sgh)
  call wrap_put_vara8 (fileid, phisid , start, count, phis)
  call wrap_put_vara8 (fileid, flandid, start, count, fland)
  call wrap_put_vara8 (fileid, landmid, start, count, landm)

  if (nf_close (fileid) == nf_noerr) then
    write(6,*) 'Successfully defined surface quantities on ', trim(filename)
  else
    write(6,*) 'ERROR CLOSING NETCDF FILE ',trim(filename)
  end if

  deallocate (nlon)
  deallocate (mlatcnts)
  deallocate (mloncnts)
  deallocate (sgh)
  deallocate (phis)
  deallocate (fland)
  deallocate (landm)

  stop 0
end program fmain

subroutine usage_exit (arg)
  implicit none
  character*(*) arg
  
  if (arg /= ' ') write (6,*) arg
  write (6,*) 'Usage: definesurf -t topofile -l landmfile [-v] [-r] [-z] filename'
  write (6,*) '       -v   verbose mode'
  write (6,*) '       -r   Do *not* extend Ross Ice Shelf as land ice'
  write (6,*) '       -z   use terrain filter (not a valid option for reduced grid)'
  stop 999
end subroutine usage_exit
