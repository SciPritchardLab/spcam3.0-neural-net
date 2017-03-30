subroutine mkrgrid (cvsid, file1, file2, cmdline)
!
! $Id: mkrgrid.f90,v 1.2 2000/05/18 19:35:41 rosinski Exp $
!
  use constants
  use control
  use gridspecs

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  character*(*) cvsid
  character*(*) file1, file2
  character*(*) cmdline
!
! Local workspace
!
  character*(nf_max_name) :: name
  character*(nf_max_name) :: attname
  character, allocatable :: history(:) ! history attribute
  character*80 :: text

  integer cmdlen, hislen, totlen       ! character array lengths
  integer ncidi, ncido
  integer ret                 ! return code
  integer ngatts              ! number of global attributes
  integer nvars               ! number of variables
  integer dimid               ! returned dimension id
  integer londimid
  integer latdimid
  integer levdimid
  integer ilevdimid
  integer nvdims
  integer plon, nlat, nlev, nlevp
  integer ntime               ! size of unlimited dimension (if present)
  integer i, irow, j, n
  integer length, size
  integer xtype, xtype_rlon
  integer ndims
  integer natts
  integer nlonid, lonid, rlonid
  integer wnummaxid
  integer unlimdimid
  integer vardids(nf_max_var_dims) ! variable dimension id's

  logical itexists

  real(r8), allocatable :: alon(:)
  real(r8), allocatable :: rlon(:,:)
!
! Namelist values cannot be dynamically allocated
!
  namelist /reduced/ nlon, wnummax

  write(6,*) 'This is mkrgrid with cvsid:', cvsid

  inquire (file=file1, exist=itexists)
  if (.not.itexists) then
    write(6,*)'Unable to find input file1 =',file1
    stop 99
  end if

  call t_initializef()
  call t_startf('total')
!
! Open input and output netcdf files
!
  call wrap_open (file1, NF_NOWRITE, ncidi)
  call wrap_create (file2, NF_CLOBBER, ncido)

  if (reverse) then
    write(6,*)'Reduced -> full'
  else
    write(6,*)'Full -> reduced'
  end if
!
! Determine interpolation type
!
  if (default_interp == 'fourier' .or. default_interp == 'linear' .or. &
      default_interp == 'cubic') then

    write(6,*) 'Using ', default_interp, ' for grid interpolation'
    do i=1,monosiz
      write(6,*)'field ',monolist(i),' will be cubic monotonic'
    end do

  else

    write(6,*)default_interp, ' is not a valid interpolation type'
    stop 99
  end if
!
! Copy dimension and attribute information from input file to output file
!
  call wrap_inq (ncidi, ndims, nvars, ngatts, unlimdimid)

  ret = nf_inq_dimlen (ncidi, unlimdimid, ntime)
  if (ret /= NF_NOERR) then
    ntime = -1
    write(6,*)'INFO: Input file has no unlimited dimension'
  end if

  do n=1,ndims
    call wrap_inq_dim (ncidi, n, name, size)
    if (n == unlimdimid) then
      call wrap_def_dim (ncido, name, NF_UNLIMITED, dimid)
    else
      call wrap_def_dim (ncido, name, size, dimid)
    end if

    if (name == 'lon') then
      londimid = dimid
      call wrap_inq_dimlen (ncidi, londimid, plon)
    else if (name == 'lat') then
      latdimid = dimid
      call wrap_inq_dimlen (ncidi, latdimid, nlat)
    else if (name == 'lev') then
      levdimid = dimid
      call wrap_inq_dimlen (ncidi, levdimid, nlev)
    else if (name == 'ilev') then
      ilevdimid = dimid
      call wrap_inq_dimlen (ncidi, ilevdimid, nlevp)
    end if

    if (dimid /= n) then
      write(6,*)'Input dimid not equal to output dimid'
      stop 999
    end if
  end do
!
! Copy global attributes
!
  call wrap_inq (ncidi, ndims, nvars, ngatts, unlimdimid)
  do n=1,ngatts
    call wrap_inq_attname (ncidi, NF_GLOBAL, n, attname)
    if (attname == 'case') then
      text = ' '
      call wrap_get_att_text (ncidi, NF_GLOBAL, attname, text)
      write(6,*)'case =',trim(text)
    else if (attname == 'title') then
      text = ' '
      call wrap_get_att_text (ncidi, NF_GLOBAL, attname, text)
      write(6,*)'title =',trim(text)
    end if
    call wrap_copy_att (ncidi, NF_GLOBAL, attname, ncido, NF_GLOBAL)
  end do
!
! Add to or define history attribute.
!
  cmdlen = len_trim (cmdline)

  if (nf_inq_attlen (ncido, nf_global, 'history', hislen) == nf_noerr) then
    totlen = cmdlen + hislen
    allocate (history(totlen))
    if (nf_get_att_text (ncido, nf_global, 'history', history) /= nf_noerr) then
      write(6,*)'MKRGRID: bad attempt to get history attribute'
      stop 999
    end if
  else
    hislen = 0
    totlen = cmdlen
    allocate (history(totlen))
  end if

  do i=1,cmdlen
    history(hislen+i) = cmdline(i:i)
  end do
  
  call wrap_put_att_text (ncido, nf_global, 'history', totlen, history)

  if (nlat.gt.maxnlat) then
    write(6,*)'maxnlat too small: recompile with this parameter > ',nlat
    stop 999
  end if
!
! Allocate arrays dependent on dimension sizes just determined
!
  allocate (rlon(plon,nlat))
  allocate (alon(plon))
!
! Ensure that special case variables are there, and define appropriate output
! variables
!
  if (reverse) then
    call wrap_inq_varid (ncidi, 'nlon', nlonid)
    call wrap_inq_varid (ncidi, 'wnummax', wnummaxid)
    call wrap_inq_varid (ncidi, 'rlon', rlonid)

    call wrap_inq_var (ncidi, nlonid, name, xtype, ndims, vardids, natts)

    if (xtype /= NF_INT .or. ndims /= 1 .or. vardids(1) /= latdimid) then
      write(6,*)'Variable nlon on input tape is not as expected'
      stop 999
    end if

    call wrap_inq_var (ncidi, wnummaxid, name, xtype, ndims, vardids, natts)
    if (xtype /= NF_INT .or. ndims /= 1 .or. vardids(1) /= latdimid) then
      write(6,*)'Variable wnummax on input tape is not as expected'
      stop 999
    end if

    call wrap_inq_var (ncidi, rlonid, name, xtype_rlon, ndims, vardids, natts)
    if ((xtype_rlon /= nf_double .and. xtype_rlon /= nf_float) .or. &
         ndims /= 2 .or. vardids(1) /= londimid .or. vardids(2) /= latdimid) then
      write(6,*)'Variable rlon on input tape is not as expected'
      stop 999
    end if
!
! Get the variables from the input file
!
    call wrap_get_var_int (ncidi, nlonid, nlon)
    call wrap_get_var_int (ncidi, wnummaxid, wnummax)
    call wrap_get_var_double (ncidi, rlonid, rlon)
!
! Define special case variables on the output file
! nlon
!
    nvdims = 1
    vardids(1) = latdimid
    call wrap_def_var (ncido, 'old_nlon', nf_int, nvdims, vardids, nlonid)

    text = 'old number of longitudes'
    length = len('old number of longitudes')
    call wrap_put_att_text (ncido, nlonid, 'long_name', length, text)
!
! wnummax
!
    nvdims = 1
    vardids(1) = latdimid
    call wrap_def_var (ncido, 'old_wnummax', nf_int, nvdims, vardids, wnummaxid)

    text = 'old cutoff Fourier wavenumber'
    length = len('old cutoff Fourier wavenumber')
    call wrap_put_att_text (ncido, wnummaxid, 'long_name', length, text)
!
! rlon
!
    nvdims = 2
    vardids(1) = londimid
    vardids(2) = latdimid
    call wrap_def_var (ncido, 'old_rlon', xtype_rlon, nvdims, vardids, rlonid)

    text = 'old reduced longitude'
    length = len('old reduced longitude')
    call wrap_put_att_text (ncido, rlonid, 'long_name', length, text)

    text = 'degrees_east'
    length = len('degrees_east')
    call wrap_put_att_text (ncido, rlonid, 'units', length, text)
!
! lon
!
    nvdims = 1
    vardids(1) = londimid
    call wrap_def_var (ncido, 'lon', xtype_rlon, nvdims, vardids, lonid)

    text = 'longitude'
    length = len('longitude')
    call wrap_put_att_text (ncido, lonid, 'long_name', length, text)

    text = 'degrees_east'
    length = len('degrees_east')
    call wrap_put_att_text (ncido, lonid, 'units', length, text)

    do i=1,plon
      alon(i) = (i-1) * 360.0 / plon
    end do
!
! End define mode, write the variables, and go back to define mode
!
    if (nf_enddef (ncido) /= nf_noerr) then
      write(6,*)'bad call to nf_enddef'
      stop 999
    end if

    call wrap_put_var_int (ncido, nlonid, nlon)
    call wrap_put_var_int (ncido, wnummaxid, wnummax)
    call wrap_put_var_double (ncido, rlonid, rlon)
    call wrap_put_var_double (ncido, lonid, alon)

    if (nf_redef (ncido) /= nf_noerr) then
      write(6,*)'bad call to nf_redef'
      stop 999
    end if

  else
!
! Full -> reduced.
!
    nlon(:) = plon
    wnummax(:) = -1
      
    read (5,reduced)

    do irow=1,nlat/2
      if (wnummax(irow).eq.-1) then
        wnummax(irow) = nlon(irow)/2
      end if
    end do

    do irow=1,nlat/2
      wnummax(nlat-irow+1) = wnummax(irow)
    end do

    do j=1,nlat
      rlon(:,j) = fillvalue
      do i=1,nlon(j)
        rlon(i,j) = (i-1) * 360.0 / nlon(j)
      end do
    end do
!
! Make output "rlon" the same type as input "lon"
!
    nvdims = 2
    vardids(1) = londimid
    vardids(2) = latdimid
    call wrap_inq_varid (ncidi, 'lon', lonid)
    if (nf_inq_vartype (ncidi, lonid, xtype) /= nf_noerr) stop 999
    call wrap_def_var (ncido, 'rlon', xtype, nvdims, vardids, rlonid)

    nvdims = 1
    vardids(1) = latdimid
    call wrap_def_var (ncido, 'nlon', nf_int, nvdims, vardids, nlonid)

    nvdims = 1
    vardids(1) = latdimid
    call wrap_def_var (ncido, 'wnummax', nf_int, nvdims, vardids, wnummaxid)
!
! End define mode, write the variables, and go back to define mode
!
    if (nf_enddef (ncido) /= nf_noerr) then
      write(6,*)'bad call to nf_enddef'
      stop 999
    end if

    call wrap_put_var_double (ncido, rlonid, rlon)
    call wrap_put_var_int (ncido, nlonid, nlon)
    call wrap_put_var_int (ncido, wnummaxid, wnummax)

    if (nf_redef (ncido) /= nf_noerr) then
      write(6,*)'bad call to nf_redef'
      stop 999
    end if

  end if

  if (.not. silent) then
    if (reverse) then
      write(6,*)'Input reduced grid:'
    else
      write(6,*)'Output reduced grid:'
    end if
    do j=1,nlat
      write(6,*)'nlon(',j,')=',nlon(j)
    end do

    if (reverse) then
      write(6,*)'Input fourier wavenumber truncation'
    else
      write(6,*)'Output fourier wavenumber truncation'
    end if
    do irow=1,nlat/2
      write(6,*)'wnummax(',irow,')=',wnummax(irow)
    end do
  end if

  call rgconvert (ncidi, ncido, plon, nlat, nlev, &
                  nlevp, londimid, latdimid, levdimid, &
                  ilevdimid, nvars, unlimdimid, ntime)

  call wrap_close (ncidi)
  call wrap_close (ncido)

  call t_stopf('total')
  call t_prf(0)

  return
end
