subroutine rgconvert (ncidi, ncido, plon, nlat, nlev, &
                      nlevp, londimid, latdimid, levdimid, &
                      ilevdimid, nvars, unlimdimid, ntime)
!
! $Id: rgconvert.f90,v 1.2.22.1 2002/06/15 13:50:07 erik Exp $
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use control

  implicit none

  include 'netcdf.inc'
!
! Input arguments
!
  integer ncidi, ncido            ! input and output netcdf file ids
  integer plon, nlat, nlev, nlevp ! spatial dimension sizes
  integer londimid, latdimid, levdimid, ilevdimid ! spatial dimension ids
  integer nvars               ! number of variables
  integer unlimdimid          ! unlimited dimension id
  integer ntime               ! size of unlimited dimension (if present)
!
! Local workspace
!
  character*8 :: shape
  character*(nf_max_name) :: name
  character*(nf_max_name) :: attname

  integer natts               ! number of attributes for a given variable
  integer nvdims              ! number of dimensions for this variable
  integer vardids(nf_max_var_dims) ! variable dimension id's
  integer id                  ! variable id
  integer j, k       ! spatial indices
  integer numlev ! number of levels
  integer dimlen ! dimension length
  integer n      ! index
  integer t      ! index over unlimited dimension
  integer v                    ! loop index over variable id
  integer vo         ! returned variable id on output file
  integer xtype                ! variable type (netcdf)
  integer tpos               ! position of unlimited dimension
  integer ncp_unlim
  integer ncp_nounlim
  integer nintp_unlim
  integer nintp_nounlim
  integer start(nf_max_var_dims)
  integer count(nf_max_var_dims)
  integer size
  integer :: indx_cp_unlim(nvars)
  integer :: indx_cp_nounlim(nvars)
  integer :: indx_intp_unlim(nvars)
  integer :: indx_intp_nounlim(nvars)

  logical xzy, xyz, xy
  logical has_unlim
  logical copy
!
! Allocatables
!
  character, allocatable :: cbuf(:)

  integer, allocatable :: ibuf(:)

  real(r8), allocatable :: buf(:)
  real(r8), allocatable :: arrxyz(:,:,:)
  real(r8), allocatable :: arrxzy(:,:,:)

  type varspecs
    character*(nf_max_name) :: name
    character*8 :: shape

    integer :: id
    integer :: size
    integer :: xtype
    integer :: numlev
    integer :: tpos
    integer :: count(nf_max_var_dims)
  end type varspecs

  type (varspecs) :: var(nvars)

  logical is_special_case
  external is_special_case
!
! Initialize indices to invalid values
!
  indx_cp_unlim(:) = -1
  indx_cp_nounlim(:) = -1
  indx_intp_unlim(:) = -1
  indx_intp_nounlim(:) = -1

  ncp_unlim = 0
  nintp_unlim = 0
  ncp_nounlim = 0
  nintp_nounlim = 0

  call t_startf('define_output')
  do v=1,nvars
    copy = .true.

    call wrap_inq_var (ncidi, v, name, xtype, nvdims, vardids, natts)
!
! Skip any special case variables: they have already been dealt with
!
    if (is_special_case (name, reverse)) then
      cycle
    end if
!
! Normal case variables
!
    call wrap_def_var (ncido, name, xtype, nvdims, vardids, vo)
!
! Copy attributes
!
    do n=1,natts
      call wrap_inq_attname (ncidi, v, n, attname)
      call wrap_copy_att (ncidi, v, attname, ncido, vo)
    end do

    xy  = .false.
    xyz = .false.
    xzy = .false.

    if (nvdims > 1) then
       xy = vardids(1) == londimid .and. &
            vardids(2) == latdimid .and. &
            (vardids(3) /= levdimid .and. vardids(3) /= ilevdimid)

       xyz = vardids(1) == londimid .and. &
            vardids(2) == latdimid .and. &
            (vardids(3) == levdimid .or. vardids(3) == ilevdimid)
      
       xzy = vardids(1) == londimid .and. &
            (vardids(2) == levdimid .or. vardids(2) == ilevdimid) .and. &
            vardids(3) == latdimid

       if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then
!
! Interpolated variables must be of a floating point type and have dimensions
! xy, xyz, or xzy.
!
          if (xy .or. xyz .or. xzy) then
             copy = .false.
             if (xy .and. vardids(3) == unlimdimid .or. &
                  (xyz .or. xzy) .and. vardids(4) == unlimdimid) then
                
                nintp_unlim = nintp_unlim + 1
                indx_intp_unlim(nintp_unlim) = v
             else
                nintp_nounlim = nintp_nounlim + 1
                indx_intp_nounlim(nintp_nounlim) = v
             end if
          end if
       end if
    end if
!
! Variables to be copied will not be interpolated.  Determine which do and
! do not have an unlimited dimension
!
    if (copy) then
      has_unlim = .false.
      do n=1,nvdims
        if (vardids(n) == unlimdimid) then
          has_unlim = .true.
          exit
        end if
      end do

      if (has_unlim) then
        ncp_unlim = ncp_unlim + 1
        indx_cp_unlim(ncp_unlim) = v
        if (verbose) then
           write(6,*)trim(name),' added to ncp_unlim=', ncp_unlim
        end if
      else
        ncp_nounlim = ncp_nounlim + 1
        indx_cp_nounlim(ncp_nounlim) = v
        if (verbose) then
           write(6,*)trim(name),' added to ncp_nounlim=', ncp_nounlim
        end if
      end if
    end if
!
! Copy useful information for copying or interpolating into the struct
!
    var(v)%name   = name
    var(v)%id     = vo
    var(v)%xtype  = xtype
    var(v)%shape  = 'unknown'
    if (xy) var(v)%shape = 'xy'
    if (xyz) var(v)%shape = 'xyz'
    if (xzy) var(v)%shape = 'xzy'
!
! Determine sizes, counts, and position of unlimited dimension for this 
! variable.  Size value does not include number of time samples.
!
    var(v)%size     = 1    ! inito to size which can be multiplied
    var(v)%numlev   = 1    ! init to single level field
    var(v)%tpos     = -1   ! init to invalid value
    var(v)%count(:) = -1   ! init to invalid value

    do n=1,nvdims
      if (vardids(n) == unlimdimid) then
        var(v)%tpos = n
      else if (vardids(n) == levdimid) then
        var(v)%numlev   = nlev
      else if (vardids(n) == ilevdimid) then
        var(v)%numlev   = nlevp
      end if

      if (vardids(n) /= unlimdimid) then
        call wrap_inq_dimlen (ncidi, vardids(n), dimlen)
        var(v)%count(n) = dimlen
        var(v)%size = var(v)%size * dimlen
      end if
    end do
  end do                           ! loop over input variables

  if (nf_enddef (ncido) /= NF_NOERR) stop 999
  call t_stopf('define_output')
!
! End out output file specification.  Now copy or interpolate data to the
! output file.  First: copy data which have no unlimited dimension
!
  do n=1,ncp_nounlim
    v     = indx_cp_nounlim(n)

    name  = var(v)%name
    id    = var(v)%id
    size  = var(v)%size
    xtype = var(v)%xtype

    if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then
      if (verbose) then
         write(6,*)'copying floating point var ',trim(name)
      end if
      allocate (buf(size))
      call wrap_get_var_double (ncidi, v, buf)
      call wrap_put_var_double (ncido, id, buf)
      deallocate (buf)
    else if (xtype == NF_INT) then
      if (verbose) then
         write(6,*)'copying integer var ',trim(name)
      end if
      allocate (ibuf(size))
      call wrap_get_var_int (ncidi, v, ibuf)
      call wrap_put_var_int (ncido, id, ibuf)
      deallocate (ibuf)
    else if (xtype == NF_CHAR) then
      if (verbose) then
         write(6,*)'copying character var ',trim(name)
      end if
      allocate (cbuf(size))
      call wrap_get_var_text (ncidi, v, cbuf)
      call wrap_put_var_text (ncido, id, cbuf)
      deallocate (cbuf)
    else
      write(6,*)'Unknown type for variable ',var(v)%name
      stop 999
    end if
  end do
!
! Interpolate data which have no unlimited dimension
!
  do n=1,nintp_nounlim
    v      = indx_intp_nounlim(n)

    name   = var(v)%name
    id     = var(v)%id
    size   = var(v)%size
    xtype  = var(v)%xtype
    shape  = var(v)%shape
    numlev = var(v)%numlev

    allocate (arrxzy(plon,numlev,nlat))

    if (shape == 'xy' .or. shape == 'xyz') then

      allocate (arrxyz(plon,nlat,numlev))
      call wrap_get_var_double (ncidi, v, arrxyz)

      do j=1,nlat
        do k=1,numlev
          arrxzy(:,k,j) = arrxyz(:,j,k)
        end do
      end do

      call interp_driver (name, arrxzy, plon, nlat, numlev)

      do j=1,nlat
        do k=1,numlev
          arrxyz(:,j,k) = arrxzy(:,k,j)
        end do
      end do

      call wrap_put_var_double (ncido, id, arrxyz)
      deallocate (arrxyz)

    else if (shape == 'xzy') then

      call wrap_get_var_double (ncidi, v, arrxzy)
      call interp_driver (name, arrxzy, plon, nlat, numlev)
      call wrap_put_var_double (ncido, id, arrxzy)

    else

      write(6,*)'Unknown shape=',shape,' for variable ',var(v)%name
      stop 999

    end if
    deallocate (arrxzy)
  end do
!
! Now loop over the unlimited dimension.  First do copies
!
  do t=1,ntime
    if (.not.silent) then
      write(6,*)'Starting time sample ',t
    end if

    do n=1,ncp_unlim
      v           = indx_cp_unlim(n)

      id          = var(v)%id
      size        = var(v)%size
      xtype       = var(v)%xtype
      tpos        = var(v)%tpos
      start(:)    = 1
      start(tpos) = t
      count(:)    = var(v)%count(:)
      count(tpos) = 1

      if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then
        allocate (buf(size))
        call wrap_get_vara_double (ncidi, v, start, count, buf)
        call wrap_put_vara_double (ncido, id, start, count, buf)
        deallocate (buf)
      else if (xtype == NF_INT) then
        allocate (ibuf(size))
        call wrap_get_vara_int (ncidi, v, start, count, ibuf)
        call wrap_put_vara_int (ncido, id, start, count, ibuf)
        deallocate (ibuf)
      else if (xtype == NF_CHAR) then
        allocate (cbuf(size))
        call wrap_get_vara_text (ncidi, v, start, count, cbuf)
        call wrap_put_vara_text (ncido, id, start, count, cbuf)
        deallocate (cbuf)
      else
        write(6,*)'Unknown type for variable ',var(v)%name
        stop 999
      end if
    end do
!
! Now the data with an unlimited dimension which need to be interpolated
!
    do n=1,nintp_unlim
      v           = indx_intp_unlim(n)

      name        = var(v)%name
      id          = var(v)%id
      size        = var(v)%size
      xtype       = var(v)%xtype
      shape       = var(v)%shape
      numlev      = var(v)%numlev
      tpos        = var(v)%tpos
      start(:)    = 1
      start(tpos) = t
      count(:)    = var(v)%count(:)
      count(tpos) = 1

      allocate (arrxzy(plon,numlev,nlat))

      if (shape == 'xy' .or. shape == 'xyz') then

        allocate (arrxyz(plon,nlat,numlev))
        call wrap_get_vara_double (ncidi, v, start, count, arrxyz)

        call t_startf('transpose_copy1')
        do j=1,nlat
          do k=1,numlev
            arrxzy(:,k,j) = arrxyz(:,j,k)
          end do
        end do
        call t_stopf('transpose_copy1')

        call interp_driver (name, arrxzy, plon, nlat, numlev)

        call t_startf('transpose_copy2')
        do j=1,nlat
          do k=1,numlev
            arrxyz(:,j,k) = arrxzy(:,k,j)
          end do
        end do
        call t_stopf('transpose_copy2')

        call wrap_put_vara_double (ncido, id, start, count, arrxyz)
        deallocate (arrxyz)

      else if (shape == 'xzy') then

        call wrap_get_vara_double (ncidi, v, start, count, arrxzy)
        call interp_driver (name, arrxzy, plon, nlat, numlev)
        call wrap_put_vara_double (ncido, id, start, count, arrxzy)

      else

        write(6,*)'Unknown shape=',shape,' for variable ',var(v)%name
        stop 999

      end if
      deallocate (arrxzy)
    end do
  end do    
  return
end subroutine

logical function is_special_case (name, reverse)

  implicit none

  character*(*) name
  logical reverse

  logical bad

  if (reverse) then

    is_special_case = name == 'rlon' .or. name == 'nlon' .or. name == 'wnummax'
    bad = name == 'lon' .or. name == 'old_rlon' .or. name == 'old_nlon' .or. &
          name == 'old_wnummax'

  else

    bad = name == 'rlon'
    is_special_case = name == 'nlon' .or. name == 'wnummax' .or. &
                      name == 'lon' .or.name == 'old_rlon' .or. &
                      name == 'old_nlon' .or. name == 'old_wnummax'
  end if

  if (bad) then
    write(6,*)'Variable ',name,' is inconsistent with conversion'
    stop 999
  end if

  return
end function is_special_case
