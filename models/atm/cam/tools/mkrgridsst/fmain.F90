  program fmain
  implicit none
  include 'netcdf.inc'
!
! Local workspace
!
  logical verbose               ! Add print statements
  data verbose /.false./

  integer fileid
  integer sstid
  integer lonid, latid, timid
  integer nlon, nlat, ntim
  integer ret
  integer nargs
  integer i, j, n
  integer start(3)
  integer count(3)
  integer num, num2
  integer s_nlon(10000)
  integer s_wnummax(10000)
  namelist /reduced/ s_nlon, s_wnummax

  character*80 arg, filename

  real*8, allocatable :: sst(:,:), sst_reduced(:,:)
  real*8, allocatable :: ice(:), ice_reduced(:)
  real*8, parameter :: tsice = -1.8d0
!
! Externals
!
  integer iargc
  external iargc

  s_nlon(:) = -999

  nargs = iargc()
  n = 1
  filename = ' '
  do while (n .le. nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1
    if (arg .eq. '-v') then
      verbose = .true.
    else
      if (filename .eq. ' ') then
        filename = arg
      else
        write (6,*) 'Argument ', arg,' is not known'
        stop 999
      end if
    end if
  end do

  if (filename .eq. ' ') call usage_exit (' ')

  ret = nf_open (filename, NF_WRITE, fileid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'file ', filename, ' cannot be opened for writing'
    stop 999
  end if

  call wrap_inq_dimid (fileid, 'lon', lonid)
  call wrap_inq_dimlen (fileid, lonid, nlon)
  
  call wrap_inq_dimid (fileid, 'lat', latid)
  call wrap_inq_dimlen (fileid, latid, nlat)
  
  ret = nf_inq_dimid (fileid, 'time', timid)
  if (ret.eq.NF_NOERR) then
    ret = nf_inq_dimlen (fileid, timid, ntim)
  end if

  call chkdims (fileid, 'SST', sstid, lonid, latid, verbose)

  allocate (sst(nlon,nlat), sst_reduced(nlon,nlat))
  allocate (ice(nlon), ice_reduced(nlon))
  read (5, reduced)

  start(1) = 1
  start(2) = 1
  count(1) = nlon
  count(2) = nlat
  count(3) = 1
  do n=1,ntim
    start(3) = n
    call wrap_get_vara8 (fileid, sstid, start, count, sst)
    do j=1,nlat
      do i=1,nlon
        ice(i) = 0.
        if (sst(i,j).lt.tsice+0.0001) ice(i) = tsice
      end do
      call lininterp (nlon, nlon, s_nlon(j), 1, ice, ice_reduced)
      call lininterp (nlon, nlon, s_nlon(j), 1, sst(1,j), sst_reduced(1,j))
      do i=1,s_nlon(j)
        if (ice_reduced(i) .lt. 0.5*tsice) sst_reduced(i,j) = tsice
      end do

      write(6,*)'Stats time n=',n
      num = 0
      do i=1,nlon
        if (sst(i,j).lt.tsice+0.0001) num = num + 1
      end do
      num2 = 0
      do i=1,s_nlon(j)
        if (sst_reduced(i,j).eq.tsice) num2 = num2 + 1
      end do
      write(6,*)num,' full ice pts. ',num2,' reduced out of ', s_nlon(j)
    end do
    call wrap_put_vara8 (fileid, sstid, start, count, sst_reduced)
  end do

  deallocate (sst, sst_reduced)
  deallocate (ice, ice_reduced)

  if (nf_close (fileid).ne.NF_NOERR) then
    stop 999
  end if

  stop 0
end

subroutine usage_exit (arg)
  implicit none
  character*(*) arg

  if (arg.ne.' ') write (6,*) arg
  write (6,*) 'Usage: mkrgridsst [-v] filename < namelist'
  stop 999
end subroutine usage_exit
