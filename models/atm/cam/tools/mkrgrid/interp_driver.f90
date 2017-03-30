subroutine interp_driver (name, arrxzy, plon, nlat, numlev)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use control
  use gridspecs
!
! Input arguments
!
  character*(*) :: name

  integer :: plon, nlat, numlev

  real(r8) :: arrxzy(plon,numlev,nlat)
!
! Local workspace
!
  integer :: lenin
  integer :: lenot
  integer :: j, irow
  integer :: endchar
  integer :: nfull
  integer :: nred

  logical :: mono

  real(r8), parameter :: tsice = -1.8d0    ! sea ice flag value
  real(r8) :: arrxz(plon,numlev)
  real(r8) :: ice(plon)
  real(r8) :: ice_reduced(plon)

  mono = .false.
  do i=1,monosiz
    if (monolist(i).eq.name) mono = .true.
  end do

  if (.not. silent) then
    endchar = 1
    do i=1,len(name)
      if (name(i:i) == ' ') then
        endchar = i
        exit
      end if
    end do

    if (mono) then
      write(6,*)'Interpolating ',name(1:endchar),' cubic monotonic'
    else 
      write(6,*)'Interpolating ',name(1:endchar), default_interp
    end if
  end if
        
  do j=1,nlat
    irow = j
    if (j > nlat/2) irow = nlat - j + 1
!
! If converting *from* reduced grid, pull number of longitudes off of
! input data record.  Compute wnummax as this number divided by 2
!
    if (reverse) then
      lenin = nlon(j)
      lenot = plon
    else
      lenin = plon
      lenot = nlon(j)
    end if
    
    call t_startf('interpolation')

    if (default_interp == 'mono') then
      call cubinterp (plon, lenin, lenot, numlev, arrxzy(1,1,j), &
                      arrxz, mono)
      
    else if (default_interp == 'fourier') then
      
      if (verbose) then
        write(6,*)'irow=',irow,': keep max wnum=',wnummax(irow), ' of ', &
             lenin/2, ' input and ',lenot/2,' output'
      end if
      call fourier (plon, lenin, lenot, numlev, arrxzy(1,1,j), &
                    arrxz, wnummax(irow), name)
      
    else if (default_interp == 'linear') then

      call lininterp (plon, lenin, lenot, numlev, arrxzy(1,1,j), &
                      arrxz)
      
    else if (default_interp == 'cubic') then
      
      call cubinterp (plon, lenin, lenot, numlev, arrxzy(1,1,j), &
                      arrxz, mono)
      
    else
      
      write(6,*)default_interp, ' is not a valid interpolation type'
      stop 999
      
    end if
!
! Handle sea ice flag values if this is an sst dataset
!
    if (iceflag) then
      where (arrxzy(:lenin,1,j) < tsice + 0.0001)
        ice(:lenin) = tsice
      elsewhere
        ice(:lenin) = 0.
      end where

      call lininterp (plon, lenin, lenot, 1, ice, &
                      ice_reduced)

      where (ice_reduced(:lenot) < 0.5*tsice) 
        arrxz(:,1) = tsice
      end where
!
! Print stats on conversion if this is an sst dataset
!
      nfull = 0
      do i=1,lenin
        if (arrxzy(i,1,j) < tsice + 0.0001) nfull = nfull + 1
      end do

      nred = 0
      do i=1,lenot
        if (arrxz(i,1) == tsice) nred = nred + 1
      end do
      write(6,*)nfull,' full ice pts. ',nred,' reduced out of ', lenot
    end if

    call t_stopf('interpolation')

    call t_startf('interp_copy')
    arrxzy(:,:,j) = arrxz(:,:)
    call t_stopf('interp_copy')
  end do

  return
  end subroutine
