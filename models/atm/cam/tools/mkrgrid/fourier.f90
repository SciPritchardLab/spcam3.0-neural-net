subroutine fourier (plon, nlonin, nlonot, numlev, arrin, &
                    arrot, wnummax, fname)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constants
  use fortfft
!
! $Id: fourier.f90,v 1.1.44.1 2002/06/15 13:50:06 erik Exp $
!
  implicit none
!
! Arguments
!
  character*(*) :: fname

  integer plon, nlonin, nlonot, numlev
  integer :: wnummax

  real(r8), dimension (plon,numlev) :: arrin, arrot
!
! Local workspace
!
  integer ifaxin(19)
  integer ifaxot(19)
  integer throwindx

  real(r8), dimension (plon+2,numlev) :: x
  real(r8) trigin(3*nlonin/2+1)
  real(r8) trigot(3*nlonot/2+1)
  real(r8) work((plon+1)*numlev)

  if (nlonin.gt.plon .or. nlonot.gt.plon) then
    write(6,*)'Number of longitudes cannot exceed dimension'
    stop 99
  end if
  if (nlonin.lt.2 .or. nlonot.lt.1) then
    write(6,*)'Must have at least 2 points for input, 1 for output'
    stop 99
  end if
!
! Call fft setup routines for input and output grids
!
  call set99 (trigin, ifaxin, nlonin)
  call set99 (trigot, ifaxot, nlonot)
!
! Need to use a temp array for transform due to length requirements of fft
!
  if (fname.eq.'PS') then
    x(:nlonin,:) = log(arrin(:nlonin,:))
  else
    x(:nlonin,:) = arrin(:nlonin,:)
  end if
  x(nlonin+1:,:) = 0.

  call fft991 (x, work, trigin, ifaxin, 1, plon+2, nlonin, numlev, -1)
!
! Zero out waves which will not be included
! The 3 in throwindx counts 2 for the mean plus 1 for fortran indexing 
! starting at 1
!
  throwindx = 2*wnummax + 3
      
  if (throwindx.lt.nlonot) then
    x(throwindx:,:) = 0.
  end if

  call fft991 (x, work, trigot, ifaxot, 1, plon+2, nlonot, numlev, +1)
!
! If interpolating to reduced grid, Fill remainder with fillvalue
!
  if (fname.eq.'PS') then
    arrot(:nlonot,:) = exp(x(:nlonot,:))
  else
    arrot(:nlonot,:) = x(:nlonot,:)
  end if
  arrot(nlonot+1:,:) = fillvalue

  return
end subroutine
