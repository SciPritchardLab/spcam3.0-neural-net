subroutine lininterp (plon, nlonin, nlonot, numlev, arrin, arrot)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constants
!
! $Id: lininterp.f90,v 1.1.44.1 2002/06/15 13:50:07 erik Exp $
!
  implicit none
!
! Arguments
!
  integer plon, nlonin, nlonot, numlev
  real(r8), dimension (0:plon-1,numlev) :: arrin, arrot
!
! Local workspace
!
  integer i,ii,iiright
  real(r8) dxin, dxot, xlocot, left, right, lfact, rfact

  if (nlonin.gt.plon .or. nlonot.gt.plon) then
    write(6,*)'Number of longitudes cannot exceed dimension'
    stop 99
  end if

  if (nlonin.lt.2 .or. nlonot.lt.1) then
    write(6,*)'Must have at least 2 points for input, 1 for output'
    stop 99
  end if
!
! 1st longitude point: *assume* overlapping
!
  arrot(0,:) = arrin(0,:)
  dxin = 1./nlonin
  dxot = 1./nlonot
!
! Loop over longitude of output grid, finding interpolation points
!
  do i=1,nlonot-1
    xlocot = i*dxot
    right = 0.
    do ii=1,nlonin
      left = right
      right = ii*dxin
      if (right.ge.xlocot) then
        lfact = (right-xlocot)/(right-left)
        rfact = (xlocot-left)/(right-left)
!
! Wrap rhs index of input grid if interpolating reduced -> full, and we are 
! at the rhs of the input grid
!
        iiright = ii
        if (iiright.eq.nlonin) iiright = 0
        if (arrin(ii-1,1) > 1.e35 .or. arrin(iiright,1) > 1.e35) then
          write(6,*)'arrin(',ii-1,')=',arrin(ii-1,1)
          write(6,*)'arrin(',iiright,')=',arrin(iiright,1)
          write(6,*)'i,ii,nlonin,nlonot=',i,ii,nlonin,nlonot
          stop
        end if
        arrot(i,:) = arrin(ii-1,:)*lfact + arrin(iiright,:)*rfact
        goto 10
      end if
    end do
    write(6,*)'Interpolation failed for i=',i
    stop 99
10  continue
  end do
!
! If interpolating to reduced grid, Fill remainder with fillvalue
!
  arrot(nlonot:,:) = fillvalue

  return
end subroutine
