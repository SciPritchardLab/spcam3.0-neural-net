      subroutine lininterp(plon, nlonin, nlonot, numlev, arrin, arrot)
!
! $Id: lininterp.F90,v 1.1 2001/07/25 00:26:28 rosinski Exp $
!
      implicit none
!
! Parameters
!
      real*8 spval
      parameter (spval = 1.e36)
!
! Arguments
!
      integer plon, nlonin, nlonot, numlev
      real*8, dimension (0:plon-1,numlev) :: arrin, arrot
!
! Local workspace
!
      integer i,ii,iiright
      real*8 dxin, dxot, xlocot, left, right, lfact, rfact

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
            arrot(i,:) = arrin(ii-1,:)*lfact + arrin(iiright,:)*rfact
            goto 10
          end if
        end do
        write(6,*)'Interpolation failed for i=',i
        stop 99
   10   continue
      end do
!
! If interpolating to reduced grid, Fill remainder with spval
!
      arrot(nlonot:,:) = spval

      return
      end
