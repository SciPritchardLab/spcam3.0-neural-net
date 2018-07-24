subroutine cubinterp (plon, nlonin, nlonot, numlev, arrin, &
                      arrot, mono)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constants
!
! $Id: cubinterp.f90,v 1.1.44.1 2002/06/15 13:50:06 erik Exp $
!
  implicit none
!
! Parameters
!
  real(r8), parameter :: epsmac = 1.e-13
!
! Arguments
!
  integer plon            ! input/output longitude dimension
  integer nlonin          ! number of longitudes input
  integer nlonot          ! number of longitudes output
  integer numlev          ! number of vertical levels
!
! input (interpolated from) and output (interpolated to) arrays
!
  real(r8), dimension (0:plon-1,numlev) :: arrin, arrot
  logical mono            ! whether or not do impose monotonicity
!
! Local workspace
!
  integer i, ii           ! longitude indices

  real(r8) dxin             ! delta-x on input grid
  real(r8) dxot             ! delta-x on output grid
  real(r8) xlocot           ! position on output grid (i * dxot)
  real(r8) left             ! left position on input grid (ii * dxin)
  real(r8) right            ! right position on input grid ((ii+1) * dxin)
  real(r8) beta             ! (Xi - Xl) / (Xr - Xl)
  real(r8) onemb            ! 1. - beta
  real(r8) onepb            ! 1. + beta
  real(r8) twomb            ! 2. - beta
  real(r8) delq(numlev)     ! delta-q (fm Williamson paper)
  real(r8) test(numlev)     ! compared with 0 and 3 (fm Williamson paper)
  real(r8) dl(numlev)       ! left derivative estimate (fm Williamson paper)
  real(r8) dr(numlev)       ! right derivative estimate (fm Williamson paper)
  real(r8) tmp(numlev)
!
! Input array with wrap points included
!
  real(r8), dimension (-1:nlonin+1,numlev) :: arrin_wr

  if (nlonin.gt.plon .or. nlonot.gt.plon) then
    write(6,*)'Number of longitudes cannot exceed dimension'
    stop 99
  end if
  if (nlonin.lt.2 .or. nlonot.lt.1) then
    write(6,*)'Must have at least 2 points for input, 1 for output'
    stop 99
  end if
!
! Build local copy of input array with wrap points
!
  arrin_wr(-1,:) = arrin(nlonin-1,:)
  arrin_wr(0:nlonin-1,:) = arrin(:,:)
  arrin_wr(nlonin,:) = arrin(0,:)
  arrin_wr(nlonin+1,:) = arrin(1,:)
!
! 1st longitude point: *assume* overlapping
!
  arrot(0,:) = arrin(0,:)
  dxin = 1.d0/nlonin
  dxot = 1.d0/nlonot
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
        beta = (xlocot-left)/(right-left)
        onepb = 1. + beta
        onemb = 1. - beta
        twomb = 2. - beta

        if (mono) then
          dl(:) = -arrin_wr(ii-2,:)/(3.*dxin) - &
                   arrin_wr(ii-1,:)/(2.*dxin) + &
                   arrin_wr(ii  ,:)/dxin -      &
                   arrin_wr(ii+1,:)/(6.*dxin)

          dr(:) =  arrin_wr(ii-2,:)/(6.*dxin) - &
                   arrin_wr(ii-1,:)/dxin +      &
                   arrin_wr(ii  ,:)/(2.*dxin) + &
                   arrin_wr(ii+1,:)/(3.*dxin)

          delq(:) = (arrin_wr(ii,:) - arrin_wr(ii-1,:))/dxin
          test(:) = abs(3.*delq(:))*(1. - epsmac)

          where (dl(:)*delq(:) .lt. 0.)
            dl(:) = 0.
          elsewhere
            tmp(:) = abs(dl(:))
            dl(:) = sign (min (tmp(:), test(:)), delq(:))
          end where

          where (dr(:)*delq(:) .lt. 0.)
            dr(:) = 0.
          elsewhere
            tmp(:) = abs(dr(:))
            dr(:) = sign (min (tmp(:), test(:)), delq(:))
          end where

          arrot(i,:) = (3. - 2.*onemb)*onemb**2*arrin_wr(ii-1,:) + &
                       dxin*beta*onemb**2*dl(:) +                  &
                       (3. - 2.*beta)*beta**2*arrin_wr(ii,:) -     &
                       dxin*beta**2*onemb*dr(:)            
        else

          arrot(i,:) = -beta/6.*onemb*twomb*arrin_wr(ii-2,:) +     &
                       0.5*onepb*onemb*twomb*arrin_wr(ii-1,:) +    &
                       0.5*onepb*beta*twomb*arrin_wr(ii,:) -       &
                       onepb/6.*beta*onemb*arrin_wr(ii+1,:)
        end if

        goto 10

      end if
    end do             ! do ii=1,nlonin
    write(6,*)'Interpolation failed for i=',i
    stop 99
10  continue
  end do               ! do i=1,nlonot-1
!
! If interpolating to reduced grid, Fill remainder with fillvalue
!
  arrot(nlonot:,:) = fillvalue

  return
end subroutine
