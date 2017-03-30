! Terrain Filter
! 
! Contributed by S.J. Lin.
!
! Added to the definesurf program by G. Grant, 30 June 2000.
! Updated with latest version from S.J. by B. Eaton, 23 August 2001
!
! Notes from S.J.:
!
! "I compute the unsmoothed mean height and the variance
! exactly the same as the standard CCM utility. The only difference
! is the grid being uniformly spaced from North pole to South pole.
! The filter is applied to the mean height and the sqaure root of 
! the variance (the standard deviation).
!
! For the 2x2.5 deg resolution 
!
!      mlon = 144
!      mlat = 91
!
! Assuming the mean height is Z(mlon,mlat), and the standard deviation
! (the sqaure root of the variance) is SD(moln,mlat), the filter 
! algorithm goes like this:
!
!       call sm2(mlon, mlat,  Z, itmax_Z,  0.25D0)
!       call sm2(mlon, mlat, SD, itmax_SD, 0.25D0)
!
! where 0.25D0 is the dimensionless filter coefficient, and 
!
!       itmax_Z  = 2*mlat
!       itmax_SD = mlon
!
! [As discussed elsewhere] the above filtering is a bit too strong.
! But it is the filter I used up to now.
! I am currently testing the following setting 
!
!       itmax_Z  = mlat/2
!       itmax_SD = mlon/4
! "


      subroutine sm2(im, jm, ht, itmax, c)
!
! Del-2 diffusion on the sphere
!
      implicit none

! Input:
      integer im               ! e-w dimension (eg, 144 for 2.5 deg resolution)
      integer jm               ! n-s doemsnion (eg, 91 for 2 deg resolution)
      integer itmax            ! iteration count
      real*8  c                ! filter coefficient

! Input/Output
      real*8 ht(im,jm)         ! array to be filtered

! Local
      real*8 dg(im,jm)   ! del2 of h
      real*8 cose(jm), cosp(jm), sinp(jm), sine(jm)
      real*8 dl
      real*8 dp
      real*8 fmin, fmax
      integer jm1
      integer mnk, mxk
      integer ndeg
      integer it, i, j
      real*8 s1, s2

      jm1 = jm-1

      call setrig(im, jm, dp, DL, cosp, cose, sinp, sine)

      call pmnx(ht, im, jm, fmin, fmax, mnk, mxk)
      write(6,*) 'hmax=', fmax,' at j= ',mxk
      write(6,*) 'hmin=', fmin,' at j= ',mnk

      ndeg = 60                    ! starting latitude for the monotonicity
                                   ! preserving polar filter

      call pmnx(ht,im,jm,fmin,fmax,mnk,mxk)
      write(6,*) 'hmax=', fmax,' at j= ',mxk
      write(6,*) 'hmin=', fmin,' at j= ',mnk

! Apply Monotonicity preserving polar filter
      call plft2d(im, jm, ht, 2, jm1, ndeg)
      call avgp2(ht, sine, im, jm)

      do it=1,itmax
	call del2(ht, im, jm, dg, cosp, cose, sine, DL, dp, ndeg)
        call plft2d(im, jm, dg, 2, jm1, ndeg)

	do j=1,jm
	   do i=1,im
             ht(i,j) = ht(i,j) + c*dg(i,j)
           enddo
        enddo
      enddo

! Final polar filter
      call plft2d(im, jm, ht, 2, jm1, ndeg)

      return
      end

      subroutine del2(h, im, jm, dg, cosp, cose, sine, dL, dp, ndeg)
      implicit none

! AE = 1  (unit radius)
! Input:
        integer im
        integer jm
        integer ndeg
! Input-output

	real*8 h(im,jm)
	real*8 dg(im,jm)              ! del2 of h
	real*8 cose(jm),cosp(jm)
	real*8 sine(jm)
	real*8 PI, ycrit, coszc, CD
	real*8 DL, dp

! Local
	real*8 fx(im,jm)   ! e-w fluxes
	real*8 fy(im,jm)   ! n-s fluxes
        integer i, j
 
	call grad(h, im, jm, fx, fy, cosp, dl, dp)

        PI = 4. * ATAN(1.)
        ycrit = float(ndeg)*PI/180.
        coszc = cos(ycrit)

 	CD = 0.25*DL*DP*coszc**2
!	CD = 0.25*DL*DP*cosp(2)**2

        do j=2,jm-1
           do i=1,im
              fx(i,j) = fx(i,j) * CD
           enddo
        enddo

        do j=2,jm
	   do i=1,im
	      fy(i,j) = fy(i,j) * CD
           enddo
        enddo

	call divg(im,jm,fx,fy,DG,cosp,cose,sine, dl, dp)

	return
	end

      subroutine divg(im, jm, fx, fy, dg, cosp, cose, sine, dl, dp)
      implicit none

      integer im
      integer jm
      real*8 fx(im,jm)   ! e-w fluxes
      real*8 fy(im,jm)   ! n-s fluxes
      real*8 DG(im,jm)   ! del2 of h
      real*8 wk(im,jm)
      real*8 cosp(jm),  cose(jm), sine(jm)
      real*8 rdx
      real*8 dl, dp, CDP, sum1, sum2
      integer i,j

	do j=2,jm-1

	rdx = 1./ (cosp(j)*DL)

           do i=1,im-1
              DG(i,j) = (fx(i+1,j) - fx(i,j)) * rdx
           enddo
              DG(im,j) = (fx(1,j) - fx(im,j)) * rdx
        enddo

	do j=2,jm
           do i=1,im
              wk(i,j) = fy(i,j) * cose(j)
           enddo
        enddo

	do j=2,jm-1
           CDP = 1./ (DP*cosp(j))
           do i=1,im
              DG(i,j) = DG(i,j) + (wk(i,j+1) - wk(i,j)) * CDP
           enddo
	enddo

! Poles;

	sum1 = wk(im, 2)
	sum2 = wk(im,jm)

	do i=1,im-1
           sum1 = sum1 + wk(i, 2)
           sum2 = sum2 + wk(i,jm)
	enddo

	sum1 =  sum1 / ( float(im)*(1.+sine(2)) )
	sum2 = -sum2 / ( float(im)*(1.+sine(2)) )

	do i=1,im
           DG(i, 1) =  sum1
           DG(i,jm) =  sum2
	enddo

	return
	end

      subroutine grad(h, im, jm, fx, fy, cosp, DL, DP)
      implicit none
      integer im
      integer jm
      real*8 h(im,jm)
      real*8 fx(im,jm)   ! e-w fluxes
      real*8 fy(im,jm)   ! n-s fluxes
      real*8 cosp(jm)
      real*8 RDP, DL, DP, rdx
      integer i, j

	RDP = 1./ DP

      do j=2,jm
         do i=1,im
            fy(i,j) = (h(i,j) - h(i,j-1)) * RDP
         enddo
      enddo

      do j=2,jm-1

	rdx = 1./ (cosp(j)*DL)
        fx(1,j) = (h(1,j) - h(im,j)) * rdx
        do i=2,im
           fx(i,j) = (h(i,j) - h(i-1,j)) * rdx
        enddo
      enddo

      return
      end

      subroutine avgp2(p, sine, im, jm)
      implicit none
      integer im, jm
      real*8 p(im,jm)
      real*8 sine(jm)
      real*8 sum1, sum2
      real*8 sum3, sum4
      real*8 rim
      integer i
      integer j
      integer jm1

      jm1 = jm-1
      rim = 1./ float(im)

      call sump2(p(1,1),p(1,jm),IM,sum1,sum2)
      sum1 = sum1*(1.+sine(2))
      sum2 = sum2*(1.+sine(2))

      call sump2(p(1,2),p(1,jm1),IM,sum3,sum4)
      sum1 = rim * ( sum1 + sum3*(sine(3)-sine(2)) ) / (1.+sine(3))
      sum2 = rim * ( sum2 + sum4*(sine(3)-sine(2)) ) / (1.+sine(3))

      do i=1,im
      P(i,  1) = sum1
      P(i,  2) = sum1
      P(i,jm1) = sum2
      P(i, jm) = sum2
      enddo
      return
      end

      subroutine sump2(p1,p2,im,s1,s2)
      implicit none
      integer im,i
      real*8 s1,s2
      real*8 p1(*),p2(*)
 
         s1 =  p1(im)
         s2 =  p2(im)
 
      do i=1,im-1
         s1 =  s1 + p1(i)
         s2 =  s2 + p2(i)
      enddo
      return
      end

      subroutine plft2d(im, jm, p, JS, JN, ndeg)
!
! This is a weak LOCAL polar filter.
! Developer: Shian-Jiann Lin

      implicit none

      integer im
      integer jm
      integer js, jn, ndeg
      real*8 p(im,jm)

      integer i, j, n, ideg, jj, jc
      real*8 cosp(jm),cose(jm)
      real*8 a(0:im/2+1)

      real*8 sine(jm),sinp(jm)
      real*8, allocatable, save :: se(:), sc(:)

      real*8 pi, dp, dl, e0, ycrit, coszc, smax, rn, rn2, esl, tmp

      data IDEG /0/

      if(IDEG .ne. ndeg) then
	IDEG = ndeg
!  (e0 = 2.6)
      e0 = 0.5 * sqrt(27.)
      PI = 4. * ATAN(1.)

      allocate( sc(jm), se(jm))

      call setrig(im, jm, DP, DL, cosp, cose, sinp, sine)

      ycrit = IDEG*PI/180.
      coszc = cos(ycrit)

      smax = (jm-1)/2
      write(6,*) 'Critical latitude in local pft = ',ndeg

      a(0) = 1.
      do n=1,im/2+1
      rn = n
      rn2 = 2*n
      a(n) = sqrt(rn2+1.) * ((rn2+1.)/rn2)**rn
      enddo

      do j=2,jm-1
      sc(j) = coszc / cosp(j)

      IF(sc(j) .gt.1. .and. sc(j) .le. 1.5 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. +  (1.-esl) / (1.+esl)
      ELSEIF(sc(j) .gt.1.5 .and. sc(j) .le. e0 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(sc(j).gt. e0) THEN
! Search
      do jj=1,im/2
      if(sc(j).le. a(jj)) then
      jc = jj
      write(*,*) 'jc=', jc
      goto 111
      endif
      enddo

      jc = im/2 + 1
111   continue
      tmp = ((sc(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      sc(j) =  jc + min(1., tmp)
!     sc(j) =  min(smax,sc(j))
      ENDIF
      enddo
! ====================================================
      do j=2,jm
      se(j) = coszc / cose(j)
      IF(se(j) .gt.1. .and. se(j) .le. 1.5 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + (1.-esl) / (1.+esl)
      ELSEIF(se(j) .gt.1.5 .and. se(j) .le. e0 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(se(j).gt. e0) THEN
! Search
      do jj=1,im/2
      if(se(j) .le. a(jj)) then
      jc = jj
      goto 222
      endif
      enddo

      jc = im/2 + 1
222   continue
      tmp = ((se(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      se(j) =  jc + min(1., tmp)
!     se(j) =  min(smax,se(j))
      ENDIF
      enddo

      do i=1,im
        se( 2) = sc(2)
        se(jm) = sc(jm-1)
      enddo

      do j=2,jm-1
         write(*,*) j,sc(j)
      enddo
      ENDIF

      if(JN .eq. (jm-1) ) then
! Cell-centered variables
         call lpft(im, jm, p, 2, jm-1, Sc)
      else
! Cell-edge variables
         call lpft(im, jm, p, 2, jm, Se)
      endif
      return
      end


      subroutine lpft(im, jm, p, j1, j2, s)
      implicit none

      integer im, jm, j1, j2
      real*8  p(im,jm)
      real*8  s(jm)

! Local
      integer i, j, n, nt

      real*8 ptmp(0:im+1)
      real*8 q(0:im+1)
      real*8 frac, rsc, bt

      do 2500 j=j1,j2
      if(s(j) .gt. 1.02) then

        NT  = INT(S(j))
        frac = S(j) - NT
        NT = NT-1

        rsc = 1. / (1.+frac)
        bt = 0.5 * frac

        do i=1,im
           ptmp(i) = p(i,j)
        enddo

        ptmp(0)    = p(im,j)
        ptmp(im+1) = p(1 ,j)

        if( NT .lt. 1 ) then
          do i=1,im
             p(i,j) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo
        else
          do i=1,im
             q(i) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo

           do 500 N=1,NT
              q(0) = q(im)
              do i=1,im
                 ptmp(i) = q(i) + q(i-1)
              enddo
                 ptmp(im+1) = ptmp(1)

             if ( n .eq. nt ) then
               do i=1,im
                  p(i,j) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             else
               do i=1,im
                  q(i) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             endif
500        continue
        endif
      endif
2500  continue

      return
      end

      subroutine pmnx(a,nx,ny,fmin,fmax,mnk,mxk)
      implicit none
      integer nx
      integer ny
      integer mnk
      integer mxk
      real*8 a(nx,*)
      real*8 fmax, fmin, temp
      integer i,j

      fmax = a(1,1)
      fmin = a(1,1)
      mnk = 1
      mxk = 1

      do j=1,ny
        do i=1,nx
           temp = a(i,j)
           if(temp.gt.fmax) then
              fmax = temp
              mxk = j
           elseif(temp .lt. fmin) then
              fmin = temp
              mnk = j
           endif
        enddo
      enddo

      return
      end

      subroutine setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

      implicit none

      integer im, jm
      integer j, jm1
      real*8 sine(jm),cosp(jm),sinp(jm),cose(jm)
      real*8 dp, dl
      double precision pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

! Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      return
      end
