	subroutine stepout()

	implicit none	
	include 'common.inc'
	include 'params.inc'
	
	integer i,j,k,ic,jc
	real div, divmax, divmin
	real rdx, rdy, rdz, coef
	real *8 tav
	integer im,jm,km
	real wmax
	character *6 nstepchar

c	if(rank.eq.0) print *,'NSTEP = ',nstep,'    NCYCLE=',ncycle


	if(mod(nstep,nprint).eq.0) then
	
	 divmin=1.e20
	 divmax=-1.e20
	 
	 rdx = 1./dx
	 rdy = 1./dy

	 wmax=0.
	 do k=1,nzm
	  coef = rho(k)*adz(k)*dz
	  rdz = 1./coef
	  if(ny.ne.1) then
	   do j=1,ny-1*YES3D
	    jc = j+1*YES3D
	    do i=1,nx-1
	     ic = i+1
	     div = (u(ic,j,k)-u(i,j,k))*rdx + (v(i,jc,k)-v(i,j,k))*rdy +
     &		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
             divmax = max(divmax,div)
             divmin = min(divmin,div)
		if(w(i,j,k).gt.wmax) then
			wmax=w(i,j,k)
			im=i
			jm=j
			km=k
		endif
	    end do
	   end do
	  else
	    j = 1
	    do i=1,nx-1
	     ic = i+1
	     div = (u(ic,j,k)-u(i,j,k))*rdx +
     &		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
             divmax = max(divmax,div)
             divmin = min(divmin,div)
		if(w(i,j,k).gt.wmax) then
			wmax=w(i,j,k)
			im=i
			jm=j
			km=k
		endif
	    end do
	  endif
	 end do

	if(rank.eq.0) then
	
c	 write(6,*)
c	 write(6,'(16f8.5)')((w(i,j,10),i=1,16),j=16,1,-1)
c	 write(6,*)
c	 write(6,'(16f8.5)')((u(i,j,10),i=1,16),j=16,1,-1)
c	 write(6,*)
c	 write(6,'(16f8.5)')((v(i,j,10),i=1,16),j=16,1,-1)
c	 write(6,*)

c--------------------------------------------------------
	    print*,'DAY = ',day	
	    write(6,*) nstep,ncycle,divmax,divmin
	    write(6,*) 'SST=',tabs_s, '  pres0=',pres0,
     &		'  H=',fluxt0*cp*rho(1),'  LE=',fluxq0*lcond*rho(1)	
c	print*,wmax,im,jm,km

	endif

  	call fminmax_print('u:',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
  	call fminmax_print('v:',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm)
	call fminmax_print('w:',w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
	call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
	call fminmax_print('t:',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
	call fminmax_print('tabs:',tabs,1,nx,1,ny,nzm)
	call fminmax_print('q:',q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
	call fminmax_print('qn:',qn,1,nx,1,ny,nzm)
	call fminmax_print('qp:',qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
	call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
	call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)
	call fminmax_print('shf:',fluxbt,1,nx,1,ny,1)
	call fminmax_print('lhf:',fluxbq,1,nx,1,ny,1)
	call fminmax_print('uw:',fluxbu,1,nx,1,ny,1)
	call fminmax_print('vw:',fluxbv,1,nx,1,ny,1)

	 end if ! (mod(nstep,nprint).eq.0)
	
	
	return
	end
	
	
	
