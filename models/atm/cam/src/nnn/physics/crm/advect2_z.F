

c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c       momentum tendency due to the 2nd-order-central vertical advection
c----------------------------------------------------------------


	subroutine advect2_z()

 	implicit none
	include 'common.inc'

        real fuz(nx,ny,nz),fvz(nx,ny,nz),fwz(nx,ny,nzm)
        integer i, j, k, kc, kb
        real dz2, dz25, www, rhoi

	dz25=1./(4.*dz)
	dz2=dz25*2.

	do j=1,ny
	 do i=1,nx
	  fuz(i,j,1) = 0.
	  fvz(i,j,1) = 0.
	  fuz(i,j,nz) = 0.
	  fvz(i,j,nz) = 0.
	  fwz(i,j,1) = 0.
	  fwz(i,j,nzm) = 0.
	 end do
	end do


	if(RUN3D) then

	do k=2,nzm
	 kb = k-1
	 rhoi = dz25 * rhow(k)
	 do j=1,ny
	  do i=1,nx
	   fuz(i,j,k) = rhoi*(w(i,j,k)+w(i-1,j,k))*(u(i,j,k)+u(i,j,kb))
	   fvz(i,j,k) = rhoi*(w(i,j,k)+w(i,j-1,k))*(v(i,j,k)+v(i,j,kb))
	  end do
	 end do
	end do

	else

	do k=2,nzm
	 kb = k-1
	 rhoi = dz25 * rhow(k)
	 do j=1,ny
	  do i=1,nx
	    www = rhoi*(w(i,j,k)+w(i-1,j,k)) 
	    fuz(i,j,k) = www*(u(i,j,k)+u(i,j,kb))
	    fvz(i,j,k) = www*(v(i,j,k)+v(i,j,kb))
	  end do
	 end do
	end do


	endif
	
	do k=1,nzm
	 kc = k+1
	 rhoi = 1./(rho(k)*adz(k))
	 do j=1,ny
	  do i=1,nx
	   dudt(i,j,k,na)=dudt(i,j,k,na)-(fuz(i,j,kc)-fuz(i,j,k))*rhoi
	   dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fvz(i,j,kc)-fvz(i,j,k))*rhoi	 
           fwz(i,j,k)=dz25*(w(i,j,kc)*rhow(kc)+w(i,j,k)*rhow(k))*
     &                     (w(i,j,kc)+w(i,j,k))
	  end do
	 end do
	end do

	do k=2,nzm
	 kb=k-1
	 rhoi = 1./(rhow(k)*adzw(k))
	 do j=1,ny
	  do i=1,nx
 	   dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fwz(i,j,k)-fwz(i,j,kb))*rhoi
	  end do
	 end do
 	end do ! k

c	if(dostatis) call stat_tke_adv()

	return
	end

