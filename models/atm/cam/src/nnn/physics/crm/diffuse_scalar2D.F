	subroutine diffuse_scalar2D(field,fluxb,fluxt,tkh,flux)

	implicit none
	include 'grid.inc'
	
! input
        real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
        real tkh(0:nxp1, 1-YES3D:nyp1, nzm)	! eddy conductivity
        real fluxb(nx,ny)		! bottom flux
        real fluxt(nx,ny)		! top flux
	real flux(nz)
    
! local        
        real flx(0:nx,1,0:nzm)
    	real dfdt(nx,ny,nzm) 
	real rdx2,rdz2,rdz,rdx5,rdz5,tmp
	real dxz,dzx,tkx,tkz,rhoi,coef
	integer i,j,k,ib,ic,kc,kb


	rdx2=1./(dx*dx)
	rdz2=1./(dz*dz)
	rdz=1./dz
	rdz5=0.5*rdz2
	dxz=dx/dz
	dzx=dz/dx
	coef = 1./float(nx)

	j=1

	do k=1,nzm	
	  do i=1,nx
  	    dfdt(i,j,k)=0.
	  end do 
	end do 



	do k=1,nzm
	
	  rdx5 = 0.5*rdx2 !  *(dx/(dz*adz(k)))**2

	  do i=0,nx
	    ic=i+1
	    tkx=rdx5*(tkh(i,j,k)+tkh(ic,j,k)) 	
	    flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))
	  end do 
	  do i=1,nx
	    ib=i-1
  	    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
	  end do 

	end do 

	tmp=1./adzw(nz)
	do i=1,nx	
	   flx(i,j,0)=fluxb(i,j)*rdz*rhow(1)
	   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
	   flux(1) = flux(1) + flx(i,j,0)
	end do

	do k=1,nzm-1
	 kc=k+1	
	 rhoi = rhow(kc)/adzw(kc)
	 do i=1,nx
	    tkz=rdz5*(tkh(i,j,k)+tkh(i,j,kc))
	    flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi
	    flux(kc) = flux(kc) + flx(i,j,k)
	 end do
	end do
	
	do k=1,nzm
	 kb=k-1
	 rhoi = 1./(adz(k)*rho(k))
	 do i=1,nx		 
  	  dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
	  field(i,j,k)=field(i,j,k) + dfdt(i,j,k)
	 end do 
	end do 

	return
	end



