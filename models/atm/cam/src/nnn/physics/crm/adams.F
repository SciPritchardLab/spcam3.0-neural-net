c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c       Adams-Bashforth scheme
c----------------------------------------------------------------

	subroutine adams()

	implicit none
	include 'common.inc'

	real dtdx, dtdy, dtdz, rhox, rhoy, rhoz	
	integer i,j,k

        dtdx = dtn/dx
        dtdy = dtn/dy
        dtdz = dtn/dz

	do k=1,nzm 
          rhox = rho(k)*dtdx  
	  rhoy = rho(k)*dtdy  
	  rhoz = rhow(k)*dtdz  

	  do j=1,ny
	   do i=1,nx
	     dudt(i,j,k,nc) = u(i,j,k) + dt3(na)*
     &	        (at*dudt(i,j,k,na)+bt*dudt(i,j,k,nb)+ct*dudt(i,j,k,nc))
	   
	     dvdt(i,j,k,nc) = v(i,j,k) + dt3(na)*
     &          (at*dvdt(i,j,k,na)+bt*dvdt(i,j,k,nb)+ct*dvdt(i,j,k,nc))
	   
	     dwdt(i,j,k,nc) = w(i,j,k)+dt3(na)*
     &          (at*dwdt(i,j,k,na)+bt*dwdt(i,j,k,nb)+ct*dwdt(i,j,k,nc))
	   end do 
	  end do 

          if(nstep.eq.1) then 	

          do j=1,ny
           do i=1,nx
             u(i,j,k) = dudt(i,j,k,nc) * rhox
             v(i,j,k) = dvdt(i,j,k,nc) * rhoy
             w(i,j,k) = dwdt(i,j,k,nc) * rhoz
           end do
          end do

	  else

          do j=1,ny
           do i=1,nx
             u(i,j,k) = 0.5*(u(i,j,k)+dudt(i,j,k,nc)) * rhox
             v(i,j,k) = 0.5*(v(i,j,k)+dvdt(i,j,k,nc)) * rhoy
             w(i,j,k) = 0.5*(w(i,j,k)+dwdt(i,j,k,nc)) * rhoz
           end do
          end do

	  end if

	end do 

	
	return
	end

	
