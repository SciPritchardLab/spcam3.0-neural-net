c  These routines exchanges overlaping information for various subdomains to
c  be able to calculate various quantities in common /com3d/
c  near the subdomain boundaries.
c
c

        subroutine boundaries(flag)

        implicit none
	include 'common.inc'

        integer flag

        if(flag.eq.0) then

         call bound_exchange
     &         (u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
         call bound_exchange
     &         (v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
         call bound_exchange
     &         (w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)

        endif

        if(flag.eq.1) then

         call bound_exchange
     &         (u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2,2,1)
         call bound_exchange
     &         (v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2,2,2,3,2)
         call bound_exchange
     &         (w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2,2,2,2,3)


        endif


        if(flag.eq.2) then

         call bound_exchange
     &         (t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,4)
         call bound_exchange
     &         (q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
         if(docloud.and.doprecip)
     &     call bound_exchange
     &         (qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,9)

        endif
        
        if(flag.eq.3) then
        
         call bound_exchange
     &         (t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
         call bound_exchange
     &         (q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
         if(docloud.and.doprecip)
     &     call bound_exchange
     &         (qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,9)
        
        endif
        
        if(dosgs.and.flag.eq.4) then
        
         call bound_exchange(tk,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7)
         call bound_exchange(tkh,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,8)
        
        endif
        
        
        return
        end

	
	

	subroutine bound_exchange(f,dimx1,dimx2,dimy1,dimy2,dimz,
     &                 i_1, i_2, j_1, j_2, id)
	
	implicit none
	include 'grid.inc'
	
	integer dimx1, dimx2, dimy1, dimy2, dimz
        integer i_1, i_2, j_1, j_2
	real f(dimx1:dimx2, dimy1:dimy2, dimz)
	
	real buffer((nx+ny)*3*nz)	! buffer for sending data
	integer id   ! id of the sent field
	
	integer i, j, k, n
	integer i1, i2, j1, j2
	
	i1 = i_1 - 1
	i2 = i_2 - 1
	j1 = j_1 - 1
	j2 = j_2 - 1

c----------------------------------------------------------------------
c  Send buffers to neighbors
c----------------------------------------------------------------------


	if(RUN3D) then

c "North" -> "South":	

	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

c "North-East" -> "South-West":	
	
	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

c "South-East" -> "North-West":

	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

c "South" -> "North":

	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=1,nx
	           n = n+1
	           buffer(n) = f(i,j,k) 
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

c "South-West" -> "North-East":
	  
	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k) 
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do


c To "North-West" -> "South-East":
	  	  
	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do
	     

	endif

c  "East" -> "West":
 	  
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

c "West" -> "East":

	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do


	     return
	     end
	     
	     
	     
	     
	
	
	subroutine bound_duvdt()

        implicit none
	include 'common.inc'
         
	integer i,j,k

	  do k=1,nzm
	   do j=1,ny
	     dudt(nxp1,j,k,na) = dudt(1,j,k,na)
	   end do
	  end do

	  if(RUN3D) then

	    do k=1,nzm
	     do i=1,nx
	      dvdt(i,nyp1,k,na) = dvdt(i,1,k,na)
	     end do
	    end do

	  endif
	
	return
	end
