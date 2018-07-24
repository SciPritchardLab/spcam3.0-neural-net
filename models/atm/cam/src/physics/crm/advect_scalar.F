c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c     positively definite monotonic advection with non-oscillatory option
c-----------------------------------------------------------------


 	subroutine advect_scalar(f,fadv,flux)
 	
 	implicit none
	include 'common.inc'

	real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
	real flux(nz), fadv(nz)

	real df(nx, ny, nzm)
	integer i,j,k

          do k=1,nzm
            do j=1,ny
             do i=1,nx
              df(i,j,k) = f(i,j,k)
             end do
            end do
          end do

	if(RUN3D) then
	  call advect_scalar3D(f, u, v, w, flux)
	else
	  call advect_scalar2D(f, u, w, flux)	  
	endif

          do k=1,nzm
            do j=1,ny
             do i=1,nx
              fadv(k)=fadv(k)+f(i,j,k)-df(i,j,k)
             end do
            end do
          end do


	return
	end

