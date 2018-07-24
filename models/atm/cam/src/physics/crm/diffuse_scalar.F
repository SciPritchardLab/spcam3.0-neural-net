	subroutine diffuse_scalar(f,fluxb,fluxt,fdiff,flux)

	implicit none
	include 'common.inc'
! input:	
        real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
        real fluxb(nx,ny)		! bottom flux
        real fluxt(nx,ny)		! top flux
	real flux(nz)
	real fdiff(nz)
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
	  call diffuse_scalar3D
     &              (f,fluxb,fluxt,tkh,flux)
	else  
	  call diffuse_scalar2D
     &              (f,fluxb,fluxt,tkh,flux)
	endif

          do k=1,nzm
            do j=1,ny
             do i=1,nx
              fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
             end do
            end do
          end do


	return
	end



