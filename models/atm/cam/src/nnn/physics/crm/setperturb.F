c
c  Set random perturbations and other stuff:
c

	subroutine setperturb()

	implicit none
	include 'common.inc'
	
	real rrr(nx,ny),ranf_,rrrm
	integer i,j,k

	call ranset_(10)

	do k=1,nzm
	 rrrm = 0.
	 do j=1,ny
	  do i=1,nx
	    if(k.le.5) then
	      rrr(i,j) = 1.-2.*ranf_()
	      rrrm = rrrm + rrr(i,j)
	    endif
	  end do
	 end do
         rrrm = rrrm/dble(nx*ny)
	 do j=1,ny
	  do i=1,nx
	    if(k.le.5) then
	      t(i,j,k)=t(i,j,k)+0.02*(rrr(i,j)-rrrm)*(6-k)
	    endif
	  end do
	 end do
	end do

	return
	end
