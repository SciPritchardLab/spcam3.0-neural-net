c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------

	subroutine buoyancy()

	implicit none
	include 'common.inc'
	integer i,j,k,kb


	do k=2,nzm	
         kb=k-1
	 do j=1,ny
	  do i=1,nx
	   dwdt(i,j,k,na)=dwdt(i,j,k,na)+0.5*(bet(k)*
     &      (tabs0(k)*(0.61*(q(i,j,k)-q0(k))-1.61*qn(i,j,k)-qp(i,j,k))
     &         +(tabs(i,j,k)-tabs0(k))*(1.+0.61*q0(k))) + bet(kb)*
     &   (tabs0(kb)*(0.61*(q(i,j,kb)-q0(kb))-1.61*qn(i,j,kb)-qp(i,j,kb))
     &         +(tabs(i,j,kb)-tabs0(kb))*(1.+0.61*q0(kb))) 
     &          - betp(k)*p(i,j,k) - betp(kb)*p(i,j,kb))
	  end do ! i
	 end do ! j
	end do ! k

	 
 
	return
	end


