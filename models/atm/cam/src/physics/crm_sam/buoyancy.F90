
subroutine buoyancy()

use vars
use params
implicit none
	
integer i,j,k,kb

if(docolumn) return


do k=2,nzm	
 kb=k-1
 do j=1,ny
  do i=1,nx
   dwdt(i,j,k,na)=dwdt(i,j,k,na)  &
       +0.5*(bet(k)* &
              (tabs0(k)*(0.61*(q(i,j,k)-q0(k))-1.61*qn(i,j,k)-qp(i,j,k)) &
             +(tabs(i,j,k)-tabs0(k))*(1.+0.61*q0(k))) &
           + bet(kb)* &
              (tabs0(kb)*(0.61*(q(i,j,kb)-q0(kb))-1.61*qn(i,j,kb)-qp(i,j,kb)) &
             +(tabs(i,j,kb)-tabs0(kb))*(1.+0.61*q0(kb)))  )
  end do ! i
 end do ! j
end do ! k

	 
 
end subroutine buoyancy


