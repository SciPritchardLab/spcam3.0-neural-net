
subroutine setperturb

!  Random noise

use vars

implicit none
	
integer i,j,k
real rrr,ranf_

call ranset_(30*rank)

do k=1,nzm
 do j=1,ny
  do i=1,nx
    rrr=1.-2.*ranf_()

    if(k.le.5) then
      t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
    endif

    if(k.le.4.and..not.dosmagor) then
      tke(i,j,k)=tke(i,j,k)+0.04*(5-k)
    endif

    if(doscalar) then
      tke(i,j,k) = q(i,j,k)
    end if
  end do
 end do
end do



end
