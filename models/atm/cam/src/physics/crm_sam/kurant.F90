
subroutine kurant

use vars

implicit none

integer i, j, k, ncycle1(1),ncycle2(1)
real wm(nz)
real um(nz)
real vm(nz)
real tkhmax(nz)
real cfl

ncycle = 1
	
w_max=0.
do k = 1,nzm
 wm(k) = 0.
 tkhmax(k) = 0.
 um(k) = 0.
 vm(k) = 0.
 do j = 1,ny
  do i = 1,nx
   wm(k) = max(wm(k),abs(w(i,j,k)))
   tkhmax(k)=max(tkhmax(k),tkh(i,j,k))
   w_max=max(w_max,w(i,j,k))
   um(k) = max(um(k),abs(u(i,j,k)))
  end do
 end do
end do

if(RUN3D) then

 do k=1,nzm
  do j = 1,ny
   do i = 1,nx
    vm(k) = max(vm(k),abs(v(i,j,k)))
   end do
  end do
 end do

endif

cfl = 0.
do k=1,nzm
 cfl = max(cfl	&
     ,wm(k)*dt/(dz*adzw(k)), 0.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2 &
     ,um(k)*dt/dx, 0.5*tkhmax(k)*grdf_x(k)*dt/dx**2 &
     ,vm(k)*dt/dy, 0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
end do
	
ncycle = max(1,ceiling(cfl/0.7))

if(dompi) then
  ncycle1(1)=ncycle
  call task_max_integer(ncycle1,ncycle2,1)
  ncycle=ncycle2(1)
end if

if(ncycle.gt.4) then
   if(masterproc) print *,'the number of cycles exceeded 4.'
   call task_abort()
end if

end subroutine kurant	
