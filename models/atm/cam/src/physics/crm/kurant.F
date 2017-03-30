
c-----------------------------------------------------------------
c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------

	subroutine kurant(stop_signal)

   	implicit none
	include 'common.inc'

        integer lchnk, icol
	integer i, j, k, ncycle1(1),ncycle2(1)
        real wm(nz)
        real um(nz)
        real vm(nz)
        real tkhmax(nz)
	real cfl
        logical stop_signal
        character(5) text

	ncycle = 1
        stop_signal = .false.
	
	do k = 1,nzm
	 wm(k) = 0.
	 um(k) = 0.
	 vm(k) = 0.
	 tkhmax(k) = 0.
	 do j = 1,ny
	  do i = 1,nx
	   wm(k) = max(wm(k),abs(w(i,j,k)))
	   um(k) = max(um(k),abs(u(i,j,k)))
	   tkhmax(k)=max(tkhmax(k),tkh(i,j,k))
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
	 cfl = max(cfl,
     &		wm(k)*dt/(dz*adzw(k)) 
     &          ,0.5*tkhmax(k)*dt/(dz*adzw(k))**2
     &		,um(k)*dt/dx, 0.5*tkhmax(k)*dt/dx**2
     &		,vm(k)*dt/dy, 0.5*tkhmax(k)*dt/dy**2 
     &	 )
	end do
	
 10	continue
 
	if(cfl .gt. 0.7) then
	 ncycle = ncycle + 1
	 if(ncycle.gt.4) then
	   nprint=1
           print*,'lchnk=',lchnk,'  icol=',icol
	   call stepout()
           print *,'cfl=',cfl,' the number of cycles exceeded 4.'
           stop_signal = .true.
           return
	 end if
	 cfl = cfl / ncycle
	 goto 10 
	end if 
	
	if(dompi) then
	  ncycle1(1)=ncycle
	  call task_max_integer(ncycle1,ncycle2,1)
	  ncycle=ncycle2(1)
	end if
	return
	end	
	
