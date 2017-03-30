c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c       pressure solver
c----------------------------------------------------------------


	subroutine pressure()
	
	implicit none
	include 'common.inc'
	
	integer nzslab
	parameter (nzslab = nzm / npressureslabs)
	integer nx2, ny2
	parameter (nx2=nx_gl+2,ny2=ny_gl+2*YES3D)
	integer n3i, n3j
	parameter (n3i=3*nx_gl/2+1,n3j=3*ny_gl/2+1)

	real f(nx2,ny2,nzslab) ! global rhs and array for FTP coefficeients
	real ff(nx+1,ny+2*YES3D,nzm)	! local (subdomain's) version of f
	real buff_slabs(nxp1,nyp2,nzslab,npressureslabs)
	real buff_subs(nxp1,nyp2,nzslab,nsubdomains) 
	real bufp_slabs(0:nx,1-YES3D:ny,nzslab,npressureslabs)  
	real bufp_subs(0:nx,1-YES3D:ny,nzslab,nsubdomains)  
	common/tmpstack/f,ff,buff_slabs,buff_subs
!tBM*   threadlocal /tmpstack/ 
!$OMP   threadprivate (/tmpstack/) 
	equivalence (buff_slabs,bufp_slabs)
	equivalence (buff_subs,bufp_subs)

	real work(nx2,ny2),trigxi(n3i),trigxj(n3j) ! FFT stuff
	integer ifaxj(100),ifaxi(100)

	real *8 a(nzm),b(nxp1,nyp2),c(nzm),e	
	real *8 xi,xj,xnx,xny,ddx2,ddy2,pii,eign(nxp1,nyp2)
	real alfa(nxp1,nyp2,nzm-1),beta(nxp1,nyp2,nzm-1)

	integer reqs_out(nsubdomains), reqs_in(nsubdomains)
	integer i, j, k, id, jd, m, n, it, jt, ii, jj, tag, rf
	integer nyp22, n_in, n_out, count
	integer iii(0:nx_gl),jjj(0:ny_gl)
        logical flag(nsubdomains)
c-----------------------------------------------------------------
	do i=1,nx_gl
	 iii(i)=i
	end do
	iii(0)=nx_gl
	do j=1,ny_gl
	 jjj(j)=j
	end do
	jjj(0)=ny_gl

	if(RUN2D) then  
	  nyp22=1
	else
	  nyp22=nyp2
	endif
c-----------------------------------------------------------------
c  Compute the r.h.s. of the Poisson equation for pressure

	call press_rhs()

c-----------------------------------------------------------------	 
c   Form the horizontal slabs of right-hand-sides of Poisson equation 
c   for the global domain. Request sending and receiving tasks.

	n_in = 0
	n_out = 0
	do m = 0,nsubdomains-1

	  if(m.lt.npressureslabs.and.m.ne.rank) then

	    n_out = n_out + 1
	    n = m*nzslab + 1
	    call task_send_float(m,p(0,1-YES3D,n),nzslab*nxp1*nyp1,
     &                                            33, reqs_out(n_out))
	  endif

	  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then

	    n_in = n_in + 1
	    call task_receive_float(bufp_subs(0,1-YES3D,1,n_in),
     &                        nzslab*nxp1*nyp1,reqs_in(n_in))
            flag(n_in) = .false.
 
	  endif

	  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

	    call task_rank_to_index(rank,it,jt)	  
	    n = rank*nzslab
	    do k = 1,nzslab
	     do j = 1,ny
	       do i = 1,nx
	         f(i+it,j+jt,k) = p(i,j,k+n)
	       end do
	     end do
	    end do
	  endif

	end do ! m

c
c Fill slabs when receive buffers are full:
c

	count = n_in
        do while (count .gt. 0)
	  do m = 1,n_in
	   if(.not.flag(m)) then
		call task_test(reqs_in(m), flag(m), rf, tag)
	        if(flag(m)) then 
		   count=count-1
	           if(tag.ne.33) then
	             print*,'MPI:Wrong message tag in poisson().'
	             print*,'    expected 33,  Received:',tag
	             call task_abort() 
	           endif
	           call task_rank_to_index(rf,it,jt)	  
	           do k = 1,nzslab
	            do j = 1,ny
	             do i = 1,nx
	               f(i+it,j+jt,k) = bufp_subs(i,j,k,m)
	             end do
	            end do
	           end do
		endif   
	   endif
	  end do
	end do

c-------------------------------------------------
c Perform Fourier transformation for a slab:

 	if(rank.lt.npressureslabs) then

	 call fftfax_crm(nx_gl,ifaxi,trigxi)
	 if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)

	 do k=1,nzslab
	
	  call fft991_crm
     &             (f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
	  if(RUN3D)
     $	  call fft991_crm
     &            (f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
	 
	 end do 

	endif

c
c Wait until all pending buffers are sent:
c
        call task_waitall(n_out,reqs_out,reqs_in,reqs_in)

c
c Synchronize all slabs:
c
        call task_barrier()

c-------------------------------------------------
c   Send Fourier coeffiecients back to subdomains:

	n_in = 0
	n_out = 0

	do m = 0, nsubdomains-1
		
	   call task_rank_to_index(m,it,jt)

	   if(rank.lt.npressureslabs.and.m.ne.rank) then

	     n_out = n_out + 1 
	     do k = 1,nzslab
	      do j = 1,nyp22
	       do i = 1,nxp1
	         buff_subs(i,j,k,n_out) = f(i+it,j+jt,k) 
	       end do
	      end do
	     end do 

	     call task_send_float(m, buff_subs(1,1,1,n_out), 
     &	                         nzslab*nxp1*nyp22,44,reqs_out(n_out))
	     
	   endif

	   if(rank.lt.npressureslabs.and.m.eq.rank) then

	     n = rank*nzslab
	     do k = 1,nzslab
	      do j = 1,nyp22
	        do i = 1,nxp1
	          ff(i,j,k+n) = f(i+it,j+jt,k) 
	        end do
	      end do
	     end do 
	
	   end if

	   if(m.lt.npressureslabs-1.or.
     &		m.eq.npressureslabs-1.and.rank.ge.npressureslabs) then

 	     n_in = n_in + 1
	     call task_receive_float(buff_slabs(1,1,1,n_in),
     &                                nzslab*nxp1*nyp22,reqs_in(n_in))
	     flag(n_in) = .false.	    
	   endif

	end do ! m

c
c Fill slabs when receive buffers are complete:
c

	count = n_in
        do while (count .gt. 0)
	  do m = 1,n_in
	   if(.not.flag(m)) then
		call task_test(reqs_in(m), flag(m), rf, tag)
	        if(flag(m)) then 
		   count=count-1
	           if(tag.ne.44) then
	             print*,'MPI:Wrong message tag in pressure().'
	             print*,'    expected 44   Received:',tag
	             call task_abort() 
	           endif
	           n = rf*nzslab           
                   do k = 1,nzslab
                    do j=1,nyp22
                     do i=1,nxp1
                        ff(i,j,k+n) = buff_slabs(i,j,k,m)
                     end do
                    end do
                   end do
		endif   
	   endif
	  end do
	end do

c-------------------------------------------------
c   Solve the tri-diagonal system for Fourier coeffiecients 
c   in the vertical for each subdomain:

	call task_rank_to_index(rank,it,jt)
	
	ddx2=dx*dx
	ddy2=dy*dy
	pii = dacos(-1.d0)
	do j=1,nyp22
	 do i=1,nxp1
	  id=(i+it-0.1)/2.
	  jd=(j+jt-0.1)/2.
	  xi=id
	  xj=jd
	  xnx=nx_gl
	  xny=ny_gl 	 
	  eign(i,j)=(2.d0*dcos(2.d0*pii/xnx*xi)-2.d0)/ddx2+
     $              (2.d0*dcos(2.d0*pii/xny*xj)-2.d0)/ddy2
	 end do 
	end do 

  	 do k=1,nzm
	    a(k)=rhow(k)/(adz(k)*adzw(k)*dz*dz)
	    c(k)=rhow(k+1)/(adz(k)*adzw(k+1)*dz*dz)	 
	 end do 

	  do j=1,nyp22
	    jd=(j-0.1)/2.
	   do i=1,nxp1	 
	     id=(i-0.1)/2.	   
	     if(id+jd.eq.0.and.rank.eq.0) then
	       b(i,j)=eign(i,j)*rho(1)-a(1)-c(1)
	       alfa(i,j,1)=-c(1)/b(i,j)
	       beta(i,j,1)=ff(i,j,1)/b(i,j)
	     else
	       b(i,j)=eign(i,j)*rho(1)-c(1)
	       alfa(i,j,1)=-c(1)/b(i,j)
	       beta(i,j,1)=ff(i,j,1)/b(i,j)
	     end if
	   end do
	  end do

	  do k=2,nzm-1
            do j=1,nyp22
	      do i=1,nxp1	 
	        e=eign(i,j)*rho(k)-a(k)-c(k)+a(k)*alfa(i,j,k-1)
	        alfa(i,j,k)=-c(k)/e
	        beta(i,j,k)=(ff(i,j,k)-a(k)*beta(i,j,k-1))/e
	      end do
	    end do
	  end do

          do j=1,nyp22
	    do i=1,nxp1	 
	     ff(i,j,nzm)=(ff(i,j,nzm)-a(nzm)*beta(i,j,nzm-1))/
     &	        (eign(i,j)*rho(nzm)-a(nzm)+a(nzm)*alfa(i,j,nzm-1))
	    end do
	  end do
	  
	  do k=nzm-1,1,-1
           do j=1,nyp22
	    do i=1,nxp1	 	  
	       ff(i,j,k)=alfa(i,j,k)*ff(i,j,k+1)+beta(i,j,k)
	    end do
	   end do  
	  end do 

c
c Wait until all pending buffers are sent:
c
        call task_waitall(n_out,reqs_out,reqs_in,reqs_in)


	call task_barrier()

c-----------------------------------------------------------------	 
c   Send the Fourier coefficient to the tasks performing
c   the inverse Fourier transformation:

        n_in = 0
	n_out = 0

	do m = 0,nsubdomains-1
	
	  if(m.lt.npressureslabs.and.m.ne.rank) then

	    n_out = n_out + 1
	    n = m*nzslab+1
	    call task_send_float(m,ff(1,1,n),nzslab*nxp1*nyp22,
     &                                          33, reqs_out(n_out))
	  endif

	  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then

	    n_in = n_in + 1
	    call task_receive_float(buff_subs(1,1,1,n_in),
     &                            nzslab*nxp1*nyp22, reqs_in(n_in))
	    flag(n_in) = .false.	    
	  endif

	  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

	    call task_rank_to_index(rank,it,jt)	  
	    n = rank*nzslab
	    do k = 1,nzslab
	     do j = 1,nyp22
	       do i = 1,nxp1
	         f(i+it,j+jt,k) = ff(i,j,k+n)
	       end do
	     end do
	    end do

	  endif

	end do ! m

c
c Fill slabs when receive buffers are full:
c

	count = n_in
        do while (count .gt. 0)
	  do m = 1,n_in
	   if(.not.flag(m)) then
		call task_test(reqs_in(m), flag(m), rf, tag)
	        if(flag(m)) then 
		   count=count-1
	           if(tag.ne.33) then
	             print*,'MPI:Wrong message tag in poisson().'
	             print*,'    expected 33,  Received:',tag
	             call task_abort() 
	           endif
	           call task_rank_to_index(rf,it,jt)	  
	           do k = 1,nzslab
	            do j = 1,nyp22
	             do i = 1,nxp1
	                f(i+it,j+jt,k) = buff_subs(i,j,k,m)
	             end do
	            end do
	           end do
		endif   
	   endif
	  end do
	end do

c-------------------------------------------------
c   Perform inverse Fourier transformation:

 	if(rank.lt.npressureslabs) then

	 do k=1,nzslab
	
	  if(RUN3D)
     $	  call fft991_crm
     &         (f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
	  call fft991_crm
     &           (f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
	 
	 end do 

	endif

c
c Wait until all pending buffers are sent:
c
        call task_waitall(n_out,reqs_out,reqs_in,reqs_in)

	call task_barrier()

c-----------------------------------------------------------------	 
c   Fill the pressure field for each subdomain: 

	n_in = 0
	n_out = 0

	do m = 0, nsubdomains-1
		
	   call task_rank_to_index(m,it,jt)

	   if(rank.lt.npressureslabs.and.m.ne.rank) then

	     n_out = n_out + 1
	     do k = 1,nzslab
	      do j = 1-YES3D,ny
	       jj=jjj(j+jt)
	       do i = 0,nx
		 ii=iii(i+it)
	         bufp_subs(i,j,k,n_out) = f(ii,jj,k) 
	       end do
	      end do
	     end do 

	     call task_send_float(m, bufp_subs(0,1-YES3D,1,n_out),
     &                             nzslab*nxp1*nyp1,44,reqs_out(n_out))

	   endif

	   if(m.lt.npressureslabs-1.or.
     &		m.eq.npressureslabs-1.and.rank.ge.npressureslabs) then

	     n_in = n_in + 1
	     call task_receive_float(bufp_slabs(0,1-YES3D,1,n_in),
     &                                nzslab*nxp1*nyp1, reqs_in(n_in))
	     flag(n_in) = .false.    

	   endif

	   if(rank.lt.npressureslabs.and.m.eq.rank) then

	     n = rank*nzslab
	     do k = 1,nzslab
	      do j = 1-YES3D,ny
	       jj=jjj(j+jt)
	        do i = 0,nx
		 ii=iii(i+it)
	          p(i,j,k+n) = f(ii,jj,k) 
	        end do
	      end do
	     end do 
	
	   end if

	end do ! m

	count = n_in
        do while (count .gt. 0)
	  do m = 1,n_in
	   if(.not.flag(m)) then
		call task_test(reqs_in(m), flag(m), rf, tag)
	        if(flag(m)) then 
		   count=count-1
	           if(tag.ne.44) then
	             print*,'MPI:Wrong message tag in pressure().'
	             print*,'    expected 44   Received:',tag
	             call task_abort() 
	           endif
	           n = rf*nzslab           
                   do k = 1,nzslab
                    do j=1-YES3D,ny
                     do i=0,nx
                         p(i,j,k+n) = bufp_slabs(i,j,k,m)
                     end do
                    end do
                   end do
		endif   
	   endif
	  end do
	end do

c
c Wait until all pending buffers are sent:
c
        call task_waitall(n_out,reqs_out,reqs_in,reqs_in)

	call task_barrier()

c
c  Add pressure gradient term to the rhs of the momentum equation:
c

	call press_grad()

	return
	end



