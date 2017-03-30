	subroutine diagnose()
	

	implicit none
	include 'common.inc'
	
	include 'params.inc'
	integer i,j,k
	real *8 coef

	coef = 1./float(nx*ny)
	
        do k=1,nzm
          u0(k)=0.
          v0(k)=0.
          t0(k)=0.
          tabs0(k)=0.
          q0(k)=0.
          p0(k)=0.
          do j=1,ny
           do i=1,nx
           
	     tabs(i,j,k) = t(i,j,k)-gamaz(k)
     &	    + (fac_cond+(1.-omegan(i,j,k))*fac_fus)*qn(i,j,k)+
     &	      (fac_cond+(1.-omegap(i,j,k))*fac_fus)*qp(i,j,k)
	     u(i,j,k) = dudt(i,j,k,nc)
	     v(i,j,k) = dvdt(i,j,k,nc)
	     w(i,j,k) = dwdt(i,j,k,nc)

             u0(k)=u0(k)+u(i,j,k)
             v0(k)=v0(k)+v(i,j,k)
             p0(k)=p0(k)+p(i,j,k)
             t0(k)=t0(k)+t(i,j,k)
             tabs0(k)=tabs0(k)+tabs(i,j,k)
             q0(k)=q0(k)+q(i,j,k)

           end do
          end do
          u0(k)=u0(k)*coef
          v0(k)=v0(k)*coef
          t0(k)=t0(k)*coef
          tabs0(k)=tabs0(k)*coef
          q0(k)=q0(k)*coef
          p0(k)=p0(k)*coef
        end do ! k

	return
	end
