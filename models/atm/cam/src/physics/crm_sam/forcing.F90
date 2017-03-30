
	subroutine forcing()
	
        use vars
	implicit none

	real coef,qneg,qpoz, factor
	integer i,j,k,nneg

	coef = 1./3600.

	do k=1,nzm

	    qpoz = 0.
	    qneg = 0.
	    nneg = 0

	    do j=1,ny
	     do i=1,nx
	      t(i,j,k)=t(i,j,k) + ttend(k) * dtn
	      q(i,j,k)=q(i,j,k) + qtend(k) * dtn
	      if(q(i,j,k).lt.0.) then
	           nneg = nneg + 1
	           qneg = qneg + q(i,j,k)
	      else
	           qpoz = qpoz + q(i,j,k)
	      end if
	      dudt(i,j,k,na)=dudt(i,j,k,na) + utend(k)
!for 2d crm have added a term that damps perturbations to the mean on
!   a 1 hour time scale (Don Dazlich)
!	      dvdt(i,j,k,na)=dvdt(i,j,k,na) + vtend(k) - (1-YES3D)*coef*(v(i,j,k)-v0(k))
              dvdt(i,j,k,na)=dvdt(i,j,k,na) + vtend(k)
	     end do
	    end do

            if(nneg.gt.0.and.qpoz+qneg.gt.0.) then
             factor = 1. + qneg/qpoz
             do j=1,ny
              do i=1,nx
               q(i,j,k) = max(0.,q(i,j,k)*factor)
              end do
             end do
            end if

	end do


	return
	end
