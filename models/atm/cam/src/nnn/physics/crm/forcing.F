
	subroutine forcing()
	
	implicit none
	include 'common.inc'

	real coef,qneg,qpoz, factor
	integer i,j,k,nneg

	coef = 1./7200.

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
!              dudt(i,j,k,na)=dudt(i,j,k,na) -(u0(k)-ug0(k))*coef
!              dvdt(i,j,k,na)=dvdt(i,j,k,na) -(v0(k)-vg0(k))*coef
	      dudt(i,j,k,na)=dudt(i,j,k,na) + utend(k)
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

!
! Nudge the upper level temperature and vapor to the observe
!

!        coef = dtn / 3600.
!        do k=nzm-2,nzm
!          do j=1,ny
!            do i=1,nx
!              t(i,j,k)=t(i,j,k)-(t(i,j,k)-tg0(k))*coef
!              q(i,j,k)=q(i,j,k)-(q(i,j,k)-qg0(k))*coef
!            end do
!          end do
!        end do


	return
	end
