subroutine stepout(nstatsteps)

use vars
use tracers
use params
!use hbuffer
!use isccp, only : isccp_write
implicit none	
	
integer i,j,k,ic,jc,nstatsteps
real div, divmax, divmin
real rdx, rdy, rdz, coef
integer im,jm,km
real wmax, qnmax(1), qnmax1(1)
real buffer(5), buffer1(5)

if(masterproc) print *,'NSTEP = ',nstep,'    NCYCLE=',ncycle

if(mod(nstep,nprint).eq.0) then
	

 divmin=1.e20
 divmax=-1.e20
	 
 rdx = 1./dx
 rdy = 1./dy

 wmax=0.
 do k=1,nzm
  coef = rho(k)*adz(k)*dz
  rdz = 1./coef
  if(ny.ne.1) then
   do j=1,ny-1*YES3D
    jc = j+1*YES3D
    do i=1,nx-1
     ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx + (v(i,jc,k)-v(i,j,k))*rdy + &
		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
   end do
  else
    j = 1
    do i=1,nx-1
    ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx +(w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
  endif
 end do

 if(dompi) then
   buffer(1) = total_water_before
   buffer(2) = total_water_after
   buffer(3) = total_water_evap
   buffer(4) = total_water_prec
   buffer(5) = total_water_ls
   call task_sum_real(buffer, buffer1,5)
   total_water_before = buffer1(1)
   total_water_after = buffer1(2)
   total_water_evap = buffer1(3)
   total_water_prec = buffer1(4)
   total_water_ls = buffer1(5)
 end if

         j=nx/32
         write(6,*)
         write(6,'(32f4.0)')((u(i,1,k),i=1,nx,j),k=nzm,1,-1)
         write(6,*)
         write(6,'(32f4.0)')((v(i,1,k),i=1,nx,j),k=nzm,1,-1)
         write(6,*)
         write(6,'(32f4.0)')((w(i,1,k),i=1,nx,j),k=nzm,1,-1)
         write(6,*)
         write(6,'(32f4.0)')((tk(i,1,k),i=1,nx,j),k=nzm,1,-1)
         write(6,*)
         write(6,'(32f4.0)')((t(i,1,k)-t0(k),i=1,nx,j),k=nzm,1,-1)
         write(6,*)
         write(6,'(32f4.0)')((q(i,1,k)*1.e3,i=1,nx,j),k=nzm,1,-1)
         write(6,*)


!--------------------------------------------------------
 if(masterproc) then
	
    print*,'DAY = ',day	
    write(6,*) 'NSTEP=',nstep
    write(6,*) 'div:',divmax,divmin
    write(6,*) 'SST=',tabs_s, '  surface pressure=',pres0

 endif

 call fminmax_print('u:',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
 call fminmax_print('v:',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm-5)
 call fminmax_print('w:',w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
 call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
 call fminmax_print('t:',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tabs:',tabs,1,nx,1,ny,nzm)
 call fminmax_print('q:',q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('qn:',qn,1,nx,1,ny,nzm)
 call fminmax_print('qp:',qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',misc*86400.,1,nx,1,ny,nzm)
 if(dotracers) then
   do k=1,ntracers
     call fminmax_print(trim(tracername(k))//':',tracer(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
   end do
 end if
 call fminmax_print('shf:',fluxbt*cp*rho(1),1,nx,1,ny,1)
 call fminmax_print('lhf:',fluxbq*lcond*rho(1),1,nx,1,ny,1)
 call fminmax_print('uw:',fluxbu,1,nx,1,ny,1)
 call fminmax_print('vw:',fluxbv,1,nx,1,ny,1)
 call fminmax_print('sst:',sstxy,1,nx,1,ny,1)

 total_water_before = total_water_before/float(nx_gl*ny_gl)
 total_water_after = total_water_after/float(nx_gl*ny_gl)
 total_water_evap = total_water_evap/float(nx_gl*ny_gl)
 total_water_prec = total_water_prec/float(nx_gl*ny_gl)
 total_water_ls = total_water_ls/float(nx_gl*ny_gl)
 
 if(masterproc) then
   
   print*,'total water budget:'
   print*,'before (mm):',total_water_before
   print*,'after (mm) :',total_water_after
   print*,'evap (mm/day):',total_water_evap/dtn*86400.
   print*,'prec (mm/day):',total_water_prec/dtn*86400.
   print*,'ls (mm/day):',total_water_ls/dtn*86400.
   print*,'disbalance (mm/day)', &
     (total_water_after-(total_water_before+total_water_evap+total_water_ls-total_water_prec))/dtn * 86400.

 end if



end if ! (mod(nstep,nprint).eq.0)
	
end
