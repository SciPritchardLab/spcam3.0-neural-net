subroutine diagnose
	
! Diagnose some useful stuff

use vars
use params
implicit none
	
integer i,j,k,kb,kc,k200,k500
double precision coef, coef1, buffer(nzm,7), buffer1(nzm,7)
real omn, omp

coef = 1./float(nx*ny)

k200 = nzm
	
do k=1,nzm
  u0(k)=0.
  v0(k)=0.
  t01(k) = tabs0(k)
  q01(k) = q0(k)
  t0(k)=0.
  tabs0(k)=0.
  q0(k)=0.
  p0(k)=0.
  tke0(k)=0.
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  if(pres(kc).le.200..and.pres(kb).gt.200.) k200=k
  coef1 = rho(k)*dz*adz(k)*dtfactor
  do j=1,ny
    do i=1,nx
     omn  = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
     omp  = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
     tabs(i,j,k) = t(i,j,k)-gamaz(k)+ &
	      (fac_cond+(1.-omn)*fac_fus)*qn(i,j,k)+ &
	      (fac_cond+(1.-omp)*fac_fus)*qp(i,j,k)
     u(i,j,k) = dudt(i,j,k,nc)
     v(i,j,k) = dvdt(i,j,k,nc)
     w(i,j,k) = dwdt(i,j,k,nc)

     u0(k)=u0(k)+u(i,j,k)
     v0(k)=v0(k)+v(i,j,k)
     p0(k)=p0(k)+p(i,j,k)
     t0(k)=t0(k)+t(i,j,k)
     tabs0(k)=tabs0(k)+tabs(i,j,k)
     q0(k)=q0(k)+q(i,j,k)
     tke0(k)=tke0(k)+tke(i,j,k)

     pw_xy(i,j) = pw_xy(i,j)+q(i,j,k)*coef1
     cw_xy(i,j) = cw_xy(i,j)+qn(i,j,k)*omn*coef1
     iw_xy(i,j) = iw_xy(i,j)+qn(i,j,k)*(1.-omn)*coef1
    end do
  end do
  u0(k)=u0(k)*coef
  v0(k)=v0(k)*coef
  t0(k)=t0(k)*coef
  tabs0(k)=tabs0(k)*coef
  q0(k)=q0(k)*coef
  p0(k)=p0(k)*coef
  tke0(k)=tke0(k)*coef
  rel0(k) = q0(k)/qsatw_crm(tabs0(k),pres(k))

end do ! k

k500 = nzm
do k = 1,nzm
   if((pres(kc).le.500.).and.(pres(k).gt.500.)) then
      if ((500.-pres(kc)).lt.(pres(k)-500.))then
         k500=kc
      else
         k500=k
      end if
   end if
end do


do j=1,ny
 do i=1,nx
  usfc_xy(i,j) = usfc_xy(i,j) + u(i,j,1)*dtfactor  
  vsfc_xy(i,j) = vsfc_xy(i,j) + v(i,j,1)*dtfactor  
  u200_xy(i,j) = u200_xy(i,j) + u(i,j,k200)*dtfactor  
  v200_xy(i,j) = v200_xy(i,j) + v(i,j,k200)*dtfactor  
  w500_xy(i,j) = w500_xy(i,j) + w(i,j,k500)*dtfactor
 end do
end do

if(dompi) then

  coef1 = 1./float(nsubdomains)
  do k=1,nzm
    buffer(k,1) = u0(k)
    buffer(k,2) = v0(k)
    buffer(k,3) = t0(k)
    buffer(k,4) = q0(k)
    buffer(k,5) = p0(k)
    buffer(k,6) = tabs0(k)
    buffer(k,7) = tke0(k)
  end do
  call task_sum_real8(buffer,buffer1,nzm*7)
  do k=1,nzm
    u0(k)=buffer1(k,1)*coef1
    v0(k)=buffer1(k,2)*coef1
    t0(k)=buffer1(k,3)*coef1
    q0(k)=buffer1(k,4)*coef1
    p0(k)=buffer1(k,5)*coef1
    tabs0(k)=buffer1(k,6)*coef1
    tke0(k)=buffer1(k,7)*coef1
  end do

end if ! dompi

!   recompute pressure levels:

!call pressz()

total_water_after = 0.
do k=1,nzm
 total_water_after = total_water_after + &
            (sum(q(1:nx,1:ny,k))+sum(qp(1:nx,1:ny,k)))*adz(k)*dz *rho(k)
end do
	
end subroutine diagnose
