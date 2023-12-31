
subroutine cloud

!  Condensation of cloud water/cloud ice.

use vars
use params

implicit none

integer i,j,k, kb, kc, kmax, kmin
real dtabs, tabs1, an, bn, ap, bp, om, ag, omp
real fac1,fac2  
real fff,dfff,qsat,dqsat
real lstarn,dlstarn,lstarp,dlstarp
integer niter

an = 1./(tbgmax-tbgmin)	
bn = tbgmin * an
ap = 1./(tprmax-tprmin)	
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
ag = 1./(tgrmax-tgrmin)	

kmax=0
kmin=nzm+1


do k = 1, nzm
 do j = 1, ny
  do i = 1, nx

    q(i,j,k)=max(0.,q(i,j,k))


! Initail guess for temperature assuming no cloud water/ice:


    tabs(i,j,k) = t(i,j,k)-gamaz(k)
    tabs1=(tabs(i,j,k)+fac1*qp(i,j,k))/(1.+fac2*qp(i,j,k))

! Warm cloud:

    if(tabs1.ge.tbgmax) then

      tabs1=tabs(i,j,k)+fac_cond*qp(i,j,k)
      qsat = qsatw_crm(tabs1,pres(k))

! Ice cloud:

    elseif(tabs1.le.tbgmin) then

      tabs1=tabs(i,j,k)+fac_sub*qp(i,j,k)
      qsat = qsati_crm(tabs1,pres(k))

! Mixed-phase cloud:

    else

      om = an*tabs1-bn
      qsat = om*qsatw_crm(tabs1,pres(k))+(1.-om)*qsati_crm(tabs1,pres(k))

    endif


!  Test if condensation is possible:


    if(q(i,j,k).gt.qsat) then

      niter=0
      dtabs = 100.
      do while(abs(dtabs).gt.0.01.and.niter.lt.10)
	if(tabs1.ge.tbgmax) then
	   om=1.
	   lstarn=fac_cond
	   dlstarn=0.
	   qsat=qsatw_crm(tabs1,pres(k))
	   dqsat=dtqsatw_crm(tabs1,pres(k))
        else if(tabs1.le.tbgmin) then
	   om=0.
	   lstarn=fac_sub
	   dlstarn=0.
	   qsat=qsati_crm(tabs1,pres(k))
	   dqsat=dtqsati_crm(tabs1,pres(k))
	else
	   om=an*tabs1-bn
	   lstarn=fac_cond+(1.-om)*fac_fus
	   dlstarn=an*fac_fus
	   qsat=om*qsatw_crm(tabs1,pres(k))+(1.-om)*qsati_crm(tabs1,pres(k))
	   dqsat=om*dtqsatw_crm(tabs1,pres(k))+(1.-om)*dtqsati_crm(tabs1,pres(k))
	endif
	if(tabs1.ge.tprmax) then
	   omp=1.
	   lstarp=fac_cond
	   dlstarp=0.
        else if(tabs1.le.tprmin) then
	   omp=0.
	   lstarp=fac_sub
	   dlstarp=0.
	else
	   omp=ap*tabs1-bp
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   dlstarp=ap*fac_fus
	endif
	fff = tabs(i,j,k)-tabs1+lstarn*(q(i,j,k)-qsat)+lstarp*qp(i,j,k)
	dfff=dlstarn*(q(i,j,k)-qsat)+dlstarp*qp(i,j,k)-lstarn*dqsat-1.
	dtabs=-fff/dfff
	niter=niter+1
	tabs1=tabs1+dtabs
      end do   

      qsat = qsat + dqsat * dtabs
      qn(i,j,k) = max(0.,q(i,j,k)-qsat)
      if(tabs1.lt.tbgmax) then
        kmin = min(kmin,k)
        kmax = max(kmax,k)
      end if

    else

      qn(i,j,k) = 0.

    endif

    tabs(i,j,k) = tabs1
    qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case

  end do
 end do
end do


end subroutine cloud

