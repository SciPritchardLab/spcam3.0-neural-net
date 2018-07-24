c
c  Condensation of cloud water/cloud ice.
c



	subroutine cloud()
	
	implicit none	
	include 'common.inc'
	include 'params.inc'

	integer i,j,k,n,nitermax,kb,kc,kmax,kmin
	real dtabs, tabs1, an, bn, ap, bp, om, ag, omp
        real fac1,fac2  
	real fff,dfff,dqsat
	real lstarn,dlstarn,lstarp,dlstarp
        real coef,dqi,lat_heat,vt_ice
        real omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
        real fz(nx,ny,nz), rhofac(nz)
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

c
c Initail guess for temperature assuming no cloud water/ice:
c

	    tabs(i,j,k) = t(i,j,k)-gamaz(k)
	    tabs1=(tabs(i,j,k)+fac1*qp(i,j,k))/(1.+fac2*qp(i,j,k))

c      Warm cloud:

	    if(tabs1.ge.tbgmax) then

	      tabs1=tabs(i,j,k)+fac_cond*qp(i,j,k)
              qsat(i,j,k) = qsatw_crm(tabs1,pres(k))
	      omegan(i,j,k) = 1.

c      Ice cloud:

	    elseif(tabs1.le.tbgmin) then

	      tabs1=tabs(i,j,k)+fac_sub*qp(i,j,k)	     
              qsat(i,j,k) = qsati_crm(tabs1,pres(k))
	      omegan(i,j,k) = 0.

c      Mixed-phase cloud:

	    else

	      om = an*tabs1-bn
              qsat(i,j,k) = om*qsatw_crm(tabs1,pres(k))+
     &			 (1.-om)*qsati_crm(tabs1,pres(k))
	      omegan(i,j,k) = om

	    endif

c
c  Test if condensation is possible:
c

	    if(q(i,j,k).gt.qsat(i,j,k)) then

	      niter=0
	      dtabs = 100.
	      do while(abs(dtabs).gt.0.01.and.niter.lt.10)
		if(tabs1.ge.tbgmax) then
		   om=1.
		   lstarn=fac_cond
		   dlstarn=0.
		   qsat(i,j,k)=qsatw_crm(tabs1,pres(k))
		   dqsat=dtqsatw_crm(tabs1,pres(k))
	        else if(tabs1.le.tbgmin) then
		   om=0.
		   lstarn=fac_sub
		   dlstarn=0.
		   qsat(i,j,k)=qsati_crm(tabs1,pres(k))
		   dqsat=dtqsati_crm(tabs1,pres(k))
		else
		   om=an*tabs1-bn
		   lstarn=fac_cond+(1.-om)*fac_fus
		   dlstarn=an
		   qsat(i,j,k)=om*qsatw_crm(tabs1,pres(k))+
     &			(1.-om)*qsati_crm(tabs1,pres(k))
		   dqsat=om*dtqsatw_crm(tabs1,pres(k))+
     &			(1.-om)*dtqsati_crm(tabs1,pres(k))
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
		   dlstarp=ap
		endif
		fff = tabs(i,j,k)-tabs1+
     &                lstarn*(q(i,j,k)-qsat(i,j,k))+lstarp*qp(i,j,k)
		dfff=dlstarn*(q(i,j,k)-qsat(i,j,k))+dlstarp*qp(i,j,k)-
     &                                 lstarn*dqsat-1.
		dtabs=-fff/dfff
		niter=niter+1
		tabs1=tabs1+dtabs
	      end do   

     	      qsat(i,j,k) = qsat(i,j,k) + dqsat * dtabs
              qn(i,j,k) = max(0.,q(i,j,k)-qsat(i,j,k))
	      omegan(i,j,k) = om
	      omegap(i,j,k) = omp
	      if(om.lt.1.) then
                 kmin = min(kmin,k)
	         kmax = max(kmax,k)
              end if
	    else

	      qn(i,j,k) = 0.
	      omegap(i,j,k) = max(0.,min(1.,ap*(tabs(i,j,k)-tprmin)))

	    endif

	    tabs(i,j,k) = tabs1
	    omegag(i,j,k) = max(0.,min(1.,ag*(tabs(i,j,k)-tgrmin)))
            qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case

	  end do
	 end do
	end do

c
c Sedimentation of ice:
c
! Compute cloud ice flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002).

        fz = 0.

	do k = max(1,kmin-1),kmax
           rhofac(k) = sqrt(1.29/rho(k))
	   ! Set up indices for x-y planes above and below current plane.
   	kc = min(nzm,k+1)
   	kb = max(1,k-1)
   	! CFL number based on grid spacing interpolated to interface i,j,k-1/2
   	coef = dtn/(0.5*(adz(kb)+adz(k))*dz)
   	do j = 1,ny
      	 do i = 1,nx
         ! Compute cloud ice density in this cell and the ones above/below.
         ! Since cloud ice is falling, the above cell is u (upwind),
         ! this cell is c (center) and the one below is d (downwind).
           omnu = max(0.,min(1.,(tabs(i,j,kc)-tbgmin)*a_bg))
           omnc = max(0.,min(1.,(tabs(i,j,k) -tbgmin)*a_bg))
           omnd = max(0.,min(1.,(tabs(i,j,kb)-tbgmin)*a_bg))

           qiu = rho(kc)*qn(i,j,kc)*(1.-omnu)
           qic = rho(k) *qn(i,j,k) *(1.-omnc)
           qid = rho(kb)*qn(i,j,kb)*(1.-omnd)

! Ice sedimentation velocity depends on ice content. The fiting is
! based on the data by Heymsfield (JAS,2003). -Marat
               vt_ice = 0.1*rhofac(k)
!                vt_ice = max(0.1,0.5*log10(qic+1.e-12)+3.)
!                vt_ice = 8.66*(qic+1.e-12)**0.24   ! Heymsfield (JAS, 2003, p.2607)

         ! Use MC flux limiter in computation of flux correction.
         ! (MC = monotonized centered difference).
           if (qic.eq.qid) then
            tmp_phi = 0.
           else
             tmp_theta = (qiu-qic)/(qic-qid)
             tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
           end if

         ! Compute limited flux.
         ! Since falling cloud ice is a 1D advection problem, this
         ! flux-limited advection scheme is monotonic.
           fz(i,j,k) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid))
      	 end do
   	end do
	end do
	fz(:,:,nz) = 0.

 	do k=max(1,kmin-2),kmax
   	coef=dtn/(dz*adz(k)*rho(k))
   	do j=1,ny
      	do i=1,nx
         ! The cloud ice increment is the difference of the fluxes.
         dqi=coef*(fz(i,j,k)-fz(i,j,k+1))
         ! Add this increment to both non-precipitating and total water.
         qn(i,j,k) = qn(i,j,k) + dqi
         q(i,j,k)  = q(i,j,k)  + dqi

         ! The latent heat flux induced by the falling cloud ice enters
         ! the liquid-ice static energy budget in the same way as the
         ! precipitation.  Note: use latent heat of sublimation.
         lat_heat  = (fac_cond+fac_fus)*dqi
         ! Add divergence of latent heat flux to liquid-ice static energy.
         t(i,j,k)  = t(i,j,k)  - lat_heat
      	end do
   	end do
	end do

        coef=dtn/dz
        do j=1,ny
        do i=1,nx
          dqi=-coef*fz(i,j,1)
	  precsfc(i,j) = precsfc(i,j)+dqi
	  precssfc(i,j) = precssfc(i,j)+dqi
	 end do
	end do


	return
	end

