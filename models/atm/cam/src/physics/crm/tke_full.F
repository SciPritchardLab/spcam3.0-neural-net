c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------
c	this subroutine solves the TKE equation
c-----------------------------------------------------------------

	subroutine tke_full()

	implicit none
	include 'common.inc'
	include 'params.inc'

        real def2(nx,ny,nzm)	
	real grid,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
	real buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss,tkk
        real lstarn, lstarp, bbb, omn, omp
	real dqsat
	integer i,j,k,kc,kb

	Cs = 0.15
	Ck=0.1
	Ce=Ck**3/Cs**4
	Ces=Ce/0.7*3.0	

 	if(RUN3D) then
	  call shear_prod3D(def2)
	else
	  call shear_prod2D(def2)
	endif

	do k=1,nzm      
	  kb=k-1
	  kc=k+1

	  grid=dz*adz(k)

	  betdz=bet(k)/dz/(adzw(kc)+adzw(k))
	  Ce1=Ce/0.7*0.19
	  Ce2=Ce/0.7*0.51
	  if(k.eq.1) then
	    kb=1
	    kc=2
	    betdz=bet(k)/dz/adzw(kc)
	    Ce1=Ces/0.7*0.19
	    Ce2=Ces/0.7*0.51
	  end if
	  if(k.eq.nzm) then
	    kb=nzm-1
	    kc=nzm
	    betdz=bet(k)/dz/adzw(k)
	    Ce1=Ces/0.7*0.19
	    Ce2=Ces/0.7*0.51
	  end if
	  do j=1,ny
	  do i=1,nx
c  SGS buoyancy flux

           omn = omegan(i,j,k)
           omp = omegap(i,j,k)
           lstarn = fac_cond+(1.-omn)*fac_fus
           lstarp = fac_cond+(1.-omp)*fac_fus

!           if(qn(i,j,kb)+qn(i,j,k)+qn(i,j,kc) .gt. 0.) then
           if(qn(i,j,k) .gt. 0.) then
     
	      dqsat = omn*dtqsatw_crm(tabs(i,j,k),pres(k))+
     &		(1.-omn)*dtqsati_crm(tabs(i,j,k),pres(k))
     	      bbb = 1. + 0.61*qsat(i,j,k)-qn(i,j,k)
     &		  -qp(i,j,k)+1.61*tabs(i,j,k)*dqsat
              bbb = bbb / (1.+lstarn*dqsat)
              buoy_sgs=betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) 
     &		+(bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))*
     &                                    (q(i,j,kc)-q(i,j,kb))
     &		+(bbb*lstarp - (1.+lstarp*dqsat)*tabs(i,j,k))*
     &                                    (qp(i,j,kc)-qp(i,j,kb)) )
	   else

     	      bbb = 1.+0.61*q(i,j,k)-qp(i,j,k)
	      buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb))
     &          +0.61*tabs(i,j,k)*(q(i,j,kc)-q(i,j,kb))
     &		+(bbb*lstarp-tabs(i,j,k))*(qp(i,j,kc)-qp(i,j,kb)) )    	     

	   end if


	   if(buoy_sgs.le.1.e-10) then
	     smix=grid
	   else
	     smix=min(grid,
     $         max(0.1*grid,sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs))))
	   end if

	   ratio=smix/grid
	   Pr=1.   !+2.*ratio
	   Cee=Ce1+Ce2*ratio

	   tk(i,j,k)=sqrt(Ck**3/Cee*max(0.,def2(i,j,k)-Pr*buoy_sgs))*
     &                smix**2
	   a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
	   a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
	   a_diss=a_prod_sh+a_prod_bu

	   tkh(i,j,k)=Pr*tk(i,j,k)
           tke(i,j,k) =(tk(i,j,k)/(Ck*smix))**2

	  end do ! i
	  end do ! j

	end do ! k

c----------------------------------------------------------------
	return
	end


