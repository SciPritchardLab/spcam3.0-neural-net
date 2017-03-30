  
       subroutine precip_init()

	implicit none
	include 'common.inc'
	include 'params.inc'
      
        real pratio, coef1, coef2,estw,esti,rrr1,rrr2
        real gamma1
        integer k 

	gam3 = 3. 
	gamr1 = 3.+b_rain
	gamr2 = (5.+b_rain)/2.
	gamr3 = 4.+b_rain
	gams1 = 3.+b_snow
	gams2 = (5.+b_snow)/2.
	gams3 = 4.+b_snow
	gamg1 = 3.+b_grau
	gamg2 = (5.+b_grau)/2.
	gamg3 = 4.+b_grau
	gam3 = gamma1(gam3) 
	gamr1 = gamma1(gamr1)
	gamr2 = gamma1(gamr2)
	gamr3 = gamma1(gamr3)
	gams1 = gamma1(gams1)
	gams2 = gamma1(gams2)
	gams3 = gamma1(gams3)
	gamg1 = gamma1(gamg1)
	gamg2 = gamma1(gamg2)
	gamg3 = gamma1(gamg3)

	do k=1,nzm
	
c	  pratio = (1000. / pres(k)) ** 0.4
	  pratio = sqrt(1.29 / rho(k))	

	  rrr1=393./(tabs0(k)+120.)*(tabs0(k)/273.)**1.5
	  rrr2=(tabs0(k)/273.)**1.94*(1000./pres(k))

	  estw = 100.*esatw_crm(tabs0(k))
	  esti = 100.*esati_crm(tabs0(k))

c	  accretion by snow:

          coef1 = 0.25 * pi * nzeros * a_snow * gams1 * pratio/
     &            (pi * rhos * nzeros/rho(k) ) ** ((3+b_snow)/4.)
          coef2 = exp(0.025*(tabs0(k) - 273.15))                        
          accrsi(k) =  coef1 * coef2 * esicoef
          accrsc(k) =  coef1 * esccoef 
          coefice(k) =  coef2
          
c         evaporation of snow:
 
          coef1  =(lsub/(tabs0(k)*rv)-1.)*lsub/(therco*rrr1*tabs0(k))
          coef2  = rv*tabs0(k)/(diffelq*rrr2*esti)
          evaps1(k)  =  0.65 * 4. * nzeros / 
     &        sqrt(pi * rhos * nzeros) / (coef1+coef2) / sqrt(rho(k)) 
          evaps2(k)  =  0.49 * 4. * nzeros * gams2 * 
     &        sqrt(a_snow/(muelq*rrr1))/ 
     &        (pi * rhos * nzeros)**((5+b_snow)/8.) / (coef1+coef2) 
     &               * rho(k)**((1+b_snow)/8.)*sqrt(pratio)  

c	  accretion by graupel:

          coef1 = 0.25 * pi * nzerog * a_grau * gamg1 * pratio/
     &            (pi * rhog * nzerog/rho(k) ) ** ((3+b_grau)/4.)
          coef2 = exp(0.025*(tabs0(k) - 273.15))                        
          accrgi(k) =  coef1 * coef2 * egicoef
          accrgc(k) =  coef1 * egccoef 
          
c         evaporation of graupel:
 
          coef1  =(lsub/(tabs0(k)*rv)-1.)*lsub/(therco*rrr1*tabs0(k))
          coef2  = rv*tabs0(k)/(diffelq*rrr2*esti)
          evapg1(k)  =  0.65 * 4. * nzerog / 
     &        sqrt(pi * rhog * nzerog) / (coef1+coef2) / sqrt(rho(k)) 
          evapg2(k)  =  0.49 * 4. * nzerog * gamg2 * 
     &        sqrt(a_grau/(muelq*rrr1))/ 
     &        (pi * rhog * nzerog)**((5+b_grau)/8.) / (coef1+coef2) 
     &               * rho(k)**((1+b_grau)/8.)*sqrt(pratio)  

 
c	  accretion by rain:

          accrrc(k)=  0.25 * pi * nzeror * a_rain * gamr1 * pratio/
     &              (pi * rhor * nzeror / rho(k)) ** ((3+b_rain)/4.)
     &		    * erccoef   

c         evaporation of rain:

          coef1  =(lcond/(tabs0(k)*rv)-1.)*lcond/(therco*rrr1*tabs0(k))
          coef2  = rv*tabs0(k)/(diffelq * rrr2 * estw)
          evapr1(k)  =  0.78 * 2. * pi * nzeror / 
     &        sqrt(pi * rhor * nzeror) / (coef1+coef2) / sqrt(rho(k)) 
          evapr2(k)  =  0.31 * 2. * pi  * nzeror * gamr2 * 
     &		0.89 * sqrt(a_rain/(muelq*rrr1))/ 
     &        (pi * rhor * nzeror)**((5+b_rain)/8.) / (coef1+coef2) 
     &              * rho(k)**((1+b_rain)/8.)*sqrt(pratio) 

c	if(rank.eq.0) then
c          print *, k,accrsi(k),accrsc(k),accrgi(k),accrgc(k),
c     &       accrrc(k),coefice(k),evaps1(k),evaps2(k),
c     &	          evapg1(k),evapg2(k),evapr1(k),evapr2(k)
c	end if        

	end do

           
!	if(rank.eq.0) then
!	 coef1=b_rain / 4.	   
!	 coef2=a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** coef1
!	 print*,'rain vel:', coef2*1.e-3**coef1*sqrt(1.2)
!	 coef1=b_snow / 4.	   
!	 coef2=a_snow * gams3 / 6. / (pi * rhos * nzeros) ** coef1
!	 print*,'snow vel:', coef2*1.e-3**coef1*sqrt(1.2)
!	 coef1=b_grau / 4.	   
!	 coef2=a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** coef1
!	 print*,'grau vel:', coef2*1.e-3**coef1*sqrt(1.2)
!	end if

      return
      end

	real function gamma1(arg1)
	
	real arg1
	real(4) gammafff, arg
	arg = arg1
	gamma1 = gammafff(arg)
	return
	
	end
