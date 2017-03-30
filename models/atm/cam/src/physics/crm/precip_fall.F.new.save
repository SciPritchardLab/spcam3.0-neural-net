	subroutine precip_fall

!     positively definite monotonic advection with non-oscillatory option
!     and gravitational sedimentation (rain and snow advection)

	implicit none
        include 'common.inc'
        include 'params.inc'

 	
	real mx(nzm),mn(nzm)
	real www(nz),fz(nz)
	real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
	real f0(nzm),df0(nzm)
	real eps
	integer i,j,k,kc,kb
	logical nonos
	real irho

	real x1, x2, a, b, a1, a2, y 
	real rhoz
	real vrain, vsnow, vgrau, crain, csnow, cgrau
	real lfac(nz), qrr, qss, qgg
	real andiff,pp,pn, lat_heat
	real wmax, omp, omg

	real wp(nzm), tmp_qp(nzm), irhoadz(nzm), iwmax(nzm), rhofac(nzm), prec_cfl
	integer nprec, iprec

	andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
	pp(y)= max(0.,y)
	pn(y)=-min(0.,y)


!--------------------------------------------------------


	eps = 1.e-10
	nonos = .true.
  
	crain = b_rain / 4.
	csnow = b_snow / 4. 
	cgrau = b_grau / 4. 
	vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain	  
	vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
	vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau


	 do k = 1,nzm
	    rhofac(k) = sqrt(1.29/rho(k)) ! Factor in precipitation velocity formula
	    irhoadz(k) = 1./(rho(k)*adz(k)) ! Useful factor
	    kb = max(1,k-1)
	    wmax       = dz*adz(kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
	    iwmax(k)   = 1./wmax
	 end do

! 	Add sedimentation of precipitation field to the vert. vel.

	do j=1,ny
	   do i=1,nx

      ! Compute precipitation velocity and flux column-by-column
      
	      prec_cfl = 0.
	      do k=1,nzm
!!$         rhoz = rhow(k)*dtn/dz
	         wp(k)=0.
	         omp = omegap(i,j,k)
	         lfac(k) = fac_cond+(1.-omp)*fac_fus
	         if(qp(i,j,k).gt.qp_threshold) then
	            if(omp.eq.1.) then
	               wp(k)= rhofac(k)*vrain*(rho(k)*qp(i,j,k))**crain
	            elseif(omp.eq.0.) then
	               omg = omegag(i,j,k)
	               qgg=omg*qp(i,j,k)
	               qss=qp(i,j,k)-qgg
	               wp(k)= rhofac(k)*(omg*vgrau*(rho(k)*qgg)**cgrau 
     &                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
	            else
	               omg = omegag(i,j,k)
	               qrr=omp*qp(i,j,k)
	               qss=qp(i,j,k)-qrr
	               qgg=omg*qss
	               qss=qss-qgg
	               wp(k)=rhofac(k)*(omp*vrain*(rho(k)*qrr)**crain 
     &	                     +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau 
     &	                          +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
	            endif
	            prec_cfl = max(prec_cfl,wp(k)*iwmax(k)) ! Keep column maximum CFL
	            wp(k) = -wp(k)*rhow(k)*dtn/dz 
	         endif
	      end do

	      fz(nz)=0.
	      www(nz)=0.
	      lfac(nz)=0.

      ! If maximum CFL due to precipitation velocity is greater than 0.9,
      ! take more than one advection step to maintain stability.
	      if (prec_cfl.gt.0.9) then
	         nprec = CEILING(prec_cfl/0.9)
	         do k = 1,nzm
	            ! wp already includes factor of dt, so reduce it by a
	            ! factor equal to the number of precipitation steps.
	            wp(k) = wp(k)/float(nprec) 
	         end do
	      else
	         nprec = 1
	      end if

	      do iprec = 1,nprec

	         do k = 1,nzm
	            tmp_qp(k) = qp(i,j,k) ! Temporary array for qp in this column
	         end do

         !-----------------------------------------

	         if(nonos) then

	            do k=1,nzm
	               kc=min(nzm,k+1)
	               kb=max(1,k-1)
	               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))
	               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))	  
	            end do

	         end if  ! nonos

         !  loop over iterations

	         do k=1,nzm
            ! Define upwind precipitation flux
	            fz(k)=tmp_qp(k)*wp(k)
	         end do

	         do k=1,nzm
	            kc=k+1
	            tmp_qp(k)=tmp_qp(k)-(fz(kc)-fz(k))*irhoadz(k) !Update temporary qp
	         end do

	         do k=1,nzm
            ! Also, compute anti-diffusive correction to previous
            ! (upwind) approximation to the flux
	            kb=max(1,k-1)
!!$            irho=1./(rho(k)*adz(k))
!!$            www(k)=andiff(tmp_qp(kb),tmp_qp(k),wp(k),irhoadz(k))
            ! The precipitation velocity is a cell-centered quantity,
            ! since it is computed from the cell-centered
            ! precipitation mass fraction.  Therefore, a reformulated
            ! anti-diffusive flux is used here which accounts for
            ! this and results in reduced numerical diffusion.
	            www(k) = 0.5*(1.+wp(k)*irhoadz(k)) 
     &	                 *(tmp_qp(kb)*wp(kb) - tmp_qp(k)*wp(k)) ! works for wp(k)<0
	         end do

         !---------- non-osscilatory option ---------------

	         if(nonos) then

	            do k=1,nzm
	               kc=min(nzm,k+1)
	               kb=max(1,k-1)
	               mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mx(k))
	               mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mn(k))	  
	            end do

	            do k=1,nzm
	               kc=min(nzm,k+1)
	               mx(k)=rho(k)*adz(k)*(mx(k)-tmp_qp(k)) 
     &	                      /(pn(www(kc)) + pp(www(k))+eps)
	               mn(k)=rho(k)*adz(k)*(tmp_qp(k)-mn(k)) 
     &	                      /(pp(www(kc)) + pn(www(k))+eps)
	            end do

	            do k=1,nzm
	               kb=max(1,k-1)
               ! Add limited flux correction to fz(k).
	               fz(k) = fz(k)                         ! Upwind flux
     &	                    + pp(www(k))*min(1.,mx(k), mn(kb)) 
     &	                    - pn(www(k))*min(1.,mx(kb),mn(k)) ! Anti-diffusive flux
	            end do

	         endif ! nonos

         ! Update precipitation mass fraction and liquid-ice static
         ! energy using precipitation fluxes computed in this column.
	         do k=1,nzm
	            kc=k+1
            ! Update precipitation mass fraction.
            ! Note that fz is the total flux, including both the
            ! upwind flux and the anti-diffusive correction.
	            qp(i,j,k)=qp(i,j,k)-(fz(kc)-fz(k))*irhoadz(k)
	            qpfall(k)=qpfall(k)-(fz(kc)-fz(k))*irhoadz(k)
	            lat_heat = -(lfac(kc)*fz(kc)-lfac(k)*fz(k))*irhoadz(k)
	            t(i,j,k)=t(i,j,k)-lat_heat
	            precflux(k) = precflux(k) - fz(k)   ! For statistics
	         end do
	         precsfc(i,j) = precsfc(i,j) - fz(1) ! For statistics
	         precssfc(i,j) = precssfc(i,j) - fz(1)*(1.-omegap(i,j,1)) ! For statistics

	         if (iprec.lt.nprec) then

            ! Re-compute precipitation velocity using new value of qp.
	            do k=1,nzm
	               wp(k)=0.
	               omp = omegap(i,j,k)
	               lfac(k) = fac_cond+(1.-omp)*fac_fus
	               if(qp(i,j,k).gt.qp_threshold) then
	                  if(omp.eq.1.) then
	                     wp(k)= rhofac(k)*vrain*(rho(k)*qp(i,j,k))**crain
	                  elseif(omp.eq.0.) then
	                     omg = omegag(i,j,k)
	                     qgg=omg*qp(i,j,k)
	                     qss=qp(i,j,k)-qgg
	                     wp(k)= rhofac(k)*(omg*vgrau*(rho(k)*qgg)**cgrau 
     &	                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
	                  else
	                     omg = omegag(i,j,k)
	                     qrr=omp*qp(i,j,k)
	                     qss=qp(i,j,k)-qrr
	                     qgg=omg*qss
	                     qss=qss-qgg
	                     wp(k)=rhofac(k)*(omp*vrain*(rho(k)*qrr)**crain 
     &	                           +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau 
     &	                                +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
	                  endif
                  ! Decrease precipitation velocity by factor of nprec
	                  wp(k) = -wp(k)*rhow(k)*dtn/dz/float(nprec)
                  ! Note: Don't bother checking CFL condition at each
                  ! substep since it's unlikely that the CFL will
                  ! increase very much between substeps when using
                  ! monotonic advection schemes.
	               endif
	            end do

	            fz(nz)=0.
	            www(nz)=0.
	            lfac(nz)=0.

	         end if

	      end do !iprec

	  end do
	end do	
	 


	end 


