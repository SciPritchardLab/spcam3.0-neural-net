! Parameterization of momentum transport by deep convection.
! Marat Khairoutdinov, May 2005


subroutine convmom(zm_ls, zi_ls, dp_ls, u_ls, cape, mc, ztop, nz, utend)
  implicit none

  real, intent(in) :: zm_ls(nz) ! grid levels (m)
  real, intent(in) :: zi_ls(nz) ! grid levels (m)
  real, intent(in) :: dp_ls(nz) ! pressure thickness (Pa)
  real, intent(in) :: u_ls(nz)   ! velocity profile
  real, intent(in) :: mc(nz)   ! convective mass-flux profile (kg/m2/s)
  real, intent(in) :: cape    ! CAPE (J/kg)
  real, intent(in) :: ztop    ! convectopn top height (m)
  integer, intent(in) :: nz ! number of levels
  real, intent(out) :: utend(nz) ! tendency of u due to convection
  real, external:: w2sim


! local space

 real wstar ! convective velocity scale
 real zstar ! convective height scale
 real taud  ! momentum relaxation timescale
 real tauc  ! convective time scale
 real z(nz) ! grid levels (m)
 real rho(nz) ! air density (kg/kg)
 real u(nz)   ! velocity profile
 real ss(nz)    ! source
 real zf(nz)
 real w2(nz)
 real mom(nz)   ! momentum flux
 real dudt(nz)
 real alpha, beta

 real, parameter :: cw = 0.089  ! coefficient in CAPE->wstar
 real, parameter :: aw = 7.  ! coefficient in third-moment approximation
 real, parameter :: ad = 14.  ! coefficient in relaxation term
 integer k,m

 wstar = cw*sqrt(cape) ! from Khairoutdinov and Randall, JAS 2002
 zstar = ztop/1.2
 tauc = zstar/wstar
 taud = tauc/ad ! initial guess

 do k=1,nz
  m = nz-k+1 
  z(k) = zm_ls(m)
  zf(k) = zi_ls(m+1)
  u(k) = u_ls(m)
  rho(k) = dp_ls(m)/9.81/(zi_ls(m)-zi_ls(m+1))
  w2(k) = 200.*(mc(m)/rho(k))**2
 end do

 dudt = 0.
 mom = 0.
 ss = 0.
 do k=2,nz-1
  if(z(k).gt.ztop) EXIT
  ss(k) = 0.5*(rho(k)+rho(k-1))*w2(k)*(u(k)-u(k-1))/(z(k)-z(k-1))
!   ss(k) = 0.5*(rho(k)+rho(k-1))*w2sim(zf(k)/zstar)*wstar**2*(u(k)-u(k-1))/(z(k)-z(k-1))
  alpha = 1./taud
  beta = aw/tauc/(zf(k)-zf(k-1))
  mom(k) = (beta*zf(k-1)*mom(k-1) - ss(k))/(alpha+beta*zf(k))
  dudt(k-1) = -(mom(k)-mom(k-1))/(zf(k)-zf(k-1))/rho(k-1)
 end do
 do k=1,nz
  m = nz-k+1 
  utend(m) = max(-1./3600.,min(1./3600.,dudt(k)))
!  write(*,'(11g13.5)')z(k),zf(k),u(k),rho(k),mom(k)/rho(k),utend(m)*3600.,ss(k)
 end do

end



real function w2sim(z)
 real z
 real w2(13)/0.,0.01,0.02,0.085,0.115,0.130,0.150,0.155,0.135,0.105,0.06,0.025,0.01/
 integer ind
 ind = z*10.+1
 if(ind.ge.13) then
   w2sim=0.
   return
 end if
 w2sim = w2(ind)+(z-(ind-1)*0.1)/0.1*(w2(ind+1)-w2(ind))
 return
end


