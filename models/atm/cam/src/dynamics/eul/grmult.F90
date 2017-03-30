#include <misc.h>
#include <params.h>

subroutine grmult(rcoslat ,d       ,qm1     ,tm1     ,um1     ,&
                  vm1     ,z       ,tm2     ,phis    ,dpsl    ,&
                  dpsm    ,omga    ,pdel    ,pbot    ,logpsm2 ,&
                  logpsm1 ,rpmid   ,rpdel   ,fu      ,fv      ,&
                  t2      ,ut      ,vt      ,drhs    ,pmid    ,&
                  etadot  ,etamid  ,engy    ,ddpn    ,vpdsn   ,&
                  dpslon  ,dpslat  ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Non-linear dynamics calculations in grid point space
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id: grmult.F90,v 1.5.8.2 2002/06/15 13:47:44 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use commap
   use physconst, only: rair, cappa, cpvir, zvir

   implicit none

#include <comhyb.h>
!
! Input arguments
!
   real(r8), intent(in) :: rcoslat               ! 1./cosine(latitude)
   real(r8), intent(in) :: d(plond,plev)         ! divergence
   real(r8), intent(in) :: qm1(plond,plev)         ! specific humidity
   real(r8), intent(in) :: tm1(plond,plev)         ! temperature
   real(r8), intent(in) :: um1(plond,plev)         ! zonal wind * cos(lat)
   real(r8), intent(in) :: vm1(plond,plev)         ! meridional wind * cos(lat)
   real(r8), intent(in) :: z(plond,plev)         ! vorticity
   real(r8), intent(in) :: phis(plond)           ! surface geopotential
   real(r8), intent(in) :: dpsl(plond)           ! longitudinal component of grad ln(ps)
   real(r8), intent(in) :: dpsm(plond)           ! latitudinal component of grad ln(ps)
   real(r8), intent(in) :: omga(plond,plev)      ! vertical pressure velocity
   real(r8), intent(in) :: pdel(plond,plev)      ! layer thicknesses (pressure)
   real(r8), intent(in) :: pbot(plond)           ! bottom interface pressure
   real(r8), intent(in) :: logpsm2(plond)        ! log(psm2)
   real(r8), intent(in) :: logpsm1(plond)          ! log(ps)
   real(r8), intent(in) :: rpmid(plond,plev)     ! 1./pmid
   real(r8), intent(in) :: rpdel(plond,plev)     ! 1./pdel
   real(r8), intent(in) :: tm2(plond,plev)       ! temperature at previous time step
   integer, intent(in) :: nlon
!
! Input/Output arguments
!
   real(r8), intent(inout) :: fu(plond,plev)        ! nonlinear term - u momentum eqn
   real(r8), intent(inout) :: fv(plond,plev)        ! nonlinear term - v momentum eqn
   real(r8), intent(inout) :: t2(plond,plev)        ! nonlinear term - temperature
   real(r8), intent(inout) :: ut(plond,plev)        ! (u*TM1) - heat flux - zonal 
   real(r8), intent(inout) :: vt(plond,plev)        ! (u*TM1) - heat flux - meridional
   real(r8), intent(inout) :: drhs(plond,plev)      ! RHS of divergence eqn (del^2 term)
   real(r8), intent(inout) :: pmid(plond,plev)      ! pressure at full levels
   real(r8), intent(inout) :: etadot(plon,plevp)    ! vertical velocity in eta coordinates
   real(r8), intent(in)    :: etamid(plev)          ! midpoint values of eta (a+b)
   real(r8), intent(inout) :: engy(plond,plev)      ! kinetic energy
!
! Output arguments
!
   real(r8), intent(out) :: ddpn(plond)           ! complete sum of d*delta p
   real(r8), intent(out) :: vpdsn(plond)          ! complete sum V dot grad(ln(ps)) delta b
   real(r8), intent(out) :: dpslat(plond,plev)    ! ln(ps) component of lon press gradient
   real(r8), intent(out) :: dpslon(plond,plev)    ! ln(ps) component of lat press gradient

!
!---------------------------Local workspace-----------------------------
!
   real(r8) tv(plond,plev)        ! virtual temperature
   real(r8) ddpk(plond)           ! partial sum of d*delta p
   real(r8) vkdp                  ! V dot grad(ln(ps))
   real(r8) vpdsk(plond)          ! partial sum  V dot grad(ln(ps)) delta b
   real(r8) tk0(plond)            ! tm1 at phony level 0
   real(r8) uk0(plond)            ! u at phony level 0
   real(r8) vk0(plond)            ! v at phone level 0
   real(r8) rtv(plond,plev)       ! rair*tv
   real(r8) pterm(plond,plev)     ! intermediate term for hydrostatic eqn
   real(r8) tterm(plond,plev)     ! intermediate term for hydrostatic eqn
   real(r8) tmp                   ! temporary workspace
   real(r8) tmpk                  ! temporary workspace
   real(r8) tmpkp1                ! temporary workspace
   real(r8) edotdpde(plond,plevp) ! etadot*dp/deta
   real(r8) udel(plond,0:plev-1)  ! vertical u difference
   real(r8) vdel(plond,0:plev-1)  ! vertical v difference
   real(r8) tdel(plond,0:plev-1)  ! vertical TM1 difference

   integer i,k,kk             ! longitude, level indices
!
! Initialize arrays which represent vertical sums (ddpk, ddpn, vpdsk,
! vpdsn).  Set upper boundary condition arrays (k=0: tk0, uk0, vk0).
!
   do i=1,nlon
      ddpk(i)  = 0.0
      ddpn(i)  = 0.0
      vpdsk(i) = 0.0
      vpdsn(i) = 0.0
      tk0(i)  = 0.0
      uk0(i)  = 0.0
      vk0(i)  = 0.0
   end do
!
! Virtual temperature
!
   call virtem(nlon,plond,plev,tm1,qm1,zvir,tv)
!
! sum(plev)(d(k)*dp(k))
!
   do k=1,plev
      do i=1,nlon
         ddpn(i) = ddpn(i) + d(i,k)*pdel(i,k)
         rtv(i,k) = rair*tv(i,k)
      end do
   end do
!
! sum(plev)(v(k)*grad(lnps)*db(k))
!
   do k=nprlev,plev
      do i=1,nlon
         vkdp = rcoslat*(um1(i,k)*dpsl(i) + vm1(i,k)*dpsm(i))*pbot(i)
         vpdsn(i) = vpdsn(i) + vkdp*hybd(k)
      end do
   end do
!
! Compute etadot (dp/deta) (k+1/2).  Note: sum(k)(d(j)*dp(j)) required in
! pressure region. sum(k)(d(j)*dp(j)) and sum(k)(v(j)*grad(ps)*db(j)) 
! required in hybrid region
!
   do i=1,nlon
      edotdpde(i,1) = 0.
      edotdpde(i,plevp) = 0.
   end do
   do k=1,nprlev-1
      do i=1,nlon
         ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
         edotdpde(i,k+1) = -ddpk(i)
      end do
   end do
   do k=nprlev,plev-1
      do i=1,nlon
         ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
         vkdp = rcoslat*(um1(i,k)*dpsl(i) + vm1(i,k)*dpsm(i))*pbot(i)
         vpdsk(i) = vpdsk(i) + vkdp*hybd(k)
         edotdpde(i,k+1) = -ddpk(i) - vpdsk(i) + hybi(k+1)*(ddpn(i)+vpdsn(i))
      end do
   end do
!
! Nonlinear advection terms.  u*tm1, v*tm1, kinetic energy first
!
   do k=1,plev
      do i=1,nlon
         ut(i,k) = um1(i,k)*tm1(i,k)
         vt(i,k) = vm1(i,k)*tm1(i,k)
         engy(i,k) = 0.5*(um1(i,k)**2 + vm1(i,k)**2)
      end do
   end do
!
! Compute workspace arrays for delta-u, delta-v, delta-tm1 (k)
!
   do i=1,nlon
      udel(i,0) = um1(i,1) - uk0(i)
      vdel(i,0) = vm1(i,1) - vk0(i)
      tdel(i,0) = tm1(i,1) - tk0(i)
   end do
   do k=1,plev-1
      do i=1,nlon
         udel(i,k) = um1(i,k+1) - um1(i,k)
         vdel(i,k) = vm1(i,k+1) - vm1(i,k)
         tdel(i,k) = tm1(i,k+1) - tm1(i,k)
      end do
   end do
!
! Horizontal advection: u*z, v*z, energy conversion term (omega/p),
! vertical advection for interface above.  Pure pressure region first.
!
   do k=1,nprlev-1
      do i=1,nlon
         dpslat(i,k) = 0.
         dpslon(i,k) = 0.
         tmpk   = 0.5*rpdel(i,k)*edotdpde(i,k  )
         tmpkp1 = 0.5*rpdel(i,k)*edotdpde(i,k+1)
         fu(i,k) = fu(i,k) + vm1(i,k)*z(i,k) - udel(i,k-1)*tmpk - udel(i,k  )*tmpkp1
         fv(i,k) = fv(i,k) - um1(i,k)*z(i,k) - vdel(i,k-1)*tmpk - vdel(i,k  )*tmpkp1
#ifdef HADVTEST
!
!jr Modify so TM1 only has horizontal advection
!
         t2(i,k) = t2(i,k) + d(i,k)*tm1(i,k)
#else
         t2(i,k) = t2(i,k) + d(i,k)*tm1(i,k) - tdel(i,k-1)*tmpk + &
            cappa*tv(i,k)/(1. + cpvir*qm1(i,k))* &
            omga(i,k)*rpmid(i,k) - tdel(i,k)*tmpkp1
#endif
      end do
   end do
!
! Hybrid region above bottom level: Computations are the same as in pure
! pressure region, except that pressure gradient terms are added to 
! momentum tendencies.
!
   do k=nprlev,plev-1
      do i=1,nlon
         tmpk   = 0.5*rpdel(i,k)*edotdpde(i,k  )
         tmpkp1 = 0.5*rpdel(i,k)*edotdpde(i,k+1)
         tmp = rtv(i,k)*hybm(k)*rpmid(i,k)*pbot(i)
         dpslon(i,k) = rcoslat*tmp*dpsl(i)
         dpslat(i,k) = rcoslat*tmp*dpsm(i)
         fu(i,k) = fu(i,k) + vm1(i,k)*z(i,k) - udel(i,k-1)*tmpk - &
            udel(i,k  )*tmpkp1 - dpslon(i,k)
         fv(i,k) = fv(i,k) - um1(i,k)*z(i,k) - vdel(i,k-1)*tmpk - &
            vdel(i,k  )*tmpkp1 - dpslat(i,k)
#ifdef HADVTEST
!
!jr Modify so TM1 only has horizontal advection
!
         t2(i,k) = t2(i,k) + d(i,k)*tm1(i,k)
#else
         t2(i,k) = t2(i,k) + d(i,k)*tm1(i,k) - tdel(i,k-1)*tmpk + &
            cappa*tv(i,k)/(1. + cpvir*qm1(i,k))* &
            omga(i,k)*rpmid(i,k) - tdel(i,k)*tmpkp1
#endif
      end do
   end do
!
! Bottom level
!
   do i=1,nlon
      tmpk = 0.5*rpdel(i,plev)*edotdpde(i,plev  )
      tmp  = rtv(i,plev)*hybm(plev)*rpmid(i,plev)*pbot(i)
      dpslon(i,plev) = rcoslat*tmp*dpsl(i)
      dpslat(i,plev) = rcoslat*tmp*dpsm(i)
      fu(i,plev) = fu(i,plev) + vm1(i,plev)*z(i,plev) - &
         udel(i,plev-1)*tmpk - dpslon(i,plev)
      fv(i,plev) = fv(i,plev) - um1(i,plev)*z(i,plev) - &
         vdel(i,plev-1)*tmpk - dpslat(i,plev)
#ifdef HADVTEST
!
!jr Modify so TM1 only has horizontal advection
!
      t2(i,plev) = t2(i,plev) + d(i,plev)*tm1(i,plev)
#else
      t2(i,plev) = t2(i,plev) + d(i,plev)*tm1(i,plev) - &
         tdel(i,plev-1)*tmpk + &
         cappa*tv(i,plev)/(1. + cpvir*qm1(i,plev))* &
         omga(i,plev)*rpmid(i,plev)
#endif
   end do
!
! Convert eta-dot(dp/deta) to eta-dot (top and bottom = 0.)
!
   do i=1,nlon
      etadot(i,1) = 0.
      etadot(i,plevp) = 0.
   end do
   do k=2,plev
      tmp = etamid(k) - etamid(k-1)
      do i=1,nlon
#ifdef HADVTEST
!
!jr Set etadot to zero for horizontal advection test
!
         etadot(i,k) = 0.
#else
         etadot(i,k) = edotdpde(i,k)*tmp/(pmid(i,k) - pmid(i,k-1))
#endif
      end do
   end do
!
!-----------------------------------------------------------------------
!
! Divergence and hydrostatic equations
! Store some temporary terms
!
   do k=2,plev-1
      do i=1,nlon
         pterm(i,k) = rtv(i,k)*rpmid(i,k)*pdel(i,k)
         tterm(i,k) = 0.5*tm2(i,k) - tm1(i,k)
      end do
   end do
   do i=1,nlon
      tterm(i,1) = 0.5*tm2(i,1) - tm1(i,1)
      tterm(i,plev) = 0.5*tm2(i,plev) - tm1(i,plev)
   end do
!
! Del squared part of RHS of divergence equation.
! Kinetic energy and diagonal term of hydrostatic equation.
! Total temperature as opposed to  perturbation temperature is acceptable
! since del-square operator will operate on this term.
!
   do k=1,plev-1
      do i=1,nlon
         drhs(i,k) = phis(i) + engy(i,k) + rtv(i,k)*0.5* &
            rpmid(i,k)*pdel(i,k) + href(k,k)*tterm(i,k) + &
            bps(k)*(0.5*logpsm2(i) - logpsm1(i))
      end do
   end do
   do i=1,nlon
      drhs(i,plev) = phis(i) + engy(i,plev) + rtv(i,plev)*0.5* &
         rpmid(i,plev)*pdel(i,plev) + &
         href(plev,plev)*tterm(i,plev) + &
         bps(plev)*(0.5*logpsm2(i) - logpsm1(i))
   end do
!
! Bottom level term of hydrostatic equation
!
   do k=1,plev-1
      do i=1,nlon
         drhs(i,k) = drhs(i,k) + rtv(i,plev)* &
            rpmid(i,plev)*pdel(i,plev) + &
            href(plev,k)*tterm(i,plev)
      end do
   end do
!
! Interior terms of hydrostatic equation
!
   do k=1,plev-2
      do kk=k+1,plev-1
         do i=1,nlon
            drhs(i,k) = drhs(i,k) + pterm(i,kk) + href(kk,k)*tterm(i,kk)
         end do
      end do
   end do
!    
   return
end subroutine grmult
