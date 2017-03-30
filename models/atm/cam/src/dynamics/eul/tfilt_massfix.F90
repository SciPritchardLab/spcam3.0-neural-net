#include <misc.h>
#include <params.h>

subroutine tfilt_massfix (ztodt,         lat,    u3m1,   u3,     &
                          v3m1,   v3,    t3m1,   t3,     q3m1,   &
                          q3,     psm1,  ps,             alpha,  &
                          etamid, qfcst, vort,   div,    vortm2, &
                          divm2,         qminus, psm2,   um2,    &
                          vm2,    tm2,   qm2,    vortm1, divm1,  &
                          omga,   dpsl,  dpsm,   beta,   nlon)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Time filter (second half of filter for vorticity and divergence only)
! 
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: outfld
   use pmgrid
   use pspect
   use commap
   use constituents, only: pcnst, pnats, qmin
   use constituents, only: tottnam, tendnam
   use time_manager, only: get_nstep
   use physconst, only: cpair, gravit

   implicit none

#include <comctl.h>
#include <comlun.h>
#include <comqfl.h>
#include <comtfc.h>
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                  ! two delta t (unless nstep<2)

   real(r8), intent(in) :: qfcst(plond,plev,pcnst)! slt moisture forecast
   real(r8), intent(in) :: vort(plond,plev)
   real(r8), intent(in) :: div(plond,plev)
   real(r8), intent(inout) :: vortm2(plond,plev)
   real(r8), intent(inout) :: divm2(plond,plev)
   real(r8), intent(in) :: qminus(plond,plev,pcnst)
   real(r8), intent(inout) :: psm2(plond)
   real(r8), intent(inout) :: um2(plond,plev)
   real(r8), intent(inout) :: vm2(plond,plev)
   real(r8), intent(inout) :: tm2(plond,plev)
   real(r8), intent(inout) :: qm2(plond,plev,pcnst+pnats)
   real(r8), intent(inout) :: omga(plond,plev)
   real(r8), intent(in) :: dpsl(plond)
   real(r8), intent(in) :: dpsm(plond)
   real(r8), intent(in) :: beta                   ! energy fixer coefficient
   real(r8), intent(in) :: alpha(pcnst)
   real(r8), intent(in) :: etamid(plev)           ! vertical coords at midpoints 
   real(r8), intent(in) :: u3(plond,plev)
   real(r8), intent(in) :: v3(plond,plev)
   real(r8), intent(inout) :: t3(plond,plev)

   integer, intent(in) :: lat
   integer, intent(in) :: nlon

! Input/Output arguments

   real(r8), intent(inout) :: q3(plond,plev,pcnst+pnats)
   real(r8), intent(inout) :: ps(plond)
   real(r8), intent(inout) :: vortm1(plond,plev)
   real(r8), intent(inout) :: psm1(plond)
   real(r8), intent(inout) :: u3m1(plond,plev)
   real(r8), intent(inout) :: v3m1(plond,plev)
   real(r8), intent(inout) :: t3m1(plond,plev)
   real(r8), intent(inout) :: divm1(plond,plev)
   real(r8), intent(inout) :: q3m1(plond,plev,pcnst+pnats)
!
! Local workspace
!
   integer ifcnt                   ! Counter
   integer :: nstep                ! current timestep number

   real(r8) qfcst1(plond,plev,pcnst) ! Workspace for qfcst
   real(r8) engycorr(plond,plev)   ! energy equivalent to T correction
   real(r8) rpmid(plond,plev)      ! 1./pmid
   real(r8) pdel(plond,plev)       ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) pint(plond,plevp)      ! pressure at model interfaces (n  )
   real(r8) pmid(plond,plev)       ! pressure at model levels (time n)
   real(r8) utend(plond,plev)      ! du/dt
   real(r8) vtend(plond,plev)      ! dv/dt
   real(r8) ttend(plond,plev)      ! dT/dt
   real(r8) qtend(plond,plev,pcnst)! dq/dt
   real(r8) pstend(plond)          ! d(ps)/dt
   real(r8) pintm1(plond,plevp)    ! pressure at model interfaces (n-1)
   real(r8) pmidm1(plond,plev)     ! pressure at model levels (time n-1)
   real(r8) pdelm1(plond,plev)     ! pdelm1(k) = pintm1(k+1)-pintm1(k)
   real(r8) om2eps
   real(r8) corm
   real(r8) wm
   real(r8) absf
   real(r8) worst
   logical lfixlim               ! flag to turn on fixer limiter

   real(r8) ta(plond,plev,pcnst)   ! total advection of constituents
   real(r8) dqfx3(plond,plev,pcnst)! q tendency due to mass adjustment
   real(r8) coslat                 ! cosine(latitude)
   real(r8) rcoslat(plond)         ! 1./cosine(latitude)
!   real(r8) engt                   ! Thermal   energy integral
!   real(r8) engk                   ! Kinetic   energy integral
!   real(r8) engp                   ! Potential energy integral

   integer i, k, m
!-----------------------------------------------------------------------
   qfcst1(1:nlon,:,:) = qfcst(i1:nlon+i1-1,:,:)

   nstep = get_nstep()

   coslat = cos(clat(lat))
   do i=1,nlon
     rcoslat(i) = 1./coslat
   enddo
   lfixlim = .true.
   corm    = 0.1
!
! Set average dry mass to specified constant preserving horizontal
! gradients of ln(ps). Proportionality factor was calculated in STEPON
! for nstep=0 or SCAN2 otherwise from integrals calculated in INIDAT
! and SCAN2 respectively.
! Set p*.
!
   do i=1,nlon
      ps(i) = ps(i)*fixmas
   end do
!
! Set current time pressure arrays for model levels etc.
!
   call plevs0(nlon    ,plond   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
   rpmid(:nlon,:plev) = 1./pmid(:nlon,:plev)
!
! Add temperature correction for energy conservation
!
   do k=1,plev
      do i=1,nlon
         engycorr(i,k) = (cpair/gravit)*beta*pdel(i,k)/ztodt
         t3      (i,k) = t3(i,k) + beta
      end do
   end do
!
! Output Energy correction term
!
   call outfld('ENGYCORR',engycorr ,plond   ,lat     )
!
! Compute q tendency due to mass adjustment
! If LFIXLIM = .T., then:
! Check to see if fixer is exceeding a desired fractional limit of the
! constituent mixing ratio ("corm").  If so, then limit the fixer to
! that specified limit.
!
   do m=1,pcnst
      do k=1,plev
         do i=1,nlon
            dqfx3(i,k,m) = alpha(m)*etamid(k)*abs(qfcst1(i,k,m) - qminus(i,k,m))
         end do
         if (lfixlim) then
            ifcnt = 0
            worst = 0.
            wm    = 0.
            do i = 1,nlon
               absf = abs(dqfx3(i,k,m))
               if (absf.gt.corm) then
                  ifcnt = ifcnt + 1
                  worst = max(absf,worst)
                  wm = wm + absf
                  dqfx3(i,k,m) = sign(corm,dqfx3(i,k,m))
               endif
            end do
            if (ifcnt.gt.0) then
               wm = wm/float(ifcnt)
#ifndef HADVTEST
! TBH:  Commented out as of CAM CRB meeting on 6/20/03
!               write (6,1000) m,corm,ifcnt,k,lat,wm,worst
#endif
            endif
         endif
         do i=1,nlon
            dqfx3(i,k,m) = qfcst1(i,k,m)*dqfx3(i,k,m)/ztodt
#ifdef HADVTEST
            q3(i,k,m) = qfcst1(i,k,m)
#else
            q3(i,k,m) = qfcst1(i,k,m) + ztodt*dqfx3(i,k,m)
#endif
            ta(i,k,m) = (q3(i,k,m) - qminus(i,k,m))/ztodt
         end do
      end do
   end do
!
! Check for and correct invalid constituents
!
   call qneg3 ('TFILT_MASSFIX',lat   ,nlon    ,plond   ,plev    , &
               pcnst+pnats,qmin ,q3(1,1,1))
!
! Send slt tendencies to the history tape
!
   do m=1,pcnst
      call outfld(tottnam(m),ta(1,1,m),plond   ,lat     )
   end do
!
! Calculate vertical motion field
!
   call omcalc (rcoslat ,div     ,u3      ,v3      ,dpsl    ,  &
                dpsm    ,pmid    ,pdel    ,rpmid   ,pint(1,plevp), &
                omga    ,nlon    )

!   write(6,*)'tfilt: lat=',lat
!   write(6,*)'omga=',omga
!
! Time filter (second half of filter for vorticity and divergence only)
!
!   if(lat.eq.2) then
!      write(6,*)'tfilt: ps=',psm2(13),psm1(13),ps(13)
!      write(6,*)'tfilt: u=',um2(13,18),u3m1(13,18),u3(13,18)
!      write(6,*)'tfilt: t=',tm2(13,18),t3m1(13,18),t3(13,18)
!      write(6,*)'tfilt: water=',qm2(13,18,1),q3m1(13,18,1),q3(13,18,1)
!      write(6,*)'tfilt: cwat=',qm2(13,18,2),q3m1(13,18,2),q3(13,18,2)
!      write(6,*)'tfilt: vort=',vortm2(13,18),vortm1(13,18),vort(13,18)
!      write(6,*)'tfilt: div=',divm2(13,18),divm1(13,18),div(13,18)
!   end if

   om2eps = 1. - 2.*eps
   if (nstep.ge.2) then
      do k=1,plev
         do i=1,nlon
            u3m1(i,k) = om2eps*u3m1(i,k) + eps*um2(i,k) + eps*u3(i,k)
            v3m1(i,k) = om2eps*v3m1(i,k) + eps*vm2(i,k) + eps*v3(i,k)
            t3m1(i,k) = om2eps*t3m1(i,k) + eps*tm2(i,k) + eps*t3(i,k)
            q3m1(i,k,1) = om2eps*q3m1(i,k,1) + eps*qm2(i,k,1) + eps*q3(i,k,1)
            vortm1(i,k) = om2eps*vortm1(i,k) + eps*vortm2(i,k) + eps*vort(i,k)
            divm1(i,k) = om2eps*divm1(i,k) + eps*divm2(i,k) + eps*div(i,k)
         end do
         do m=2,pcnst+pnats
            do i=1,nlon
               q3m1(i,k,m) = om2eps*q3m1(i,k,m) + eps*qm2(i,k,m) + eps*q3(i,k,m)
            end do
         end do
      end do
      do i=1,nlon
         psm1(i) = om2eps*psm1(i) + eps*psm2(i) + eps*ps(i)
      end do
   end if

   call plevs0 (nlon    ,plond   ,plev    ,psm1    ,pintm1  ,pmidm1  ,pdelm1)
!
! Compute time tendencies:comment out since currently not on h-t
!
   do k=1,plev
      do i=1,nlon
         ttend(i,k) = (t3(i,k)-tm2(i,k))/ztodt
         utend(i,k) = (u3(i,k)-um2(i,k))/ztodt
         vtend(i,k) = (v3(i,k)-vm2(i,k))/ztodt
      end do
   end do

   do m=1,pcnst
      do k=1,plev
         do i=1,nlon
            qtend(i,k,m) = (q3(i,k,m) - qm2(i,k,m))/ztodt
         end do
      end do
   end do

   do i=1,nlon
      pstend(i) = (ps(i) - psm2(i))/ztodt
   end do
!
   do m=1,pcnst
      call outfld (tendnam(m),qtend(1,1,m),plond,lat)
   end do

   call outfld ('UTEND   ',utend,plond,lat)
   call outfld ('VTEND   ',vtend,plond,lat)
   call outfld ('TTEND   ',ttend,plond,lat)
   call outfld ('LPSTEN  ',pstend,plond,lat)

   return
1000 format(' TIMEFILTER: WARNING: fixer for tracer ',i3,' exceeded ', &
      f8.5,' for ',i5,' points at k,lat = ',2i4, &
      ' Avg/Worst = ',1p2e10.2)
end subroutine  tfilt_massfix
