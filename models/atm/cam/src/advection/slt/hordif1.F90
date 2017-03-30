#include <misc.h>
#include <params.h>

subroutine hordif1(rearth,phi)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Horizontal diffusion of z,d,t,q
! 
! Method: 
! 1. implicit del**2 form above level kmnhd4
! 2. implicit del**4 form at level kmnhd4 and below
! 3. courant number based truncation at level kmxhdc and above
! 4. increased del**2 coefficient at level kmxhd2 and above
!
! Computational note: this routine is multitasked by level, hence it 
! is called once for each k
!
! Note: Most of the "ifdef" constructs employed in this routine are related
! to the fact that storage order for spectral coefficients is different
! depending upon whether the target architecture is PVP or not.
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! $Id: hordif1.F90,v 1.1.2.2 2003/10/23 21:52:20 jmccaa Exp $
! $Author: jmccaa $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in)    :: rearth     ! radius of earth
  real(r8), intent(inout) :: phi(psp)   ! used in spectral truncation of phis
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
#if ( defined PVP )
  real(r8) tqfac(2*pnmax)  ! time-split implicit diffusion factors (t,q)
  integer isp              ! spectral pointer
  integer ne               ! spectral pointer
#else
  integer ir,ii            ! spectral indices       
  integer mr,mc            ! spectral indices
#endif
  real(r8) k42             ! Nominal  Del^4 diffusion coeff at T42
  real(r8) k63             ! Nominal  Del^4 diffusion coeff at T63
  real(r8) knn             ! Computed Del^4 diffusion coeff at TNN
  real(r8) tmp             ! temp space
  real(r8) hdfst4(pnmax)
  integer  expon
  integer  m,n             ! spectral indices
!-----------------------------------------------------------------------
!
! Compute Del^4 diffusion coefficient
!
  k42   = 1.e+16
  k63   = 5.e+15
  expon = 25

  if(pmax-1 <= 42) then
     knn = k42
  elseif(pmax-1 == 63) then
     knn = k63
  else
     if(pmax-1 < 63) then
        tmp = log(k42/k63)/log(63._r8*64._r8/42._r8/43._r8)
     else
        tmp = 2.
     endif
     knn = k63*(63._r8*64._r8/float(pmax)/float(pmax-1))**tmp
  endif
!
! Set the Del^4 diffusion coefficients for each wavenumber
!
  hdfst4(1) = 0.
  do n=2,pnmax
     hdfst4(n) = knn * (n*(n-1)*n*(n-1)  ) / rearth**4
  end do
!
! Set the horizontal diffusion factors for each wavenumer at this level
! del^4 diffusion is to be applied and compute time-split implicit
! factors.
!
#if ( defined PVP )

  do n = 1, pnmax
     tqfac(2*n)   = 1./(1. + 3600.*hdfst4(n))
     tqfac(2*n-1) = 1./(1. + 3600.*hdfst4(n))
  end do
!
! Apply the horizontal diffusion
!
  do n=1,pmax
     isp = nco2(n) - 2
     ne = 2*(n-1)
     do m=1,2*nm(n)
        phi(isp+m)  = phi(isp+m)*tqfac(ne+m)**expon
     end do
  end do

#else
  do m=1,pmmax
     mr = nstart(m)
     mc = 2*mr
     do n=1,nlen(m)
        ir = mc + 2*n - 1
        ii = ir + 1
        phi(ir)  = phi(ir)/(1. + 3600.*hdfst4(n+m-1))**expon
        phi(ii)  = phi(ii)/(1. + 3600.*hdfst4(n+m-1))**expon
     end do
  end do

#endif

  return
end subroutine hordif1
