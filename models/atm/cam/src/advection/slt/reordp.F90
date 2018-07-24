#include <misc.h>
#include <params.h>

subroutine reordp(irow    ,iy      ,zalp    ,zdalp   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Renormalize associated Legendre polynomials and their derivatives.
! 
! Method: 
! Reorder associated Legendre polynomials and their derivatives from
! column rectangular storage to diagonal pentagonal storage. The
! reordered polynomials and derivatives are returned via common/comspe/
! Note: Most of the "ifdef" constructs employed in this routine are related
! to the fact that storage order for spectral coefficients is different
! depending up whether the target architecture is PVP or not.  The token
! SPMD has an "ifdef" test associated with it since the message-passing 
! implementation of CCM3 distributes subregions of Fourier wavenumber space
! to individual processes.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: reordp.F90,v 1.1.2.1 2002/06/15 13:46:59 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: irow            ! latitude pair index
  integer , intent(in)  :: iy              ! dimension of input polynomials
  real(r8), intent(in)  :: zalp(iy)        ! Legendre polynomial
  real(r8), intent(in)  :: zdalp(iy)       ! Legendre polynomial derivative
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
#if ( defined PVP )
  integer im              ! length of current diagonal
  integer itrn            ! number of rows
  integer ik2             ! K+2
  integer is              ! off-set to current diagonal
  integer in              ! index of row
#else
  integer mr              ! spectral index
#ifdef QVORTDAMP
integer :: counter
integer :: mn_distance (iy),nvalue(iy)
#endif
#endif
  integer m               ! index along diagonal and row
  integer n               ! index of diagonal
  real(r8) sqrt2              ! sqrt(2)
!-----------------------------------------------------------------------
!
! Multiply ALP and DALP by SQRT(2.) in order to get proper
! normalization. DALP is multiplied by -1 to correct for - sign
! in Copenhagen definition.
!
  sqrt2 = sqrt(2.)
#if ( defined PVP )
  itrn = ptrn + 1
  ik2 = ptrk + 2
  is = 0
  do n=1,itrn
     im = min0(pmmax,ik2-n)
     in = n - itrn
!CDIR$ IVDEP
     do m=1,im
        alp(is+m,irow) = zalp(m*itrn+in)*sqrt2
        dalp(is+m,irow) = -zdalp(m*itrn+in)*sqrt2
     end do
     is = is + im
  end do
#else
#ifdef QVORTDAMP
! pritch, I think this is the right m,n vectorization corresponding to zalp above. 
! (based on browsing phcs.F90 )
counter = 0
do m=1,pmmax
do n=1,(ptrn+1)
    counter=counter+1
    mn_distance (counter) = abs( (m-1) - (n-1) )
    nvalue(counter) = n-1
  end do
end do
#endif
  do m=1,pmmax
     mr = nstart(m)
     do n=1,nlen(m)
        alp(mr+n,irow) = zalp((m-1)*pmax + n)*sqrt2
        dalp(mr+n,irow) = -zdalp((m-1)*pmax + n)*sqrt2
#ifdef QVORTDAMP
!        Pritch, save the absolute (m-n) difference in the maximal-vector
!         ordering for use in equatorial spectral filtering in grcalc.F90
        mn_distance_reordered (mr+n,irow) = mn_distance((m-1)*pmax + n)
        nvalue_reordered (mr+n,irow) = nvalue((m-1)*pmax+n)
! nb, mn_distance_reordered is a new global added to comspe.F90
#endif
     end do
  end do
#endif

  return
end subroutine reordp

