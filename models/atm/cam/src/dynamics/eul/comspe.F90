#include <misc.h>
#include <params.h>

module comspe

!----------------------------------------------------------------------- 
! 
! Purpose: Spectral space arrays
! 
! Method: 
! 
! Author: CCM Core Group
! $Author: pworley $
! $Id: comspe.F90,v 1.2.6.3 2003/08/13 19:16:11 pworley Exp $
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use infnan
  use pmgrid, only: plev, plat
  use pspect

  implicit none
!
! $Id: comspe.F90,v 1.2.6.3 2003/08/13 19:16:11 pworley Exp $
! $Author: pworley $
!
! Spectral space arrays
!
  real(r8), dimension(:,:), allocatable :: vz   ! Vorticity spectral coefficients
  real(r8), dimension(:,:), allocatable :: d    ! Divergence spectral coefficients
  real(r8), dimension(:,:), allocatable :: t    ! Temperature spectral coefficients
  real(r8), dimension(:), allocatable :: alps ! Log-pressure spectral coefficients

#if ( defined SPMD )
  integer :: numm(0:plat-1) = bigint  ! number of Fourier wavenumbers owned per task
  integer :: maxm           = bigint  ! max number of Fourier wavenumbers per MPI task
  integer :: lpspt          = bigint  ! number of local spectral coefficients
  integer, dimension(:,:), allocatable :: locm 
                                      ! assignment of wavenumbers to MPI tasks
  integer, dimension(:), allocatable   :: lnstart 
                                      ! Starting indices for local spectral arrays (real)
#else
  integer :: numm(0:0)      = pmmax
  integer :: maxm           = pmmax
  integer :: lpspt          = pspt
  integer :: locm(1:pmmax, 0:0) = bigint
  integer :: lnstart(1:pmmax) = bigint
#endif

#if ( defined PVP )
  integer :: ncoefi(pmaxp) = bigint   ! Pointer to start of coefficient diagonals
  integer :: nm(pmax)      = bigint   ! Number of coeffs stored on a given diagonal
  integer :: nco2(pmax)    = bigint   ! Complex form of ncoefi
  integer :: nalp(pmax)    = bigint   ! Pointer into polynomial arrays
  integer :: ncutoff       = bigint   ! Break-even point for vector lengths in GRCALC
#else
  integer :: nstart(pmmax) = bigint   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)   = bigint   ! Length vectors for spectral arrays
#endif

  real(r8), dimension(:,:), allocatable :: alp  ! Legendre polynomials (pspt,plat/2)
  real(r8), dimension(:,:), allocatable :: dalp ! Legendre polynomial derivatives (pspt,plat/2)
!
  real(r8), dimension(:,:), allocatable :: lalp  ! local Legendre polynomials
  real(r8), dimension(:,:), allocatable :: ldalp ! local Legendre polynomial derivatives

end module comspe
