#include <misc.h>
#include <params.h>

module comslt

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialization of Semi-Lagrangian transport variables.
! 
! Author: 
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plond, plev, plat, beglat, endlat
  use constituents, only: pcnst
  use infnan

  implicit none

  real(r8) hw1(pcnst)   ! Pre-SLT global integral of constituent
  real(r8) hw2(pcnst)   ! Post-SLT global integral of const.
  real(r8) hw3(pcnst)   ! Global integral for denom. of expr. for alpha
  real(r8) alpha(pcnst) ! alpha(m) = ( hw1(m) - hw2(m) )/hw3(m)
  real(r8) hw1lat(pcnst,plat) ! lat contribution to const. mass integral
  real(r8) engy1lat(plat)          ! lat contribution to total energy integral

  real(r8), allocatable :: lammp(:,:,:)
  real(r8), allocatable :: phimp(:,:,:)
  real(r8), allocatable :: sigmp(:,:,:)
  real(r8), allocatable :: qfcst(:,:,:,:)

contains

  subroutine initialize_comslt

    allocate (lammp(plon,plev,beglat:endlat))
    allocate (phimp(plon,plev,beglat:endlat))
    allocate (sigmp(plon,plev,beglat:endlat))
    allocate (qfcst(plond,plev,pcnst,beglat:endlat))

    lammp (:,:,:)   = inf
    phimp (:,:,:)   = inf
    sigmp (:,:,:)   = inf
    qfcst (:,:,:,:) = inf

    return
  end subroutine initialize_comslt

end module comslt

