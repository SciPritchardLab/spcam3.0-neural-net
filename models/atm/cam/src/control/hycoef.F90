#include <misc.h>
#include <params.h>

subroutine hycoef

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Defines the locations of model interfaces from input data in the
! hybrid coordinate scheme.  Actual pressure values of model level
! interfaces are determined elsewhere from the fields set here.
! 
! Method: 
! the following fields are set:
! hyai     fraction of reference pressure used for interface pressures
! hyam     fraction of reference pressure used for midpoint pressures
! hybi     fraction of surface pressure used for interface pressures
! hybm     fraction of surface pressure used for midpoint pressures
! hybd     difference of hybi's
! hypi     reference state interface pressures
! hypm     reference state midpoint pressures
! hypd     reference state layer thicknesses
! hypdln   reference state layer thicknesses (log p)
! hyalph   distance from interface to level (used in integrals)
! prsfac   log pressure extrapolation factor (used to compute psl)
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------
!
! $Id: hycoef.F90,v 1.2.4.1 2002/06/15 13:47:27 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none
#include <comhyb.h>

!---------------------------Local workspace-----------------------------
  integer k                 ! Level index
  real(r8) amean,bmean,atest,btest,eps
!-----------------------------------------------------------------------
!
  eps    = 1.e-05
  nprlev = 0
  ps0    = 1.0e5            ! Base state surface pressure (pascals)
  psr    = 1.0e5            ! Reference surface pressure (pascals)
!
! Set layer locations
!
  do k=1,plev
!
! Interfaces. Set nprlev to the interface above, the first time a 
! nonzero surface pressure contribution is found. "nprlev" 
! identifies the lowest pure pressure interface.
!
     if (nprlev==0 .and. hybi(k).ne.0.0) nprlev = k - 1
  end do
!
! Set nprlev if no nonzero b's have been found. All interfaces are 
! pure pressure. A pure pressure model requires other changes as well. 
!
  if (nprlev==0) nprlev = plev + 2
!
! Set delta sigma part of layer thickness and reference state midpoint 
! pressures
!
  do k=1,plev
     hybd(k) = hybi(k+1) - hybi(k)
     hypm(k) = hyam(k)*ps0 + hybm(k)*psr
  end do
!
! Reference state interface pressures
!
  do k=1,plevp
     hypi(k) = hyai(k)*ps0 + hybi(k)*psr
  end do
!
! Reference state layer thicknesses
!
  do k=1,plev
     hypd(k) = hypi(k+1) - hypi(k)
  end do
!
! Calculate the log pressure extrapolation factor
!
  prsfac = log( hyam(plev  ) + hybm(plev)) / &
           log((hyam(plev  ) + hybm(plev)) / (hyam(plev-1) + hybm(plev-1)))
!
! Test that A's and B's at full levels are arithmetic means of A's and
! B's at interfaces
!
  do k = 1,plev
     amean = ( hyai(k+1) + hyai(k) )*0.5
     bmean = ( hybi(k+1) + hybi(k) )*0.5
     if(amean == 0. .and. hyam(k) == 0.) then
        atest = 0.
     else
        atest = abs( amean - hyam(k) )/ ( 0.5*( abs(amean + hyam(k)) ) )
     endif
     if(bmean == 0. .and. hybm(k) == 0.) then
        btest = 0.
     else
        btest = abs( bmean - hybm(k) )/ ( 0.5*( abs(bmean + hybm(k)) ) )
     endif
     if (atest > eps) then
        if (masterproc) then
           write(6,9850)
           write(6,*)'k,atest,eps=',k,atest,eps
        end if
!        call endrun
     endif

     if (btest > eps) then
        if (masterproc) then
           write(6,9850)
           write(6,*)'k,btest,eps=',k,btest,eps
        end if
!        call endrun
     endif
  end do

  if (masterproc) then
     write(6,'(a)')'1 Layer Locations (*1000) '
     do k=1,plev
        write(6,9800)k,hyai(k),hybi(k),hyai(k)+hybi(k)
        write(6,9810) hyam(k), hybm(k), hyam(k)+hybm(k)
     end do

     write(6,9800)plevp,hyai(plevp),hybi(plevp),hyai(plevp)+hybi(plevp)
     write(6,9820)
     do k=1,plev
        write(6,9830) k, hypi(k)
        write(6,9840) hypm(k), hypd(k)
     end do
     write(6,9830) plevp, hypi(plevp)
  end if

9800 format( 1x, i3, 3p, 3(f10.4,10x) )
9810 format( 1x, 3x, 3p, 3(10x,f10.4) )
9820 format(1x,'reference pressures (Pa)')
9830 format(1x,i3,f15.4)
9840 format(1x,3x,15x,2f15.4)
9850 format('HYCOEF: A and/or B vertical level coefficients at full',/, &
            ' levels are not the arithmetic mean of half-level values')

  return
end subroutine hycoef

