
subroutine advect_scalar (f,fadv,flux,f2leadv,f2legrad,fwleadv,doit)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
use vars, only: u, v, w, rho, rhow

implicit none

real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real flux(nz), fadv(nz)
real f2leadv(nzm),f2legrad(nzm),fwleadv(nzm)
logical doit

real coef
integer i,j,k

if(docolumn) flux = 0.

if(docolumn) return


if(RUN3D) then
  call advect_scalar3D(f, u, v, w, rho, rhow, flux)
else
  call advect_scalar2D(f, u, w, rho, rhow, flux)	  
endif



end subroutine advect_scalar

