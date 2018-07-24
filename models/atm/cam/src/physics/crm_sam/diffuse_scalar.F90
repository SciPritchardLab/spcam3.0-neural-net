subroutine diffuse_scalar (f,fluxb,fluxt, &
                          fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

use grid
use vars, only: tkh, rho, rhow
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real flux(nz)
real fdiff(nz)
real f2lediff(nzm)
real f2lediss(nzm)
real fwlediff(nzm)
logical doit
! Local
real r2dx,r2dy,r2dx0,r2dy0,r2dz
integer i,j,k,kb,kc,jb,jc


if(RUN3D) then
  call diffuse_scalar3D (f,fluxb,fluxt,tkh,rho,rhow,flux)
else  
  call diffuse_scalar2D (f,fluxb,fluxt,tkh,rho,rhow,flux)
endif


end subroutine diffuse_scalar 
