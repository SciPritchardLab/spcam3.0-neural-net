subroutine diffuse_mom

!  Interface to the diffusion routines

use vars
implicit none
integer i,j,k

if(RUN3D) then
   call diffuse_mom3D()
else
   call diffuse_mom2D()
endif


end subroutine diffuse_mom

