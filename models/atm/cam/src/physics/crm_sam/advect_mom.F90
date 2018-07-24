subroutine advect_mom

use vars

implicit none
integer i,j,k

if(docolumn) return


call advect2_mom_xy()
call advect2_mom_z()



end subroutine advect_mom

