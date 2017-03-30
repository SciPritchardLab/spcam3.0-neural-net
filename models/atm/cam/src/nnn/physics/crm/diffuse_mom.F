	subroutine diffuse_mom()

	implicit none
	include 'grid.inc'
	integer i,j,k

        if(RUN3D) then
          call diffuse_mom3D()
        else
          call diffuse_mom2D()
        endif


	return
	end

