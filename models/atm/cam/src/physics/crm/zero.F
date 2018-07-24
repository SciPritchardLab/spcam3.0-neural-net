c-----------------------------------------------------------------
c-----------------------------------------------------------------
c
c       CIMMS/CSU LES - Version X 
c
c       (C) 1994-1998, Marat Khairoutdinov 
c-----------------------------------------------------------------

	subroutine zero()
	
	implicit none
	include 'common.inc'
	
	integer i,j,k
	
	dudt(:,:,:,na) = 0.
	dvdt(:,:,:,na) = 0.
	dwdt(:,:,:,na) = 0.

	return
	end
