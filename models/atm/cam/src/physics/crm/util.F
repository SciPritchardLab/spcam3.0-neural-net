      integer function lenstr (string)
      
      ! returns string's length ignoring the rightmost blank and null characters
      implicit none
      character *(*) string
      integer k
      lenstr = 0
      do k = 1,len(string)
       if (string(k:k).ne.' '.and.string(k:k).ne.char(0)) then
            lenstr = lenstr+1
       end if
      end do
 111  return
      end



		
	
	
	subroutine fminmax_print(name,f,dimx1,dimx2,dimy1,dimy2,dimz)

	implicit none
	include 'grid.inc'
	integer dimx1, dimx2, dimy1, dimy2, dimz
	real f(dimx1:dimx2, dimy1:dimy2, dimz),fmn(nz),fmx(nz)
	character *(*) name
	real fmin(1),fmax(1),fff(1)
	integer i,j,k,m
	
	do k=1,dimz
	 fmn(k) = 1.e30
	 fmx(k) =-1.e30
	 do j=1,ny
	  do i=1,nx
	   fmn(k) = min(fmn(k),f(i,j,k))
	   fmx(k) = max(fmx(k),f(i,j,k))
	  end do
	 enddo
	enddo
	fmin(1) = 1.e30
	fmax(1) =-1.e30
	do k=1,dimz
	 fmin(1) = min(fmin(1),fmn(k))
	 fmax(1) = max(fmax(1),fmx(k))
	end do
		
	if(dompi) then
	  fff(1)=fmax(1)
	  call task_max_real(fff(1),fmax(1),1)
	  fff(1)=fmin(1)
	  call task_min_real(fff(1),fmin(1),1)
	end if
	if(rank.eq.0) print *,name,fmin,fmax
	return
	end



	
	subroutine setvalue(f,n,f0)
	real f(n), f0
	integer k
	do k=1,n
	 f(k)=f0
	end do
	return
	end
