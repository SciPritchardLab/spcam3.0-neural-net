#include <misc.h>

subroutine print_memusage (string)
!
! Purpose: interface function to C routine print_memusage
!
! Method: Handle all character strings locally to avoid machine-dependent
!         inter-language communication of character variables. 
!         print_memusage will return variables in native types.
!
   use pmgrid, only: iam
!
! Arguments
!
   character(len=*), intent(in) :: string
!
! Local workspace
!
   integer :: size         ! process size
   integer :: rss          ! process resident set size
   integer :: share        ! process shared memory size
   integer :: text         ! process text size
   integer :: datastack    ! data + stack memory
   integer :: ret          ! return code from get_memusage
!
! Externals
!
   integer, external :: get_memusage

   ret = get_memusage (size, rss, share, text, datastack)
   write(0,'(a,i3," ",a,a)')'print_memusage iam ', iam, string, '. -1 in the next line means unavailable'
   if (ret == 0) then
      write(0,'(a,5i8)')'print_memusage: size, rss, share, text, datastack=', &
                                         size, rss, share, text, datastack
   else
      write(0,*)'print_memusage: get_memusage returns -1'
   end if

   return
end subroutine print_memusage
