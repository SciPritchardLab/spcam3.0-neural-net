#include <misc.h>
#include <params.h>

subroutine cldinti ()

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,  only: plev, plevp, masterproc
   use cldconst

#include <comhyb.h>
!
! Find vertical level nearest 700 mb
!
   k700 = 1
   do k=1,plev-1
      if (hypm(k) < 7.e4 .and. hypm(k+1) >= 7.e4) then
         if (7.e4-hypm(k) < hypm(k+1)-7.e4) then
            k700 = k
         else
            k700 = k + 1
         end if
         goto 20
      end if
   end do

   write(6,*)'INTI: model levels bracketing 700 mb not found'
   call endrun

20 continue

   if (masterproc) then
      write(6,*)'INTI: model level nearest 700 mb is',k700,'which is',hypm(k700),'pascals'
   end if

   return
end subroutine cldinti
