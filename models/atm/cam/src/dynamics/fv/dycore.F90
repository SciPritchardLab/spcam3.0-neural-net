module dycore
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: dycore
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC     dycore_is,get_resolution

!
! !DESCRIPTION:
!
! Utility routines related to the dycore
!
!      \begin{tabular}{|l|l|} \hline \hline
!        dycore         & dycore being used \\ \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.01   Rosinski    Creation
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dycore_is --- determine dynamical core in use
!
! !INTERFACE: 
   logical function dycore_is (name)

      implicit none

! !INPUT PARAMETERS:
      character(len=*) :: name

! !DESCRIPTION:
!   Determine the dynamical core in use.  
!   True if the dycore name is "lr" or "LR"
!
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC


      if (name == 'lr' .or. name == 'LR') then
         dycore_is = .true.
      else
         dycore_is = .false.
      end if

      return
!EOC
   end function dycore_is

   character(len=7) function get_resolution()
     use pmgrid, only: plat
!
! Input arguments
!
     
     select case ( plat )
     case ( 91 )
        get_resolution = '2x2.5'
     case default
        get_resolution = 'UNKNOWN'
     end select

     return
   end function get_resolution

end module dycore
