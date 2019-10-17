#define BOMBDEBUG
module heat_bombs

   use shr_kind_mod, only: r8 => shr_kind_r8

implicit none

integer, parameter :: max_bombs = 5
integer, parameter :: bomb_maxsteps = 480 !TODO make stochastic
integer, parameter :: taper_steps = 48

integer :: nactive

real(r8) :: bomb_lon(max_bombs),bomb_lat(max_bombs)
real(r8) :: bomb_taperval(max_bombs)
real(r8) :: bomb_pmid(max_bombs)
integer :: bomb_bornstep(max_bombs)

contains

subroutine heating_bombs (lchnk   ,ncol   )

   use ppgrid
   use phys_grid,     only: get_rlat_p, get_rlon_p
   use pmgrid, only: iam, masterproc
   use mpishorthand
   use time_manager,    only: get_nstep

   integer, intent(in) :: lchnk, ncol
   integer :: nstep,step,i,intaux
   real :: xrand 
   real :: lon,lat
   nstep = get_nstep()
   if (nstep .le. 8) then
     nactive = 0
   endif

! Manage the bomb portfolio and make sure all MPI tasks are aware.
   if (masterproc) then
     write (6,*) 'Masterproc made it nstep=',nstep
     if (nstep .le. 8) then
       nactive = 0
     else ! far enough past init to start dropping bombs?
       if (nactive .lt. max_bombs) then  ! need to drop another?
         call random_number(xrand)
         if (xrand .le. 0.005 ) then ! % chance of new bomb drop
            if (nactive .eq. 0 .or. (nactive .gt. 0 .and. nstep .gt. bomb_bornstep(nactive))) then ! Don't make two on same step
            nactive = nactive+1
            i = nactive
            call random_number(xrand)
            bomb_lon(i) = xrand*360
            call random_number(xrand)
            bomb_lat(i) = (xrand*160. - 80.)
            call random_number(xrand)
            bomb_pmid(i) = 100000. - xrand*80000. ! 1000-200 hPa
            bomb_bornstep(i) = nstep
#ifdef BOMBDEBUG
            write (6,*) 'HEY master just dropped bomb number',i,':lon,lat,pres=',bomb_lon(i),bomb_lat(i),bomb_pmid(i)
#endif
         endif
         endif
       endif ! need to drop another bomb?
     
       ! Need to kill any bombs?
       if ( (nstep-bomb_bornstep(1)+1) .eq. bomb_maxsteps) then
       ! Note that by constructions i=1 is always the oldest bomb in the
       ! portfolio, so will be the only slated to die.
         write (6,*) 'HEY master decided to kill a bomb...'
         nactive = nactive - 1
         do i = 1,nactive-1
           bomb_lon(i) = bomb_lon(i+1)
           bomb_lat(i) = bomb_lat(i+1)
           bomb_pmid(i) = bomb_pmid(i+1)
           bomb_bornstep(i) = bomb_bornstep(i+1)
         end do 
       end if

       ! Ensure gentle births and deaths, i.e.
       ! Calculate bomb taper based on current time step:
       do i= 1,nactive
         step = nstep-bomb_bornstep(i) + 1
         if (step .le. taper_steps) then
           bomb_taperval(i) = step/taper_steps ! 0-->1
         elseif(step .gt. bomb_maxsteps-taper_steps) then
           bomb_taperval(i) = (bomb_maxsteps-step)/taper_steps ! 1-->0 
         else
           bomb_taperval(i) = 1.
         endif
       end do
     endif ! far enough past init to do anything?
   endif ! masterproc's first assigned grid cell? (need to only do this once or trouble!)
   
    ! Broadcast the info to the other MPI tasks:
   call mpibcast (nactive,1,mpiint,0,mpicom)
   call mpibcast (bomb_bornstep(:),max_bombs,mpiint,0,mpicom)          

#ifdef BOMBDEBUG
   if (nactive .gt. 0 .and. iam .eq. 4) then
     do i=1,nactive 
       write (6,*) 'HEY iam=',iam,' got bomb(',i,') born on step ',bomb_bornstep(i)
     end do
   end if
#endif 
end subroutine heating_bombs


end module heat_bombs
