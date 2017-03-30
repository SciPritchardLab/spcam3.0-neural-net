#include <misc.h>
#include <params.h>

subroutine courlim (vmax2d, vmax2dt, vcour)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Find out whether Courant limiter needs to be applied
! 
! Method: 
! 
! Author: 
! Original version:  CCM2
! 
!-----------------------------------------------------------------------
!
! $Id: courlim.F90,v 1.5.2.3 2003/12/15 18:52:50 hender Exp $
! $Author: hender $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use physconst, only: rga
  use time_manager, only: get_nstep, is_first_step
  use comhd
#ifdef TIMING_BARRIERS
  use mpishorthand
#endif

  implicit none

#include <comsta.h>

!
! Arguments
!
  real(r8), intent(inout) :: vmax2d(plev,plat)   ! Max. wind at each level, latitude
  real(r8), intent(inout) :: vmax2dt(plev,plat)  ! Max. truncated wind at each lvl,lat
  real(r8), intent(inout) :: vcour(plev,plat)    ! Maximum Courant number in slice
!
!--------------------------Local Variables------------------------------
!
  integer k,lat                ! Indices
  integer latarr(1)            ! Output from maxloc (needs to be array for conformability)
  integer :: nstep             ! Current timestep number

  real(r8) vcourmax            ! Max courant number in the vertical wind field
  real(r8) vmax1d(plev)        ! Sqrt of max wind speed 
  real(r8) vmax1dt(plev)       ! Sqrt of max wind speed
  real(r8) cn                  ! Estimate of truncated Courant number 
  real(r8) cnmax               ! Max. courant no. horiz. wind field
  real(r8) psurfsum            ! Summing variable - global mass
  real(r8) stqsum              ! Summing variable - global moisture
  real(r8) rmszsum             ! Summing variable - global vorticity
  real(r8) rmsdsum             ! Summing variable - global divergence
  real(r8) rmstsum             ! Summing variable - global temperature
  real(r8) stps                ! Global Mass integral
  real(r8) stqf                ! Global Moisture integral
  real(r8) rmszf               ! Global RMS Vorticity
  real(r8) rmsdf               ! Global RMS Divergence
  real(r8) rmstf               ! Global RMS Temperature
!
!-----------------------------------------------------------------------
!
#if ( defined SPMD )
#if ( defined TIMING_BARRIERS )
  call t_startf ('sync_realloc7')
  call mpibarrier (mpicom)
  call t_stopf ('sync_realloc7')
#endif
  call t_startf ('realloc7')
  call realloc7 (vmax2d, vmax2dt, vcour)
  call t_stopf ('realloc7')
#endif

  nstep = get_nstep()
!
! Compute maximum wind speed for each level
!
  do k=1,plev
     vmax1d(k)  = sqrt (maxval (vmax2d(k,:)))
     vmax1dt(k) = sqrt (maxval (vmax2dt(k,:)))
  end do
!
! Compute max. vertical Courant number (k is index to Model interfaces)
!
  vcourmax = maxval (vcour(2:,:))
!
! Determine whether the CFL limit has been exceeded for each level
! within the specified range (k<=kmxhdc). Set the truncation wave number
! (for each level independently) so that the CFL limit will not be
! violated and print a message (information only). The trunc wavenumber
! is used later in "hordif" to adjust the diffusion coefficients for
! waves beyond the limit. Store the maximum Courant number for printing
! on the stats line. Note that the max Courant number is not computed
! for the entire vertical domain, just the portion for which the limiter
! is actually applied.
!
  cnmax = 0.
  do k=1,kmxhdc
     cn = vmax1dt(k)*cnfac  ! estimate of truncated Courant number 
     cnmax = max(cnmax,cn)
     if (cn .gt. cnlim) then
        nindex(k) = int(nmaxhd*cnlim/cn + 1.)
        latarr = maxloc (vmax2dt(k,:))
        if (masterproc) write(6,800)k,latarr,cn,nindex(k)-1
     else
        nindex(k) = 2*nmaxhd
     endif
  end do
!
! Write out estimate of original Courant number if limit is exceeded
!
  do k=1,kmxhdc
     cn = vmax1d(k)*cnfac   ! estimate of original Courant number 
     if (cn .gt. cnlim) then
        latarr = maxloc (vmax2d(k,:))
        if (masterproc) write(6,805) k,latarr,cn
     end if
  end do
!
! Compute Max Courant # for whole atmosphere for diagnostic output
!
  cnmax = 0.
  do k=1,plev-1
     cn = vmax1dt(k)*cnfac  ! estimate of Courant number
     cnmax = max(cnmax,cn)
  end do
!
! Write out statisitics to standard output
!
  psurfsum = 0.
  stqsum = 0.
  rmszsum = 0.
  rmsdsum = 0.
  rmstsum = 0.

  do lat=1,plat
     psurfsum = psurfsum + psurf(lat)
     stqsum = stqsum + stq(lat)
     rmszsum = rmszsum + rmsz(lat)
     rmsdsum = rmsdsum + rmsd(lat)
     rmstsum = rmstsum + rmst(lat)
  end do

  stps = 0.5*psurfsum
  stqf = 0.5*rga*stqsum
  rmszf = sqrt(0.5*rmszsum)
  rmsdf = sqrt(0.5*rmsdsum)
  rmstf = sqrt(0.5*rmstsum)
  if (masterproc) then
     if (is_first_step()) write(6,810)
     write (6,820) nstep, rmszf, rmsdf, rmstf, stps, stqf, cnmax, vcourmax
  end if
!
  return
!
! Formats
!
800 format('COURLIM: *** Courant limit exceeded at k,lat=',2i3,  &
       ' (estimate = ',f6.3, '), solution has been truncated to ',  &
       'wavenumber ',i3,' ***')
805 format(' *** Original  Courant limit exceeded at k,lat=',2i3,  &
       ' (estimate = ',f6.3,')',' ***')
810 format(/109x,'COURANT'/10x,'NSTEP',4x,'RMSZ',19x,'RMSD',19x,  &
       'RMST',4x,'STPS',9x,'STQ',19x,'HOR   VERT')
820 format(' NSTEP =',i8,1x,1p2e23.15,0p1f8.3,1p1e13.5,e23.15,  &
       0p1f5.2,f6.2)
end subroutine courlim

