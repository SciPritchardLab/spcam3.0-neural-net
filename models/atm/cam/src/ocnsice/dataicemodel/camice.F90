#include <misc.h>
#include <params.h>

subroutine camice(srf_state,srfflx)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! CAM sea ice surface fluxes.
!
! Method: 
! 
! Author:
! 
!-----------------------------------------------------------------------
!
! $Id: camice.F90,v 1.2.2.1 2002/06/15 13:48:48 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use pspect
  use comsrf, only: surface_state,srfflx_parm,icefrac,snowhice
  use phys_grid, only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use sst_data,  only: sstan
  use time_manager, only: get_nstep, get_step_size, get_curr_calday

  implicit none
!
! Input/Output arguments
!
   type(surface_state), intent(inout), dimension(begchunk:endchunk) :: srf_state
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx

#include <comctl.h>

!---------------------------Local variables-----------------------------
  integer :: nstep          ! current timestep number
  integer :: dtime          ! timestep size [seconds]
  real(r8) cdaynext         ! calendar day for next timestep
  real(r8) clat(pcols)      ! current latitudes(radians)
  real(r8) clon(pcols)      ! current longitudes(radians)
  real(r8) cosznext(pcols)  ! cosine solar zenith angle next timestep
  integer ncol              ! number of columns in chunk
  integer lchnk             ! chunk index
  integer idum1,idum2,idum3,idum4 ! temporary variables
!-----------------------------------------------------------------------
!
! Calendar day for next time step
!
  nstep = get_nstep()
  dtime = get_step_size()
  cdaynext = get_curr_calday(offset=dtime)


  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
!
! Sea ice surface fluxes and temperatures
!
     call srfsice(lchnk     , ncol    ,icefrac(1,lchnk),&
                    snowhice(1,lchnk) , srf_state(lchnk)%ubot, &
                  srf_state(lchnk)%vbot ,srf_state(lchnk)%tbot, &
                    srf_state(lchnk)%qbot ,srf_state(lchnk)%thbot, &
                    srf_state(lchnk)%zbot , &
                  srf_state(lchnk)%pbot, srf_state(lchnk)%srfrad, &
                    srf_state(lchnk)%tssub, srfflx(lchnk)%cflx, &
                    srfflx(lchnk)%wsx , &
                  srfflx(lchnk)%wsy ,srfflx(lchnk)%ts  , &
                    srfflx(lchnk)%shf ,srfflx(lchnk)%lhf , &
                    srfflx(lchnk)%lwup, &
                  srfflx(lchnk)%tref   )
!
! Albedos for next time step 
!
     call get_rlat_all_p(lchnk, ncol, clat)
     call get_rlon_all_p(lchnk, ncol, clon)
     call zenith (cdaynext, clat, clon, cosznext, ncol)
     call albice(lchnk,ncol, &
                 snowhice(1,lchnk), cosznext, &
                 srfflx(lchnk)%asdir, srfflx(lchnk)%aldir, &
                 srfflx(lchnk)%asdif, srfflx(lchnk)%aldif)
  end do

  return
end subroutine camice
