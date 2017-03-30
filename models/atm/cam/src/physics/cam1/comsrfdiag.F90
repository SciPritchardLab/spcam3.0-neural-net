#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! !MODULE: comsrfdiag
!
! !DESCRIPTION:	Module to handle surface fluxes for the subcomponents of cam/csm
!
! Public interfaces:
!
!	update       
!	init           
!	zero           
!
!-----------------------------------------------------------------------
module comsrfdiag
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, pnats
  use ppgrid, only: pcols, pvermx
  use phys_grid, only: get_ncols_p
  use infnan
  use pmgrid, only: plon, plat
  use comsrf
  use history, only: outfld
  implicit none

!----------------------------------------------------------------------- 
! PUBLIC: Make default data and interfaces private
!----------------------------------------------------------------------- 
!
! ! PUBLIC MEMBER FUNCTIONS:
!
! Public interfaces

  public output_flns_fsns_fluxes
  public output_shf_lhf_fluxes
  public output_shfoi_lhfoi_fluxes

!===============================================================================
CONTAINS
!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================
!----------------------------------------------------------------------- 

  subroutine output_flns_fsns_fluxes(surface_state1d,c)
!-----------------------------------------------------------------------
!
! Purpose:
! output flsn/fsns fluxes for component models and som off line code
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------
    use time_manager,  only: get_nstep, is_first_step
    implicit none
    type(surface_state), intent(in) :: surface_state1d
    integer, intent(in) :: c

    real(r8) :: ocnwt(pcols)
    real(r8) :: icewt(pcols)
    real(r8) :: tmpocn(pcols)
    real(r8) :: tmpocnice(pcols)
    real(r8) :: tmplnd(pcols)
    real(r8) :: tmpice(pcols)
    integer :: ncol,i
    integer itim, ifld

    ncol = get_ncols_p(c)
    
    tmplnd(:ncol) = lwuplnd(:ncol,c) - surface_state1d%flwds(:ncol)
    tmpice(:ncol) = lwupice(:ncol,c) - surface_state1d%flwds(:ncol)
    tmpocn(:ncol) = lwupocn(:ncol,c) - surface_state1d%flwds(:ncol)

    do i = 1,ncol
       if (1-landfrac(i,c) .ne. 0.) then
          ocnwt(i) = ocnfrac(i,c)/(1-landfrac(i,c))
          icewt(i) = icefrac(i,c)/(1-landfrac(i,c))
       else
          ocnwt(i) = 0.
          icewt(i) = 0.
       endif
    end do
! flux over ocn+ice needed for som, weight each 

    tmpocnice(:ncol) = tmpice(:ncol)*icewt(:ncol) +tmpocn(:ncol)*ocnwt(:ncol)

    tmplnd(:ncol)=tmplnd(:ncol)*landfrac(:ncol,c)
    tmpocn(:ncol)=tmpocn(:ncol)*ocnfrac(:ncol,c)
    tmpice(:ncol)=tmpice(:ncol)*icefrac(:ncol,c)
    call outfld ('FLNSOI', tmpocnice, pcols, c)
    call outfld ('FLNSLND', tmplnd, pcols, c)
    call outfld ('FLNSICE', tmpice, pcols, c)
    call outfld ('FLNSOCN', tmpocn, pcols, c)

    tmpocn(:ncol) = (surface_state1d%sols (:ncol)*(1. - asdirocn(:ncol,c)) + &
         surface_state1d%solsd(:ncol)*(1. - asdifocn(:ncol,c)) + &
         surface_state1d%soll (:ncol)*(1. - aldirocn(:ncol,c)) + &
         surface_state1d%solld(:ncol)*(1. - aldifocn(:ncol,c)))
    tmpice(:ncol) = (surface_state1d%sols (:ncol)*(1. - asdirice(:ncol,c)) + &
         surface_state1d%solsd(:ncol)*(1. - asdifice(:ncol,c)) + &
         surface_state1d%soll (:ncol)*(1. - aldirice(:ncol,c)) + &
         surface_state1d%solld(:ncol)*(1. - aldifice(:ncol,c)))
    tmplnd(:ncol) = (surface_state1d%sols (:ncol)*(1. - asdirlnd(:ncol,c)) + &
         surface_state1d%solsd(:ncol)*(1. - asdiflnd(:ncol,c)) + &
         surface_state1d%soll (:ncol)*(1. - aldirlnd(:ncol,c)) + &
         surface_state1d%solld(:ncol)*(1. - aldiflnd(:ncol,c)))

    tmpocnice(:ncol)= tmpice(:ncol)*icewt(:ncol)+tmpocn(:ncol)*ocnwt(:ncol)

    tmplnd(:ncol)=tmplnd(:ncol)*landfrac(:ncol,c)
    tmpice(:ncol)=tmpice(:ncol)*icefrac(:ncol,c)
    tmpocn(:ncol)=tmpocn(:ncol)*ocnfrac(:ncol,c)

    call outfld ('FSNSOI',  tmpocnice, pcols, c)
    call outfld ('FSNSLND', tmplnd, pcols, c)
    call outfld ('FSNSICE', tmpice, pcols, c)
    call outfld ('FSNSOCN', tmpocn, pcols, c)
    return
  end subroutine output_flns_fsns_fluxes

  subroutine output_shf_lhf_fluxes (srfflx, c, ncol, sfcfrac, sfctype)
!-----------------------------------------------------------------------
!
! Purpose:
! output shf and lhf fluxes for each surface types
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------
     use ppgrid, only: pcols
     use phys_grid, only: get_ncols_p
     use history, only: outfld

     implicit none !

! Input/Output arguments

     type(srfflx_parm), intent(in) :: srfflx
     integer, intent(in) :: c
     integer, intent(in) :: ncol
     character(len=*),intent(in) :: sfctype !  ! Local workspace !
     real(r8), intent(in) :: sfcfrac(pcols)
     
     integer :: i
     real(r8) :: lhflx(pcols)
     real(r8) :: shflx(pcols)
     
     do i=1,ncol
        if (sfcfrac(i) > 0.) then
           lhflx(i) = srfflx%lhf(i)*sfcfrac(i)
           shflx(i) = srfflx%shf(i)*sfcfrac(i)
        else
           lhflx(i) = 0.
           shflx(i) = 0.
        end if
     end do
     call outfld ('LHFLX'//sfctype, lhflx, pcols, c)
     call outfld ('SHFLX'//sfctype, shflx, pcols, c)
     
     return
  end subroutine output_shf_lhf_fluxes

  subroutine output_shfoi_lhfoi_fluxes (srfflx_ocn, srfflx_ice, c)
!-----------------------------------------------------------------------
!
! Purpose:
! output shf and lhf fluxes for each surface types
!
! Method:
!
! Author: John Truesdale
!
!-----------------------------------------------------------------------

    use ppgrid, only: pcols
    use phys_grid, only: get_ncols_p
    use history, only: outfld
    
    implicit none !  ! Input/Output arguments !
    type(srfflx_parm), intent(in) :: srfflx_ocn,srfflx_ice
    integer, intent(in) :: c
    
    integer :: i
    integer :: ncol
    real(r8) :: ocnwt(pcols)
    real(r8) :: icewt(pcols)
    real(r8) :: lhflxoi(pcols)
    real(r8) :: shflxoi(pcols)

    ncol = get_ncols_p(c)
    do i=1,ncol
       if (1-landfrac(i,c) .ne. 0.) then
          ocnwt(i) = ocnfrac(i,c)/(1-landfrac(i,c))
          icewt(i) = icefrac(i,c)/(1-landfrac(i,c))
       else
          ocnwt(i) = 0.
          icewt(i) = 0.
       endif
       if (icefrac(i,c) > 0. .or. ocnfrac(i,c) > 0.) then
          lhflxoi(i) = srfflx_ice%lhf(i)*icewt(i)+ srfflx_ocn%lhf(i)*ocnwt(i)
          shflxoi(i) = srfflx_ice%shf(i)*icewt(i)+ srfflx_ocn%shf(i)*ocnwt(i)
       else
          lhflxoi(i) = 0.
          shflxoi(i) = 0.
       end if
    end do
    call outfld ('LHFLXOI', lhflxoi, pcols, c)
    call outfld ('SHFLXOI', shflxoi, pcols, c)

    return
  end subroutine output_shfoi_lhfoi_fluxes
  
end module comsrfdiag
