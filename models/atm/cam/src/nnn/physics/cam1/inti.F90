#include <misc.h>
#include <params.h>

subroutine inti ()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set constants and call initialization procedures for time independent
! physics routines
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,             only: plev, plevp          ! Needed for hypm passed to vd_inti
   use chemistry,          only: trace_gas, chem_init
   use aerosol_intr,       only: prognostic_aerosol_initialize
   use gw_drag,            only: gw_inti
   use vertical_diffusion, only: vd_inti
   use moistconvection,    only: mfinti
   use cloud_fraction,     only: cldfrc_init
   use cldcond,            only: cldcond_init
   use param_cldoptics,    only: param_cldoptics_init
   use zm_conv,            only: zm_convi
   use shr_const_mod,      only: shr_const_zvir, shr_const_cpwv, shr_const_rwv
   use physconst,          only: rair, cpair, cpwv, gravit, stebol, epsilo, tmelt, &
                                 latvap, latice, rh2o, zvir, cpvir, rhoh2o, pstd,  &
                                 karman, rhodair
   use aerosols,           only: aerosol_initialize
   use aer_optics,         only: aer_optics_initialize
   use check_energy,       only: check_energy_init

   implicit none

#include <comctl.h>
#include <comhyb.h>

!
! Initialize radiation data for aerosol forcing calculation
! Initialize aerosol fields from files
!
   call aer_optics_initialize()
   call aerosol_initialize()
   
   call prognostic_aerosol_initialize()

!
!-----------------------------------------------------------------------
!
! Initialize physconst variables
! In adiabatic case, set zvir and cpvir explicitly to zero instead of 
! computing as (rh2o/rair - 1.) and (cpwv/cpair - 1.) respectively, in order 
! to guarantee an identical zero.
!
   if (adiabatic .or. ideal_phys) then
      rh2o  = rair
      zvir  = 0.
      cpwv  = cpair
      cpvir = 0.
   else
      rh2o  = shr_const_rwv
      zvir  = shr_const_zvir
      cpwv  = shr_const_cpwv
      cpvir = cpwv/cpair - 1.
   end if
!
! Call time independent initialization routines for parameterizations.
!
   if (trace_gas) call chem_init
   call gw_inti (cpair   ,cpwv    ,gravit  ,rair    ,hypi    )
   call vd_inti (cpair   ,cpwv    ,gravit  ,rair    ,zvir   , &
                 hypm    ,karman    )
   call tsinti  (tmelt   ,latvap  ,rair    ,stebol  ,latice  )
   call radini  (gravit  ,cpair   ,epsilo  ,stebol  ,pstd*10.0 )
   call esinti  (epsilo  ,latvap  ,latice  ,rh2o    ,cpair  , &
                 tmelt   )
   call mfinti  (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )
   call cldfrc_init
   call zm_convi( tmelt, epsilo, latvap, cpair )
   call cldinti ()

   call cldcond_init
   call param_cldoptics_init
   call check_energy_init

   return
end subroutine inti
