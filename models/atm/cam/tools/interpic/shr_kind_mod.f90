!===============================================================================
! CVS: $Id: shr_kind_mod.f90,v 1.1.2.1 2002/06/15 13:50:05 erik Exp $
! CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/tools/atm_tools/ccm_src/interpic/Attic/shr_kind_mod.f90,v $
! CVS: $Name:  $
!===============================================================================

MODULE shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer

END MODULE shr_kind_mod
