#include <params.h>
#define CAER_NBR 4


module caer

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Simple carbon aerosol model for MATCH.
   !
   ! Author: B. Eaton
   !-----------------------------------------------------------------------

   use pmgrid
#ifdef MATCH
   use precision
#else
  use shr_kind_mod,only: r8 => shr_kind_r8
#endif
  use ppgrid,      only: pcols, pver

   implicit none
   save
   private
   public ::  caerini      ! initialization
   public ::  caersf       ! return surface fluxes
   public ::  caercv       ! 1st order chemistry
   public ::  docaer       ! returns .true. when carbon aerosol species are present
   public ::  caer_ncomp   ! returns number of carbon aerosol components
   public ::  caer_idx1    ! returns constituent index of 1st carbon aerosol component
   public ::  caer_idx_set ! sets constituent index of 1st carbon aerosol component

   real(r8), parameter ::&
      freq = 1./(1.2*24.*3600.),  &! 1 / (e-folding time in seconds)
      mwrat = 1.                   ! ratio of mol. wts of cphil to cphob

   integer, parameter ::  ncomp = CAER_NBR    ! number of carbon aerosol components

   integer idx1  ! start constituent index for first carbon aerosol component

!##############################################################################
contains
!##############################################################################

   subroutine caerini()

      !----------------------------------------------------------------------- 
      ! Purpose: Initialize carbon aerosol module.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      !-----------------------------------------------------------------------

      if ( docaer() ) then
         write(6,*)'Carbon aerosol configured with ',ncomp,' components.'
         write(6,*)'  Hydrophobic components converted to hydrophilic with 1.2 day e-folding time.'
      else
         write(6,*)'Carbon aerosols disabled.'
      end if

   end subroutine caerini

   subroutine caer_idx_set(m)

      !----------------------------------------------------------------------- 
      ! Purpose: set the start index for the carbon species
      ! 
      !-----------------------------------------------------------------------

      implicit none
      integer m

      idx1 = m
      write (6,*) ' start index for carbon species is ', idx1

      !-----------------------------------------------------------------------


   end subroutine caer_idx_set

!#######################################################################

   subroutine caersf( ncol, nemis, lat, lon, sflx )

      ! Set carbon aerosol surface fluxes.

      use caerbnd, only: caerbndget

      implicit none

      integer, intent(in) ::       ncol              ! number of columns to produce emissions for
      integer, intent(in) ::       nemis             ! number of carbon aerosol species with emissions
      integer, intent(in) ::       lat(pcols)         ! model latitude index
      integer, intent(in) ::       lon(pcols)         ! model longitude index

      real(r8), intent(out) ::&
         sflx(pcols,nemis)   ! surface flux in kg/m^2/s

      ! Local variables:
      integer i, m
      real(r8) ::&
         tmp(pcols,5)
      !------------------------------------------------------------------------------

      ! Get a chunk of carbon aerosol emissions data.
      call caerbndget( lat, lon, ncol, tmp )

      if ( nemis .eq. 1 ) then
         do i = 1, ncol
            sflx(i,1) = tmp(i,1)+tmp(i,2)+tmp(i,3)+tmp(i,4)+tmp(i,5)  ! total carbon
         end do
      else if ( nemis .eq. 2 ) then
         do i = 1, ncol
            sflx(i,1) = tmp(i,1) + tmp(i,3) + tmp(i,5)  ! organic carbon
            sflx(i,2) = tmp(i,2) + tmp(i,4)             ! black carbon
         end do
      else if ( nemis .eq. 5 ) then
         do m = 1, 5
            do i = 1, ncol
               sflx(i,m) = tmp(i,m)
            end do
         end do
      end if

   end subroutine caersf

!##############################################################################

   subroutine caercv( ncol, A, deltat, Asrc, Bsrc )

      ! A undergoes 1st order decay to B.

      implicit none

      integer ncol

      real(r8), intent(in) ::&
         A(pcols,pver),         &! mixing ratio (kg A/(kg moist air))
         deltat                  ! time step in seconds

      real(r8), intent(out) ::&
         Asrc(pcols,pver),      &! conversion rate (kg A /(s kg moist air))
         Bsrc(pcols,pver)        ! conversion rate (kg B /(s kg moist air))

      ! Local variables:
      integer i, k
      !-----------------------------------------------------------------------

      ! calculate tendencies using Backward Euler

      do k = 1,pver
         do i = 1,ncol
            Asrc(i,k) = -A(i,k)*freq / (1. + freq*deltat)
            Bsrc(i,k) = -Asrc(i,k)*mwrat
         end do
      end do

   end subroutine caercv

!#######################################################################

   logical function docaer ()

      !----------------------------------------------------------------------- 
      ! Purpose: Return .true. when carbon aerosols are active.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      docaer = ( ncomp > 0 )

   end function docaer

!#######################################################################

   integer function caer_ncomp()
      implicit none
      caer_ncomp = ncomp
   end function caer_ncomp

!#######################################################################

   integer function caer_idx1()
      implicit none
      caer_idx1 = idx1
   end function caer_idx1

!#######################################################################

end module caer
