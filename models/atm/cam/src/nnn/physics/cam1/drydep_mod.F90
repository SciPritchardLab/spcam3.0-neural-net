module drydep_mod
#include <params.h>

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use pmgrid,  only: plat
  use constituents, only: pcnst, pnats

      ! Shared Data for dry deposition calculation.

      real(r8) rair                ! Gas constant for dry air (J/K/kg)
      real(r8) gravit              ! Gravitational acceleration
      real(r8) phi(plat)           ! grid latitudes (radians)11

contains

      subroutine drydepdr( lchnk, ncol, month, dtime, landfrac, icefrac, ocnfrac, &
                           sn, prect, &
                           t2m, wvflx, shflx, pblh, ustar, obklen, &
                           pmbot, rpdel, tvmbot, zi, um1, vm1, tm1, &
                           swradsf, tr, tend, obuf )

! Driver for dry deposition parameterization.  The default version of this
! routine does nothing.  It is provided only as a template.
! N.B. The input arguments may not all be required for some parameterizations
!      but they are included in an attempt to create a general interface that
!      is meets the needs of a range of parameterizations.

#ifdef MATCH
      use histout, only: outfld
#else 
      use history, only: outfld
#endif
      use scyc, only: doscyc, scyc_idx1
#ifdef CAER
!++caer
      use caer, only: docaer, caer_ncomp, caer_idx1
!--caer
#endif

!-----------------------------------------------------------------------
#ifdef MATCH
#include <tracnam.h>
#else
  use constituents,only: cnst_name
#endif

!-----------------------------------------------------------------------
      implicit none

      ! Input arguments
      integer, intent(in)  ::   lchnk                 ! chunk
      integer, intent(in)  ::   ncol                  ! number of cols
      integer, intent(in)  ::   month                 ! month of year [1-12]
      real(r8), intent(in) ::   dtime                 ! timestep size
      real(r8), intent(in) ::   landfrac(pcols)       ! land fraction
      real(r8), intent(in) ::   icefrac(pcols)        ! ice fraction
      real(r8), intent(in) ::   ocnfrac(pcols)        ! ocean fraction
      real(r8), intent(in) ::   sn(pcols)             ! snow; equivalent water amount (m)
      real(r8), intent(in) ::   prect(pcols)          ! total precipitation (units arbitrary as long
                                                      ! the value is > 1.e-12 when it's raining)
      real(r8), intent(in) ::   t2m(pcols)            ! temperature at 2 m
      real(r8), intent(in) ::   wvflx(pcols)          ! water vapor flux (kg/m^2/s)
      real(r8), intent(in) ::   shflx(pcols)          ! sensible heat flux (w/m^2)
      real(r8), intent(in) ::   pblh(pcols)           ! planetary boundary layer height
      real(r8), intent(in) ::   ustar(pcols)          ! friction velocity
      real(r8), intent(in) ::   obklen(pcols)         ! Obukhov length
      real(r8), intent(in) ::   pmbot(pcols)          ! midpoint pressure in bottom layer
      real(r8), intent(in) ::   rpdel(pcols,pver)     ! 1./(pint(k+1)-pint(k))
      real(r8), intent(in) ::   tvmbot(pcols)         ! midpoint virtual temperature in bottom layer
      real(r8), intent(in) ::   zi(pcols,pverp)       ! height at interface of levels
      real(r8), intent(in) ::   um1(pcols,pver)       ! u horizontal wind component
      real(r8), intent(in) ::   vm1(pcols,pver)       ! v horizontal wind component
      real(r8), intent(in) ::   tm1(pcols,pver)       ! temperature
      real(r8), intent(in) ::   swradsf(pcols)        ! solar radiation at the surface (W/m2)
      real(r8), intent(in) :: tr(pcols,pver,pcnst)    ! tracers

! Output arguments.
      real(r8) obuf(*)               ! output buffer
      real(r8), intent(out) :: tend(pcols,pver,pcnst)    ! tracer tends

! Local variables
      integer        i, m, mm
      real(r8)  dvel(pcols)            ! deposition velocity
      real(r8)  flx(pcols,pcnst)       ! dry deposition fluxes
!----------------------------------------------------------------------c


      do m = 1, pcnst
         do i = 1, ncol
            flx(i,m) = 0.0
         end do
      end do


!     Update tracers with dry deposition contribution.
      tend = 0
      do m = 1, pcnst
         do i = 1, ncol
            tend(i,pver,m) =  flx(i,m)*gravit*rpdel(i,pver)
         end do
      end do

      return
      end subroutine drydepdr

!##############################################################################

! $Id: drydep_mod.F90,v 1.1.6.2 2003/08/06 17:48:51 rosinski Exp $
#include <params.h>

      subroutine inidrydep( xrair, xgravit, xphi )

! Initialize dry deposition parameterization.

      implicit none

! Input arguments:
      real(r8), intent(in) :: xrair                ! Gas constant for dry air
      real(r8), intent(in) :: xgravit              ! Gravitational acceleration
      real(r8), intent(in) :: xphi(plat)           ! grid latitudes (radians)

! Local variables:
      integer i, j, ncid, vid
!-----------------------------------------------------------------------

      rair = xrair
      gravit = xgravit
      do j = 1, plat
         phi(j) = xphi(j)
      end do

      return
      end subroutine inidrydep

!##############################################################################

      subroutine setdvel( ncol, landfrac, icefrac, ocnfrac, vgl, vgo, vgsi, vg )

! Set the deposition velocity depending on whether we are over
! land, ocean, and snow/ice


      implicit none

! Input arguments:

      integer, intent(in) :: ncol
      real (r8), intent(in) :: landfrac(pcols)       ! land fraction
      real (r8), intent(in) :: icefrac(pcols)       ! ice fraction
      real (r8), intent(in) :: ocnfrac(pcols)       ! ocean fraction

      real(r8), intent(in) :: vgl                  ! dry deposition velocity in m/s (land)
      real(r8), intent(in) :: vgo                  ! dry deposition velocity in m/s (ocean)
      real(r8), intent(in) :: vgsi                 ! dry deposition velocity in m/s (snow/ice)

! Output arguments:
      real(r8), intent(out) ::  vg(pcols) ! dry deposition velocity in m/s

! Local variables:

      integer i
      real(r8) a


      do i = 1, ncol
         vg(i) = landfrac(i)*vgl + ocnfrac(i)*vgo + icefrac(i)*vgsi
!         if (ioro(i).eq.0) then
!            vg(i) = vgo
!         else if (ioro(i).eq.1) then
!            vg(i) = vgl
!         else
!            vg(i) = vgsi
!         endif
      end do

      return
      end subroutine setdvel

!##############################################################################

      subroutine ddflux( ncol, vg, q, p, tv, flux )

! Compute surface flux due to dry deposition processes.


      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vg(pcols)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    p(pcols)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    tv(pcols)  ! midpoint virtual temperature in surface layer (K)

! Output arguments:

      real(r8), intent(out) ::    flux(pcols) ! flux due to dry deposition in kg/m^s/sec

! Local variables:

      integer i

      do i = 1, ncol
         flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
      end do

      return
      end subroutine ddflux

!##############################################################################
end module drydep_mod
