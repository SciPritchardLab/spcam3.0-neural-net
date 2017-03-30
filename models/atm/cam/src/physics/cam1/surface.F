      module surface

      use shr_kind_mod, only: r8 => shr_kind_r8

      real(r8)
     $  gravit       ! Gravitational acceleration
     $, rair         ! Gas constant for dry air

      real(r8) :: zero = 0.

      contains

c $Id: surface.F,v 1.1.6.1 2003/07/29 00:13:54 hender Exp $
      subroutine inisflx( xrair, xgravit, xrearth, xlam, xphi
     $,                   xgw )

c----------------------------------------------------------------------------
      use pmgrid
c----------------------------------------------------------------------------
      implicit none

c Input arguments:
      real(r8) 
     $  xrair             ! Gas constant for dry air
     $, xgravit           ! Gravitational acceleration
     $, xrearth           ! radius of earth (m)
     $, xlam(plon)        ! longitudes in grid (radians)
     $, xphi(plat)        ! Gauss latitudes in grid (radians)
     $, xgw(plat)         ! Gauss weights (sum to 2.0)
c----------------------------------------------------------------------------

c     Set constants for surface flux routines.
      rair = xrair
      gravit = xgravit

      return
      end subroutine inisflx

c##############################################################################


c##############################################################################

      end module surface
