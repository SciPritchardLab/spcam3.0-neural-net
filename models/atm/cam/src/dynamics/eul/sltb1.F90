#include <misc.h>
#include <params.h>

subroutine sltb1(pmap    ,jcen    ,jgc     ,dt      ,ra      , &
                 iterdp  ,uxl     ,uxr     ,vxl     ,vxr     , &
                 wb      ,fxl     ,fxr     ,lam     ,phib    , &
                 dphib   ,sig     ,sigh    ,dsig    ,dsigh   , &
                 lbasdy  ,lbasdz  ,lbassd  ,lbasiy  ,kdpmpf  , &
                 kdpmph  ,lammp   ,phimp   ,sigmp   ,fbout   , &
                 u3      ,v3      ,qminus  ,n3m1    , &
                 nlon    ,nlonex  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Drive the slt algorithm on a given latitude slice in the extended
! data arrays using information from the entire latitudinal extent
! of the arrays.
! 
! Method: 
! Compute departure points and corresponding indices.
! Poleward of latitude phigs (radians), perform the computation in
! local geodesic coordinates.
! Equatorward of latitude phigs, perform the computation in global
! spherical coordinates
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: 
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use prognostics, only: ptimelevels
  use constituents,only: pcnst

  implicit none

#include <parslt.h>

!------------------------------Parameters-------------------------------
  real(r8), parameter :: phigs = 1.221730 ! cut-off latitude: about 70 degrees
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                ! longitude dimension
  integer , intent(in) :: nlonex(platd)       ! extended longitude dimension
  integer , intent(in) :: pmap                ! artificial vert grid dim.
  integer , intent(in) :: jcen                ! index of lat slice(extend)
  integer , intent(in) :: jgc                 ! index of lat slice (model)
  real(r8), intent(in) :: dt                  ! time step (seconds)
  real(r8), intent(in) :: ra                  ! 1./(radius of earth)
  integer , intent(in) :: iterdp              ! iteration count
  real(r8), intent(in) :: uxl(plond,plev,beglatex:endlatex) ! left  x-deriv of ub
  real(r8), intent(in) :: uxr(plond,plev,beglatex:endlatex) ! right x-deriv of ub
  real(r8), intent(in) :: vxl(plond,plev,beglatex:endlatex) ! left  x-deriv of vb
  real(r8), intent(in) :: vxr(plond,plev,beglatex:endlatex) ! right x-deriv of vb
  real(r8), intent(in) :: wb(plon,plevp)                    ! eta-dot
  real(r8), intent(in) :: fxl(plond,plev,  pcnst,beglatex:endlatex) ! left  fb x-deriv
  real(r8), intent(in) :: fxr(plond,plev,  pcnst,beglatex:endlatex) ! right fb x-deriv
  real(r8), intent(in) :: lam  (plond,platd)  ! long. coord of model grid
  real(r8), intent(in) :: phib (platd)        ! lat.  coord of model grid
  real(r8), intent(in) :: dphib(platd)        ! increment between lats.
  real(r8), intent(in) :: sig  (plev)         ! vertical full levels
  real(r8), intent(in) :: sigh (plevp)        ! vertical half levels
  real(r8), intent(in) :: dsig (plev)         ! inc. between full levs
  real(r8), intent(in) :: dsigh(plevp)        ! inc. between half levs
  real(r8), intent(in) :: lbasdy(4,2,platd)   ! lat deriv weights
  real(r8), intent(in) :: lbasdz(4,2,plev)    ! vert full level deriv wts
  real(r8), intent(in) :: lbassd(4,2,plevp)   ! vert half level deriv wts
  real(r8), intent(in) :: lbasiy(4,2,platd)   ! lat interp wts(lagrng)
  integer , intent(in) :: kdpmpf(pmap)        ! artificial vert grid index
  integer , intent(in) :: kdpmph(pmap)        ! artificial vert grid index
  real(r8), intent(inout) ::  u3(plond, plev, beglatex:endlatex,ptimelevels)     ! u wind vel
  real(r8), intent(inout) ::  v3(plond, plev, beglatex:endlatex,ptimelevels)     ! v wind vel
  real(r8), intent(inout) ::  qminus(plond, plev, pcnst, beglatex:endlatex) !moist
  integer , intent(inout) ::  n3m1               ! time indicies
  real(r8), intent(inout) ::  lammp(plon,plev)       ! long coord of mid-point
  real(r8), intent(inout) ::  phimp(plon,plev)       ! lat  coord of mid-point
  real(r8), intent(inout) ::  sigmp(plon,plev)       ! vert coord of mid-point
  real(r8), intent(out) :: fbout(plond,plev,pcnst)   ! advected constituents
!
!  pmap    Dimension of kdpmpX arrays
!  jcen    Latitude index in extended grid corresponding to lat slice
!          being forecasted.
!  jgc     Latitude index in model    grid corresponding to lat slice 
!          being forecasted.
!  dt      Time interval that parameterizes the parcel trajectory.
!  ra      Reciprocal of radius of earth.
!  iterdp  Number of iterations used for departure point calculation.
!  uxl     x-derivatives of u at the left  (west) edge of given interval
!  vxl     x-derivatives of v at the left  (west) edge of given interval
!  uxr     x-derivatives of u at the right (east) edge of given interval
!  vxr     x-derivatives of v at the right (east) edge of given interval
!  wb      z-velocity component (eta-dot).
!  fxl     x-derivatives at the left  edge of each interval containing 
!          the departure point.
!  fxr     x-derivatives at the right edge of each interval containing 
!          the departure point.
!  lam     Longitude values for the extended grid.
!  phib    Latitude  values for the extended grid.
!  dphib   Interval between latitudes in the extended grid.
!  sig     Hybrid eta values at the "full-index" levels.
!  sigh    Half-index eta-levels including sigh(i,1) = eta(1/2) = 0.0
!          and sigh(i,plev+1) = eta(plev+1/2) = 1.  Note that in general
!          sigh(i,k) .lt. sig(i,k)  where sig(i,k) is the hybrid value
!          at the k_th full-index level.
!  dsig    Interval lengths in full-index hybrid level grid.
!  dsigh   Interval lengths in half-index hybrid level grid.
!  lbasdy  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced latitude grid.
!  lbasdz  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (full levels).
!  lbassd  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (half levels).
!  lbasiy  Weights for Lagrange cubic interpolation on the unequally
!          spaced latitude grid.
!  kdpmpf  indices of artificial grid mapped into the full level grid
!  kdpmph  indices of artificial grid mapped into the half level grid
!  lammp   Longitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecasted.  On entry lammp
!          is an initial guess.
!  phimp   Latitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecasted.  On entry phimp
!          is an initial guess.
!  sigmp   Hybrid value at the trajectory midpoint for each gridpoint
!          in a vertical slice from the global grid.  On entry sigmp is
!          an initial guess.
!  fbout   Extended array only one latitude of which, however, is filled
!          with forecasted (transported) values.  This routine must be
!          called multiple times to fill the entire array.  This is
!          done to facilitate multi-tasking.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer m                            ! constituent index
  integer idp(plon,plev,4)             ! zonal      dep point index
  integer jdp(plon,plev)               ! meridional dep point index
  integer kdp(plon,plev)               ! vertical   dep point index
  real(r8) fhr(plon,plev,pcnst)        ! horizontal interpolants
  real(r8) lamdp(plon,plev)            ! zonal      departure pt. coord.
  real(r8) phidp(plon,plev)            ! meridional departure pt. coord.
  real(r8) sigdp(plon,plev)            ! vertical   departure pt. coord.
  real(r8) fhst(plon,plev,pcnst)       ! derivative at top of interval
  real(r8) fhsb(plon,plev,pcnst)       ! derivative at bot of interval
  real(r8) wst(plon,plevp)             ! w derivative at top of interval
  real(r8) wsb(plon,plevp)             ! w derivative at bot of interval
  real(r8) fint(plon,plev,ppdy,pcnst)  ! work space
  real(r8) fyb(plon,plev,pcnst)        ! work space
  real(r8) fyt(plon,plev,pcnst)        ! work space
  real(r8) fbout1(plond,plev,pcnst)    ! work space to please lf95 compiler
  logical locgeo                       ! flag indicating coordinate sys
!-----------------------------------------------------------------------
!
! Horizontal interpolation
!
  locgeo = abs(phib(jcen))>=phigs
!
  call sphdep(jcen    ,jgc     ,dt      ,ra      ,iterdp  ,                     &
              locgeo  ,u3(1,1,beglatex,n3m1)       ,uxl     ,uxr     ,lam     ,   &
              phib    ,lbasiy  ,lammp   ,phimp   ,lamdp   ,                     &
              phidp   ,idp     ,jdp     ,v3(1,1,beglatex,n3m1),                   &
              vxl     ,vxr     ,nlon    ,nlonex  )
!
! Interpolate scalar fields to the departure points.
!
  call hrintp(pcnst   ,pcnst   ,qminus(1,1,1,beglatex), fxl     ,fxr     , &
              lam     ,phib    ,dphib   ,lbasdy  ,lamdp   ,                    &
              phidp   ,idp     ,jdp     ,jcen    ,plimdr  ,                    &
              fint    ,fyb     ,fyt     ,fhr     ,nlon    ,                    &   
              nlonex  )
!
! Vertical interpolation.
! Compute vertical derivatives of vertical wind
!
  call cubzdr(nlon    ,plevp   ,wb      ,lbassd  ,wst     , &
              wsb     )
!
! Compute departure points and corresponding indices.
!
  call vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
              wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , &
              kdpmph  ,sigmp   ,sigdp   ,kdp     ,nlon    )
!
! Vertical derivatives of scalar fields.
! Loop over constituents.
!
  do m = 1,pcnst
     call cubzdr(nlon    ,plev    ,fhr(1,1,m), lbasdz  ,fhst(1,1,m), &
                 fhsb(1,1,m) )
  end do
  if( plimdr )then
     call limdz(fhr     ,dsig    ,fhst    ,fhsb    ,nlon    )
  end if
!
! Vertical interpolation of scalar fields.
!
  call herzin(plev    ,pcnst   ,fhr     ,fhst    ,fhsb    , &
              sig     ,dsig    ,sigdp   ,kdp     ,fbout1  , &
              nlon    )

  fbout(i1:nlon+i1-1,:,:) = fbout1(1:nlon,:,:)

  return
end subroutine sltb1

!============================================================================================

subroutine vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
                  wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , & 
                  kdpmph  ,sigmp   ,sigdp   ,kdp     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute vertical departure point and departure point index.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                ! longitude dimension
  integer , intent(in) :: pmap                ! dimension of artificial vert grid
  real(r8), intent(in) :: dt                  ! time step (seconds)
  integer , intent(in) :: iterdp              ! number of iterations
  real(r8), intent(in) :: wb (plon,plevp)     ! vertical velocity
  real(r8), intent(in) :: wst(plon,plevp)     ! z-derivative of wb at top of interval
  real(r8), intent(in) :: wsb(plon,plevp)     ! z-derivative of wb at bot of interval
  real(r8), intent(in) :: sig  (plev )        ! sigma values of model full levels
  real(r8), intent(in) :: sigh (plevp)        ! sigma values of model half levels
  real(r8), intent(in) :: dsigh(plevp)        ! increment between half levels
  integer , intent(in) :: kdpmpf(pmap)        ! artificial grid indices
  integer , intent(in) :: kdpmph(pmap)        ! artificial grid indices
  real(r8), intent(inout) :: sigmp(plon,plev) ! vert coords of traj mid-points
  real(r8), intent(out) :: sigdp(plon,plev)   ! vert coords of traj departure points
  integer , intent(out) :: kdp(plon,plev)     ! vertical departure point indices
!
!  pmap    Dimension of kdpmap arrays
!  dt      Time interval that parameterizes the parcel trajectory.
!  iterdp  Number of iterations used for departure point calculation.
!  wb      Vertical velocity component (sigma dot).
!  wst     z-derivs at the top edge of each interval contained in wb
!  wsb     z-derivs at the bot edge of each interval contained in wb
!  sig     Sigma values at the full-index levels.
!  sigh    Half-index sigma levels including sigh(1) = sigma(1/2) = 0.0
!          sigh(plev+1) = sigma(plev+1/2) = 1.0 .  Note that in general
!          sigh(k) .lt. sig(k)  where sig(k) is the sigma value at the
!          k_th full-index level.
!  dsigh   Increment in half-index sigma levels.
!  kdpmpf  Array of indices of the model full levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  kdpmph  Array of indices of the model half levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  sigmp   Sigma value at the trajectory midpoint for each gridpoint
!          in a vertical slice from the global grid.  On entry sigmp is
!          an initial guess.
!  sigdp   Sigma value at the trajectory endpoint for each gridpoint
!          in a vertical slice from the global grid.
!  kdp     Vertical index for each gridpoint.  This index points into a
!          vertical slice array whose vertical grid is given by sig.
!          E.g.,   sig(kdp(i,k)) .le. sigdp(i,k) .lt. sig(kdp(i,k)+1).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! |
  integer iter              ! |-- indices
  integer k                 ! |
  real(r8) wmp(plond,plev)  ! vert vel. at midpoint
!-----------------------------------------------------------------------
!
! Loop over departure point iterates.
!
  do iter = 1,iterdp
!
! Compute midpoint indices in half-index sigma-level arrays (use kdp
! as temporary storage).
!
     call kdpfnd(plevp   ,pmap    ,sigh    ,sigmp   ,kdpmph  , &
                 kdp     ,nlon    )
!
! Interpolate sigma dot field to trajectory midpoints using Hermite
! cubic interpolant.
!
     call herzin(plevp   ,1       ,wb      ,wst     ,wsb     , &
                 sigh    ,dsigh   ,sigmp   ,kdp     ,wmp     , &
                 nlon    )
!
! Update estimate of trajectory midpoint.
!
     do k = 1,plev
        do i = 1,nlon
           sigmp(i,k) = sig(k) - .5*dt*wmp(i,k)
        end do
     end do
!
! Restrict vertical midpoints to be between the top and bottom half-
! index sigma levels.
!
     call vdplim(plevp   ,sigh    ,sigmp   ,nlon)
  end do
!
! Compute trajectory endpoints.
!
  do k = 1,plev
     do i = 1,nlon
        sigdp(i,k) = sig(k) - dt*wmp(i,k)
     end do
  end do
!
! Restrict vertical departure points to be between the top and bottom
! full-index sigma levels.
!
  call vdplim(plev    ,sig     ,sigdp   ,nlon)
!
! Vertical indices for trajectory endpoints that point into full-index
! sigma level arrays.
!
  call kdpfnd(plev    ,pmap    ,sig     ,sigdp   ,kdpmpf  , &
              kdp     ,nlon    )
!
  return
end subroutine vrtdep


!============================================================================================

subroutine vdplim(pkdim   ,sig     ,sigdp   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Restrict vertical departure points to be between the top and bottom
! sigma levels of the "full-" or "half-" level grid
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: sltb1.F90,v 1.6.4.3 2002/06/15 13:47:52 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none

!-----------------------------------------------------------------------
  integer , intent(in)    :: nlon               ! longitude dimension
  integer , intent(in)    :: pkdim              ! vertical dimension
  real(r8), intent(in)    :: sig(pkdim)         ! vertical coordinate of model grid
  real(r8), intent(inout) :: sigdp(plon,plev)   ! vertical coords. of departure points.
! pkdim   Vertical dimension of "sig"
! sig     Sigma values at the "full" or "half" model levels
! sigdp   Sigma value at the trajectory endpoint or midpoint for each
!         gridpoint in a vertical slice from the global grid.  This
!         routine restricts those departure points to within the
!         model's vertical grid.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k                 ! index
!-----------------------------------------------------------------------
!
  do k=1,plev
     do i = 1,nlon
        if (sigdp(i,k) < sig(1)) then
           sigdp(i,k) = sig(1)
        end if
        if (sigdp(i,k) >= sig(pkdim)) then
           sigdp(i,k) = sig(pkdim)*(1. - 10.*epsilon(sigdp))
        end if
     end do
  end do

  return
end subroutine vdplim
