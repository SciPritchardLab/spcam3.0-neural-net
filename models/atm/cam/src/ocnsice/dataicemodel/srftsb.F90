#include <misc.h>
#include <params.h>

subroutine srftsb(isrfty  ,indx    ,npts    ,fnt     ,dfntdt  , &
                  snowh   ,tsbsf   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute surface and subsurface temperatures over sea-ice surfaces.
!
! Method: 
! Sea ice temperatures are specified in 'plevmx' layers of fixed
! thickness and thermal properties.  The forecast temperatures are
! determined from a backward/implicit diffusion calculation using
! linearized sensible/latent heat fluxes. The bottom ocean temperature
! is fixed at -2C, allowing heat flux exchange with underlying ocean.
! 
! Sub-surface layers are indexed 1 at the surface, increasing downwards
! to plevmx.  Layers have mid-points and interfaces between layers.
!
! Temperatures are defined at mid-points, while fluxes between layers
! and the top/bottom media are defined at layer interfaces.
! 
! Author: B. Briegleb, CCM2
! 
!-----------------------------------------------------------------------
!
! $Id: srftsb.F90,v 1.2.2.1 2002/06/15 13:48:50 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use comsrf, only: plevmx
  use time_manager, only: get_step_size
  implicit none

#include <comtsc.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: isrfty(pcols)          ! Surface type index (1..7)
  integer , intent(in) :: indx(pcols)            ! Index of points to be computed
  integer , intent(in) :: npts                   ! Number of points to be computed
  real(r8), intent(in) :: fnt(pcols)             ! Top surface/atmosphere net energy flux
  real(r8), intent(in) :: dfntdt(pcols)          ! ts partial derivative of net sfc flux
  real(r8), intent(in) :: snowh(pcols)           ! Snow depth (liquid water equivalent)
  real(r8), intent(inout) :: tsbsf(pcols,plevmx) ! Initial/final surface/sub-surface tmps
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer, parameter :: psrfty = 7  ! Number of surface types
  integer i,ii              ! Column indices
  integer ier               ! Error return on tridiagonal solver
  integer jndx              ! Surface type index
  integer n                 ! Sub-surface layer index
! 
  real(r8) cmass(pcols,plevmx)  ! Specific heat of soil (J/kg/K)
  real(r8) rho(pcols,plevmx)    ! Mass densty of sub-sfc mat (kg/m3)
  real(r8) tk(pcols,plevmx)     ! Thermal conductivity (watts/m/K)
  real(r8) z(pcols,0:plevmx)    ! Interface geometrical depth (m)
! 
  real(r8) diag(pcols,plevmx)   ! Diagonal matrix elements
  real(r8) htsrc(pcols,plevmx)  ! External heat source (W/m3)
  real(r8) rhs(pcols,plevmx)    ! Rhs of tri-diagonal matrix equation
  real(r8) sbdiag(pcols,plevmx) ! Sub-diagonal matrix elements
  real(r8) spdiag(pcols,plevmx) ! Super-diagonal matrix elements
  real(r8) tin(pcols,plevmx)    ! Initial sub-surface temperatures
  real(r8) ws(pcols*plevmx)     ! Working storage for mtdlss
!
  real(r8) cmsnow               ! Snow mass heat capacity
  real(r8) cmty(pcols)          ! Layer mass heat capacity
  real(r8) crt                  ! cmass*rho*rdtime
  real(r8) delz                 ! Layer thickness
  real(r8) delzmn               ! Thick from mid-point to lyr above mid-point
  real(r8) delzpl               ! Thick from mid-point to lyr below mid-point
  real(r8) fbt(pcols)           ! Ocean heat flux into sea-ice
  real(r8) fmns                 ! 1/(delz*delzmn)
  real(r8) fpls                 ! 1/(delz*delzpl)
  real(r8) msnow                ! Mass path of snow
  real(r8) msoil                ! Mass path of soil
  real(r8) rdtime               ! Inverse model time step
  real(r8) rhsnow               ! Snow mass density
  real(r8) rhty(pcols)          ! Layer mass density
  real(r8) rztop                ! 1/ztop
  real(r8) thck(pcols)          ! Layer thickness
  real(r8) tkbot                ! Bottom layer top interf thermal conduct
  real(r8) tkmns                ! Layer bottom interface thermal conduct
  real(r8) tkpls                ! Layer top interface thermal conductivity
  real(r8) tksnow               ! Snow thermal conducitivity
  real(r8) tktop                ! Top layer bottom interface thermal conduct
  real(r8) tkty(pcols)          ! Layer thermal conductivity
  real(r8) tmp                  ! crt - dfntdt(i)*rztop
  real(r8) zbot                 ! Bottom layer thickness
  real(r8) zm                   ! Present layer mid-point depth
  real(r8) zmmn                 ! Layer above mid-point depth
  real(r8) zmpl                 ! Layer below mid-point depth
  real(r8) zsnow                ! Snow geometric depth
  real(r8) ztop                 ! Top layer thickness
!
  logical scvr(pcols)       ! True if surface snow covered
!-----------------------------------------------------------------------

!--------------------------Data Statements------------------------------
!
! Specified (and invariant) thermal properties for surface types
!
  real(r8) cmtype(psrfty,plevmx)   ! Mass specific heat (J/kg/K)
  real(r8) rhtype(psrfty,plevmx)   ! Mass density (kg/m3)
  real(r8) thckly(psrfty,plevmx)   ! Layer thicknesses (m)
  real(r8) tktype(psrfty,plevmx)   ! Thermal conductivity (J/m/s)
  save cmtype,rhtype,thckly,tktype
!
  data cmtype /4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2, &
               4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2, &
               4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2, &
               4.20e3,2.07e3,2.07e3,1.04e3,7.20e2,5.60e2,4.16e2/
!
  data rhtype /1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3, &
               1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3, &
               1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3, &
               1.00e3,9.20e2,9.20e2,2.50e3,2.50e3,2.50e3,2.50e3/
!
  data thckly / 2., .500, .250, .050, .090, .080, .120, &
                5., .500, .500, .366, .390, .435, .492, &
               10., .500, .500,1.369,1.459,1.628,1.841, &
               33.,.500,8.500,6.990,7.450,8.310,9.400/
!
  data tktype /15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 , &
               15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 , &
               15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 , &
               15.0 ,2.200 ,2.200 ,1.408 ,1.104 ,1.071 ,1.019 /
!
! Properties of ice and air used to determine snow thermal properties
!
  real(r8) cmair     ! Mass specific heat of air
  real(r8) cmice     ! Mass specific heat of ice
  real(r8) frcair    ! Fraction of air assumed in mix of ice and air
  real(r8) rhair     ! Mass density of surface air
  real(r8) rhice     ! Mass density of ice
  real(r8) tkair     ! Thermal conductivity of air
  real(r8) tkice     ! Thermal conductivity of ice
  save cmair,cmice,frcair,rhair,rhice,tkair,tkice
!
  data cmair  /1.00e3/
  data cmice  /2.07e3/
  data frcair /0.90/
  data rhair  /1.25 /
  data rhice  /9.20e2/
  data tkair  /0.025/
  data tkice  /2.200/
!
!-----------------------------------------------------------------------
!
  rdtime = 1._r8/get_step_size()
!
! Calculate snow properties
!
  cmsnow = (1.-frcair)*cmice + frcair*cmair
  rhsnow = (1.-frcair)*rhice + frcair*rhair
  tksnow = (1.-frcair)*tkice + frcair*tkair
!
! No external heat source
!
  do n=1,plevmx
     do ii=1,npts
        i = indx(ii)
        htsrc(i,n) = 0.0
     end do
  end do
!
! Define logical for snow covered surfaces:
!
  do ii=1,npts
     i = indx(ii)
     scvr(i) = snowh(i)>0.
  end do
!
! Define thermal properities for each sub/surface layer, starting
! with the top layer
!
  do ii=1,npts
     i = indx(ii)
     jndx    = isrfty(i)
     thck(i) = thckly(jndx,1)
     cmty(i) = cmtype(jndx,1)
     rhty(i) = rhtype(jndx,1)
     tkty(i) = tktype(jndx,1)
  end do
!
!CDIR$ IVDEP
  do ii=1,npts
     i = indx(ii)
!
! Initialize fields for no snow cover
!
     z(i,0)     = 0.0
     z(i,1)     = thck(i)
     cmass(i,1) = cmty(i)
     rho(i,1)   = rhty(i)
     tk(i,1)    = tkty(i)
!
! Modify layer 1 fields for snow cover if present
!
     if (scvr(i)) then
!
! Snow equivlnt depth times snow liquid water depth gives the physical
! depth of snow for thermal conduction computation; snow is mixed
! uniformly by mass with the top surface layer
!
        zsnow     = snowh(i)*snwedp
        msnow     = rhsnow*zsnow
        msoil     = rhty(i)*thck(i)
        rho(i,1)  = (msnow*rhsnow + msoil*rhty(i))/(msnow+msoil)
        cmass(i,1)= (msnow*cmsnow + msoil*cmty(i))/(msnow+msoil)
        tk(i,1)   = (msnow*tksnow + msoil*tkty(i))/(msnow+msoil)
        z(i,1)    = (msnow+msoil) / rho(i,1)
     end if
!
  end do
!
! Set surface thermal properties for the lower sub/surface layers:
!
  do n=2,plevmx
     do ii=1,npts
        i = indx(ii)
        jndx       = isrfty(i)
        thck(i)    = thckly(jndx,n)
        cmass(i,n) = cmtype(jndx,n)
        rho(i,n)   = rhtype(jndx,n)
        tk(i,n)    = tktype(jndx,n)
     end do
     do ii=1,npts
        i = indx(ii)
        z(i,n)  = z(i,n-1) + thck(i)
     end do
  end do
!
! Define set of linear equations for temperature
!
  do n=1,plevmx
     do ii=1,npts
        i = indx(ii)
        tin(i,n) = tsbsf(i,n)
     end do
  end do
!
! If sea ice, compute heat flux from underlying ocean, assumed to be at
! the temperature of -2C
!
!CDIR$ IVDEP
  do ii=1,npts
     i = indx(ii)
     fbt(i) = 0.0
     if(isrfty(i)==2) then
        zbot    = 0.5*(z(i,plevmx) - z(i,plevmx-1))
        fbt(i)  = -tk(i,plevmx)*(271.16 - tin(i,plevmx))/zbot
     end if
  end do
!
! Set up linear equations
!
  do ii=1,npts
     i = indx(ii)
     sbdiag(i,1)      = 0.
     spdiag(i,plevmx) = 0.
  end do

  if (plevmx==1) then
!
! Single layer
!
     do ii=1,npts
        i = indx(ii)
        rztop     = 1./(z(i,1) - z(i,0))
        crt       = (cmass(i,1)*rho(i,1)*rdtime)
        diag(i,1) = crt - dfntdt(i)*rztop
        rhs(i,1)  = diag(i,1)*tin(i,1) + fnt(i)*rztop - fbt(i)*rztop + htsrc(i,1)
     end do
  else if (plevmx>1) then
!
! More than one layer: top layer first
!
     do ii=1,npts
        i = indx(ii)
        crt         = cmass(i,1)*rho(i,1)*rdtime
        ztop        = z(i,1) - z(i,0)
        rztop       = 1./ztop
        tktop       = 0.5*(tk(i,1) + tk(i,2))
        zmpl        = 0.5*(z(i,2) + z(i,1))
        zm          = 0.5*(z(i,1) + z(i,0))
        delzpl      = zmpl - zm
        fpls        = 1./(ztop*delzpl)
        tmp         = crt - dfntdt(i)*rztop
        diag(i,1)   = tmp + tktop*fpls
        spdiag(i,1) = -tktop*fpls
        rhs(i,1)    = tmp*tin(i,1) + fnt(i)*rztop + htsrc(i,1)
     end do
!
! Intermediate layers
!
     do n=2,plevmx-1
        do ii=1,npts
           i = indx(ii)
           crt         = cmass(i,n)*rho(i,n)*rdtime
           delz        = z(i,n) - z(i,n-1)
           zmpl        = 0.5*(z(i,n+1) + z(i,n))
           zm          = 0.5*(z(i,n)   + z(i,n-1))
           zmmn        = 0.5*(z(i,n-1) + z(i,n-2))
           delzpl      = zmpl - zm
           delzmn      = zm - zmmn
           fpls        = 1./(delz*delzpl)
           fmns        = 1./(delz*delzmn)
           tkpls       = 0.5*(tk(i,n+1)+tk(i,n))
           tkmns       = 0.5*(tk(i,n)+tk(i,n-1))
           sbdiag(i,n) = -tkmns*fmns
           diag(i,n)   = crt + (tkpls*fpls + tkmns*fmns)
           spdiag(i,n) = -tkpls*fpls
           rhs(i,n)    = crt*tin(i,n) + htsrc(i,n)
        end do
     end do
!
! Bottom layer
!
     do ii=1,npts
        i = indx(ii)
        crt        = cmass(i,plevmx)*rho(i,plevmx)*rdtime
        zbot       = z(i,plevmx) - z(i,plevmx-1)
        zm         = 0.5*(z(i,plevmx)   + z(i,plevmx-1))
        zmmn       = 0.5*(z(i,plevmx-1) + z(i,plevmx-2))
        delzmn     = zm - zmmn
        tkbot      = 0.5*(tk(i,plevmx-1) + tk(i,plevmx))
        fmns       = 1./(zbot*delzmn)
        sbdiag(i,plevmx) = -tkbot*fmns
        diag(i,plevmx)   = crt + (tkbot*fmns)
        rhs(i,plevmx)    = crt*tin(i,plevmx) - fbt(i)/zbot + htsrc(i,plevmx)
     end do
  end if
!
! For the linear equation ax = b,  a and b arrays are now set;
! solve for temperatures (x):
!
  if (plevmx==1) then
     do ii=1,npts
        i = indx(ii)
        tsbsf(i,1) = rhs(i,1)/diag(i,1)
     end do
  else
     call mtdlss(sbdiag  ,diag   ,spdiag   ,rhs     ,tsbsf , &
                 plevmx  ,pcols  ,npts     ,indx    ,ws    , &
                 pcols*plevmx    ,ier     )
     if(ier/=0) then
        write(6,*) 'SRFTSB: Error returned from mtdlss'
        call endrun
     end if
  end if
!
  return
end subroutine srftsb

