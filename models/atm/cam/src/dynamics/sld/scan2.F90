#include <misc.h>
#include <params.h>

subroutine scan2 (ztodt, cwava, etamid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Second gaussian latitude scan, converts from spectral coefficients to 
! grid point values, from poles to equator, with read/calculate/write cycle.
! 
! Method: 
! The latitude pair loop in this routine is multitasked.
!
! The grid point values of ps, t, u, v, z (vorticity), and d (divergence)
! are calculated and stored for each latitude from the spectral coefficients.
! In addition, the pressure-surface corrections to the horizontal diffusion
! are applied and the global integrals of the constituent fields are 
! computed for the mass fixer.
!
! Author: 
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: scan2.F90,v 1.11.2.6 2003/12/18 16:22:34 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use comslt
   use prognostics
   use comspe, only: maxm
   use rgrid
   use mpishorthand
   use physconst, only: cpair
   use history, only: write_inithist, inithist
   use time_manager, only: get_step_size, get_curr_date, get_nstep
   use infnan
#ifdef CRM
   use runtime_opts, only: crminitsave, crmsavechunks
#endif
#if ( defined QVORTDAMP )
   use runtime_opts, only: qvortdampfac, qvortnstepspinup, qvort_dailymean_interference
#endif
#ifdef QVORTDAMP
    use vorticity_anncycle, only: accumulate_daily_mean_vz, ref_vorticity_int
#endif

!-----------------------------------------------------------------------
   implicit none
!------------------------------Commons----------------------------------
#include <comqfl.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                ! twice the timestep unless nstep = 0
   real(r8), intent(in) :: cwava(plat)          ! weight applied to global integrals
   real(r8), intent(in) :: etamid(plev)         ! vertical coords at midpoints 
!
!---------------------------Local workspace-----------------------------
!
   real(r8) engy1         ! component of global energy integral (for time step n)
   real(r8) engy2         ! component of global energy integral (for time step n+1)
   real(r8) engy2a        ! component of global energy integral (for time step n+1)
   real(r8) engy2b        ! component of global energy integral (for time step n+1)
   real(r8) difft         ! component of global delta-temp integral ( (n+1) - n )
   real(r8) diffta        ! component of global delta-temp integral ( (n+1) - n )
   real(r8) difftb        ! component of global delta-temp integral ( (n+1) - n )
   real(r8) hw2a(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "A" portion of the hybrid grid)
   real(r8) hw2b(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "B" portion of the hybrid grid)
   real(r8) hw3a(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "A" portion of the hybrid grid)
   real(r8) hw3b(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "B" portion of the hybrid grid)
   real(r8) hwxa(pcnst,4)
   real(r8) hwxb(pcnst,4)
   real(r8) engy2alat(plat)     ! lat contribution to total energy integral
   real(r8) engy2blat(plat)     ! lat contribution to total energy integral
   real(r8) difftalat(plat)     ! lat contribution to delta-temperature integral
   real(r8) difftblat(plat)     ! lat contribution to delta-temperature integral
   real(r8) hw2al(pcnst,plat)   ! |------------------------------------
   real(r8) hw2bl(pcnst,plat)   ! |
   real(r8) hw3al(pcnst,plat)   ! | latitudinal contributions to the
   real(r8) hw3bl(pcnst,plat)   ! | components of global mass integrals
   real(r8) hwxal(pcnst,4,plat) ! |
   real(r8) hwxbl(pcnst,4,plat) ! |-----------------------------------
!                                
! Symmetric fourier coefficient arrays for all variables transformed 
! from spherical harmonics (see subroutine grcalc)
!                                
   real(r8) grdpss(2*maxm,plat/2)      ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grpss(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grpls(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grpms(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*H(n,m)

   real(r8) grds(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruhs(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvhs(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grths(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grus(2*maxm,plev,plat/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvs(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
#if ( defined QVORTDAMP )
   real(r8) auxgruhs(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrvhs(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrus(2*maxm,plev,plat/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrvs(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
#endif
   real(r8) grts(2*maxm,plev,plat/2)   ! sum(n) of t(n,m)*P(n,m)
   real(r8) grqs(2*maxm,plev,plat/2)
   real(r8) grtms(2*maxm,plev,plat/2)
   real(r8) grtls(2*maxm,plev,plat/2)
   real(r8) grqms(2*maxm,plev,plat/2)
   real(r8) grqls(2*maxm,plev,plat/2)
!
! Antisymmetric fourier coefficient arrays for all variables transformed
! from spherical harmonics (see grcalc)
!
   real(r8) grdpsa(2*maxm,plat/2)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grpsa(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grpla(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grpma(2*maxm,plat/2)       ! sum(n) of lnps(n,m)*H(n,m)
   real(r8) grda(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruha(2*maxm,plev,plat/2)  ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvha(2*maxm,plev,plat/2)  ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grtha(2*maxm,plev,plat/2)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grua(2*maxm,plev,plat/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grva(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
#if ( defined QVORTDAMP )
   real(r8) auxgruha(2*maxm,plev,plat/2)  ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrvha(2*maxm,plev,plat/2)  ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrua(2*maxm,plev,plat/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) auxgrva(2*maxm,plev,plat/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
#endif
   real(r8) grta(2*maxm,plev,plat/2)   ! sum(n) of t(n,m)*P(n,m)
   real(r8) grqa(2*maxm,plev,plat/2)
   real(r8) grtma(2*maxm,plev,plat/2)
   real(r8) grtla(2*maxm,plev,plat/2)
   real(r8) grqma(2*maxm,plev,plat/2)
   real(r8) grqla(2*maxm,plev,plat/2)
   real(r8) residual                           ! residual energy integral
   real(r8) beta                               ! energy fixer coefficient
!
! FFT buffers
!
   real(r8), allocatable:: fftbuf_in(:,:,:,:)  ! fftbuf_in(2*maxm,11,plevp,plat) 
   real(r8), allocatable:: fftbuf_out(:,:,:,:) ! fftbuf_out(plond,11,plevp,beglat:endlat) 
!
   integer m,n, irow                           ! indices
   integer lat,j                               ! latitude indices
   integer endi                                ! index
   integer :: yr, mon, day    ! year, month, and day components of a date
   integer :: ncsec           ! current time of day [seconds]
   logical l_write_inithist
   integer  :: dtime               !timestep size [seconds]

   integer :: nstep

   nstep = get_nstep()
!
!-----------------------------------------------------------------------
!
#if ( defined SPMD )
#if ( defined QVORTDAMP )
   allocate(fftbuf_in(2*maxm,13,plevp,plat))
#else
   allocate(fftbuf_in(2*maxm,11,plevp,plat))
#endif
#else
   allocate(fftbuf_in(1,1,1,1))
#endif

#if ( defined QVORTDAMP )
   allocate(fftbuf_out(plond,13,plevp,beglat:endlat))
#else
   allocate(fftbuf_out(plond,11,plevp,beglat:endlat))
#endif
!
fftbuf_in(:,:,:,:) = inf
fftbuf_out(:,:,:,:) = inf
grdpss(:,:)       = inf
grpss(:,:)        = inf
grpls(:,:)        = inf
grpms(:,:)        = inf
grds(:,:,:)    = inf
gruhs(:,:,:)   = inf
grvhs(:,:,:)   = inf
grths(:,:,:)   = inf
grus(:,:,:)    = inf
grvs(:,:,:)    = inf
#if ( defined QVORTDAMP )
auxgruhs(:,:,:)   = inf
auxgrvhs(:,:,:)   = inf
auxgrus(:,:,:)    = inf
auxgrvs(:,:,:)    = inf
#endif

grts(:,:,:)    = inf
grqs(:,:,:) = inf
grtms(:,:,:) = inf
grtls(:,:,:) = inf
grqms(:,:,:) = inf
grqls(:,:,:) = inf
grdpsa(:,:)       = inf
grpsa(:,:)        = inf
grpla(:,:)        = inf
grpma(:,:)        = inf
grda(:,:,:)    = inf
gruha(:,:,:)   = inf
grvha(:,:,:)   = inf
grtha(:,:,:)   = inf
grua(:,:,:)    = inf
grva(:,:,:)    = inf
#if ( defined QVORTDAMP )
auxgruha(:,:,:)   = inf
auxgrvha(:,:,:)   = inf
auxgrua(:,:,:)    = inf
auxgrva(:,:,:)    = inf
#endif
grta(:,:,:)    = inf
grqa(:,:,:) = inf
grtma(:,:,:) = inf
grtla(:,:,:) = inf
grqma(:,:,:) = inf
grqla(:,:,:) = inf
   call t_startf ('grcalc')

#ifdef QVORTDAMP
  ! Update the internal daily mean vorticity
  call accumulate_daily_mean_vz
  ! Interpoalte the offline vorticity annual cycle to current model time
  call ref_vorticity_int
#endif

#if ( defined SPMD )

!$OMP PARALLEL DO PRIVATE (J)

   do j=1,plat/2
      call grcalcs (j, ztodt, grts(1,1,j), grqs(1,1,j), grths(1,1,j), &
                    grds(1,1,j), grus(1,1,j), gruhs(1,1,j), grvs(1,1,j), grvhs(1,1,j), &
#if ( defined QVORTDAMP )
           auxgrus(1,1,j) ,auxgruhs(1,1,j),auxgrvs(1,1,j) ,auxgrvhs(1,1,j), qvortdampfac, qvortnstepspinup,  nstep, qvort_dailymean_interference, &
#endif
                    grpss(1,j), grdpss(1,j), grpms(1,j), grpls(1,j), grtms(1,1,j), &
                    grtls(1,1,j), grqms(1,1,j), grqls(1,1,j))

      call grcalca (j, ztodt, grta(1,1,j), grqa(1,1,j), grtha(1,1,j), &
                    grda(1,1,j), grua(1,1,j), gruha(1,1,j), grva(1,1,j), grvha(1,1,j), &
#if ( defined QVORTDAMP )
                    auxgrua(1,1,j), auxgruha(1,1,j), auxgrva(1,1,j), auxgrvha(1,1,j), &
                    qvortdampfac, qvortnstepspinup, nstep, qvort_dailymean_interference, &
#endif
                    grpsa(1,j), grdpsa(1,j), grpma(1,j), grpla(1,j), grtma(1,1,j), &
                    grtla(1,1,j), grqma(1,1,j), grqla(1,1,j))
   end do

#else

!$OMP PARALLEL DO PRIVATE (LAT, J)

   do lat=beglat,endlat
      if (lat > plat/2) then
         j = plat - lat + 1
         call grcalcs (j, ztodt, grts(1,1,j), grqs(1,1,j), grths(1,1,j), &
                       grds(1,1,j), grus(1,1,j), gruhs(1,1,j), grvs(1,1,j), grvhs(1,1,j), &
#if ( defined QVORTDAMP )
					auxgrus(1,1,j),auxgruhs(1,1,j),auxgrvs(1,1,j) ,auxgrvhs(1,1,j), qvortdampfac, qvortnstepspinup, nstep, qvort_dailymean_interference, &
#endif
                       grpss(1,j), grdpss(1,j), grpms(1,j), grpls(1,j), grtms(1,1,j), &
                       grtls(1,1,j), grqms(1,1,j), grqls(1,1,j))

      else
         j = lat
         call grcalca (j, ztodt, grta(1,1,j), grqa(1,1,j), grtha(1,1,j), &
                       grda(1,1,j), grua(1,1,j), gruha(1,1,j), grva(1,1,j), grvha(1,1,j), &
#if ( defined QVORTDAMP )
                    auxgrua(1,1,j), auxgruha(1,1,j), auxgrva(1,1,j), auxgrvha(1,1,j), &
                    qvortdampfac, qvortnstepspinup, nstep, qvort_dailymean_interference, &
#endif
                   grpsa(1,j), grdpsa(1,j), grpma(1,j), grpla(1,j), grtma(1,1,j), &
                       grtla(1,1,j), grqma(1,1,j), grqla(1,1,j))
      end if
   end do

#endif

   call t_stopf ('grcalc')

   call t_startf('spegrd')

!$OMP PARALLEL DO PRIVATE (LAT, IROW)

   do lat=1,plat
      irow = lat
      if (lat > plat/2) irow = plat - lat + 1
#if ( defined SPMD )
      call spegrd_bft (lat, &
                       grts(1,1,irow), grqs(1,1,irow), grths(1,1,irow), &
                       grds(1,1,irow), grus(1,1,irow), gruhs(1,1,irow), grvs(1,1,irow), grvhs(1,1,irow), &
                       grpss(1,irow), grdpss(1,irow), grpms(1,irow), grpls(1,irow), grtms(1,1,irow), &
                       grtls(1,1,irow), grqms(1,1,irow), grqls(1,1,irow), grta(1,1,irow), grqa(1,1,irow), &
                       grtha(1,1,irow), grda(1,1,irow), grua(1,1,irow), gruha(1,1,irow), grva(1,1,irow), &
                       grvha(1,1,irow), grpsa(1,irow), grdpsa(1,irow), grpma(1,irow), grpla(1,irow), &
#if ( defined QVORTDAMP )
                       auxgrua(1,1,irow), auxgruha(1,1,irow), auxgrva(1,1,irow), auxgrvha(1,1,irow), &
					   auxgrus(1,1,irow),auxgruhs(1,1,irow),auxgrvs(1,1,irow),auxgrvhs(1,1,irow), &
#endif
                       
                       grtma(1,1,irow), grtla(1,1,irow), grqma(1,1,irow), grqla(1,1,irow), fftbuf_in(1,1,1,lat))
#else
      call spegrd_bft (lat, &
                       grts(1,1,irow), grqs(1,1,irow), grths(1,1,irow), &
                       grds(1,1,irow), grus(1,1,irow), gruhs(1,1,irow), grvs(1,1,irow), grvhs(1,1,irow), &
                       grpss(1,irow), grdpss(1,irow), grpms(1,irow), grpls(1,irow), grtms(1,1,irow), &
                       grtls(1,1,irow), grqms(1,1,irow), grqls(1,1,irow), grta(1,1,irow), grqa(1,1,irow), &
                       grtha(1,1,irow), grda(1,1,irow), grua(1,1,irow), gruha(1,1,irow), grva(1,1,irow), &
                       grvha(1,1,irow), grpsa(1,irow), grdpsa(1,irow), grpma(1,irow), grpla(1,irow), &
#if ( defined QVORTDAMP )
                       auxgrua(1,1,irow), auxgruha(1,1,irow), auxgrva(1,1,irow), auxgrvha(1,1,irow), &
					   auxgrus(1,1,irow),auxgruhs(1,1,irow),auxgrvs(1,1,irow),auxgrvhs(1,1,irow), &
#endif
                       grtma(1,1,irow), grtla(1,1,irow), grqma(1,1,irow), grqla(1,1,irow), fftbuf_out(1,1,1,lat))
#endif
   end do

   call spegrd_ift ( fftbuf_in, fftbuf_out )
                   
!$OMP PARALLEL DO PRIVATE (LAT, J)

   do lat=beglat,endlat
      j = j1 - 1 + lat
      call spegrd_aft (ztodt, lat, nlon(lat), cwava(lat), qfcst(1,1,1,lat), q3(i1,1,1,j,n3), &
                       etamid, ps(1,lat,n3m1), u3(i1,1,j,n3m1), v3(i1,1,j,n3m1), t3(i1,1,j,n3m1), &
#if ( defined QVORTDAMP )
                       u3aux(i1,1,j,n3m1), v3aux(i1,1,j,n3m1), &
#endif
                       div(1,1,lat,n3m1), hw2al(1,lat), hw2bl(1,lat), hw3al(1,lat), hw3bl(1,lat), &
                       hwxal(1,1,lat), hwxbl(1,1,lat), dps(1,lat), &
                       dpsl(1,lat), dpsm(1,lat), tl(1,1,lat), tm(1,1,lat), ql(1,1,lat), &
                       qm(1,1,lat), t3(i1,1,j,n3) ,engy2alat(lat), engy2blat(lat), &
                       difftalat(lat), difftblat(lat), phis(1,lat), fftbuf_out(1,1,1,lat))
   end do

   call t_stopf('spegrd')
!
   deallocate(fftbuf_in)
   deallocate(fftbuf_out)

#ifdef SPMD
#ifdef TIMING_BARRIERS
   call t_startf ('sync_realloc5')
   call mpibarrier (mpicom)
   call t_stopf ('sync_realloc5')
#endif
   call t_startf('realloc5')
   call realloc5 (hw2al   ,hw2bl   ,hw3al   ,hw3bl   ,tmass    , &
                  hw1lat  ,hwxal   ,hwxbl   ,engy1lat,engy2alat, &
                  engy2blat, difftalat, difftblat)
   call t_stopf('realloc5')
#endif

!
! Accumulate and normalize global integrals for mass fixer (dry mass of
! atmosphere is held constant).
!
   call t_startf ('scan2_single')
   tmassf = 0.
   do lat=1,plat
      tmassf = tmassf + tmass(lat)
   end do
   tmassf = tmassf*.5
!
! Initialize moisture, mass, energy, and temperature integrals
!
   hw1(1) = 0.
   engy1  = 0.
   engy2a = 0.
   engy2b = 0.
   diffta = 0.
   difftb = 0.
   do m=1,pcnst
      hw2a(m) = 0.
      hw2b(m) = 0.
      hw3a(m) = 0.
      hw3b(m) = 0.
      do n=1,4
         hwxa(m,n) = 0.
         hwxb(m,n) = 0.
      end do
   end do
!
! Sum water and energy integrals over latitudes
!
   do lat=1,plat
      engy1   = engy1   + engy1lat (lat)
      engy2a  = engy2a  + engy2alat(lat)
      engy2b  = engy2b  + engy2blat(lat)
      diffta  = diffta  + difftalat(lat)
      difftb  = difftb  + difftblat(lat)
      hw1(1)  = hw1(1)  + hw1lat(1,lat)
      hw2a(1) = hw2a(1) + hw2al(1,lat)
      hw2b(1) = hw2b(1) + hw2bl(1,lat)
      hw3a(1) = hw3a(1) + hw3al(1,lat)
      hw3b(1) = hw3b(1) + hw3bl(1,lat)
   end do
!
! Compute atmospheric mass fixer coefficient
!
   qmassf     = hw1(1)
   if (adiabatic .or. ideal_phys) then
      fixmas = tmass0/tmassf
   else
      fixmas = (tmass0 + qmassf)/tmassf
   end if
!
! Compute alpha for water ONLY
!
   hw2(1)    = hw2a(1) + fixmas*hw2b(1)
   hw3(1)    = hw3a(1) + fixmas*hw3b(1)
   if(hw3(1) .ne. 0.) then
      alpha(1)  = ( hw1(1) - hw2(1) )/hw3(1)
   else
      alpha(1)  = 1.
   endif
!
! Compute beta for energy
!
   engy2    = engy2a + fixmas*engy2b
   difft    = diffta + fixmas*difftb
   residual = (engy2 - engy1)/ztodt
   if(difft .ne. 0.) then
     beta = -residual*ztodt/(cpair*difft)
   else
     beta = 0.
   endif
!SR Add output here to check whether this is happening.
 ! call outfld('BETA', ,plond   ,lat     )  
 if (masterproc) then
  write (6,*) 'SR: beta = ', beta
  write (6,*) 'SR: residual = ', residual
  write (6,*) 'SR: fixmas = ', fixmas
  write (6,*) 'SR: engy1 = ', engy1
  write (6,*) 'SR: engy2 = ', engy2
  write (6,*) 'SR: difft = ', difft
 endif

  
!!   write(6,125) residual,beta
!!125 format('      resid, beta      = ',25x,2f25.15)
!
! Compute alpha for non-water constituents
!
   do m = 2,pcnst
      hw1(m) = 0.
      do lat=1,plat
         hw1(m) = hw1(m) + hw1lat(m,lat)
      end do
      do n = 1,4
         do lat=1,plat
            hwxa(m,n) = hwxa(m,n) + hwxal(m,n,lat)
            hwxb(m,n) = hwxb(m,n) + hwxbl(m,n,lat)
         end do
      end do
      hw2a(m) = hwxa(m,1) - alpha(1)*hwxa(m,2)
      hw2b(m) = hwxb(m,1) - alpha(1)*hwxb(m,2)
      hw3a(m) = hwxa(m,3) - alpha(1)*hwxa(m,4)
      hw3b(m) = hwxb(m,3) - alpha(1)*hwxb(m,4)
      hw2 (m) = hw2a(m) + fixmas*hw2b(m)
      hw3 (m) = hw3a(m) + fixmas*hw3b(m)
      if(hw3(m) .ne. 0.) then
         alpha(m)  = ( hw1(m) - hw2(m) )/hw3(m)
      else
         alpha(m)  = 1.
      end if
   end do

   do j=beglatex,endlatex
      endi = nlonex(j) + i1 - 1
      do m=pcnst+1,pcnst+pnats
         q3(i1:endi,:,m,j,n3m1) = q3(i1:endi,:,m,j,n3)
      end do
   end do

   call t_stopf ('scan2_single')

   call t_startf ('tfilt_massfix')

!$OMP PARALLEL DO PRIVATE (LAT,J)

   do lat=beglat,endlat
      j = j1 - 1 + lat
      call tfilt_massfix (ztodt, lat, u3(i1,1,j,n3m1), v3(i1,1,j,n3m1), t3(i1,1,j,n3m1), &
                          q3(i1,1,1,j,n3), q3(i1,1,1,j,n3m1), ps(1,lat,n3m1), cwava(lat), alpha, &
                          etamid, qfcst(1,1,1,lat), div(1,1,lat,n3m1), phis(1,lat), omga(1,1,lat), &
                          dpsl(1,lat), dpsm(1,lat), nlon(lat), t3(i1,1,j,n3),beta)
   end do

   call t_stopf ('tfilt_massfix')
!
! Dump Initial file
!
   dtime = get_step_size()
   call get_curr_date(yr, mon, day, ncsec, offset=dtime)
      
! Write initial file if requested.

   l_write_inithist = .false.
   if (inithist == '6-HOURLY'  ) then
      if (mod(ncsec,21600) == 0 ) l_write_inithist = .true.
   end if
   if (ncsec == 0                               .and. inithist == 'DAILY'  ) l_write_inithist = .true.
   if (ncsec == 0 .and. day == 1                .and. inithist == 'MONTHLY') l_write_inithist = .true.
   if (ncsec == 0 .and. day == 1 .and. mon == 1 .and. inithist == 'YEARLY' ) l_write_inithist = .true.
#ifdef CRM
   if (l_write_inithist) call write_inithist(crminitsave,crmsavechunks)
#else
   if (l_write_inithist) call write_inithist
#endif
!
! Shift time pointers
!
   call shift_time_indices ()

   return
end subroutine scan2

#ifdef SPMD
subroutine realloc5 (hw2al   ,hw2bl   ,hw3al   ,hw3bl   ,tmass    , &
                     hw1lat  ,hwxal   ,hwxbl   ,engy1lat,engy2alat, &
                     engy2blat,difftalat,difftblat      )
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for slt variables.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: scan2.F90,v 1.11.2.6 2003/12/18 16:22:34 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use spmd_dyn
  use spmd_utils, only: pair, ceil2
  use prognostics
  use mpishorthand
  implicit none
#include <comsta.h>

!---------------------------------Parameters----------------------------------
!
  integer, parameter :: msgtype = 5000      ! message passing id
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: hw2al (pcnst,plat)   ! -
  real(r8), intent(in)   :: hw2bl (pcnst,plat)   !  | lat contributions to components
  real(r8), intent(in)   :: hw3al (pcnst,plat)   !  | of tracer global mass integrals 
  real(r8), intent(in)   :: hw3bl (pcnst,plat)   ! -
  real(r8), intent(in)   :: tmass (plat)         !    global atmospheric mass integral
  real(r8), intent(in)   :: hw1lat(pcnst,plat)   ! -
  real(r8), intent(in)   :: hwxal (pcnst,4,plat) !  | lat contributions to components
  real(r8), intent(in)   :: hwxbl (pcnst,4,plat) !  | of tracer global mass integrals 
!                                                ! -
  real(r8), intent(in)   :: engy1lat (plat)      ! lat contribution to total energy (n)
  real(r8), intent(in)   :: engy2alat(plat)      ! lat contribution to total energy (n+1)
  real(r8), intent(in)   :: engy2blat(plat)      ! lat contribution to total energy (n+1)
  real(r8), intent(in)   :: difftalat(plat)      ! lat contribution to delta-T integral
  real(r8), intent(in)   :: difftblat(plat)      ! lat contribution to delta-T integral
!
!---------------------------Local workspace-----------------------------
!
  integer len
  integer procid                                 ! Processor id
  integer stat(MPI_STATUS_SIZE)
  integer bpos
  integer procj,maxcount
  integer len_p,beglat_p,numlats_p
!
!-----------------------------------------------------------------------
!
! gather global data
!
  len = numlats*pcnst
  do procj=1,ceil2(npes)-1
     procid = pair(npes,procj,iam)
     if (procid.ge.0) then
        bpos = 0
        call mpipack (len               ,1      ,mpiint,buf1,bsiz,bpos,mpicom)
        call mpipack (beglat            ,1      ,mpiint,buf1,bsiz,bpos,mpicom)
        call mpipack (numlats           ,1      ,mpiint,buf1,bsiz,bpos,mpicom)

        call mpipack (tmass    (beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (engy1lat (beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (engy2alat(beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (engy2blat(beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (difftalat(beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (difftblat(beglat) ,numlats,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hw1lat(1  ,beglat),len    ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hw2al (1  ,beglat),len    ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hw2bl (1  ,beglat),len    ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hw3al (1  ,beglat),len    ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hw3bl (1  ,beglat),len    ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hwxal (1,1,beglat),len*4  ,mpir8 ,buf1,bsiz,bpos,mpicom)
        call mpipack (hwxbl (1,1,beglat),len*4  ,mpir8 ,buf1,bsiz,bpos,mpicom)

        call mpisendrecv (buf1,bpos,mpipk,procid,msgtype, &
                          buf2,bsiz,mpipk,procid,msgtype,mpicom)

        bpos = 0
        call mpiunpack (buf2,bsiz,bpos,len_p               ,1        ,mpiint,mpicom)
        call mpiunpack (buf2,bsiz,bpos,beglat_p            ,1        ,mpiint,mpicom)
        call mpiunpack (buf2,bsiz,bpos,numlats_p           ,1        ,mpiint,mpicom)

        call mpiunpack (buf2,bsiz,bpos,tmass    (beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,engy1lat (beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,engy2alat(beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,engy2blat(beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,difftalat(beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,difftblat(beglat_p) ,numlats_p,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hw1lat(1  ,beglat_p),len_p    ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hw2al (1  ,beglat_p),len_p    ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hw2bl (1  ,beglat_p),len_p    ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hw3al (1  ,beglat_p),len_p    ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hw3bl (1  ,beglat_p),len_p    ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hwxal (1,1,beglat_p),len_p*4  ,mpir8 ,mpicom)
        call mpiunpack (buf2,bsiz,bpos,hwxbl (1,1,beglat_p),len_p*4  ,mpir8 ,mpicom)
     end if
!JR         call mpibarrier(mpicom)
  end do

  return
end subroutine realloc5
#endif
