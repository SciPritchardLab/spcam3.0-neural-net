	include 'grid.inc'
c-----------------------------------------------------------------
c       Commons for the dynamical framework
c --------------------------------------------------------------------
c       prognostic variables:

	real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
	real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
	real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
        real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! moist static energy
        real q   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! total water
	real qp	 (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! rain+snow

        common /com3dvel/ u, v, w
	common /com3dthm/ t, q, qp

        real u_save   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
        real v_save   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
        real w_save   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
        real t_save   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! moist static energy
        real q_save   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! total water
        real qp_save (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! rain+snow

        common /com3dvelsave/ u_save, v_save, w_save
        common /com3dthmsave/ t_save, q_save, qp_save


c	diagnostic variables:

	real p      (0:nx, (1-YES3D):ny, nzm)	! pressure
	real tabs   (nx, ny, nzm)	! temperature
	real qn     (nx, ny, nzm)	! cloud water+cloud ice 
        real qsat   (nx, ny, nzm)       ! saturation mixing ratio
        real omegan (nx, ny, nzm)       ! ice/water partition coef.
        real omegap (nx, ny, nzm)       ! snow+grouple/rain partition coef.
        real omegag (nx, ny, nzm)       ! snow/graupel partition coef.
        real tke    (nx, ny, nzm)       ! SGS tke

        real tk     (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddyviscosity
        real tkh    (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddy conductivity
        
	common /com3dtke/ tke, tk, tkh 
	common /comdiathm/ p, tabs, qn
	common /comdiasat/ qsat, omegan, omegap, omegag

c-----------------------------------------------------------------
c    time-tendencies for prognostic variables

	real dudt   (nxp1, ny, nzm, 3)
	real dvdt   (nx, nyp1, nzm, 3)
	real dwdt   (nx, ny, nz,  3)

        common /comdu/ dudt
        common /comdv/ dvdt
        common /comdw/ dwdt
c------------------------------------------------------------------
c       fluxes at the top and bottom of the domain:

	real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
        real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
        real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny)
	real precsfc(nx,ny) ! surface precip. rate
	real precssfc(nx,ny) ! surface ice precip. rate
                
        common /comfluxbt/ fluxbu, fluxbv, fluxbt, fluxbq,  
     &    fluxtu, fluxtv, fluxtt, fluxtq, fzero, precsfc, precssfc
c-----------------------------------------------------------------
c	profiles of horizontally averaged fields

	real t0 (nzm), q0 (nzm), qc0(nzm), qv0(nzm),tabs0(nzm),
     &       qi0(nzm),tv0(nzm), rel0(nzm), u0(nzm), v0(nzm), 
     &       tg0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), p0(nzm),
     &       tke0(nzm)

c       reference vertical profiles:
c 
        real   z(nz)	  ! height of the pressure levels above surface,m
        real   pres(nzm)  ! pressure,mb at scalar levels
        real   prespot(nzm)  ! (1000./pres)**R/cp
        real   bet(nzm)	  ! = ggr/tv0
        real   betp(nzm)  ! = ggr/p0
        real   gamaz(nzm) ! ggr/cp*z
        real   wsub(nz)   ! Large-scale subsidence velocity,m/s
        real   qtend(nzm) ! Large-scale tendency for total water
        real   ttend(nzm) ! Large-scale tendency for pot. temp.
        real   utend(nzm) ! Large-scale tendency for u
        real   vtend(nzm) ! Large-scale tendency for v


        common /comhor0/ t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0, 
     &          p0, tg0, qg0, ug0, vg0, z, pres,  bet, gamaz,
     &          wsub, qtend, ttend, utend, vtend, prespot, tke0, betp
c---------------------------------------------------------------------


        logical dodamping, doupperbound, docloud, doprecip, 
     $          dolongwave, doshortwave, dosgs, dosubsidence,
     $          docoriolis, dosurface, dolargescale, doradforcing,
     $		dosfcforcing, doradsimple, donudging_uv, donudging_tq,
     $		dosmagor, doscalar, doensemble
     
        common /comswitches/ dodamping, doupperbound, docloud, doprecip,
     $          dolongwave, doshortwave, dosgs, dosubsidence,
     $          docoriolis, dosurface, dolargescale,doradforcing,
     $		dosfcforcing, doradsimple, donudging_uv, donudging_tq,
     $		dosmagor, doscalar, doensemble
c----------------------------------------------------------------------
c       Vertical profiles of quantities sampled for statitistics purposes:

        real
     & qwle(nz), qwsb(nz),  qpwle(nz), qpwsb(nz), 
     & precflux(nz), qadv(nz), qdiff(nz), qpadv(nz), qpdiff(nz),
     & qpsrc(nz), qpfall(nz), qpevp(nz)
	common /comevol/ qwle, qwsb, qpwle, qpwsb, 
     &                   precflux, qadv,qdiff,qpadv,
     &                   qpdiff, qpsrc, qpfall, qpevp
c----------------------------------------------------------------------
c
c register functions:
c


        integer lenstr
        real esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
        real qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
	external lenstr
        external esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
        external qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm

!tBM*   threadlocal /com3dvel/, /com3dthm/
!tBM*   threadlocal /com3dtke/,/comdiathm/,/comdiasat/
!tBM*   threadlocal /comdu/,/comdv/,/comdw/
!tBM*   threadlocal /comfluxbt/
!tBM*   threadlocal /comhor0/
!tBM*   threadlocal /comswitches/
!tBM*   threadlocal /comevol/
!$OMP   threadprivate (/com3dvel/,/com3dthm/)
!$OMP   threadprivate (/com3dvelsave/,/com3dthmsave/)
!$OMP   threadprivate (/com3dtke/,/comdiathm/,/comdiasat/)
!$OMP   threadprivate (/comdu/,/comdv/,/comdw/)
!$OMP   threadprivate (/comfluxbt/)
!$OMP   threadprivate (/comhor0/)
!$OMP   threadprivate (/comswitches/)
!$OMP   threadprivate (/comevol/)

