#undef DEBUG
#include <misc.h>

module cldcond

!---------------------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the prognostic cloud water and ice parameterization
!
! Author: Byron Boville  Sept 04, 2002
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver
  use physconst,     only: gravit, latvap, latice

  implicit none
  private
  save

  public :: cldcond_register, cldcond_init_cnst, cldcond_implements_cnst
  public :: cldcond_init, cldcond_tend
  public :: cldcond_zmconv_detrain
  public :: cldcond_sediment

! Private module data

  integer, parameter :: ncnst=2                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'CLDLIQ', 'CLDICE'/)

  integer :: &
       ixcldice,     &! cloud ice water index
       ixcldliq,     &! cloud liquid water index
       qcwat_idx,    &! qcwat index in physics buffer
       tcwat_idx,    &! tcwat index in physics buffer
       cld_idx,      &! cld index in physics buffer
       lcwat_idx      ! lcwat index in physics buffer

contains

!===============================================================================
  subroutine cldcond_sediment(state, ptend, dtime, &
       cloud, icefrac, landfrac, ocnfrac, prec, snow, landm, snowh)
!
! Interface to sedimentation of cloud liquid and ice particles
!
! NOTE: initial implementation by B.A. Boville from earlier code by P.J. Rasch
!
! B. A. Boville, Sept 20, 2002
!
!-----------------------------------------------------------------------
    use physics_types,    only: physics_state, physics_ptend
    use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
    use history,          only: outfld

    implicit none

! Arguments
    type(physics_state), intent(in)    :: state   ! state variables
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in)  :: cloud(pcols,pver)    ! cloud fraction
    real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
    real(r8), intent(in)  :: dtime                ! timestep

    real(r8), intent(out) :: prec(pcols)          ! surface flux of total cloud water
    real(r8), intent(out) :: snow(pcols)          ! surface flux of cloud ice
    real(r8), intent(in) :: landm(pcols)          ! land fraction ramped over water
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

! Local variables
    integer  :: i,k                               ! loop indexes
    integer  :: ncol                              ! number of atmospheric columns in chunk
    integer  :: lchnk                             ! chunk index
    real(r8) :: rain(pcols)                       ! surface flux of cloud liquid
    real(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
    real(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)
!-----------------------------------------------------------------------
    ncol = state%ncol

    ptend%name         = 'pcwsediment'
    ptend%ls           = .TRUE.
    ptend%lq(1)        = .TRUE.
    ptend%lq(ixcldice) = .TRUE.
    ptend%lq(ixcldliq) = .TRUE.

    call cld_sediment_vel (ncol,                                    &
         icefrac, landfrac, ocnfrac, state%pmid, state%pdel, state%t, &
         cloud, state%q(:,:,ixcldliq), state%q(:,:,ixcldice), pvliq, pvice, landm, snowh)

    call cld_sediment_tend (ncol, dtime ,                                       &
         state%pint, state%pmid           , state%pdel            , state%t     , &
         cloud  ,    state%q(:,:,ixcldliq), state%q(:,:,ixcldice) , pvliq, pvice, &
         ptend%q(:,:,ixcldliq), ptend%q(:,:,ixcldice), ptend%q(:,:,1), ptend%s  , &
         rain   , snow   )

! convert rain and snow from kg/m2 to m/s
    snow(:ncol) = snow(:ncol)/1000.
    rain(:ncol) = rain(:ncol)/1000.
! compute total precip (m/s)
    prec(:ncol) = rain(:ncol) + snow(:ncol)

! record history variables
    lchnk = state%lchnk
    call outfld('DQSED'   ,ptend%q(:,:,1)       , pcols,lchnk)
    call outfld('DISED'   ,ptend%q(:,:,ixcldice), pcols,lchnk)
    call outfld('DLSED'   ,ptend%q(:,:,ixcldliq), pcols,lchnk)
    call outfld('HSED'    ,ptend%s              , pcols,lchnk)
    call outfld('PRECSED' ,prec                 , pcols,lchnk)
    call outfld('SNOWSED' ,snow                 , pcols,lchnk)
    call outfld('RAINSED' ,rain                 , pcols,lchnk)

    return
  end subroutine cldcond_sediment

!===============================================================================
  subroutine cldcond_zmconv_detrain(dlf, cld, state, ptend)
!
! Partition the detrained condensed water from the ZM convection scheme.
!
! The ZM scheme does not have an ice phase, so the detrained water is paritioned
! soley between cloud liquid water and the environment. The ice/liquid partitioning 
! happens in cldcond_tend.
!
! NOTE: initial implementation by B.A. Boville just moves code here from TPHYSBC.
!
! B. A. Boville, Sept 09, 2002
!
!-----------------------------------------------------------------------
    use physics_types, only: physics_state, physics_ptend
    use history,       only: outfld

    implicit none

! Arguments
    type(physics_state), intent(in  )  :: state   ! state variables
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in) :: dlf(pcols,pver)       ! detrained water from ZM
    real(r8), intent(in) :: cld(pcols,pver)       ! cloud fraction

! Local variables
    integer :: i,k                            ! loop indexes
!-----------------------------------------------------------------------

    ptend%name         = 'pcwdetrain'
!!$    ptend%ls           = .TRUE.
!!$    ptend%lq(1)        = .TRUE.
!!$    ptend%lq(ixcldice) = .TRUE.
    ptend%lq(ixcldliq) = .TRUE.
!
! Put all of the detraining cloud water from convection into the large scale cloud.
! It all goes in liquid for the moment.
    do k = 1,pver
       do i = 1,state%ncol
!!$          ptend%q(i,k,1)        = dlf(i,k) * (1.-cld(i,k))
!!$          ptend%s(i,k)          =-dlf(i,k) * (1.-cld(i,k))*latvap
!!$          ptend%q(i,k,ixcldice) = 0.
!!$          ptend%q(i,k,ixcldliq) = dlf(i,k) * cld(i,k)
          ptend%q(i,k,ixcldliq) = dlf(i,k)
       end do
    end do
    call outfld('ZMDLF' ,dlf  , pcols,state%lchnk)

    return
  end subroutine cldcond_zmconv_detrain

!===============================================================================
  subroutine cldcond_register()
!
! Register the constituents (cloud liquid and cloud ice) and the fields
! in the physics buffer.
! 
!-----------------------------------------------------------------------
    use constituents, only: cnst_add, advected, nonadvec
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add

    implicit none

!    logical, parameter :: cldw_adv=.false.  ! true => cloud water is treated as advected tracer
    logical, parameter :: cldw_adv=.true.  ! true => cloud water is treated as advected tracer

    integer flag
!-----------------------------------------------------------------------

! Register cloud water and determine index (either advected or non-adv).
    if (cldw_adv) then
       flag = advected
    else
       flag = nonadvec
    endif
    call cnst_add(cnst_names(1), flag, mwdry, cpair, 0._r8, ixcldliq, &
         longname='Grid box averaged liquid condensate amount')
    call cnst_add(cnst_names(2), flag, mwdry, cpair, 0._r8, ixcldice, &
         longname='Grid box averaged ice condensate amount')

! Request physics buffer space for fields that persist across timesteps.
    call pbuf_add('QCWAT', 'global', 1,pver,pbuf_times, qcwat_idx)
    call pbuf_add('TCWAT', 'global', 1,pver,pbuf_times, tcwat_idx)
    call pbuf_add('CLD',   'global', 1,pver,pbuf_times, cld_idx)
    call pbuf_add('LCWAT', 'global', 1,pver,pbuf_times, lcwat_idx)

    call pbuf_add('QINI' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('TINI' , 'physpkg', 1,pver,      1, flag)

  end subroutine cldcond_register

!===============================================================================

  function cldcond_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: cldcond_implements_cnst     ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     cldcond_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           cldcond_implements_cnst = .true.
           return
        end if
     end do
  end function cldcond_implements_cnst

!===============================================================================
  subroutine cldcond_init_cnst(name, q)
!
! Initialize the cloud water mixing ratios (liquid and ice), if they are
! not read from the initial file
! 
!-----------------------------------------------------------------------
    use pmgrid,        only: plon, plev, plat

    implicit none

! Arguments
    character(len=*), intent(in)  :: name                ! constituent name
    real(r8),         intent(out) :: q(plon,plev,plat)   ! mass mixing ratio
!-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0
       return
    end if

  end subroutine cldcond_init_cnst

!===============================================================================
  subroutine cldcond_init()
!
! Initialize the cloud water parameterization
! 
!-----------------------------------------------------------------------
    use cldwat,        only: inimc
    use history,       only: addfld, add_default, phys_decomp
    use physconst,     only: tmelt, rh2o, rhodair

    implicit none
!-----------------------------------------------------------------------

! initialization routine for prognostic cloud water
    call inimc (tmelt, rhodair/1000.0, gravit, rh2o )

! register history variables
    call addfld ('FWAUT   ','fraction',pver, 'A','Relative importance of liquid autoconversion' ,phys_decomp)
    call addfld ('FSAUT   ','fraction',pver, 'A','Relative importance of ice autoconversion'    ,phys_decomp)
    call addfld ('FRACW   ','fraction',pver, 'A','Relative  importance of rain accreting liquid',phys_decomp)
    call addfld ('FSACW   ','fraction',pver, 'A','Relative  importance of snow accreting liquid',phys_decomp)
    call addfld ('FSACI   ','fraction',pver, 'A','Relative  importance of snow accreting ice'   ,phys_decomp)
    call addfld ('CME     ','kg/kg/s ',pver, 'A','Rate of cond-evap within the cloud'           ,phys_decomp)
    call addfld ('ZMDLF   ','kg/kg/s ',pver, 'A','Detrained liquid water from ZM convection'    ,phys_decomp)
    call addfld ('PRODPREC','kg/kg/s ',pver, 'A','Rate of conversion of condensate to precip'   ,phys_decomp)
    call addfld ('EVAPPREC','kg/kg/s ',pver, 'A','Rate of evaporation of falling precip'        ,phys_decomp)
    call addfld ('EVAPSNOW','kg/kg/s ',pver, 'A','Rate of evaporation of falling snow'          ,phys_decomp)
    call addfld ('HPROGCLD','W/kg'    ,pver, 'A','Heating from prognostic clouds'               ,phys_decomp)
    call addfld ('HEVAP   ','W/kg'    ,pver, 'A','Heating from evaporation of falling precip'   ,phys_decomp)
    call addfld ('HMELT   ','W/kg'    ,pver, 'A','Heating from snow melt'                       ,phys_decomp)
    call addfld ('HREPART ','W/kg'    ,pver, 'A','Heating from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('FICE    ','fraction',pver, 'A','Fractional ice content within cloud'          ,phys_decomp)
    call addfld ('ICWMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud water mixing ratio'       ,phys_decomp)
    call addfld ('ICIMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud ice mixing ratio'         ,phys_decomp)
    call addfld ('PCSNOW  ','m/s'     ,1   , 'A','Snow fall from prognostic clouds'             ,phys_decomp)

    call addfld ('DQSED   ','kg/kg/s' ,pver, 'A','Water vapor tendency from cloud sedimentation',phys_decomp)
    call addfld ('DLSED   ','kg/kg/s' ,pver, 'A','Cloud liquid tendency from sedimentation'     ,phys_decomp)
    call addfld ('DISED   ','kg/kg/s' ,pver, 'A','Cloud ice tendency from sedimentation'        ,phys_decomp)
    call addfld ('HSED    ','W/kg'    ,pver, 'A','Heating from cloud sediment evaporation'      ,phys_decomp)
    call addfld ('SNOWSED ','m/s'     ,1   , 'A','Snow from cloud ice sedimentation'            ,phys_decomp)
    call addfld ('RAINSED ','m/s'     ,1   , 'A','Rain from cloud liquid sedimentation'         ,phys_decomp)
    call addfld ('PRECSED ','m/s'     ,1   , 'A','Precipitation from cloud sedimentation'       ,phys_decomp)

    call add_default ('FICE    ', 1, ' ')
!  These fields removed 10/30/2003 per CRB decision.
!    call add_default ('CME     ', 1, ' ')
!    call add_default ('ZMDLF   ', 1, ' ')
!    call add_default ('PRODPREC', 1, ' ')
!    call add_default ('EVAPPREC', 1, ' ')
!    call add_default ('EVAPSNOW', 1, ' ')
!    call add_default ('HPROGCLD', 1, ' ')
!    call add_default ('HEVAP   ', 1, ' ')
!    call add_default ('HMELT   ', 1, ' ')
!    call add_default ('HREPART ', 1, ' ')
!    call add_default ('ICWMR   ', 1, ' ')
!    call add_default ('ICIMR   ', 1, ' ')
!    call add_default ('PCSNOW  ', 1, ' ')
!    call add_default ('DQSED   ', 1, ' ')
!    call add_default ('DLSED   ', 1, ' ')
!    call add_default ('DISED   ', 1, ' ')
!    call add_default ('HSED    ', 1, ' ')
!    call add_default ('SNOWSED ', 1, ' ')
!    call add_default ('RAINSED ', 1, ' ')
!    call add_default ('PRECSED ', 1, ' ')

    return
  end subroutine cldcond_init

!===============================================================================
!!$  subroutine cldcond_tend(state, ptend, dt, pbuf)
  subroutine cldcond_tend(state, ptend, dt, &
       tcwato, qcwato, lcwato, precip, snow, icefrac, rhdfda, rhu00, cldn, evapprec, prodprec, cme, snowh)
!
! Compute the tendency of cloud water mixing ratios (liquid and ice) from
! microphysical processes and sedimenation of cloud particles.
! 
! Author: Byron Boville  Sept 04, 2002
!  modified pjr: march 2003 to repartition snow/rain
!
!-----------------------------------------------------------------------

    use physics_types, only: physics_state, physics_ptend
!!$    use phys_buffer,   only: pbuf_size_max, pbuf_fld
    use history,       only: outfld
    use cldwat,        only: pcond, cldwat_fice
    use physconst,     only: tmelt

    implicit none

! Arguments
    type(physics_state), intent(inout) :: state          ! state variables
    type(physics_ptend), intent(inout) :: ptend          ! package tendencies
    real(r8),            intent(in)    :: dt             ! timestep
!!$    type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
    real(r8), intent(in)  :: tcwato(pcols,pver)        !cloud water old temperature
    real(r8), intent(in)  :: qcwato(pcols,pver)        ! cloud water old q
    real(r8), intent(in)  :: lcwato(pcols,pver)        ! cloud liquid water old q
    real(r8), intent(in)  :: icefrac(pcols)            ! sea ice fraction (fraction)
    real(r8), intent(in)  :: rhdfda(pcols,pver)        ! dRh/dcloud, old 
    real(r8), intent(in)  :: rhu00 (pcols,pver)        ! Rh threshold for cloud, old
    real(r8), intent(in)  :: cldn(pcols,pver)          !new cloud fraction
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

    real(r8), intent(out) :: precip(pcols)             ! sfc flux of precip (m/s)
    real(r8), intent(out) :: snow (pcols)              ! sfc flux of snow   (m/s)
    real(r8), intent(out) :: evapprec(pcols,pver)          ! local evaporation of precipitation
    real(r8), intent(out) :: prodprec(pcols,pver)          ! local production of precipitation
    real(r8), intent(out) :: cme     (pcols,pver)          ! local condensation - evaporation of cloud water

! Local variables
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! number of atmospheric columns in chunk
    integer :: i,k                            ! loop indexes
!!$    real(r8), pointer, dimension(:)   :: buffld1  ! physics buffer field1
!!$    real(r8), pointer, dimension(:,:) :: buffld2  ! physics buffer field2
    real(r8) :: rdt                          ! 1./dt
    real(r8) :: qtend (pcols,pver)            ! moisture tendencies
    real(r8) :: ttend (pcols,pver)            ! temperature tendencies
    real(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
!    real(r8) :: cme     (pcols,pver)          ! local condensation - evaporation of cloud water
    real(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
!    real(r8) :: evapprec(pcols,pver)          ! local evaporation of precipitation
    real(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
    real(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
    real(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
!    real(r8) :: prodprec(pcols,pver)          ! local production of precipitation
    real(r8) :: prodsnow(pcols,pver)          ! local production of snow
    real(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
    real(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
    real(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
    real(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
    real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
    real(r8) fwaut(pcols,pver)              
    real(r8) fsaut(pcols,pver)              
    real(r8) fracw(pcols,pver)              
    real(r8) fsacw(pcols,pver)              
    real(r8) fsaci(pcols,pver)              
    real(r8) ice2pr(pcols,pver)   ! rate of conversion of ice to precip
    real(r8) liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
    real(r8) liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
    real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
    real(r8) temp(pcols), w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

!-----------------------------------------------------------------------

! Set output flags
    ptend%name         = 'cldwat'
    ptend%ls           = .true.
    ptend%lq(1)        = .true.
    ptend%lq(ixcldice) = .true.
    ptend%lq(ixcldliq) = .true.

! Initialize chunk id and size
    lchnk = state%lchnk
    ncol  = state%ncol

! associate local pointers with fields in the physics buffer
!!$    buffld1 => pbuf(ixbuffld1)%fld_ptr(1,1:pcols,1,     lchnk,1)
!!$    buffld2 => pbuf(ixbuffld2)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

! Define fractional amount of cloud condensate in ice phase
    call cldwat_fice(ncol, state%t, fice, fsnow)

! compute total cloud water
    totcw(:ncol,:pver) = state%q(:ncol,:pver,ixcldice) + state%q(:ncol,:pver,ixcldliq)

! save input cloud ice
    repartht(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

! Repartition ice and liquid
    state%q(:ncol,:pver,ixcldice) = totcw(:ncol,:pver) * fice(:ncol,:pver)
    state%q(:ncol,:pver,ixcldliq) = totcw(:ncol,:pver) * (1.0_r8 - fice(:ncol,:pver))

! Determine heating from change in cloud ice
    repartht(:ncol,:pver) = latice/dt * (state%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver))

! calculate the tendencies for moisture, temperature and cloud fraction
    rdt = 1./dt
    qtend(:ncol,:pver) = (state%q(:ncol,:pver,1) - qcwato(:ncol,:pver))*rdt
    ttend(:ncol,:pver) = (state%t(:ncol,:pver)   - tcwato(:ncol,:pver))*rdt
    ltend(:ncol,:pver) = (totcw  (:ncol,:pver)   - lcwato(:ncol,:pver))*rdt

! call microphysics package to calculate tendencies
    call t_startf('pcond')
    call pcond (lchnk,   ncol, &
         state%t  , ttend      , state%q(1,1,1), qtend   , state%omega, &
         totcw    , state%pmid , state%pdel , cldn       , fice       , fsnow, &
         cme      , prodprec   , prodsnow, evapprec   , evapsnow   , evapheat   , prfzheat, &
         meltheat , precip     , snow       , dt         , fwaut      , &
         fsaut    , fracw      , fsacw      , fsaci      , ltend      , &
         rhdfda   , rhu00      , icefrac    , state%zi   , ice2pr, liq2pr, liq2snow, snowh)
    call t_stopf('pcond')


! make it interactive
    do i = 1,ncol
       pr1 = 0
       hs1 = 0
       qv1 = 0
       ql1 = 0
       qi1 = 0
       qs1 = 0
       qr1 = 0
       w1 = 0
       wl = 0
       wv = 0
       wi = 0
       wlf = 0
       wvf = 0
       wif = 0
       do k = 1,pver

          ptend%s(i,k)          = cme(i,k) * (latvap + latice*fice(i,k)) &
               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)

          ptend%q(i,k,1)        =-cme(i,k) + evapprec(i,k)

          ptend%q(i,k,ixcldice) =cme(i,k)*fice(i,k)      - ice2pr(i,k)
          ptend%q(i,k,ixcldliq) =cme(i,k)*(1.-fice(i,k)) - liq2pr(i,k)

!         santity checks to be removed after 5 year run
          res = state%q(i,k,1) + ptend%q(i,k,1)*dt
          if (res.lt.0.) then
             write (6,*) ' predicted neg vapor ', i,k,lchnk, res
             write (6,*) ' oldq, cme(i,k)*dt, evapprec(i,k)*dt ', state%q(i,k,1), &
                  -cme(i,k)*dt, evapprec(i,k)*dt, ptend%q(i,k,1)*dt
             call endrun()
          endif
          res = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dt
          if (res.lt.0.) then
             write (6,*) ' cldcond_tend: predicted neg ice ', i,k,lchnk, res, fice(i,k)
             write (6,*) ' oldice, cme(i,k)*dt*(fice(i,k)), ice2pr*dt, ptend_ice*dt', &
                  state%q(i,k,ixcldice), &
                  cme(i,k)*dt*(fice(i,k)), ice2pr(i,k)*dt, ptend%q(i,k,ixcldice)*dt
             write (6,*) ' oldcond, cme(i,k)*dt, prodprec*dt, ptend_cond*dt', &
                  state%q(i,k,ixcldice)+state%q(i,k,ixcldliq), &
                  cme(i,k)*dt, prodprec(i,k)*dt, &
                  (ptend%q(i,k,ixcldice)+ptend%q(i,k,ixcldliq))*dt
             call endrun()
          endif
!         end sanity check

#ifdef DEBUG
          if (lchnk.eq.248.and.i.eq.12) then

             write (6,*) 
             write (6,*) ' input state, t, q, l, i ', k, state%t(i,k), state%q(i,k,1), state%q(i,k,ixcldliq),  state%q(i,k,ixcldice)
             write (6,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
             write (6,*) ' total precip before accumulation                      ', k, pr1

             wv = wv + state%q(i,k,1       )*state%pdel(i,k)/gravit
             wl = wl + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi = wi + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit

             qvf = state%q(i,k,1) + ptend%q(i,k,1)*dt
             qlf = state%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dt
             qif = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dt

             if (qvf.lt.0.) then
                write (6,*) ' qvf is negative *******', qvf
             endif
             if (qlf.lt.0.) then
                write (6,*) ' qlf is negative *******', qlf
             endif
             if (qif.lt.0.) then
                write (6,*) ' qif is negative *******', qif
             endif
             write (6,*) ' qvf, qlf, qif ', qvf, qlf, qif

             wvf = wvf + qvf*state%pdel(i,k)/gravit
             wlf = wlf + qlf*state%pdel(i,k)/gravit
             wif = wif + qif*state%pdel(i,k)/gravit

             hs1 = hs1 + ptend%s(i,k)*state%pdel(i,k)/gravit
             pr1 = pr1 + state%pdel(i,k)/gravit*(prodprec(i,k)-evapprec(i,k))
             qv1 = qv1 - (cme(i,k)-evapprec(i,k))*state%pdel(i,k)/gravit    ! vdot
             w1  = w1  + (cme(i,k)-prodprec(i,k))*state%pdel(i,k)/gravit    ! cdot
             qi1 = qi1 + ((cme(i,k))*fice(i,k)        -ice2pr(i,k) )*state%pdel(i,k)/gravit   ! idot
             ql1 = ql1 + ((cme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state%pdel(i,k)/gravit   ! ldot

             qr1 = qr1 &
                  + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
                   -(evapprec(i,k)-evapsnow(i,k)) &     ! rain evaporation
                    )*state%pdel(i,k)/gravit
             qs1 = qs1 &
                  + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
                     -evapsnow(i,k)               &     ! snow evaporation
                    )*state%pdel(i,k)/gravit

             if (state%t(i,k).gt.tmelt) then
                qr1 = qr1 + qs1
                qs1 = 0.
             endif
             write (6,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
             write (6,*) ' total precip after accumulation      ', k, pr1
             write (6,*)
             write (6,*) ' layer prodprec, evapprec, pdel ', prodprec(i,k), evapprec(i,k), state%pdel(i,k)
             write (6,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
             write (6,*) ' layer prodprec-prodsnow, liq2pr-liq2snow ', prodprec(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
             write (6,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), evapprec(i,k)-evapsnow(i,k)
             write (6,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
             write (6,*) ' layer ice2pr+liq2pr, prodprec ', ice2pr(i,k)+liq2pr(i,k), prodprec(i,k)
             write (6,*)
             write (6,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
             write (6,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
             write (6,*) ' condensate produce after accum                   ', k, w1
             write (6,*) ' liq+ice tends accum                              ', k, ql1+qi1
             write (6,*) ' change in total water after accum                ', k, qv1+ql1+qi1
             write (6,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
             write (6,*) ' fice at this lev ', fice(i,k)
             write (6,*)

             res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36))
             write (6,*) ' relative residual in column method 1             ', k, res

             write (6,*) ' relative residual in column method 2             ', k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36))
!            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
             if (res.gt.1.e-14) then
                call endrun()
             endif

!             w3  = cme(i,k) * (latvap + latice*fice(i,k)) &
!               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

             res = qs1+qr1-pr1
             w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
             if (w4.gt.0.)  then
                if (res/w4.gt.1.e-14) then
                   write (6,*) ' imbalance in precips calculated two ways '
                   write (6,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                        res/w4, pr1, qr1, qs1, qr1+qs1
!                   call endrun()
                endif
             endif
             if (k.eq.pver) then
                write (6,*) ' pcond returned precip, rain and snow rates ', precip(i), precip(i)-snow(i), snow(i)
                write (6,*) ' I calculate ', pr1, qr1, qs1
!               call endrun
                write (6,*) ' byrons water check ', wv+wl+wi-pr1*dt, wvf+wlf+wif
             endif
          write (6,*)
          endif
#endif
       end do
    end do


#ifdef DEBUG
    if (.true.) then
    do i = 1,ncol
       if (snow(i).gt.0.01/8.64e4.and.state%t(i,pver).gt.273.16) then
          write (6,*) ' cldcond: snow, temp, ', i, lchnk, &
               snow(i), state%t(i,pver)
          write (6,*) ' t ', state%t(i,:)
          write (6,*) ' fsaut ', fsaut(i,:)
          write (6,*) ' fsacw ', fsacw(i,:)
          write (6,*) ' fsaci ', fsaci(i,:)
          write (6,*) ' meltheat ', meltheat(i,:)
          call endrun()
       endif

       if (snow(i)*8.64e4.lt.-1.e-5) then
          write (6,*) ' neg snow ', snow(i)*8.64e4
          write (6,*) ' cldcond: snow, temp, ', i, lchnk, &
               snow(i), state%t(i,pver)
          write (6,*) ' t ', state%t(i,:)
          write (6,*) ' fsaut ', fsaut(i,:)
          write (6,*) ' fsacw ', fsacw(i,:)
          write (6,*) ' fsaci ', fsaci(i,:)
          write (6,*) ' meltheat ', meltheat(i,:)
          call endrun()
       endif
    end do
    endif
#endif

! Compute in cloud ice and liquid mixing ratios
    do k=1,pver
       do i = 1,ncol
          icimr(i,k) = (state%q(i,k,ixcldice) + dt*ptend%q(i,k,ixcldice)) / max(0.01_r8,cldn(i,k))
          icwmr(i,k) = (state%q(i,k,ixcldliq) + dt*ptend%q(i,k,ixcldliq)) / max(0.01_r8,cldn(i,k))
       end do
    end do


! convert precipitation from kg/m2 to m/s
    snow  (:ncol) = snow  (:ncol)/1000.
    precip(:ncol) = precip(:ncol)/1000.

! record history variables
    call outfld('FWAUT',fwaut,    pcols,lchnk)
    call outfld('FSAUT',fsaut,    pcols,lchnk)
    call outfld('FRACW',fracw,    pcols,lchnk)
    call outfld('FSACW',fsacw,    pcols,lchnk)
    call outfld('FSACI',fsaci,    pcols,lchnk)
    call outfld('ICIMR',icimr,    pcols,lchnk)
    call outfld('ICWMR',icwmr,    pcols,lchnk)

    call outfld('PCSNOW'  ,snow    , pcols,lchnk)
    call outfld('FICE'    ,fice    , pcols,lchnk)
    call outfld('CME'     ,cme     , pcols,lchnk)
    call outfld('PRODPREC',prodprec, pcols,lchnk)
    call outfld('EVAPPREC',evapprec, pcols,lchnk)
    call outfld('EVAPSNOW',evapsnow, pcols,lchnk)
    call outfld('HPROGCLD',ptend%s , pcols,lchnk)
    call outfld('HEVAP   ',evapheat, pcols,lchnk)
    call outfld('HMELT'   ,meltheat, pcols,lchnk)
    call outfld('HREPART' ,repartht, pcols,lchnk)

! update boundary quantities
!!$    ptend%hflx_srf = 0.
!!$    ptend%hflx_top = 0.
!!$    ptend%cflx_srf = 0.
!!$    ptend%cflx_top = 0.
  end subroutine cldcond_tend
end module cldcond

