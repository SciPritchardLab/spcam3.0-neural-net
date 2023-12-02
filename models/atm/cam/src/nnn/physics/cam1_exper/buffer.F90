#include <misc.h>
#include <params.h>

module buffer

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Definition and initialization of time-cycling physics arrays
!
! Author: 
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, pnats
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use prognostics, only: ptimelevels
  use infnan
#ifdef CRM
  use crmdims, only: crm_nx, crm_ny, crm_nz
#endif

  real(r8), allocatable :: qrs(:,:,:)      ! shortwave radiative heating rate 
  real(r8), allocatable :: qrl(:,:,:)      ! longwave  radiative heating rate 

  real(r8), allocatable :: pblht(:,:)      ! PBL height
  real(r8), allocatable :: tpert(:,:)      ! temperature perturbation (PBL)
  real(r8), allocatable :: qpert(:,:,:)    ! moisture/constituent perturb.(PBL)

#ifdef CRM

! CRM (Superparameterization) prognostic arrays held in core (M.Khairoutdinov):

   real(r8), allocatable :: u_crm  (:,:,:,:,:)
   real(r8), allocatable :: v_crm  (:,:,:,:,:)
   real(r8), allocatable :: w_crm  (:,:,:,:,:)
   real(r8), allocatable :: t_crm  (:,:,:,:,:)
   real(r8), allocatable :: q_crm  (:,:,:,:,:)
   real(r8), allocatable :: qn_crm (:,:,:,:,:)
   real(r8), allocatable :: qp_crm (:,:,:,:,:)
   real(r8), allocatable :: qrs_crm(:,:,:,:,:)
   real(r8), allocatable :: qrl_crm(:,:,:,:,:)

  integer, parameter :: nrad_buffer = 23
  real(r8), allocatable :: rad_buffer(:,:,:)
  real(r8), allocatable :: qrs1(:,:,:)      ! shortwave radiative heating rate
  real(r8), allocatable :: qrl1(:,:,:)      ! longwave  radiative heating rate
  real(r8), allocatable :: fsds_crm(:,:,:,:) ! Flux Shortwave Downwelling Surface
  real(r8), allocatable :: fsns_crm(:,:,:,:) ! Surface solar absorbed flux
  real(r8), allocatable :: fsdsc_crm(:,:,:,:) ! clearsky downward solar flux at surface
  real(r8), allocatable :: fsntoa_crm(:,:,:,:) ! Net solar flux at TOA
  real(r8), allocatable :: fsntoac_crm(:,:,:,:) ! clearsky net downward solar flux at TOA
  real(r8), allocatable :: fsutoa_crm(:,:,:,:) ! Flux Shortwave Upwelling TOA
  real(r8), allocatable :: flwds_crm(:,:,:,:) ! Surface longwave down flux
  real(r8), allocatable :: flns_crm(:,:,:,:) ! Srf longwave cooling (up-down) flux
  real(r8), allocatable :: flnsc_crm(:,:,:,:) ! clearsky Srf longwave cooling (up-down) flux
  real(r8), allocatable :: flut_crm(:,:,:,:) ! Upward lw flux at top of model
  real(r8), allocatable :: flutc_crm(:,:,:,:) ! clearksy outgoing lw flux at model top
#endif

#ifdef LES

! CRM-LES (Superparameterization) prognostic arrays held in core (M.Khairoutdinov):

   real(r8), allocatable :: u_les  (:,:,:,:,:)
   real(r8), allocatable :: v_les  (:,:,:,:,:)
   real(r8), allocatable :: w_les  (:,:,:,:,:)
   real(r8), allocatable :: t_les  (:,:,:,:,:)
   real(r8), allocatable :: q_les  (:,:,:,:,:)
   real(r8), allocatable :: qn_les (:,:,:,:,:)
   real(r8), allocatable :: qp_les (:,:,:,:,:)
   real(r8), allocatable :: qrs_les(:,:,:,:,:)
   real(r8), allocatable :: qrl_les(:,:,:,:,:)

#endif

CONTAINS

  subroutine initialize_buffer
!
! Allocate memory
!
    allocate (qrs   (pcols,pver,begchunk:endchunk))     
    allocate (qrl   (pcols,pver,begchunk:endchunk))     
    allocate (pblht (pcols,begchunk:endchunk))       
    allocate (tpert (pcols,begchunk:endchunk))       
    allocate (qpert (pcols,pcnst+pnats,begchunk:endchunk)) 
!
! Initialize to NaN or Inf
!
    qrs   (:,:,:) = inf
    qrl   (:,:,:) = inf
    pblht (:,:) = inf
    tpert (:,:) = inf
    qpert (:,:,:) = inf

#ifdef CRM

    allocate (u_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (v_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (w_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (t_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (q_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qn_crm  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qp_crm (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qrs_crm(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qrl_crm(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))

    u_crm  (:,:,:,:,:) = inf
    v_crm  (:,:,:,:,:) = inf
    w_crm  (:,:,:,:,:) = inf
    t_crm  (:,:,:,:,:) = inf
    q_crm  (:,:,:,:,:) = inf
    qn_crm  (:,:,:,:,:) = inf
    qp_crm (:,:,:,:,:) = inf
    qrs_crm(:,:,:,:,:) = inf
    qrl_crm(:,:,:,:,:) = inf

    allocate (rad_buffer(pcols,nrad_buffer,begchunk:endchunk))
    allocate (qrs1   (pcols,pver,begchunk:endchunk))
    allocate (qrl1   (pcols,pver,begchunk:endchunk))
    rad_buffer(:,:,:) = inf
    qrs1   (:,:,:) = inf
    qrl1   (:,:,:) = inf

    allocate (fsds_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (fsns_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (fsntoa_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (fsdsc_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (fsntoac_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (fsutoa_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (flwds_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (flns_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (flnsc_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (flut_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))
    allocate (flutc_crm(pcols,crm_nx,crm_ny,begchunk:endchunk))

    fsds_crm(:,:,:,:) = inf
    fsns_crm(:,:,:,:) = inf
    fsntoa_crm(:,:,:,:) = inf
    fsdsc_crm(:,:,:,:) = inf
    fsntoac_crm(:,:,:,:) = inf
    fsutoa_crm(:,:,:,:) = inf
    flwds_crm(:,:,:,:) = inf
    flns_crm(:,:,:,:) = inf
    flnsc_crm(:,:,:,:) = inf
    flut_crm(:,:,:,:) = inf
    flutc_crm(:,:,:,:) = inf

#endif
    
#ifdef LES

    allocate (u_les  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (v_les  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (w_les  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (t_les  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (q_les  (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qn_les (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qp_les (pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qrs_les(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))
    allocate (qrl_les(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk))

    u_les  (:,:,:,:,:) = inf
    v_les  (:,:,:,:,:) = inf
    w_les  (:,:,:,:,:) = inf
    t_les  (:,:,:,:,:) = inf
    q_les  (:,:,:,:,:) = inf
    qn_les (:,:,:,:,:) = inf
    qp_les (:,:,:,:,:) = inf
    qrs_les(:,:,:,:,:) = inf
    qrl_les(:,:,:,:,:) = inf

#endif
    
    return
  end subroutine initialize_buffer

end module buffer
