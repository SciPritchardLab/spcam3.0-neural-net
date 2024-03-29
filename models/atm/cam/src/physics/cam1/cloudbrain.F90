#ifdef CLOUDBRAIN

#include <misc.h>
#include <params.h>
!#define BRAINDEBUG

module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use pmgrid,          only: masterproc
use string_utils,    only: to_upper

! -------- NEURAL-FORTRAN --------
! imports
use mod_kinds, only: ik, rk
use mod_network , only: network_type
use mod_ensemble, only: ensemble_type
! --------------------------------

  implicit none
  save 

  private
  !!!! CAM namelist variables !!!!
  ! these are CAM namelist variables, but they are initalized with the following values
  ! for a backward compatibility (Jerry's original NN for spreadtesting). That is,
  ! these values are used when these variables are not specified in atm_in.
  integer :: nstepNN = 48  ! a model time step at which NN turns on
  character(len=256) :: nn_in_out_vars = 'in_TBP_QBP_PS_SOLIN_SHF_LHF_out_TPHYSTND_PHQ'
                        ! A userdefined name for a specific input/output combination
  integer :: inputlength = 64  ! length of input vector for NN
  integer :: outputlength = 60 ! length of output vector for NN
  !!!!

#ifdef ENSEMBLE
  real(rk) :: noise = 0.0
  type(ensemble_type) :: cloudbrain_ensemble
#else
  type(network_type) :: cloudbrain_net
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n
#endif

  real(rk), allocatable :: inp_sub(:)
  real(rk), allocatable :: inp_div(:)
  real(rk), allocatable :: out_scale(:)
  real(rk), allocatable :: input(:)
  real(r8), allocatable :: output(:)

  real(r8), allocatable :: dtdt_m1(:,:,:), dqdt_m1(:,:,:)

  type nn_in_t
    real(r8),dimension(pver) :: &
      tbp,   &
      qbp,   &
      vbp,   &
      o3vmr, &
      dtdtm1,& ! tphystnd tendency at timestep n-1
      dqdtm1   ! phq      tendency at timestep n-1
    real(r8) :: &
      ps,     &
      lhf,    &
      shf,    &
      solin,  &
      coszrs
  end type nn_in_t

  type nn_out_t
    real(r8),dimension(pver) :: &
      phq,      &
      tphystnd
  end type nn_out_t

  public neural_net, init_keras_matrices, init_keras_norm, nstepNN, &
         nn_in_t, nn_out_t, nn_in_out_vars, inputlength, outputlength, &
         init_nn_vectors, &
         dtdt_m1, dqdt_m1 ! buffer variables for previous timestep tendencies

  contains

  subroutine neural_net (nn_in, nn_out, icol)
    type(nn_in_t), intent(in) :: nn_in
    type(nn_out_t), intent(out) :: nn_out
    integer, intent(in) :: icol
    
    ! local variables
    integer :: k, nlev, n

    ! 1. Concatenate input vector to neural network
    nlev=30
    select case (to_upper(trim(nn_in_out_vars)))
      case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_TPHYSTND_PHQ')
        input(1:nlev) = nn_in%tbp(:nlev)
        input((nlev+1):2*nlev) = nn_in%qbp(:nlev)
        input(2*nlev+1) = nn_in%ps 
        input(2*nlev+2) = nn_in%solin
        input(2*nlev+3) = nn_in%shf
        input(2*nlev+4) = nn_in%lhf
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_OUT_TPHYSTND_PHQ')
        input(1:nlev)            = nn_in%tbp(:nlev)
        input((1*nlev+1):2*nlev) = nn_in%qbp(:nlev)
        input((2*nlev+1):3*nlev) = nn_in%dtdtm1(:nlev)
        input((3*nlev+1):4*nlev) = nn_in%dqdtm1(:nlev)
        input(4*nlev+1) = nn_in%ps
        input(4*nlev+2) = nn_in%solin
        input(4*nlev+3) = nn_in%shf
        input(4*nlev+4) = nn_in%lhf
      case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ')
        input(1:nlev) = nn_in%tbp(:nlev)
        input((nlev+1):2*nlev) = nn_in%qbp(:nlev)
        input(2*nlev+1) = nn_in%ps
        input(2*nlev+2) = nn_in%solin
        input(2*nlev+3) = nn_in%shf
        input(2*nlev+4) = nn_in%lhf
        input((2*nlev+5):(2*nlev+5+nlev-1)) = nn_in%vbp(:nlev)
        input((2*nlev+5+nlev):2*nlev+5+2*nlev-1) = nn_in%o3vmr(:nlev)
        input(2*nlev+5+2*nlev) = nn_in%coszrs
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ')
        input(1:nlev) = nn_in%tbp(:nlev)
        input((1*nlev+1):2*nlev) = nn_in%qbp(:nlev)
        input((2*nlev+1):3*nlev) = nn_in%dtdtm1(:nlev)
        input((3*nlev+1):4*nlev) = nn_in%dqdtm1(:nlev)
        input(4*nlev+1) = nn_in%ps
        input(4*nlev+2) = nn_in%solin
        input(4*nlev+3) = nn_in%shf
        input(4*nlev+4) = nn_in%lhf
        input((4*nlev+5):(4*nlev+5+nlev-1)) = nn_in%vbp(:nlev)
        input((4*nlev+5+nlev):4*nlev+5+2*nlev-1) = nn_in%o3vmr(:nlev)
        input(4*nlev+5+2*nlev) = nn_in%coszrs
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ_ABLATED')
        input(1:nlev) = nn_in%tbp(:nlev)
        input(nlev+1) = 0
        input(nlev+2) = 0
        input(nlev+3) = 0
        input(nlev+4) = 0
        input(nlev+5) = 0
        input((nlev+6):2*nlev-5) = nn_in%qbp(:nlev)
        input(2*nlev+1-5) = nn_in%ps
        input(2*nlev+2-5) = nn_in%solin
        input(2*nlev+3-5) = nn_in%shf
        input(2*nlev+4-5) = nn_in%lhf
        input((2*nlev+5-5):(2*nlev+5-5+nlev-1)) = nn_in%vbp(:nlev)
        input((2*nlev+5-5+nlev):2*nlev+5-5+2*nlev-1) = nn_in%o3vmr(:nlev)
        input(2*nlev+5-5+2*nlev) = nn_in%coszrs
    end select

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input pre norm=',input
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(k) = (input(k) - inp_sub(k))/inp_div(k)
    end do

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input post norm=',input
      endif
#endif

! 3. Neural network matrix multiplications and activations
#ifdef ENSEMBLE
    output = cloudbrain_ensemble % average(input)
#else
    ! use neural fortran library
    output = cloudbrain_net % output(input)
#endif

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG output = ',output
      endif
#endif

    ! 4. Unnormalize output
    do k=1,outputlength
      output(k) = output(k) / out_scale(k)
    end do

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG out post scale = ',output
      endif
#endif

    ! 5. Split output into components
    select case (to_upper(trim(nn_in_out_vars)))
      case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_TPHYSTND_PHQ')
        nn_out%tphystnd(:nlev) = output(1:nlev)  
        nn_out%phq(:nlev) = output((nlev+1):2*nlev) ! This is still the wrong unit, needs to be converted to W/m^2
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_OUT_TPHYSTND_PHQ')
        nn_out%tphystnd(:nlev) = output(1:nlev)
        nn_out%phq(:nlev) = output((nlev+1):2*nlev) ! This is still the wrong unit, needs to be converted to W/m^2
      case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ')
        nn_out%tphystnd(:nlev) = output(1:nlev)
        nn_out%phq(:nlev) = output((nlev+1):2*nlev) ! This is still the wrong unit, needs to be converted to W/m^2
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ')
        nn_out%tphystnd(:nlev) = output(1:nlev)
        nn_out%phq(:nlev) = output((nlev+1):2*nlev) ! This is still the wrong unit, needs to be converted to W/m^2
      case('IN_TBP_QBP_TPHYSTND_PHQ_PS_SOLIN_SHF_LHF_VBP_O3VMR_COSZRS_OUT_TPHYSTND_PHQ_ABLATED')
        nn_out%tphystnd(:nlev) = output(1:nlev)
        nn_out%phq(1) = 0
        nn_out%phq(2) = 0
        nn_out%phq(3) = 0
        nn_out%phq(4) = 0
        nn_out%phq(5) = 0
        nn_out%phq(6:nlev) = output((nlev+1):2*nlev-5) ! This is still the wrong unit, needs to be converted to W/m^2     
    end select

  end subroutine neural_net


  subroutine init_keras_matrices()
#ifdef ENSEMBLE
    write (6,*) '------- NEURAL-FORTRAN: ensemble loading -------'
    cloudbrain_ensemble = ensemble_type('./Models/', noise)
    write (6,*) '------- NEURAL-FORTRAN: ensemble loaded -------'
#else
    call cloudbrain_net % load('./keras_matrices/model.txt')
    if (masterproc) then
      write (6,*) '------- NEURAL-FORTRAN: loaded network from txt file -------'
    end if
#endif
  end subroutine init_keras_matrices
    

  subroutine init_keras_norm()
    allocate(inp_sub (inputlength))
    allocate(inp_div (inputlength))
    allocate(out_scale (outputlength))
  
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: FKB is configured with ', rk, ' real number and', ik, ' integer.'
    endif
  
    ! 1. Read sub
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading inp_sub'
    endif
    open (unit=555,file='./keras_matrices/inp_sub.txt',status='old',action='read')
    read(555,*) inp_sub(:)
    close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_sub = ',inp_sub
    endif
#endif

    ! 2. Read div
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading inp_div'
    endif
    open (unit=555,file='./keras_matrices/inp_div.txt',status='old',action='read')
    read(555,*) inp_div(:)
    close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_div = ',inp_div
    endif
#endif

    ! 3. Read out_scale
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading out_scale'
    endif
    open (unit=555,file='./keras_matrices/out_scale.txt',status='old',action='read')
    read(555,*) out_scale(:)
    close (555)
#ifdef BRAINDEBUG
      if (masterproc) then
        write (6,*) 'BRAINDEBUG out_scale = ',out_scale
      endif
#endif
  end subroutine init_keras_norm

  subroutine init_nn_vectors()
    allocate(input (inputlength))
    allocate(output (outputlength))
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: nn_in_out_vars: ', trim(nn_in_out_vars)
      write (6,*) 'CLOUDBRAIN: allocate input vector: ', inputlength
      write (6,*) 'CLOUDBRAIN: allocate output vector: ', outputlength
    end if
    ! for previous tendencies
    allocate(dtdt_m1 (begchunk:endchunk,pcols,pver))
    allocate(dqdt_m1 (begchunk:endchunk,pcols,pver))
  end subroutine init_nn_vectors

end module cloudbrain
#endif
