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
  character(len=256) :: nn_in_out_vars1 = 'in_TBP_QBP_PS_SOLIN_SHF_LHF_out_TPHYSTND', &
                        nn_in_out_vars2 = 'in_TBP_QBP_PS_SOLIN_SHF_LHF_out_PHQ'
                        ! A userdefined name for a specific input/output combination
  integer :: inputlength1  = 64 ! length of input vector for NN
  integer :: inputlength2  = 64
  integer :: outputlength1 = 30 ! length of output vector for NN
  integer :: outputlength2 = 30
  !!!!

#ifdef ENSEMBLE
  real(rk) :: noise = 0.0
  type(ensemble_type) :: cloudbrain_ensemble
#else
  type(network_type) :: cloudbrain_net1, cloudbrain_net2
  integer(ik) :: n
#endif

  real(rk), allocatable :: inp_sub1(:),   inp_sub2(:)
  real(rk), allocatable :: inp_div1(:),   inp_div2(:)
  real(rk), allocatable :: out_scale1(:), out_scale2(:)
  real(rk), allocatable :: input1(:),     input2(:)
  real(r8), allocatable :: output1(:),    output2(:)

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
         nn_in_t, nn_out_t, nn_in_out_vars1, nn_in_out_vars2, &
         inputlength1, outputlength1, inputlength2, outputlength2, &
         init_nn_vectors, &
         dtdt_m1, dqdt_m1 ! buffer variables for previous timestep tendencies

  contains

  subroutine neural_net (nn_in, nn_out, icol, which_nn)
    type(nn_in_t), intent(in) :: nn_in
    type(nn_out_t), intent(out) :: nn_out
    integer, intent(in) :: icol
    integer, intent(in) :: which_nn
    
    ! local variables
    integer :: k, nlev, n

    ! 1. Concatenate input vector to neural network
    nlev=30
    if (which_nn .eq. 1) then
      select case (to_upper(trim(nn_in_out_vars1)))
        case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_TPHYSTND')
          input1(1:nlev) = nn_in%tbp(:nlev)
          input1((nlev+1):2*nlev) = nn_in%qbp(:nlev)
          input1(2*nlev+1) = nn_in%ps
          input1(2*nlev+2) = nn_in%solin
          input1(2*nlev+3) = nn_in%shf
          input1(2*nlev+4) = nn_in%lhf
      end select
    elseif (which_nn .eq. 2) then
      select case (to_upper(trim(nn_in_out_vars2)))
        case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_PHQ')
          input2(1:nlev) = nn_in%tbp(:nlev)
          input2((nlev+1):2*nlev) = nn_in%qbp(:nlev)
          input2(2*nlev+1) = nn_in%ps
          input2(2*nlev+2) = nn_in%solin
          input2(2*nlev+3) = nn_in%shf
          input2(2*nlev+4) = nn_in%lhf
      end select
    endif

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
    if (which_nn .eq. 1) then
      write (6,*) 'BRAINDEBUG input1 pre norm = ',input1
    elseif (which_nn .eq. 2) then
      write (6,*) 'BRAINDEBUG input2 pre norm = ',input2
    endif
    endif
#endif

    ! 2. Normalize input
    if (which_nn .eq. 1) then
      do k=1,inputlength1
        input1(k) = (input1(k) - inp_sub1(k))/inp_div1(k)
      end do
    elseif (which_nn .eq. 2) then
      do k=1,inputlength2
        input2(k) = (input2(k) - inp_sub2(k))/inp_div2(k)
      end do
    endif

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
    if (which_nn .eq. 1) then
      write (6,*) 'BRAINDEBUG input1 post norm = ',input1
    elseif (which_nn .eq. 2) then
      write (6,*) 'BRAINDEBUG input2 post norm = ',input2
    endif
    endif
#endif

! 3. Neural network matrix multiplications and activations
#ifdef ENSEMBLE
    output1 = cloudbrain_ensemble % average(input1)
#else
    ! use neural fortran library
    if (which_nn .eq. 1) then
      output1 = cloudbrain_net1 % output(input1)
    elseif (which_nn .eq. 2) then
      output2 = cloudbrain_net2 % output(input2)
    endif
#endif

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
    if (which_nn .eq. 1) then
      write (6,*) 'BRAINDEBUG output1 = ',output1
    elseif (which_nn .eq. 2) then
      write (6,*) 'BRAINDEBUG output2 = ',output2
    endif
    endif
#endif

    ! 4. Unnormalize output
    if (which_nn .eq. 1) then
      do k=1,outputlength1
        output1(k) = output1(k) / out_scale1(k)
      end do
    elseif (which_nn .eq. 2) then
      do k=1,outputlength2
        output2(k) = output2(k) / out_scale2(k)
      end do
    endif

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
    if (which_nn .eq. 1) then
      write (6,*) 'BRAINDEBUG output1 post scale = ',output1
    elseif (which_nn .eq. 2) then
      write (6,*) 'BRAINDEBUG output2 post scale = ',output2
    endif
    endif
#endif

    ! 5. Split output into components
    if (which_nn .eq. 1) then
      select case (to_upper(trim(nn_in_out_vars1)))
        case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_TPHYSTND')
          nn_out%tphystnd(:nlev) = output1(1:nlev)
      end select
    elseif (which_nn .eq. 2) then
      select case (to_upper(trim(nn_in_out_vars2)))
        case('IN_TBP_QBP_PS_SOLIN_SHF_LHF_OUT_PHQ')
          nn_out%phq(:nlev) = output2(1:nlev)
      end select
    endif

  end subroutine neural_net


  subroutine init_keras_matrices()
#ifdef ENSEMBLE
    write (6,*) '------- NEURAL-FORTRAN: ensemble loading -------'
    cloudbrain_ensemble = ensemble_type('./Models/', noise)
    write (6,*) '------- NEURAL-FORTRAN: ensemble loaded -------'
#else
    call cloudbrain_net1 % load('./keras_matrices/model1.txt')
    call cloudbrain_net2 % load('./keras_matrices/model2.txt')
    if (masterproc) then
      write (6,*) '------- NEURAL-FORTRAN: loaded network from txt file -------'
    end if
#endif
  end subroutine init_keras_matrices
    

  subroutine init_keras_norm()
    allocate(inp_sub1 (inputlength1))
    allocate(inp_div1 (inputlength1))
    allocate(out_scale1 (outputlength1))
    allocate(inp_sub2 (inputlength2))
    allocate(inp_div2 (inputlength2))
    allocate(out_scale2 (outputlength2))
  
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: FKB is configured with ', rk, ' real number and', ik, ' integer.'
    endif
  
    ! 1. Read sub
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading inp_sub1 and inp_sub2'
    endif
    open (unit=555,file='./keras_matrices/inp_sub1.txt',status='old',action='read')
    read(555,*) inp_sub1(:)
    close (555)
    open (unit=555,file='./keras_matrices/inp_sub2.txt',status='old',action='read')
    read(555,*) inp_sub2(:)
    close (555)
#ifdef BRAINDEBUG
    if (masterproc) then 
      write (6,*) 'BRAINDEBUG inp_sub1 = ',inp_sub1
      write (6,*) 'BRAINDEBUG inp_sub2 = ',inp_sub2
    endif
#endif

    ! 2. Read div
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading inp_div1 and inp_div2'
    endif
    open (unit=555,file='./keras_matrices/inp_div1.txt',status='old',action='read')
    read(555,*) inp_div1(:)
    close (555)
    open (unit=555,file='./keras_matrices/inp_div2.txt',status='old',action='read')
    read(555,*) inp_div2(:)
    close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_div1 = ',inp_div1
      write (6,*) 'BRAINDEBUG inp_div2 = ',inp_div2
    endif
#endif

    ! 3. Read out_scale
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: reading out_scale1 and out_scale2'
    endif
    open (unit=555,file='./keras_matrices/out_scale1.txt',status='old',action='read')
    read(555,*) out_scale1(:)
    close (555)
    open (unit=555,file='./keras_matrices/out_scale2.txt',status='old',action='read')
    read(555,*) out_scale2(:)
    close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG out_scale1 = ',out_scale1
      write (6,*) 'BRAINDEBUG out_scale2 = ',out_scale2
    endif
#endif
  end subroutine init_keras_norm

  subroutine init_nn_vectors()
    allocate(input1 (inputlength1))
    allocate(output1 (outputlength1))
    allocate(input2 (inputlength2))
    allocate(output2 (outputlength2))
    if (masterproc) then
      write (6,*) 'CLOUDBRAIN: nn_in_out_vars1: ', trim(nn_in_out_vars1)
      write (6,*) 'CLOUDBRAIN: nn_in_out_vars2: ', trim(nn_in_out_vars2)
      write (6,*) 'CLOUDBRAIN: allocate input vector1,2:  ', inputlength1, inputlength2
      write (6,*) 'CLOUDBRAIN: allocate output vector1,2: ', outputlength1, outputlength2
    end if
    ! for previous tendencies
    allocate(dtdt_m1 (begchunk:endchunk,pcols,pver))
    allocate(dqdt_m1 (begchunk:endchunk,pcols,pver))
  end subroutine init_nn_vectors

end module cloudbrain
#endif
