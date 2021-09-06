#ifdef CLOUDBRAIN

#include <misc.h>
#include <params.h>
!#define BRAINDEBUG

module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use pmgrid, only: masterproc
!use runtime_opts, only: nn_nint, inputlength, outputlength, activation_type, width

! -------- NEURAL-FORTRAN --------
! imports
use mod_kinds, only: ik, rk
use mod_network , only: network_type
use mod_ensemble, only: ensemble_type
! --------------------------------

  implicit none
  save 

  private
  ! Define variables for this entire module
  integer, parameter :: nn_nint = 8
  integer, parameter :: inputlength = 94
  integer, parameter :: outputlength = 65
  integer, parameter :: activation_type = 1
  integer, parameter :: width = 256

#ifdef ENSEMBLE
  real(rk) :: noise = 0.0
  type(ensemble_type) :: cloudbrain_ensemble
#else
  type(network_type) :: cloudbrain_net
#ifdef NNBIASCORRECTOR
  type(network_type) :: corrector_net
  integer, parameter :: corrector_inputlength = 155
  integer, parameter :: corrector_outputlength = 60
#endif
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n
#endif

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_scale(outputlength)

#ifdef NNBIASCORRECTOR
  real :: corrector_inp_sub(inputlength)
  real :: corrector_inp_div(inputlength)
  public corrector_neural_net
#endif
  public neural_net, init_keras_matrices, init_keras_norm
  contains

  subroutine neural_net (QBP, TBP, VBP, PS, SOLIN, SHFLX, LHFLX, &
                         PHQ, TPHYSTND, FSNT, FSNS, FLNT, FLNS, PRECT, &
                         icol)
    ! PNAS version: First row = inputs, second row = outputs
    ! icol is used for debugging to only output one colum
    ! Allocate inputs
    real(r8), intent(in) :: QBP(:)
    real(r8), intent(in) :: TBP(:)   
    real(r8), intent(in) :: VBP(:)
    real(r8), intent(in) :: PS
    real(r8), intent(in) :: SOLIN
    real(r8), intent(in) :: SHFLX
    real(r8), intent(in) :: LHFLX
    ! Allocate outputs
    real(r8), intent(out) :: PHQ(:)
    real(r8), intent(out) :: TPHYSTND(:)
    real(r8), intent(out) :: FSNT
    real(r8), intent(out) :: FSNS
    real(r8), intent(out) :: FLNT
    real(r8), intent(out) :: FLNS
    real(r8), intent(out) :: PRECT
    ! Allocate utilities
    real(rk) :: input(inputlength),x1(width), x2(width)
    real(rk) :: output (outputlength)
    integer :: k, nlev, n
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
    input(1:nlev) = TBP(:) 
    input((nlev+1):2*nlev) = QBP(:) 
    input((2*nlev+1):3*nlev)=VBP(:)
    input(3*nlev+1) = PS
    input(3*nlev+2) = SOLIN
    input(3*nlev+3) = SHFLX
    input(3*nlev+4) = LHFLX
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


    TPHYSTND(:) =      output(1:nlev) * cpair! JORDAN SWAPPED PHQ(:)
    PHQ(:) = output((nlev+1):2*nlev)! This is still the wrong unit, needs to be converted to W/m^2
    FSNT =        output(2*nlev+1)
    FSNS =        output(2*nlev+2)
    FLNT =        output(2*nlev+3)
    FLNS =        output(2*nlev+4)
    PRECT =       output(2*nlev+5)

  end subroutine neural_net

subroutine corrector_neural_net (QBP, TBP, VBP, PS, SOLIN, SHFLX, LHFLX, &
                         PHQ, TPHYSTND, icol)
                         
    ! icol is used for debugging to only output one colum
    ! Allocate inputs
    real(r8), intent(in) :: QBP(:)
    real(r8), intent(in) :: TBP(:)   
    real(r8), intent(in) :: VBP(:)
    real(r8), intent(in) :: PS
    real(r8), intent(in) :: SOLIN
    real(r8), intent(in) :: SHFLX
    real(r8), intent(in) :: LHFLX
    ! Allocate outputs
    real(r8), intent(inout) :: PHQ(:)  ! note these are also inputs for corrector containing NNDT and NNDQ from NN#1
    real(r8), intent(inout) :: TPHYSTND(:)
    ! Allocate utilities
    real(rk) :: input(corrector_inputlength)
    real(rk) :: output (outputlength)
    integer :: k, nlev, n
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
    ! Jerry says input order is (Slack, 9/3/21):
    !NNDT (30) NNDQ (30) NNLHF (1) NNPS (1) NNQBP (30) NNSHF (1) NNSOLIN (1) NNTBP (30) NNVBP (30)
    input(1:nlev) = TPHYSTND(:)/cpair ! NNDT
    input((nlev+1):2*nlev) = PHQ(:) ! NNDQ
    input(2*nlev + 1) = LHFLX
    input(2*nlev + 2) = PS
    n = 2*nlev+2
    input((n+1):(n+nlev)) = QBP(:)
    input(n+nlev+1) = SHFLX
    input(n+nlev+2) = SOLIN
    n = n+nlev+2
    input((n+1):(n+nlev)) = TBP(:)
    input((n+nlev+1):(n+2*nlev)) = VBP(:)
    input(corrector_inputlength) = 1. ! Jerry says imoportant to pad by 1 for y-intercept (need to understand, was this importantly missed in prior work??)
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG CORRECTOR input pre norm=',input
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(k) = (input(k) - corrector_inp_sub(k))/corrector_inp_div(k)
    end do
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG CORRECTOR input post norm=',input
      endif
#endif

! 3. Neural network matrix multiplications and activations

    output = corrector_net % output(input) ! FKB.

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG CORRECTOR output = ',output
      endif
#endif

    ! 4. Unnormalize output (not needed for corrector NN, outputs are in W/kg from training env)

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG CORRECTOR output (W/kg) = ',output
      endif
#endif

! Now overwrite the input heating rates, NNDT and NNDQ, with the bias-corrected equivalents, 
! which will be received by the GCM that called this subroutine:

    TPHYSTND(:) =      output(1:nlev) ! no cpair multiplication here, already in W/kg
    PHQ(:) = output((nlev+1):2*nlev)/latvap ! W/kg --> kg/kg/s

  end subroutine corrector_neural_net

  subroutine init_keras_matrices()    
#ifdef ENSEMBLE
    write (6,*) '------- NEURAL-FORTRAN: ensemble loading -------'
    cloudbrain_ensemble = ensemble_type('./Models/', noise)
    write (6,*) '------- NEURAL-FORTRAN: ensemble loaded -------'
#else
    call cloudbrain_net % load('./keras_matrices/model.txt')
    write (6,*) '------- NEURAL-FORTRAN: loaded network from txt file -------'
#endif

#ifdef NNBIASCORRECTOR
    call corrector_net % load('./keras_matrices/bias_corrector/model.txt')  
#endif ! ENSEMBLE

  end subroutine init_keras_matrices
    
subroutine init_keras_norm()

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

#ifdef NNBIASCORRECTOR
 ! 1. Read sub
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_sub'
  endif
  open (unit=555,file='./keras_matrices/bias_corrector/inp_sub.txt',status='old',action='read')
  read(555,*) corrector_inp_sub(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG CORRECTOR inp_sub = ',corrector_inp_sub
    endif
#endif

  ! 2. Read div
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_div'
  endif
  open (unit=555,file='./keras_matrices/bias_corrector/inp_div.txt',status='old',action='read')
  read(555,*) corrector_inp_div(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG CORRECTOR inp_div = ',corrector_inp_div
    endif
#endif
#endif ! NNBIASCORRECTOR
  end subroutine init_keras_norm

end module cloudbrain
#endif
