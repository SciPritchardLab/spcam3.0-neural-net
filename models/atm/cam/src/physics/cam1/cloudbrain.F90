#ifdef CLOUDBRAIN

#include <misc.h>
#include <params.h>
!#define BRAINDEBUG

module cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use pmgrid,          only: masterproc
!use runtime_opts, only: nn_nint, inputlength, outputlength, activation_type, width

! -------- PYTORCH BINDING --------
use torch_ftn
use iso_fortran_env
! ---------------------------------

  implicit none
  save 

  private
  ! time step at which to couple to NN (to allow SP to spin up), this is atm_in namelist variable
  integer :: nstepNN = 48  ! 48 for turning on NN from Day 3 (when dtime=1800).
                           ! nstep starts from 0.

  ! Define variables for this entire module
  integer, parameter :: inputlength = 64
  integer, parameter :: outputlength = 60

  real(r8) :: inp_sub(inputlength)
  real(r8) :: inp_div(inputlength)
  real(r8) :: out_scale(outputlength)

  ! -------- PYTORCH BINDING --------
  type(torch_module) :: torch_mod
  ! ---------------------------------

  public neural_net, init_ml_pytorch_model, init_ml_norm, nstepNN 

  contains

  subroutine neural_net (TBP, QBP, PS, SOLIN, SHFLX, LHFLX, &
                         PHQ, TPHYSTND, &                       ! only 2 output vars
                         ! FSNT, FSNS, FLNT, FLNS, PRECT, &     ! 5 extra output vars (for old PNAS version)
                         icol)
    ! PNAS version: First row = inputs, second row = outputs
    ! icol is used for debugging to only output one colum
    ! Allocate inputs
    real(r8), intent(in) :: TBP(:) ! note QBP could be RH #ifdef RHNN (see calling in tphysbc_internallythreaded)
    real(r8), intent(in) :: QBP(:)   
    real(r8), intent(in) :: LHFLX
    real(r8), intent(in) :: SHFLX
    real(r8), intent(in) :: PS
    real(r8), intent(in) :: SOLIN
    ! Allocate outputs
    real(r8), intent(out) :: PHQ(:)
    real(r8), intent(out) :: TPHYSTND(:)
    ! real(r8), intent(out) :: FSNT
    ! real(r8), intent(out) :: FSNS
    ! real(r8), intent(out) :: FLNT
    ! real(r8), intent(out) :: FLNS
    ! real(r8), intent(out) :: PRECT
    ! Allocate utilities

    integer :: k, nlev, n
    integer, intent(in) :: icol

    ! -------- PYTORCH BINDING --------
    type(torch_tensor_wrap) :: input_tensors
    type(torch_tensor) :: out_tensor
    real(real32) :: input(inputlength,1)     ! real32 is the precision used in torch_ftn 
    real(r8), pointer :: output (:,:)   ! double precision (e.g. real64, r8) is not supported yet
    ! ---------------------------------

    ! 1. Concatenate input vector to neural network
    nlev=30
    input(1:nlev,1) = TBP(:) 
    input((nlev+1):2*nlev,1) = QBP(:) 
    input(2*nlev+1,1) = PS
    input(2*nlev+2,1) = SOLIN
    input(2*nlev+3,1) = SHFLX
    input(2*nlev+4,1) = LHFLX
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input pre norm=',input
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(k,1) = (input(k,1) - inp_sub(k))/inp_div(k)
    end do
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG input post norm=',input
      endif
#endif

! 3. Neural network matrix multiplications and activations
    ! -------- PYTORCH BINDING --------
    call input_tensors%create
    call input_tensors%add_array(input)
    call torch_mod%forward(input_tensors, out_tensor)
    call out_tensor%to_array(output)
    ! ---------------------------------
#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG output = ',output
      endif
#endif

    ! 4. Unnormalize output
    do k=1,outputlength
      output(k,1) = output(k,1) / out_scale(k)
    end do

#ifdef BRAINDEBUG
      if (masterproc .and. icol .eq. 1) then
        write (6,*) 'BRAINDEBUG out post scale = ',output
      endif
#endif

    ! 5. Split output into components
    TPHYSTND(:) =      output(1:nlev,1) * cpair! JORDAN SWAPPED PHQ(:)
    PHQ(:) = output((nlev+1):2*nlev,1)! This is still the wrong unit, needs to be converted to W/m^2
    ! FSNT =        output(2*nlev+1)
    ! FSNS =        output(2*nlev+2)
    ! FLNT =        output(2*nlev+3)
    ! FLNS =        output(2*nlev+4)
    ! PRECT =       output(2*nlev+5)

  end subroutine neural_net


! -------- PYTORCH BINDING --------
  subroutine init_ml_pytorch_model()
    call torch_mod%load('./pytorch_files/model.pt') ! must be a *traced* model
    write (6,*) 'CLOUDBRAIN: loading a pytorch model -------'
  end subroutine init_ml_pytorch_model
! ---------------------------------
    

subroutine init_ml_norm()

  ! 1. Read sub
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_sub'
  endif
  open (unit=555,file='./pytorch_files/inp_sub.txt',status='old',action='read')
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
  open (unit=555,file='./pytorch_files/inp_div.txt',status='old',action='read')
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
  open (unit=555,file='./pytorch_files/out_scale.txt',status='old',action='read')
  read(555,*) out_scale(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG out_scale = ',out_scale
    endif
#endif

  end subroutine init_ml_norm

end module cloudbrain
#endif
