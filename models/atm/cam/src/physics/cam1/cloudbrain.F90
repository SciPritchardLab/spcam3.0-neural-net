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
  integer, parameter :: inputlength = 64
  integer, parameter :: outputlength = 60
  integer, parameter :: activation_type = 1
  integer, parameter :: width = 256

#ifdef ENSEMBLE
  real(rk) :: noise = 0.0
  type(ensemble_type) :: cloudbrain_ensemble
#else
  type(network_type) :: cloudbrain_net
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n
#endif

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_scale(outputlength)

  public neural_net, init_keras_matrices, init_keras_norm
  contains

  subroutine neural_net (TBP, QBP, PS, SOLIN, SHFLX, LHFLX, &
                         PHQ, TPHYSTND, FSNT, FSNS, FLNT, FLNS, PRECT, &
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
    real(r8), intent(out) :: FSNT
    real(r8), intent(out) :: FSNS
    real(r8), intent(out) :: FLNT
    real(r8), intent(out) :: FLNS
    real(r8), intent(out) :: PRECT
    ! Allocate utilities

    real(rk) :: input(inputlength),x1(width), x2(width)
    real(r8) :: output (outputlength)
    integer :: k, nlev, n
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
    input(1:nlev) = TBP(:) 
    input((nlev+1):2*nlev) = QBP(:) 
    input(2*nlev+1) = PS
    input(2*nlev+2) = SOLIN
    input(2*nlev+3) = SHFLX
    input(2*nlev+4) = LHFLX
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
    TPHYSTND(:) =      output(1:nlev) * cpair! JORDAN SWAPPED PHQ(:)
    PHQ(:) = output((nlev+1):2*nlev)! This is still the wrong unit, needs to be converted to W/m^2
    ! FSNT =        output(2*nlev+1)
    ! FSNS =        output(2*nlev+2)
    ! FLNT =        output(2*nlev+3)
    ! FLNS =        output(2*nlev+4)
    ! PRECT =       output(2*nlev+5)

  end subroutine neural_net



  subroutine init_keras_matrices()    
#ifdef ENSEMBLE
    write (6,*) '------- NEURAL-FORTRAN: ensemble loading -------'
    cloudbrain_ensemble = ensemble_type('./Models/', noise)
    write (6,*) '------- NEURAL-FORTRAN: ensemble loaded -------'
#else
    call cloudbrain_net % load('./keras_matrices/model.txt')
    write (6,*) '------- NEURAL-FORTRAN: loaded network from txt file -------'
#endif

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

  end subroutine init_keras_norm

end module cloudbrain
#endif
