#include <misc.h>
#include <params.h>
#define BRAINDEBUG
#define NEURALLIB
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
  integer, parameter :: inputlength = 94
  integer, parameter :: outputlength = 65

#ifdef NEURALLIB
!  type(network_type) :: cloudbrain_net
  type(network_type) :: cloudbrain_net(65)
  integer(ik) :: fileunit, num_layers
  integer(ik) :: n
#endif

  real :: inp_sub(inputlength)
  real :: inp_div(inputlength)
  real :: out_scale(outputlength) ! HEY needs updating.

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
#ifdef NEURALLIB
    real(rk) :: input(inputlength)
#endif
    real(r8) :: output (outputlength)
    integer :: k, nlev, n,count
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
! HEY we are not sure of the input order
! HEY the input normalization files look suspicious.
#ifdef NEURALLIB
    input(1:nlev) = TBP(:)
    input((nlev+1):2*nlev) = QBP(:)
#endif
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
#ifdef NEURALLIB
    ! use neural fortran library
   do k=1,outputlength
     output = cloudbrain_net(k) % output(input) ! note coupling to many NNs here, one for each output
   end do
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
#ifdef NEURALLIB
    TPHYSTND(:) =      output(1:nlev) * cpair! JORDAN SWAPPED PHQ(:)
    PHQ(:) = output((nlev+1):2*nlev)! This is still the wrong unit, needs to be converted to W/m^2
#else
    PHQ(:) = output(1:nlev) 
    TPHYSTND(:) = output((nlev+1):2*nlev) * cpair
#endif
    FSNT =        output(2*nlev+1)
    FSNS =        output(2*nlev+2)
    FLNT =        output(2*nlev+3)
    FLNS =        output(2*nlev+4)
    PRECT =       output(2*nlev+5)

  end subroutine neural_net

  subroutine init_keras_matrices()    
#ifdef NEURALLIB
! HEY manage the enumerated filenames here:
! Upgrade into a loop over all the laods
 ! LOGIC FOR LOADING ALL THE SEPARATE NNS IS HERE  
  integer :: k,count
  integer :: kvar,klev
  character (256) :: tmpstr
  character (256) :: kvarstr,klevstr
  count = 0
  ! 2x 1x30 output variables first...
  do kvar=1,2 
    do klev = 1,30
      count = count + 1
      write (kvarstr,"(I0)") kvar-1 ! -1 because 0 based counting
      write (klevstr,"(I0)") klev-1
      tmpstr = trim(kvarstr)//'_'//trim(klevstr)//'_model.txt'
      write (6,*) 'Attempting to load NN for: ',trim(tmpstr)
      call cloudbrain_net(count) % load('/home1/08098/tg873976/usmile/causality_convection/dummy_singleNNs_FKB_renamed/'//trim(tmpstr))
    end do
  end do
  ! 5 x scalar output variables next.
  do kvar=3,7
    write (kvarstr,"(I0)") kvar-1 ! -1 because 0 based counting
    write (klevstr,"(I0)") 0  ! 0 is the convention used here for the scalars
    tmpstr = trim(kvarstr)//'_'//trim(klevstr)//'_model.txt'
    write (6,*) 'Attempting to load NN for: ',trim(tmpstr)
    call cloudbrain_net(count+kvar) % load('/home1/08098/tg873976/usmile/causality_convection/dummy_singleNNs_FKB_renamed/'//trim(tmpstr))
  end do
#endif
! TODO: Construct the input causal 1/0 matrix (noutput)
  end subroutine init_keras_matrices
    
subroutine init_keras_norm()
  ! 1. Read sub
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading inp_sub'
  endif
  open(unit=555,file='/home1/08098/tg873976/usmile/causality_convection/dummy_singleNNs_FKB_renamed/0_24_inp_sub.txt',status='old',action='read')
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
  open(unit=555,file='/home1/08098/tg873976/usmile/causality_convection/dummy_singleNNs_FKB_renamed/0_24_inp_div.txt',status='old',action='read')
  read(555,*) inp_div(:)
  close (555)
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG inp_div = ',inp_div
    endif
#endif

  ! HEY still need to think about output scaling.
  ! 3. Read out_scale
  if (masterproc) then
    write (6,*) 'CLOUDBRAIN: reading out_scale'
  endif
!  open (unit=555,file='./keras_matrices/out_scale.txt',status='old',action='read')
!  read(555,*) out_scale(:)
!  close (555)
   out_scale(:) = 1.
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'BRAINDEBUG out_scale = ',out_scale
    endif
#endif

  end subroutine init_keras_norm

end module cloudbrain
