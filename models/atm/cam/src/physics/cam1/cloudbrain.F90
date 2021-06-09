#include <misc.h>
#include <params.h>
#define BRAINDEBUG
#define NEURALLIB
MODULE cloudbrain
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair,latvap,latice
use pmgrid, only: masterproc

! -------- NEURAL-FORTRAN --------
! imports
use mod_kinds, only: ik, rk
use mod_network , only: network_type
!use mod_ensemble, only: ensemble_type
! --------------------------------

IMPLICIT NONE
SAVE 

PRIVATE
! Define variables for this entire module
INTEGER, PARAMETER :: inputlength  = 94
INTEGER, PARAMETER :: outputlength = 65

#ifdef NEURALLIB
!  type(network_type) :: cloudbrain_net
! 210609 FIS. Hey this needs work here. It runs but will not work for the
!             CausalSingleNNs. It only "picks up" the input list from the
!             first SingleNN (i.e., 0_0_input_list.txt), which will be
!             applied to the rest.
!             Solutions?
!               - Elegant. Create new type in FKB? e.g., similar to ensemble_type
!                 that will allow the proper initialization for each
!                 singlenn.
!               - Ugly (but quick). Hard code each singlenn, i.e., 
!                 type(network_type) :: 0_0_cloudbrain_net
!                 ...
TYPE(network_type) :: cloudbrain_net(outputlength)
! 210609
INTEGER(ik) :: fileunit, num_layers
INTEGER(ik) :: n
#endif

INTEGER :: inputmask_causal_relevance(outputlength,inputlength)
REAL :: inp_sub(inputlength)
REAL :: inp_div(inputlength)
REAL :: out_scale(outputlength) ! HEY needs updating; Fine already?

PUBLIC neural_net, init_keras_matrices, init_keras_norm
CONTAINS

  SUBROUTINE neural_net (QBP, TBP, VBP, PS, SOLIN, SHFLX, LHFLX, &
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
    real(rk) :: input(inputlength),input_trimmed(inputlength)
#endif
    real(r8) :: output (outputlength),tmpoutput(outputlength)
    integer :: j,k, nlev, n,nrelevant
    integer, intent(in) :: icol

    ! 1. Concatenate input vector to neural network
    nlev=30
    input(1:nlev)            = QBP(:)
    input((nlev+1):2*nlev)   = TBP(:)
    input((2*nlev+1):3*nlev) = VBP(:)
    input(3*nlev+1)          = PS
    input(3*nlev+2)          = SOLIN
    input(3*nlev+3)          = SHFLX
    input(3*nlev+4)          = LHFLX
#ifdef BRAINDEBUG
    IF (masterproc .AND. icol .EQ. 1) THEN
      WRITE (6,*) 'BRAINDEBUG input pre norm=',input
    ENDIF
#endif

    ! 2. Normalize input
    DO k=1,inputlength
      input(k) = (input(k) - inp_sub(k))/inp_div(k)
    END DO
#ifdef BRAINDEBUG
    IF (masterproc .AND. icol .EQ. 1) THEN
      WRITE (6,*) 'BRAINDEBUG input post norm=',input
    ENDIF
#endif

    ! 3. Neural network matrix multiplications and activations
#ifdef NEURALLIB
    ! use neural fortran library
    DO k=1,outputlength
      nrelevant = 0
      DO j=1,inputlength
        IF (inputmask_causal_relevance(k,j) .EQ. 1) THEN ! is relevant input?
          nrelevant = nrelevant + 1
          input_trimmed(nrelevant) = input(j)
        ENDIF
      END DO
! Debug point: print after the trimming of input_trimmed(1:nrelevant)
#ifdef BRAINDEBUG
      IF (masterproc .AND. icol .EQ. 1 .AND. k .EQ. 1) THEN
        WRITE (6,*) 'Length of the input vector=',inputlength
        WRITE (6,*) 'input=',input(:)
        WRITE (6,*) 'Length of the trimmed input vector=',nrelevant
        WRITE (6,*) 'input_trimmed=',input_trimmed(1:nrelevant)
      ENDIF
#endif
      ! note coupling to many NNs here, one for each output
      tmpoutput = cloudbrain_net(k) % output(input_trimmed(1:nrelevant))
!      tmpoutput = cloudbrain_net(k) % output(input)
      output(k) = tmpoutput(1)
    END DO
#endif

#ifdef BRAINDEBUG
    IF (masterproc .AND. icol .EQ. 1) THEN
      WRITE (6,*) 'BRAINDEBUG output = ',output
    ENDIF
#endif

    ! 4. Unnormalize output
    DO k=1,outputlength
      output(k) = output(k) / out_scale(k)
    END DO

#ifdef BRAINDEBUG
    IF (masterproc .AND. icol .EQ. 1) THEN
      WRITE (6,*) 'BRAINDEBUG out post scale = ',output
    ENDIF
#endif

    ! 5. Split output into components
    PHQ(:)      = output(1:nlev) 
    TPHYSTND(:) = output((nlev+1):2*nlev) * cpair
    FSNT        = output(2*nlev+1)
    FSNS        = output(2*nlev+2)
    FLNT        = output(2*nlev+3)
    FLNS        = output(2*nlev+4)
    PRECT       = output(2*nlev+5)
    
  END SUBROUTINE neural_net

  SUBROUTINE init_keras_matrices()    
#ifdef NEURALLIB
! HEY manage the enumerated filenames here:
! Upgrade into a loop over all the loads
! LOGIC FOR LOADING ALL THE SEPARATE NNS IS HERE  
    INTEGER :: k,count
    INTEGER :: kvar,klev
    CHARACTER (256) :: tmpstr
    CHARACTER (256) :: kvarstr,klevstr
    count = 0
    ! 2x 1x30 output variables first...
    DO kvar=1,2 
      DO klev = 1,30
        count = count + 1
        WRITE (kvarstr,"(I0)") kvar-1 ! -1 because 0 based counting
        WRITE (klevstr,"(I0)") klev-1
        tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_model.txt'
        WRITE (6,*) 'Attempting to load NN for: ',TRIM(tmpstr)
        CALL cloudbrain_net(count) % load('./models/'//TRIM(tmpstr))
        ! INSERT TODO: read in the 1/0 files here, leveraging same string
        ! handling, then store in a outputlength x
        ! inputlength array to be added to module private contents. THEN TODO --
        ! subset the input vector in the call to the cloudbrain from the main NN
        ! subroutine. 
        tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_input_list.txt'
        OPEN(unit=555,file='./models/'//TRIM(tmpstr),status='old',action='read')
        READ(555,*) inputmask_causal_relevance(count,:)
        tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_out_scale.txt'
        OPEN(unit=555,file='./models/'//TRIM(tmpstr),status='old',action='read')
        READ(555,*) out_scale(count)
      END DO
    END DO
    ! 5 x scalar output variables next.
    DO kvar=3,7
      count = count + 1
      WRITE (kvarstr,"(I0)") kvar-1 ! -1 because 0 based counting
      WRITE (klevstr,"(I0)") 0  ! 0 is the convention used here for the scalars
      tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_model.txt'
      WRITE (6,*) 'Attempting to load NN for: ',TRIM(tmpstr)
      CALL cloudbrain_net(count) % load('./models/'//TRIM(tmpstr))
      tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_input_list.txt'
      OPEN(unit=555,file='./models/'//TRIM(tmpstr),status='old',action='read')
      READ(555,*) inputmask_causal_relevance(count,:)
      tmpstr = TRIM(kvarstr)//'_'//TRIM(klevstr)//'_out_scale.txt'
      OPEN(unit=555,file='./models/'//TRIM(tmpstr),status='old',action='read')
      READ(555,*) out_scale(count)
    END DO
#endif
  END SUBROUTINE init_keras_matrices

  SUBROUTINE init_keras_norm()
    ! 1. Read sub
    IF (masterproc) THEN
      WRITE (6,*) 'CLOUDBRAIN: reading inp_sub'
    ENDIF
    OPEN(unit=555,file='./models/0_0_inp_sub.txt',status='old',action='read')
    READ(555,*) inp_sub(:)
    CLOSE (555)
#ifdef BRAINDEBUG
    IF (masterproc) THEN
      WRITE (6,*) 'BRAINDEBUG inp_sub = ',inp_sub
    ENDIF
#endif
    
    ! 2. Read div
    IF (masterproc) THEN
      WRITE (6,*) 'CLOUDBRAIN: reading inp_div'
    ENDIF
    OPEN(unit=555,file='./models/0_0_inp_div.txt',status='old',action='read')
    READ(555,*) inp_div(:)
    CLOSE (555)
#ifdef BRAINDEBUG
    IF (masterproc) THEN
      WRITE (6,*) 'BRAINDEBUG inp_div = ',inp_div
    ENDIF
#endif
  END SUBROUTINE init_keras_norm
  
END MODULE cloudbrain
