#include <misc.h>
#include <params.h>
#define BRAINDEBUG
! Limit output to min-max
!#define LIMITOUTP
! Limit input to min-max
!#define INPLIMITER  
!#define NOADIAB
!#define DEEP

module cloudbrain_keras_dense
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair
use pmgrid, only: masterproc

  implicit none
  
  save 

  private
  ! Define the network architectures
  integer, parameter :: nlev = 30
  integer, parameter :: width1 = 512
#ifdef NOADIAB
  integer, parameter :: inputlength = 94  ! For no adiab option. 152 for original purecrm
#else
  integer, parameter :: inputlength = 152
#endif
  integer, parameter :: outputlength = 120
  integer, parameter :: nchunk = 64
  ! 1st: BASE
#ifndef DEEP
  integer, parameter :: width = 512
  !Weight and bias matrices
  real :: bias1(width)
  real :: weights1(width,inputlength)
  real :: bias2(outputlength)
  real :: weights2(outputlength,width)
#else
  ! 2nd: Deep
  integer, parameter :: width = 256
  real :: bias1(width)
  real :: weights1(width,inputlength)
  real :: bias2(width)
  real :: weights2(width,width)
  real :: bias3(width)
  real :: weights3(width,width)
  real :: bias4(width)
  real :: weights4(width,width)
  real :: bias5(width)
  real :: weights5(width,width)
  real :: bias6(width)
  real :: weights6(width,width)
  real :: bias7(width)
  real :: weights7(width,width)
  real :: bias8(width)
  real :: weights8(width,width)
  real :: bias9(width)
  real :: weights9(width,width)
  real :: bias10(outputlength)
  real :: weights10(outputlength,width)
#endif
  real :: input_norm_mean(inputlength)
  real :: input_norm_std(inputlength)
  real :: output_norm_min(outputlength)
  real :: output_norm_max(outputlength)
  real :: input_norm_min(inputlength)
  real :: input_norm_max(inputlength)
  real :: input_norm_max_rs(inputlength)

  public init_keras_norm
#ifdef DEEP
  public init_keras_matrices_deep, cloudbrain_purecrm_deep
#else
  public init_keras_matrices_base, cloudbrain_purecrm_base
#endif

  contains

#ifndef DEEP
  subroutine cloudbrain_purecrm_base (TC, QC, VC, dTdt_adiabatic, dQdt_adiabatic, PS, SOLIN, &
#ifdef NOADIAB
                                      SHFLX, LHFLX, &
#endif
                                      SPDT, SPDQ, QRL, QRS, icol)
    ! NN inputs
    real(r8), intent(in) :: TC(pver)   ! CRM-equivalent T = TAP[t-1] - DTV[t-1]*dt
    real(r8), intent(in) :: QC(pver)   ! QAP[t-1] - VD01[t-1]*dt
    real(r8), intent(in) :: VC(pver)   ! VAP[t-1]
    real(r8), intent(in) :: dTdt_adiabatic(pver) ! TBP[t]/dt - TC/dt
    real(r8), intent(in) :: dQdt_adiabatic(pver) ! QBP[t]/dt - QC/dt
#ifdef NOADIAB
    real(r8), intent(in) :: SHFLX ! Sensible heat flux from previous time step
    real(r8), intent(in) :: LHFLX ! Latent heat flux
#endif
    real(r8), intent(in) :: PS ! From t-1
    real(r8), intent(in) :: SOLIN ! From t
    ! NN outputs
    real(r8), intent(out) :: SPDT(pver) ! W/kg
    real(r8), intent(out) :: SPDQ(pver) ! W/kg
    real(r8), intent(inout) :: QRL(pver) ! W/kg
    real(r8), intent(inout) :: QRS(pver) ! W/kg

    real(r8) :: input(inputlength),x1(width)
    real(r8) :: output (outputlength)
    integer :: k,j,k1,k2
    integer, intent(in) :: icol

!    call init_keras_matrices
#ifdef NOADIAB
    ! Stacking for noadiab option
    k1=pver-nlev+1
    k2=pver
    input(1:nlev)=TC(k1:k2) 
    input((nlev+1):2*nlev)=QC(k1:k2)
    input((2*nlev+1):3*nlev)=VC(k1:k2)
    input(3*nlev+1) = SHFLX
    input(3*nlev+2) = LHFLX
    input(3*nlev+3) = PS
    input(3*nlev+4) = SOLIN
#else
    ! Stack the input variables
    k1=pver-nlev+1
    k2=pver
    input(1:nlev)=TC(k1:k2) 
    input((nlev+1):2*nlev)=QC(k1:k2)
    input((2*nlev+1):3*nlev)=VC(k1:k2)
    input((3*nlev+1):4*nlev) = dTdt_adiabatic(k1:k2)
    input((4*nlev+1):5*nlev) = dQdt_adiabatic(k1:k2)
    input(5*nlev+1) = PS
    input(5*nlev+2) = SOLIN
#endif



#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
      write (6,*) 'HEY input pre norm=',input
    endif
#endif

! Limit inputs to their min max values
#ifdef INPLIMITER
do k=1,inputlength
  input(k) = min(input(k), input_norm_max(k))
  input(k) = max(input(k), input_norm_min(k))
end do
#endif

    ! normalize input:
    do k=1,inputlength
      input(k) = (input(k) - input_norm_mean(k))/input_norm_max_rs(k)
    end do
#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY normalized = ',input
    endif
#endif
#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY weights1 = ',weights1
     write (6,*) 'HEY bias1 = ',bias1
     write (6,*) 'HEY weights2 = ',weights2
     write (6,*) 'HEY bias2 = ',bias2
    endif
#endif
! 1st layer: input length-->512.
    x1(1:width) = 0.
    do k=1,width
      do j=1,inputlength
        x1(k) = x1(k) + weights1(k,j)*input(j)
      end do
      x1(k) = x1(k) + bias1(k)
      !x1(k) = max(0.,x1(k)) ! relu activation.
      x1(k) = max(0.3 * x1(k), x1(k))  ! Leaky ReLU
    end do
! output layer: 512->output length
   output(1:outputlength) = 0.
   do k=1,outputlength
     do j=1,width
       output(k) = output(k) + weights2(k,j)*x1(j) 
     end do
     output(k) = output(k) + bias2(k)
     ! no activation for output.
   end do

#ifdef BRAINDEBUG
   if (masterproc .and. icol .eq. 1) then
    write (6,*) 'HEY output = ',output
   endif
#endif

! SR: Limit outputs to external mins and maxs
#ifdef LIMITOUTP
   do k=1,outputlength
     output(k) = min(output(k), output_norm_max(k))
     output(k) = max(output(k), output_norm_min(k))
   end do
#endif

#ifdef BRAINDEBUG
   if (masterproc .and. icol .eq. 1) then
    write (6,*) 'HEY output limit = ',output
   endif
#endif

! Unstack the output variables and unit convert them back
! ATTENTION: I confusingly, placed SPDQ before SPDT
   SPDQ(:) = 0.   ! If we are predicting all 30 levels, this should be irrelevant, right?
   SPDQ(k1:k2) = output(1:nlev)/2.5e6 ! W/kg --> kg/kg/s
   SPDT(:) = 0.
   SPDT(k1:k2) = output((nlev+1):2*nlev) ! W/kg, is this what CAM wants?
!   QRL(:) = 0. ! retain SP or upstream solution above neural net top. SR: Again, this should be irrelevant now...
   QRL(k1:k2) = output ((2*nlev+1):3*nlev) ! W/kg 
!   QRS(:) = 0. ! retain SP or upstream solution above neural net top.
   QRS(k1:k2) = output ((3*nlev+1):4*nlev) ! W/kg

  end subroutine cloudbrain_purecrm_base


#else


  subroutine cloudbrain_purecrm_deep (TC, QC, VC, dTdt_adiabatic, dQdt_adiabatic, PS, SOLIN, SPDT, SPDQ, QRL, QRS, icol)
  !subroutine cloudbrain_purecrm_base (TC, QC, VC, dTdt_adiabatic, dQdt_adiabatic, PS, SOLIN, SPDT, SPDQ, icol)
    ! NN inputs
    real(r8), intent(in) :: TC(pver)   ! CRM-equivalent T = TAP[t-1] - DTV[t-1]*dt
    real(r8), intent(in) :: QC(pver)   ! QAP[t-1] - VD01[t-1]*dt
    real(r8), intent(in) :: VC(pver)   ! VAP[t-1]
    real(r8), intent(in) :: dTdt_adiabatic(pver) ! TBP[t]/dt - TC/dt
    real(r8), intent(in) :: dQdt_adiabatic(pver) ! QBP[t]/dt - QC/dt
    real(r8), intent(in) :: PS ! From t-1
    real(r8), intent(in) :: SOLIN ! From t
    ! NN outputs
    real(r8), intent(out) :: SPDT(pver) ! W/kg
    real(r8), intent(out) :: SPDQ(pver) ! W/kg
    real(r8), intent(inout) :: QRL(pver) ! W/kg
    real(r8), intent(inout) :: QRS(pver) ! W/kg

    real(r8) :: input(inputlength),x1(width), x2(width)
    real(r8) :: output (outputlength)
    integer :: k,j,k1,k2
    integer, intent(in) :: icol

!    call init_keras_matrices


    ! Stack the input variables
    k1=pver-nlev+1
    k2=pver
    input(1:nlev)=TC(k1:k2) 
    input((nlev+1):2*nlev)=QC(k1:k2)
    input((2*nlev+1):3*nlev)=VC(k1:k2)
    input((3*nlev+1):4*nlev) = dTdt_adiabatic(k1:k2)
    input((4*nlev+1):5*nlev) = dQdt_adiabatic(k1:k2)
    input(5*nlev+1) = PS
    input(5*nlev+2) = SOLIN

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
      write (6,*) 'HEY input pre norm=',input
    endif
#endif

! Limit inputs to their min max values
#ifdef INPLIMITER
do k=1,inputlength
  input(k) = min(input(k), input_norm_max(k))
  input(k) = max(input(k), input_norm_min(k))
end do
#endif

    ! normalize input:
    do k=1,inputlength
      input(k) = (input(k) - input_norm_mean(k))/input_norm_max_rs(k)
    end do
#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY normalized = ',input
    endif
#endif


! 1st layer: input length-->256.
    x1(1:width) = 0.
    do k=1,width
      do j=1,inputlength
        x1(k) = x1(k) + weights1(k,j)*input(j)
      end do
      x1(k) = x1(k) + bias1(k)
      x1(k) = max(0.3 * x1(k), x1(k))  ! Leaky ReLU
    end do

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 1 = ',x1
    endif
#endif

! 2nd layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights2(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias2(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 2 = ',x1
    endif
#endif

! 3rd layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights3(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias3(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 3 = ',x1
    endif
#endif

! 4th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights4(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias4(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 4 = ',x1
    endif
#endif

! 5th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights5(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias5(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 5 = ',x1
    endif
#endif

! 6th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights6(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias6(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 6 = ',x1
    endif
#endif

! 7th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights7(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias7(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 7 = ',x1
    endif
#endif

! 8th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights8(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias8(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 8 = ',x1
    endif
#endif

! 9th layer: 256-->256.
    x2(1:width) = 0.
    do k=1,width
      do j=1,width
        x2(k) = x2(k) + weights9(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias9(k)
      x2(k) = max(0.3 * x2(k), x2(k))  ! Leaky ReLU
    end do
    x1(1:width) = x2(1:width)

#ifdef BRAINDEBUG
    if (masterproc .and. icol .eq. 1) then
     write (6,*) 'HEY Activations 9 = ',x1
    endif
#endif

! output layer: 256->output length
   output(1:outputlength) = 0.
   do k=1,outputlength
     do j=1,width
       output(k) = output(k) + weights10(k,j)*x1(j) 
     end do
     output(k) = output(k) + bias10(k)
     ! no activation for output.
   end do

#ifdef BRAINDEBUG
   if (masterproc .and. icol .eq. 1) then
    write (6,*) 'HEY output = ',output
   endif
#endif

! SR: Limit outputs to external mins and maxs
#ifdef LIMITOUTP
   do k=1,outputlength
     output(k) = min(output(k), output_norm_max(k))
     output(k) = max(output(k), output_norm_min(k))
   end do
#endif

#ifdef BRAINDEBUG
   if (masterproc .and. icol .eq. 1) then
    write (6,*) 'HEY output limit = ',output
   endif
#endif

! Unstack the output variables and unit convert them back
! ATTENTION: I confusingly, placed SPDQ before SPDT
   SPDQ(:) = 0.   ! If we are predicting all 30 levels, this should be irrelevant, right?
   SPDQ(k1:k2) = output(1:nlev)/2.5e6 ! W/kg --> kg/kg/s
   SPDT(:) = 0.
   SPDT(k1:k2) = output((nlev+1):2*nlev) ! W/kg, is this what CAM wants?
!   QRL(:) = 0. ! retain SP or upstream solution above neural net top. SR: Again, this should be irrelevant now...
   QRL(k1:k2) = output ((2*nlev+1):3*nlev) ! W/kg 
!   QRS(:) = 0. ! retain SP or upstream solution above neural net top.
   QRS(k1:k2) = output ((3*nlev+1):4*nlev) ! W/kg

  end subroutine cloudbrain_purecrm_deep

#endif



#ifdef DEEP
  subroutine init_keras_matrices_deep()

  integer :: nline,k,ipt1,ipt2,j,n,ios
! input layer:
 write (6,*) 'SR: reading layer1_bias'
  open (unit=555,file='./keras_matrices/layer1_bias.txt',status='old',action='read',iostat=ios)
  if (ios .ne. 0) then
    write (6,*) 'HEY keras matrices unable to load, abort.'
    stop
  endif
  read(555,*) bias1(1:width)
  close (555)
   write (6,*) 'SR: reading layer1_kernel'
  open (unit=555,file='./keras_matrices/layer1_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights1(k,1:inputlength)
  end do
  close (555)

! ----- layer 2:
 write (6,*) 'SR: reading layer2_bias'
  open (unit=555,file='./keras_matrices/layer2_bias.txt',status='old',action='read')
  read(555,*) bias2(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer2_bias'

 write (6,*) 'SR: reading layer2_kernel'
  open (unit=555,file='./keras_matrices/layer2_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights2(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer2_kernel'

 ! ----- layer 3:
 write (6,*) 'SR: reading layer3_bias'
  open (unit=555,file='./keras_matrices/layer3_bias.txt',status='old',action='read')
  read(555,*) bias3(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer3_bias'

 write (6,*) 'SR: reading layer3_kernel'
  open (unit=555,file='./keras_matrices/layer3_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights3(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer2_kernel'

 ! ----- layer 4:
 write (6,*) 'SR: reading layer4_bias'
  open (unit=555,file='./keras_matrices/layer4_bias.txt',status='old',action='read')
  read(555,*) bias4(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer4_bias'

 write (6,*) 'SR: reading layer4_kernel'
  open (unit=555,file='./keras_matrices/layer4_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights4(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer4_kernel'

 ! ----- layer 5:
 write (6,*) 'SR: reading layer5_bias'
  open (unit=555,file='./keras_matrices/layer5_bias.txt',status='old',action='read')
  read(555,*) bias5(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer5_bias'

 write (6,*) 'SR: reading layer5_kernel'
  open (unit=555,file='./keras_matrices/layer5_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights5(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer5_kernel'

 ! ----- layer :
 write (6,*) 'SR: reading layer6_bias'
  open (unit=555,file='./keras_matrices/layer6_bias.txt',status='old',action='read')
  read(555,*) bias6(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer6_bias'

 write (6,*) 'SR: reading layer6_kernel'
  open (unit=555,file='./keras_matrices/layer6_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights6(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer6_kernel'

 ! ----- layer 7:
 write (6,*) 'SR: reading layer7_bias'
  open (unit=555,file='./keras_matrices/layer7_bias.txt',status='old',action='read')
  read(555,*) bias7(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer7_bias'

 write (6,*) 'SR: reading layer7_kernel'
  open (unit=555,file='./keras_matrices/layer7_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights7(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer7_kernel'

 ! ----- layer8:
 write (6,*) 'SR: reading layer8_bias'
  open (unit=555,file='./keras_matrices/layer8_bias.txt',status='old',action='read')
  read(555,*) bias8(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer8_bias'

 write (6,*) 'SR: reading layer8_kernel'
  open (unit=555,file='./keras_matrices/layer8_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights8(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer8_kernel'

 ! ----- layer 9:
 write (6,*) 'SR: reading layer9_bias'
  open (unit=555,file='./keras_matrices/layer9_bias.txt',status='old',action='read')
  read(555,*) bias9(1:width)
  close (555)
 write (6,*) 'SR: finished reading layer9_bias'

 write (6,*) 'SR: reading layer9_kernel'
  open (unit=555,file='./keras_matrices/layer9_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights9(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer9_kernel'

! ----- layer 10:
 write (6,*) 'SR: reading layer10_bias'
  open (unit=555,file='./keras_matrices/layer10_bias.txt',status='old',action='read')
  read(555,*) bias10(1:outputlength)
  close (555)
 write (6,*) 'SR: finished reading layer10_bias'

 write (6,*) 'SR: reading layer10_kernel'
  open (unit=555,file='./keras_matrices/layer10_kernel.txt',status='old',action='read') 
  do k=1,outputlength
    read(555,*) weights10(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer10_kernel'

#ifdef BRAINDEBUG
    if (masterproc) then
     write (6,*) 'HEY weights1 = ',weights1
     write (6,*) 'HEY bias1 = ',bias1
     write (6,*) 'HEY weights2 = ',weights2
     write (6,*) 'HEY bias2 = ',bias2
     write (6,*) 'HEY weights3 = ',weights3
     write (6,*) 'HEY bias3 = ',bias3
     write (6,*) 'HEY weights4 = ',weights4
     write (6,*) 'HEY bias4 = ',bias4
     write (6,*) 'HEY weights5 = ',weights5
     write (6,*) 'HEY bias5 = ',bias5
     write (6,*) 'HEY weights6 = ',weights6
     write (6,*) 'HEY bias6 = ',bias6
     write (6,*) 'HEY weights7 = ',weights7
     write (6,*) 'HEY bias7 = ',bias7
     write (6,*) 'HEY weights8 = ',weights8
     write (6,*) 'HEY bias8 = ',bias8
     write (6,*) 'HEY weights9 = ',weights9
     write (6,*) 'HEY bias9 = ',bias9
     write (6,*) 'HEY weights10 = ',weights10
     write (6,*) 'HEY bias10 = ',bias10

     write (6,*) 'HEY weights1(5, 10) = ',weights1(5, 10)
     write (6,*) 'HEY weights2(5, 10) = ',weights2(5, 10)
    endif
#endif

end subroutine init_keras_matrices_deep


#else



  subroutine init_keras_matrices_base()

  integer :: nline,k,ipt1,ipt2,j,n,ios
! input layer:
 write (6,*) 'SR: reading layer1_bias'
  open (unit=555,file='./keras_matrices/layer1_bias.txt',status='old',action='read',iostat=ios)
  if (ios .ne. 0) then
    write (6,*) 'HEY keras matrices unable to load, abort.'
    stop
  endif
  read(555,*) bias1(1:width)
  close (555)
   write (6,*) 'SR: reading layer1_kernel'
  open (unit=555,file='./keras_matrices/layer1_kernel.txt',status='old',action='read') 
  do k=1,width
    read(555,*) weights1(k,1:inputlength)
  end do
  close (555)

! ----- layer 2:
 write (6,*) 'SR: reading layer2_bias'
  open (unit=555,file='./keras_matrices/layer2_bias.txt',status='old',action='read')
  read(555,*) bias2(1:outputlength)
  close (555)
 write (6,*) 'SR: finished reading layer2_bias'

 write (6,*) 'SR: reading layer2_kernel'
  open (unit=555,file='./keras_matrices/layer2_kernel.txt',status='old',action='read') 
  do k=1,outputlength
    read(555,*) weights2(k,1:width)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer2_kernel'

end subroutine init_keras_matrices_base
#endif


subroutine init_keras_norm()

  integer :: nline,k,ipt1,ipt2,j,n,ios

! SR: Also read mean and std 
 write (6,*) 'SR: reading means'
  open (unit=555,file='./keras_matrices/inp_means.txt',status='old',action='read')
  read(555,*) input_norm_mean(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading means'
write (6,*) 'SR: reading stds'
  open (unit=555,file='./keras_matrices/inp_stds_by_var.txt',status='old',action='read')
  read(555,*) input_norm_std(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading means'
if (masterproc) then
    write (6,*) 'SR: means = ',input_norm_mean
    write (6,*) 'SR: stds = ',input_norm_std
endif

! SR: Additionally, read min and max arrays
! OUTPUT
 write (6,*) 'SR: reading mins'
  open (unit=555,file='./keras_matrices/outp_mins.txt',status='old',action='read')
  read(555,*) output_norm_min(1:outputlength)
  close (555)
 write (6,*) 'SR: finished reading mins'
write (6,*) 'SR: reading maxs'
  open (unit=555,file='./keras_matrices/outp_maxs.txt',status='old',action='read')
  read(555,*) output_norm_max(1:outputlength)
  close (555)
 write (6,*) 'SR: finished reading maxs'
if (masterproc) then
    write (6,*) 'SR: output mins = ',output_norm_min
    write (6,*) 'SR: output maxs = ',output_norm_max
endif

! INPUT
 write (6,*) 'SR: reading mins'
  open (unit=555,file='./keras_matrices/inp_mins.txt',status='old',action='read')
  read(555,*) input_norm_min(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading mins'
write (6,*) 'SR: reading maxs'
  open (unit=555,file='./keras_matrices/inp_maxs.txt',status='old',action='read')
  read(555,*) input_norm_max(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading maxs'
write (6,*) 'SR: reading max_rs'
  open (unit=555,file='./keras_matrices/inp_max_rs.txt',status='old',action='read')
  read(555,*) input_norm_max_rs(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading max_rs'
if (masterproc) then
    write (6,*) 'SR: input mins = ',input_norm_min
    write (6,*) 'SR: input maxs = ',input_norm_max
    write (6,*) 'SR: input max_rs = ',input_norm_max_rs
endif

  end subroutine init_keras_norm

  

end module cloudbrain_keras_dense

