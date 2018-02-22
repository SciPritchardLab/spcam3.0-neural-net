#include <misc.h>
#include <params.h>
#define BRAINDEBUG

module cloudbrain_keras_dense
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair
use pmgrid, only: masterproc

  implicit none
  
  save 

  private
  ! Define the network architecture. For now just one hidden layer.
  integer, parameter :: nlev = 30
  integer, parameter :: inputlength = 152
  integer, parameter :: width1 = 512
  integer, parameter :: outputlength = 120
  integer, parameter :: nchunk = 64
  real :: bias1(width1)
  real :: weights1(width1,inputlength)
  real :: bias2(outputlength)
  real :: weights2(outputlength,width1)
  real :: input_norm_mean(inputlength)
  real :: input_norm_std(inputlength)
  real :: SPDT_max_percentile(nlev)
  real :: SPDT_min_percentile(nlev)
  real :: SPDQ_max_percentile(nlev)
  real :: SPDQ_min_percentile(nlev)
  real :: QRL_max_percentile(nlev)
  real :: QRL_min_percentile(nlev)
  real :: QRS_max_percentile(nlev)
  real :: QRS_min_percentile(nlev)

  public init_keras_matrices, cloudbrain_purecrm_base

  contains

  subroutine cloudbrain_purecrm_base (TC, QC, VC, dTdt_adiabatic, dQdt_adiabatic, PS, SOLIN, SPDT, SPDQ, QRL, QRS, icol)
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

    real(r8) :: input(inputlength),x1(width1)
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
    ! normalize input:
    do k=1,inputlength
      input(k) = (input(k) - input_norm_mean(k))/input_norm_std(k)
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
    x1(1:width1) = 0.
    do k=1,width1
      do j=1,inputlength
        x1(k) = x1(k) + weights1(k,j)*input(j)
      end do
      x1(k) = x1(k) + bias1(k)
      x1(k) = max(0.,x1(k)) ! relu activation.
    end do
! output layer: 512->output length
   output(1:outputlength) = 0.
   do k=1,outputlength
     do j=1,width1
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

! SR, I will comment this out for now. Except for QRS which must be strictly positive.
! ---- filter for outlier values ---
   do j=1,nlev
     k=j+k1-1
!      ! percentiles dfined from raw history file variables, so convert K/s -->
!      ! W/kg for SPDT, QRS, QRL to match the units predicted by the net.
!      SPDT(k) = min(SPDT(k),1.e3*SPDT_max_percentile(j)) 
!      SPDT(k) = max(SPDT(k),1.e3*SPDT_min_percentile(j))
!      SPDQ(k) = min(SPDQ(k),SPDQ_max_percentile(j)) ! no unit conversion needed.
!      SPDQ(k) = max(SPDQ(k),SPDQ_min_percentile(j))
!      QRL(k) = min(QRL(k),1.e3*QRL_max_percentile(j)) ! unit conversion since
! !history file was K/s but brain output is W/kg
!      QRL(k) = max(QRL(k),1.e3*QRL_min_percentile(j))
    !  QRS(k) = min(QRS(k),1.e3*QRS_max_percentile(j))
    !  QRS(k) = max(QRS(k),1.e3*QRS_min_percentile(j))
     QRS(k) = max(QRS(k),0.)
   end do

  end subroutine cloudbrain_purecrm_base

  subroutine init_keras_matrices()

  integer :: nline,k,ipt1,ipt2,j,n,ios
! input layer:
 write (6,*) 'SR: reading layer1_bias'
  open (unit=555,file='./keras_matrices/layer1_bias.txt',status='old',action='read',iostat=ios)
  if (ios .ne. 0) then
    write (6,*) 'HEY keras matrices unable to load, abort.'
    stop
  endif
  read(555,*) bias1(1:width1)
  close (555)
   write (6,*) 'SR: reading layer1_kernel'
  open (unit=555,file='./keras_matrices/layer1_kernel.txt',status='old',action='read') 
  do k=1,width1
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
    read(555,*) weights2(k,1:width1)
  end do
  close (555)
 write (6,*) 'SR: finished reading layer2_kernel'


! SR: Also read mean and std 
 write (6,*) 'SR: reading means'
  open (unit=555,file='./keras_matrices/means.txt',status='old',action='read')
  read(555,*) input_norm_mean(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading means'
write (6,*) 'SR: reading stds'
  open (unit=555,file='./keras_matrices/stds.txt',status='old',action='read')
  read(555,*) input_norm_std(1:inputlength)
  close (555)
 write (6,*) 'SR: finished reading means'
if (masterproc) then
    write (6,*) 'SR: means = ',input_norm_mean
    write (6,*) 'SR: stds = ',input_norm_std
endif

  end subroutine init_keras_matrices

end module cloudbrain_keras_dense

