#include <misc.h>
#include <params.h>
!#define BRAINDEBUG

module cloudbrain_keras_dense
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use history,         only: outfld, addfld, add_default, phys_decomp
use physconst,       only: gravit,cpair
use pmgrid, only: masterproc

  implicit none
  
  save 

  private

  integer, parameter :: nlev = 21
  integer, parameter :: inputlength = 87
  integer, parameter :: width1 = 1024
  integer, parameter :: width2 = 1024
  integer, parameter :: width3 = 512
  integer, parameter :: width4 = 512
  integer, parameter :: outputlength = 86
  integer, parameter :: nchunk = 64
  real :: bias1(width1)
  real :: weights1(width1,inputlength)
  real :: bias2(width2)
  real :: weights2(width2,width1)
  real :: bias3(width3)
  real :: weights3(width3,width2) 
  real :: bias4(width4)
  real :: weights4(width4,width3)
  real :: bias5(outputlength)
  real :: weights5(outputlength,width4)
  real :: input_norm_mean(inputlength)
  real :: input_norm_std(inputlength)
  real :: SPDT_99th_percentile(pver)
  real :: SPDT_1st_percentile(pver)
  real :: SPDQ_99th_percentile(pver)
  real :: SPDQ_1st_percentile(pver)
  real :: QRL_99th_percentile(pver)
  real :: QRL_1st_percentile(pver)
  real :: QRS_99th_percentile(pver)
  real :: QRS_1st_percentile(pver)
  real :: PRECT_99th_percentile, PRECT_1st_percentile
  real :: FLUT_99th_percentile, FLUT_1st_percentile

  public init_keras_matrices,cloudbrain_dense4_stephan

  contains

  subroutine cloudbrain_dense4_stephan (TAP, QAP, dTdt_adiabatic,dQdt_adiabatic,SHFLX, LHFLX, SOLIN, SPDT,SPDQ,QRL,QRS,PRECT,FLUT)

    real(r8), intent(in) :: TAP(pver)
    real(r8), intent(in) :: QAP(pver)
    real(r8), intent(in) :: dTdt_adiabatic(pver)
    real(r8), intent(in) :: dQdt_adiabatic(pver)
    real(r8), intent(in) :: SHFLX
    real(r8), intent(in) :: LHFLX
    real(r8), intent(in) :: SOLIN
    real(r8), intent(out) :: SPDT(pver) ! W/kg
    real(r8), intent(out) :: SPDQ(pver) ! W/kg
    real(r8), intent(inout) :: QRL(pver) ! W/kg
    real(r8), intent(inout) :: QRS(pver) ! W/kg
    real(r8), intent(out) :: PRECT !? 
    real(r8), intent(out) :: FLUT !?

    real(r8) :: input(inputlength),x1(width1),x2(width2),x3(width3),x4(width4)
    real(r8) :: output (outputlength)
    integer :: k,j,k1,k2

!    call init_keras_matrices

    k1=pver-nlev+1
    k2=pver
    input(1:nlev)=TAP(k1:k2) 
    input((nlev+1):2*nlev)=QAP(k1:k2)
    input((2*nlev+1):3*nlev) = dTdt_adiabatic(k1:k2)
    input((3*nlev+1):4*nlev) = dQdt_adiabatic(k1:k2)
    input(4*nlev+1) = SHFLX
    input(4*nlev+2) = LHFLX
    input(4*nlev+3) = SOLIN
#ifdef BRAINDEBUG
    if (masterproc) then
      write (6,*) 'HEY input=',input
    endif
#endif
    ! normalize input:
    do k=1,inputlength
      input(k) = (input(k) - input_norm_mean(k))/input_norm_std(k)
    end do
#ifdef BRAINDEBUG
   if (masterproc) then
    write (6,*) 'HEY normalized = ',input
    write (6,*) 'HEY bias1=',bias1
    stop
   endif
#endif
! 1st layer: input length-->1024.
    x1(1:width1) = 0.
    do k=1,width1
      do j=1,inputlength
        x1(k) = x1(k) + weights1(k,j)*input(j)
      end do
      x1(k) = x1(k) + bias1(k)
      x1(k) = max(0.,x1(k)) ! relu activation.
    end do
! 2nd layer: 1024->1024
    x2(1:width2) = 0.
    do k=1,width2
      do j=1,width1
        x2(k) = x2(k) + weights2(k,j)*x1(j)
      end do
      x2(k) = x2(k) + bias2(k)
      x2(k) = max(0.,x2(k)) ! relu activation.
    end do
    
! 3rd layer: 1024->512
    x3(1:width3) = 0.
    do k=1,width3
      do j=1,width2
        x3(k) = x3(k) + weights3(k,j)*x2(j)
      end do
      x3(k) = x3(k) + bias3(k)
      x3(k) = max(0.,x3(k)) ! relu activation.
    end do
    
! 4th layer: 512->512
    x4(1:width3) = 0.
    do k=1,width4
      do j=1,width3
        x4(k) = x4(k) + weights4(k,j)*x3(j)
      end do
      x4(k) = x4(k) + bias4(k)
      x4(k) = max(0.,x4(k)) ! relu activation.
    end do

! output layer: 512->output length
   output(1:outputlength) = 0.
   do k=1,outputlength
     do j=1,width4
       output(k) = output(k) + weights5(k,j)*x4(j) 
     end do
     output(k) = output(k) + bias5(k)
     ! no activation for output.
   end do

   SPDT(:) = 0.
   SPDT(k1:k2) = output(1:nlev) ! W/kg, check w Rasp
   SPDQ(:) = 0.
   SPDQ(k1:k2) = output((nlev+1):2*nlev)/2.5e6 ! W/kg --> kg/kg/s
!   QRL(:) = 0. ! retain SP or upstream solution above neural net top.
   QRL(k1:k2) = output ((2*nlev+1):3*nlev) ! W/kg, check w Rasp
!   QRS(:) = 0. ! retain SP or upstream solution above neural net top.
   QRS(k1:k2) = output ((3*nlev+1):4*nlev) ! W/kg, check w Rasp
   PRECT = 1.e2*output(4*nlev+1) ! mm/day? INSERT output normalization
   FLUT = 1.e4*output(4*nlev+2) ! W/m2?

! ---- filter for outlier values ---
   do j=1,nlev
     k=j+k1-1
     SPDT(k) = min(SPDT(k),SPDT_99th_percentile(k))
     SPDT(k) = max(SPDT(k),SPDT_1st_percentile(k))
     SPDQ(k) = min(SPDQ(k),SPDQ_99th_percentile(k))
     SPDQ(k) = max(SPDQ(k),SPDQ_1st_percentile(k))
     QRL(k) = min(QRL(k),QRL_99th_percentile(k))
     QRL(k) = max(QRL(k),QRL_1st_percentile(k))
     QRS(k) = min(QRS(k),QRS_99th_percentile(k))
     QRS(k) = max(QRS(k),QRS_1st_percentile(k))
   end do
     
   PRECT = min(PRECT,PRECT_99th_percentile)
   PRECT = max(PRECT,PRECT_1st_percentile)
   FLUT = min(FLUT,FLUT_99th_percentile)
   FLUT = max(FLUT,FLUT_1st_percentile)

  end subroutine cloudbrain_dense4_stephan

  subroutine init_keras_matrices()

  integer :: nline,k,ipt1,ipt2,j,n,ios
! input layer:
  open (unit=555,file='./keras_matrices/layer1_bias.txt',status='old',action='read',iostat=ios) 
  if (ios .ne. 0) then
    write (6,*) 'HEY keras matrices unable to load, abort.'
    stop
  endif
  nline = width1/nchunk
  do k=1,nline
    ipt1=(k-1)*nchunk+1
    ipt2 = (k-1)*nchunk+nchunk
    read(555,*) bias1(ipt1:ipt2) 
  end do
  close (555)
!  write (6,*) bias1(1000)
  open (unit=555,file='./keras_matrices/layer1_kernel.txt',status='old',action='read') 
  do k=1,width1
    read(555,*) weights1(k,1:inputlength)
  end do
  close (555)
!  write (6,*) weights1(55,44) !-- note this matches kernel (55,44) from the h5
! ----- layer 2:
  open (unit=555,file='./keras_matrices/layer2_bias.txt',status='old',action='read')
  nline=width2/nchunk
  do k=1,nline
    ipt1=(k-1)*nchunk+1
    ipt2 = (k-1)*nchunk+nchunk
    read(555,*) bias2(ipt1:ipt2)
  end do
!  write (6,*) bias2(1000)
  open (unit=555,file='./keras_matrices/layer2_kernel.txt',status='old',action='read')
  nline = width1/nchunk
  do k=1,width2
    do j=1,nline
      ipt1 = (j-1)*nchunk+1
      ipt2 = (j-1)*nchunk+nchunk
      read(555,*) weights2(k,ipt1:ipt2)
    end do
  end do
!  write (6,*) weights2(55,44) ! this matches untransposed kernel (55,44) from the h5
  close (555)

! ----- layer 3:
  open (unit=555,file='./keras_matrices/layer3_bias.txt',status='old',action='read')
  nline=width3/nchunk
  do k=1,nline
    ipt1=(k-1)*nchunk+1
    ipt2 = (k-1)*nchunk+nchunk
    read(555,*) bias3(ipt1:ipt2)
  end do
!  write (6,*) bias3(500)
  open (unit=555,file='./keras_matrices/layer3_kernel.txt',status='old',action='read')
  nline = width2/nchunk
  do k=1,width3
    do j=1,nline
      ipt1 = (j-1)*nchunk+1
      ipt2 = (j-1)*nchunk+nchunk
      read(555,*) weights3(k,ipt1:ipt2)
    end do
  end do
!  write (6,*) weights3(55,44) ! this matches untransposed kernel (55,44) from the h5
  close (555)

! ----- layer 4:
  open (unit=555,file='./keras_matrices/layer4_bias.txt',status='old',action='read')
  nline=width4/nchunk
  do k=1,nline
    ipt1=(k-1)*nchunk+1
    ipt2 = (k-1)*nchunk+nchunk
    read(555,*) bias4(ipt1:ipt2)
  end do
!  write (6,*) bias4(500)
  open (unit=555,file='./keras_matrices/layer4_kernel.txt',status='old',action='read')
  nline = width3/nchunk
  do k=1,width4
    do j=1,nline
      ipt1 = (j-1)*nchunk+1
      ipt2 = (j-1)*nchunk+nchunk
      read(555,*) weights4(k,ipt1:ipt2)
    end do
  end do
!  write (6,*) weights4(55,44) ! this matches untransposed kernel (55,44) from the h5
  close (555)

! ----- layer 5, output:
  open (unit=555,file='./keras_matrices/layer5_bias.txt',status='old',action='read')
  read(555,*) bias5(1:outputlength)
!  write (6,*) bias5(60)
  close (555)
  open (unit=555,file='./keras_matrices/layer5_kernel.txt',status='old',action='read')
  do k=1,width4
    read(555,*) weights5(:,k)
  end do
!  write (6,*) weights5(55,44) ! this matches untransposed kernel (55,44) from the h5
  close (555)

! define input normalization vectors:
! TAP
  n=nlev
  input_norm_mean(1:n) = (/ 2.081741e+02, 2.096764e+02 , 2.112493e+02 , &
2.137626e+02 , 2.175255e+02 , 2.229660e+02 , 2.296972e+02 , 2.373011e+02 , &
2.453548e+02 , 2.534855e+02 , 2.607477e+02 , 2.664352e+02 , 2.703428e+02 , &
2.729633e+02 , 2.750452e+02 , 2.766279e+02 , 2.782079e+02 , 2.797391e+02 , &
2.812025e+02 , 2.824956e+02 , 2.838659e+02 /)
  input_norm_std(1:n) = (/ 1.068176e+01, 8.557691e+00 , 7.126962e+00 , &
7.801861e+00 , 9.709693e+00 , 1.163176e+01 , 1.305740e+01 , 1.377499e+01 , &
1.381810e+01 , 1.342406e+01 , 1.291402e+01 , 1.236711e+01 , 1.150727e+01 , &
1.076127e+01 , 1.062279e+01 , 1.053711e+01 , 1.047981e+01 , 1.047697e+01 , &
1.051171e+01 , 1.057561e+01 , 1.058655e+01 /)
! QAP
  input_norm_mean((n+1):2*n) = (/ 2.036627e-06, 3.745755e-06 , 8.677882e-06 , &
2.006991e-05 , 4.299121e-05 , 8.704409e-05 , 1.665838e-04 , 3.049243e-04 , &
5.466265e-04 , 9.488122e-04 , 1.551427e-03 , 2.408429e-03 , 3.659426e-03 , &
5.004307e-03 , 5.784920e-03 , 6.360915e-03 , 6.886293e-03 , 7.316025e-03 , &
7.632065e-03 , 7.813137e-03 , 8.620515e-03 /)
  input_norm_std((n+1):2*n) = (/ 9.205131e-07, 3.305859e-06 , 1.066414e-05 ,&
2.863349e-05 , 6.660249e-05 , 1.417360e-04 , 2.740884e-04 , 4.825847e-04 ,&
8.077080e-04 , 1.277974e-03 , 1.884112e-03 , 2.591188e-03 , 3.381740e-03 ,&
4.085043e-03 , 4.471513e-03 , 4.755924e-03 , 5.040067e-03 , 5.287095e-03 ,&
5.454785e-03 , 5.548473e-03 , 6.002052e-03 /) 

! dTdt_adiab
  input_norm_mean((2*n+1):3*n) = (/ 8.536206e-07, 1.600277e-06 , 2.957204e-06 ,&
3.625355e-06 , 4.329734e-06 , 4.518082e-06 , 4.731569e-06 , 4.288723e-06 ,&
3.090620e-06 , 1.228350e-06 , 2.375154e-07 , 9.856326e-07 , 3.322685e-06 ,&
3.107050e-06 , -2.164754e-06 , -4.569586e-06 , -6.973220e-06 , -7.999588e-06 ,&
-8.109147e-06 , -6.780746e-06 , -4.220139e-06 /)
  input_norm_std((2*n+1):3*n) = (/ 4.442613e-05, 5.104051e-05 , 5.712672e-05 , &
5.752289e-05 , 5.278342e-05 , 5.389225e-05 , 6.355475e-05 , 7.590572e-05 , &
8.294105e-05 , 8.341112e-05 , 8.007428e-05 , 7.622122e-05 , 7.618488e-05 , &
7.383446e-05 , 7.296872e-05 , 6.615106e-05 , 5.856330e-05 , 5.045315e-05 , &
4.393260e-05 , 4.114061e-05 , 4.291857e-05 /) 

! dQdt_adiab
  input_norm_mean((3*n+1):4*n) = (/ 4.399783e-12, 2.026636e-11 , 7.435193e-11 ,&
2.240202e-10 , 5.608590e-10 , 1.100511e-09 , 1.712260e-09 , 2.380716e-09 ,&
2.866350e-09 , 3.205748e-09 , 3.226877e-09 , 2.182715e-09 , -1.291032e-09 ,&
-4.318887e-09 , -3.510408e-09 , -4.246641e-09 , -4.829205e-09 , -4.755784e-09 ,&
-5.785572e-09 , -5.379017e-09 , -5.684466e-09 /)
  input_norm_std ((3*n+1):4*n) = (/ 4.649366e-11, 1.753216e-10 , 5.676271e-10 ,&
1.608429e-09 , 3.964362e-09 , 7.938207e-09 , 1.268376e-08 , 1.795537e-08 , &
2.151038e-08 , 2.507808e-08 , 3.000117e-08 , 3.576282e-08 , 4.091234e-08 , &
4.101056e-08 , 3.923303e-08 , 3.745079e-08 , 3.503659e-08 , 3.188926e-08 , &
2.852305e-08 , 2.773234e-08 , 2.552654e-08 /) 

  input_norm_mean(4*n+1) = 1.118392e+01 ! SHFLX
  input_norm_std(4*n+1) = 1.530522e+01
  input_norm_mean(4*n+2) = 7.276066e+01 !LHFLX
  input_norm_std(4*n+2) = 8.038508e+01
  input_norm_mean(4*n+3) = 3.260061e+02 !SOLIN
  input_norm_std(4*n+3) = 4.229145e+02

FLUT_1st_percentile = 1.050357e+02
FLUT_99th_percentile = 3.006259e+02
PRECT_1st_percentile = 0.000000e+00
PRECT_99th_percentile = 5.468393e-07

QRL_1st_percentile = (/ -6.216213e-05, -4.786214e-05 , -1.256669e-05 , &
-3.463747e-05 , -2.075084e-05 , -1.487543e-05 , -1.190000e-05 , -1.010048e-05 ,&
-8.858710e-06 , -9.239218e-06 , -3.533971e-05 , -5.625170e-05 , -6.792715e-05 ,&
-7.292182e-05 , -6.816280e-05 , -6.093415e-05 , -5.520021e-05 , -5.255231e-05 ,&
-4.910039e-05 , -4.995287e-05 , -5.905762e-05 , -8.697672e-05 , -1.440127e-04 ,&
-2.254230e-04 , -1.765118e-04 , -1.305296e-04 , -9.227055e-05 , -6.312769e-05 ,&
-4.791830e-05 , -4.394767e-05  /)
QRL_99th_percentile = (/ -2.754649e-06, -3.494669e-06 , 6.443672e-06 ,&
-3.217267e-06 , -1.438849e-06 , 4.984639e-06 , 5.586671e-06 , 8.278785e-06 ,&
6.251957e-06 , 1.054514e-05 , 1.633342e-05 , 2.360345e-05 , 3.079865e-05 ,&
3.016389e-05 , 2.756721e-05 , 1.868325e-05 , 1.010371e-05 , 8.113140e-06 ,&
7.440296e-06 , 4.227221e-06 , 2.190772e-06 , 1.971648e-06 , 4.138415e-06 ,&
1.364484e-05 , 2.147574e-05 , 1.927053e-05 , 1.717588e-05 , 1.542207e-05 ,&
1.175662e-05 , 4.321834e-06  /)
QRS_1st_percentile(:) = 0.
QRS_99th_percentile = (/ 1.498364e-04, 1.111121e-04 , 7.985876e-05 ,&
4.868936e-05 , 2.913275e-05 , 1.973089e-05 , 1.402391e-05 , 1.086071e-05 ,&
8.813042e-06 , 1.012987e-05 , 2.464520e-05 , 3.939078e-05 , 4.780931e-05 ,&
5.244340e-05 , 5.080475e-05 , 4.627056e-05 , 4.470625e-05 , 4.368851e-05 ,&
4.233806e-05 , 4.243507e-05 , 4.617311e-05 , 5.793697e-05 , 7.428529e-05 ,&
8.173476e-05 , 6.425010e-05 , 4.969213e-05 , 3.755981e-05 , 3.325900e-05 ,&
3.117134e-05 , 3.105591e-05  /)  
SPDQ_1st_percentile = (/ 0.000000e+00, 0.000000e+00 , -4.046073e-12 ,&
-1.632948e-12 , -6.370369e-13 , -9.658690e-13 , -2.571320e-12 , -1.171232e-11 ,&
-6.215761e-11 , -2.558772e-10 , -8.807576e-10 , -2.654580e-09 , -7.218700e-09 ,&
-1.776687e-08 , -3.413185e-08 , -5.163957e-08 , -7.193884e-08 , -8.884352e-08 ,&
-1.044166e-07 , -1.189258e-07 , -1.378456e-07 , -1.620004e-07 , -1.893262e-07 ,&
-2.513909e-07 , -2.631886e-07 , -2.673291e-07 , -2.535809e-07 , -2.239956e-07 ,&
-2.286600e-07 , -9.556879e-07  /)
SPDQ_99th_percentile = (/ 0.000000e+00, 0.000000e+00 , 1.442080e-12 ,&
3.227176e-12 , 6.178259e-13 , 7.964458e-13 , 2.072204e-12 , 1.019333e-11 ,&
4.774195e-11 , 1.701812e-10 , 4.811421e-10 , 1.120881e-09 , 2.222386e-09 ,&
3.699655e-09 , 5.659162e-09 , 8.834332e-09 , 1.435378e-08 , 2.283158e-08 ,&
3.362205e-08 , 5.394483e-08 , 8.879985e-08 , 1.588675e-07 , 2.193465e-07 ,&
2.513876e-07 , 2.469180e-07 , 2.417496e-07 , 2.229990e-07 , 1.955417e-07 ,&
1.973605e-07 , 1.994683e-08  /) 
SPDT_1st_percentile = (/ 0.000000e+00, 0.000000e+00 , -1.244579e-04 ,&
1.359525e-07 , -2.295945e-06 , -4.367248e-06 , -1.714213e-05 , -2.709405e-05 ,&
-2.526674e-05 , -2.578611e-05 , -2.522827e-05 , -2.596968e-05 , -2.799722e-05 ,&
-2.654915e-05 , -2.700086e-05 , -2.669265e-05 , -3.055040e-05 , -3.652241e-05 ,&
-4.535123e-05 , -6.417721e-05 , -8.799575e-05 , -1.309851e-04 , -1.510721e-04 ,&
-1.442822e-04 , -1.326312e-04 , -1.191652e-04 , -1.073594e-04 , -1.076496e-04 ,&
-1.308139e-04 , -5.065778e-04  /)
SPDT_99th_percentile = (/ 0.000000e+00, 0.000000e+00 , -9.157017e-08 ,&
8.838763e-05 , 8.283755e-07 , 1.059550e-06 , 2.006552e-06 , 1.477947e-05 ,&
2.515600e-05 , 2.518402e-05 , 2.777238e-05 , 3.307331e-05 , 4.652006e-05 ,&
7.865911e-05 , 1.217936e-04 , 1.731028e-04 , 2.333037e-04 , 2.831942e-04 ,&
3.152643e-04 , 3.288059e-04 , 3.250533e-04 , 3.329234e-04 , 3.275821e-04 ,&
3.292384e-04 , 2.779141e-04 , 2.323042e-04 , 1.875224e-04 , 1.584954e-04 ,&
1.523564e-04 , 7.777996e-05  /)






  end subroutine init_keras_matrices

end module cloudbrain_keras_dense

