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
  real :: SPDT_max_percentile(nlev)
  real :: SPDT_min_percentile(nlev)
  real :: SPDQ_max_percentile(nlev)
  real :: SPDQ_min_percentile(nlev)
  real :: QRL_max_percentile(nlev)
  real :: QRL_min_percentile(nlev)
  real :: QRS_max_percentile(nlev)
  real :: QRS_min_percentile(nlev)
  real :: PRECT_max_percentile, PRECT_min_percentile
  real :: FLUT_max_percentile, FLUT_min_percentile

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
   SPDT(k1:k2) = output(1:nlev) ! W/kg, checked w Rasp
   SPDQ(:) = 0.
   SPDQ(k1:k2) = output((nlev+1):2*nlev)/2.5e6 ! W/kg --> kg/kg/s
!   QRL(:) = 0. ! retain SP or upstream solution above neural net top.
   QRL(k1:k2) = output ((2*nlev+1):3*nlev) ! W/kg 
!   QRS(:) = 0. ! retain SP or upstream solution above neural net top.
   QRS(k1:k2) = output ((3*nlev+1):4*nlev) ! W/kg
   PRECT = output(4*nlev+1)/24./3600. ! m/s
   FLUT = 1.e5*output(4*nlev+2) ! W/m2

! ---- filter for outlier values ---
   do j=1,nlev
     k=j+k1-1
     ! percentiles dfined from raw history file variables, so convert K/s -->
     ! W/kg for SPDT, QRS, QRL to match the units predicted by the net.
     SPDT(k) = min(SPDT(k),1.e3*SPDT_max_percentile(j)) 
     SPDT(k) = max(SPDT(k),1.e3*SPDT_min_percentile(j))
     SPDQ(k) = min(SPDQ(k),SPDQ_max_percentile(j)) ! no unit conversion needed.
     SPDQ(k) = max(SPDQ(k),SPDQ_min_percentile(j))
     QRL(k) = min(QRL(k),1.e3*QRL_max_percentile(j)) ! unit conversion since
!history file was K/s but brain output is W/kg
     QRL(k) = max(QRL(k),1.e3*QRL_min_percentile(j))
     QRS(k) = min(QRS(k),1.e3*QRS_max_percentile(j))
     QRS(k) = max(QRS(k),1.e3*QRS_min_percentile(j))
   end do

! avoid limiting these since jsut diagnostic anyways.
! actually precip is now used in the energy checker so needs to be suitably
! valued.
   PRECT = min(PRECT,PRECT_max_percentile)
   PRECT = max(PRECT,PRECT_min_percentile)
!   FLUT = min(FLUT,FLUT_max_percentile)
!   FLUT = max(FLUT,FLUT_min_percentile)

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

FLUT_max_percentile = 3.094235e+02
FLUT_min_percentile = 7.548374e+01
PRECT_max_percentile = 1.050470e-06
PRECT_min_percentile = 0.000000e+00
QRL_max_percentile = (/ 3.090728e-05, 4.284809e-05 , 5.996321e-05 , 7.136957e-05 , 6.484223e-05 , 5.164710e-05 , 3.828702e-05 , 2.957482e-05 , 2.845253e-05 , 1.779430e-05 , 1.109738e-05 , 8.510120e-06 , 7.348488e-06 , 7.840575e-06 , 1.349668e-05 , 4.123941e-05 , 3.043655e-05 , 3.125438e-05 , 2.569086e-05 , 1.374181e-05 , 7.312753e-06  /)
QRL_min_percentile = (/ -1.898068e-04, -1.713869e-04 , -1.776459e-04 , -1.836122e-04 , -1.720742e-04 , -1.472925e-04 , -1.214079e-04 , -9.899111e-05 , -8.221790e-05 , -7.629859e-05 , -8.143139e-05 , -1.017169e-04 , -1.330241e-04 , -1.890170e-04 , -3.525729e-04 , -3.404294e-04 , -3.561454e-04 , -3.745385e-04 , -3.936331e-04 , -3.276335e-04 , -1.045685e-04  /)
QRS_max_percentile = (/ 9.687352e-05, 1.161653e-04 , 1.174952e-04 , 1.131290e-04 , 1.084309e-04 , 9.884732e-05 , 8.629871e-05 , 7.336954e-05 , 6.304269e-05 , 5.736889e-05 , 5.580080e-05 , 5.902935e-05 , 7.132516e-05 , 1.009557e-04 , 1.490823e-04 , 1.050209e-04 , 8.642577e-05 , 6.579969e-05 , 4.876339e-05 , 3.651136e-05 , 3.365544e-05  /)
QRS_min_percentile = (/ 0.000000e+00, 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00 , 0.000000e+00  /)
SPDQ_max_percentile = (/ 3.652213e-10, 9.945122e-10 , 2.392437e-09 , 5.231546e-09 , 9.734324e-09 , 1.861035e-08 , 3.582969e-08 , 5.714112e-08 , 8.358769e-08 , 1.064119e-07 , 1.193410e-07 , 1.348140e-07 , 1.603937e-07 , 2.265732e-07 , 2.644881e-07 , 2.183017e-07 , 1.645231e-07 , 1.349110e-07 , 1.010198e-07 , 5.974822e-08 , 8.714578e-08  /)
SPDQ_min_percentile = (/ -8.194176e-10, -2.872304e-09 , -8.956697e-09 , -2.421039e-08 , -5.577062e-08 , -1.062700e-07 , -1.648740e-07 , -2.304835e-07 , -2.748349e-07 , -3.067612e-07 , -3.629529e-07 , -4.178967e-07 , -4.451102e-07 , -4.369381e-07 , -4.831560e-07 , -4.720646e-07 , -4.547246e-07 , -4.372495e-07 , -4.363778e-07 , -4.722500e-07 , -3.866966e-07  /)
SPDT_max_percentile = (/ 7.444144e-05, 8.750352e-05 , 1.101834e-04 , 1.557328e-04 , 2.091027e-04 , 3.151416e-04 , 4.353381e-04 , 5.416538e-04 , 6.063663e-04 , 5.978945e-04 , 5.608783e-04 , 5.125014e-04 , 4.689668e-04 , 3.980124e-04 , 3.935995e-04 , 3.410105e-04 , 2.904154e-04 , 2.430378e-04 , 1.880349e-04 , 1.261280e-04 , 3.189476e-05  /)
SPDT_min_percentile = (/ -1.078989e-04, -1.042879e-04 , -1.204434e-04 , -1.204263e-04 , -1.143236e-04 , -1.064950e-04 , -1.071540e-04 , -1.434145e-04 , -1.945311e-04 , -2.329269e-04 , -2.604201e-04 , -2.765352e-04 , -2.872671e-04 , -2.840540e-04 , -2.803827e-04 , -2.439642e-04 , -2.191743e-04 , -2.002433e-04 , -2.200884e-04 , -3.644076e-04 , -5.644033e-04  /)

  end subroutine init_keras_matrices

end module cloudbrain_keras_dense

