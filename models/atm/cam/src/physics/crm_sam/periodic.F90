
subroutine periodic(flag)

use vars
use tracers
implicit none

integer flag, ntr

if(flag.eq.0) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)

endif

if(flag.eq.1) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2,2,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2,2,2,3,2)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2,2,2,2,3)

endif

if(flag.eq.2) then

 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,4)
 call bound_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,6)
 if(docloud.and.doprecip) &
     call bound_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,9)
 if(dotracers) then
   do ntr=1,ntracers
     call bound_exchange(tracer(:,:,:,ntr),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,9+ntr)
   end do
 end if

endif
        
if(flag.eq.3) then
        
 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
 call bound_exchange(q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,6)
 if(docloud.and.doprecip) &
     call bound_exchange(qp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,9)
 if(dotracers) then
   do ntr=1,ntracers
     call bound_exchange(tracer(:,:,:,ntr),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,9+ntr)
   end do
 end if
        
endif
        
if(dosgs.and.flag.eq.4) then
        
 call bound_exchange(tk,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7)
 call bound_exchange(tkh,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,8)
        
endif
        
        
end subroutine periodic

