#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: trac2d --- Remap Lagrangian to fixed coordinates
!
! !INTERFACE:
 subroutine trac2d( dp1,        q,     nc,    nq, cx,           &
                    cy,       mfx,    mfy,  iord, jord,         &
                    ng,      fill,     im,    jm, km,           &
                    jfirst, jlast, kfirst, klast, va,           &
                    flx )
 
! !USES:

   use shr_kind_mod, only: r8 => shr_kind_r8
   use tp_core
   use fill_module
   use dynamics_vars, only : sine, cosp, acosp, acap, rcap

#if defined( SPMD )
   use pmgrid, only : npr_z, myid_y, npr_y, iam
   use spmd_dyn, only : comm_y, comm_z
   use parutilitiesmodule, only: maxop, parcollective
   use mod_comm, only : mp_send4d_ns, mp_recv4d_ns,  &
                        mp_send3d_2, mp_recv3d_2
#endif

   implicit none

! !INPUT PARAMETERS:

   integer, intent(in):: im, jm, km, jfirst, jlast, kfirst, klast
   integer, intent(in):: ng       ! Max number of ghost latitudes
   integer, intent(in):: nc       ! total # of tracers in q
   integer, intent(in):: nq       ! # of advected tracers
   integer, intent(in):: iord,  jord

   logical, intent(in):: fill

! !INPUT/OUTPUT PARAMETERS:
   real(r8), intent(inout):: dp1(im,jfirst:jlast,kfirst:klast)
   real(r8), intent(inout):: cx(im,jfirst-ng:jlast+ng,kfirst:klast)
   real(r8), intent(inout):: cy(im,jfirst:jlast+1,kfirst:klast)
   real(r8), intent(inout):: mfx(im,jfirst:jlast,kfirst:klast)
   real(r8), intent(inout):: mfy(im,jfirst:jlast+1,kfirst:klast)
   real(r8), intent(inout):: q(im,jfirst-ng:jlast+ng,kfirst:klast,nc)

! Input work arrays
   real(r8), intent(out):: va(im,jfirst:jlast,kfirst:klast)
   real(r8), intent(out):: flx(im,jfirst:jlast,kfirst:klast)


! !DESCRIPTION:
!
!  Perform large-time-step tracer transport using accumulated Courant
!  numbers (cx, cy) and the mass fluxes (mfx, mfy) within the Lagrangian
!  layers.  This routine is 100\% parallel in the vertical direction
!  (with SMP).  Merdional Courant number will be further split, if
!  necessary, to ensure stability.  Cy <= 1 away from poles; Cy $\le$
!  1/2 at the latitudes closest to the poles.
!
! !CALLED FROM:
!    dynpkg 
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.05.26:  Added jfirst:jlast concept; im, jm, km as parameters
!                  replaced IMR, JMR, JNP, NL with IM, JM-1, JM and KM
!   WS  99.09.27:  Documentation; indentation; jfirst:jlast 
!   WS  99.09.30:  Ghosting; loop limits; full parallelization; tested
!   SJL 99.10.15:  nsplt migrated to outermost loop to remove bug
!   SJL 99.12.19:  Local 2D arrays trimmed!
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!   WS  00.07.13:  Changed PILGRIM API
!   AAM 00.08.29:  Added kfirst, klast
!   AAM 01.06.27:  Added y communicators
!   SJL 30.07.01:  MPI optimization/simplification
!   WS  02.04.24:  New mod_comm interfaces
!   WS  02.07.04:  Fixed 2D decomposition bug dest/src for mp_send3d
!
!EOP
!---------------------------------------------------------------------
!BOC

   real(r8)  tiny
   parameter ( tiny = 1.e-10 )

! Local variables:
! 2d arrays
   real(r8)  a2(im,jfirst:jlast)
   real(r8)  fx(im,jfirst:jlast)
   real(r8)  fy(im,jfirst:jlast+1)
   real(r8) cymax(kfirst:klast), cytmp(kfirst:klast)

   real(r8) dp2(im,jfirst:jlast,kfirst:klast)
   logical ffsl(jm,kfirst:klast)

   integer i, j, k
   integer it, iq, kq, nsplt
   integer ktot
   integer js1gd, js2g0, js2gd, jn2g0,jn2gd,jn1g1,jn1gd
#if defined( SPMD )
      integer dest, src
#endif

   real(r8) cy_global
   real(r8) frac
   real(r8) cmax
   real(r8) sum1, sum2

   ktot  = klast - kfirst + 1
   js2g0 = max(2,jfirst)
   jn2g0 = min(jm-1,jlast)
   jn1g1 = min(jm,jlast+1)
   js1gd = max(1,jfirst-ng)     ! NG latitudes on S (starting at 1)
   js2gd = max(2,jfirst-ng)     ! NG latitudes on S (starting at 2)
   jn2gd = min(jm-1,jlast+ng)   ! NG latitudes on S (ending at jm-1)
   jn1gd = min(jm,jlast+ng)     ! NG latitudes on N (ending at jm)

#if defined( SPMD )
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast, kfirst, klast, &
                         ng, ng, cx )
! Send one latitude of both cy and mfy to the south
      dest = iam-1
      src  = iam+1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d_2( dest, src, im, jm, km,                     &
                        1, im, jfirst, jlast+1, kfirst, klast,        &
                        1, im, jfirst, jfirst, kfirst, klast, cy, mfy)
#endif

!$omp parallel do default(shared) private(i,j,k,cmax)
   do k=kfirst,klast
        cymax(k) = 0.
       do j=js2g0,jlast
            cmax = 0.
          do i=1,im
            cmax = max( abs(cy(i,j,k)), cmax)
          enddo
            cymax(k) = max(cymax(k), cmax*(1. + sine(j)**16) )
       enddo
   enddo

#if defined( SPMD )
   call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast, kfirst, klast, &
                      ng, ng, cx )
   call mp_recv3d_2( src, im, jm, km,                            &
                     1, im, jfirst, jlast+1, kfirst, klast,        &
                     1, im, jlast+1, jlast+1, kfirst, klast, cy, mfy)

   call parcollective( comm_y, MAXOP, ktot, cymax )
#endif

! find global max cymax

      cy_global = cymax(kfirst)
   do k=kfirst+1,klast
      cy_global = max(cymax(k), cy_global)
   enddo

#if defined( SPMD )
    if (npr_z > 1) then
        call ParCollective(comm_z, MAXOP, cy_global)
    endif

! Send q for first (and in low resolution cases the only) time split
    call mp_send4d_ns( im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                       ng, ng, q)
#endif

    nsplt = int(1. + cy_global)
    frac  = 1. / float(nsplt)

!$omp  parallel do default(shared) private(i,j,k)

 do 4000 k=kfirst,klast

    if( nsplt .ne. 1 ) then
       do j=js2gd,jn2gd                  
          do i=1,im
            cx(i,j,k) =  cx(i,j,k) * frac      ! cx ghosted on N*ng S*ng
          enddo
       enddo

       do j=js2g0,jn2g0
          do i=1,im
            mfx(i,j,k) = mfx(i,j,k) * frac
          enddo
       enddo

       do j=js2g0,jn1g1                     
          do i=1,im
             cy(i,j,k) =  cy(i,j,k) * frac    ! cy ghosted on N
            mfy(i,j,k) = mfy(i,j,k) * frac    ! mfy ghosted on N
          enddo
       enddo
    endif

       do j=js2g0,jn2g0
          do i=1,im
             if(cy(i,j,k)*cy(i,j+1,k) > 0.) then
                if( cy(i,j,k) > 0.) then
                   va(i,j,k) = cy(i,j,k)
                else
                   va(i,j,k) = cy(i,j+1,k)      ! cy ghosted on N
                endif
             else
                va(i,j,k) = 0.
             endif
          enddo
       enddo

! Check if FFSL extension is needed.

       do 2222 j=js2gd,jn2gd             ! flux needed on N*ng S*ng
          ffsl(j,k) = .false.
          do i=1,im
            if( abs(cx(i,j,k)) > 1. ) then  ! cx ghosted on N*ng S*ng
              ffsl(j,k) = .true.
              go to 2222
            endif
          enddo
2222    continue

! Scale E-W mass fluxes by CX if FFSL
       do j=js2g0,jn2g0
          if( ffsl(j,k) ) then
            do i=1,im
              flx(i,j,k) = mfx(i,j,k) / sign( max(abs(cx(i,j,k)), tiny), &
                                        cx(i,j,k) )
            enddo
          else
            do i=1,im
              flx(i,j,k) = mfx(i,j,k)
            enddo
          endif
       enddo
4000  continue

 do 6000 it=1, nsplt
 
#if defined( SPMD )
    if( it /= 1 ) then
       call mp_send4d_ns( im, jm, km, nq, jfirst, jlast,      &
                          kfirst, klast, ng, ng, q)
    endif
#endif

!$omp parallel do default(shared) private(i,j,k,sum1,sum2)

  do 3000 k=kfirst,klast
     do j=js2g0,jn2g0
        do i=1,im-1
           dp2(i,j,k) =  dp1(i,j,k) + mfx(i,j,k) - mfx(i+1,j,k) +  &
                        (mfy(i,j,k) - mfy(i,j+1,k)) * acosp(j)
        enddo
           dp2(im,j,k) = dp1(im,j,k) + mfx(im,j,k) - mfx(1,j,k) +  &
                         (mfy(im,j,k) - mfy(im,j+1,k)) * acosp(j)
     enddo

      if ( jfirst == 1  ) then
           sum1 = 0.
           do i=1,im
              sum1 = sum1 + mfy(i,2,k)
           end do

           sum1 = - sum1 * rcap
           do i=1,im
              dp2(i,1,k) = dp1(i,1,k) +  sum1
           enddo
      endif

      if ( jlast == jm ) then
           sum2 = 0.
           do i=1,im
              sum2 = sum2 + mfy(i,jm,k)
           end do

              sum2 = sum2 * rcap
           do i=1,im
              dp2(i,jm,k) = dp1(i,jm,k) +  sum2
           enddo
      endif
3000  continue

#if defined( SPMD )
   call mp_recv4d_ns( im, jm, km, nq, jfirst, jlast,        &
                      kfirst, klast, ng, ng, q)
#endif

!$omp parallel do default(shared)    &
!$omp private(i, j, k, kq, iq, fx, fy, a2)

      do 5000 kq=1,ktot*nq
         iq = 1 + (kq-1)/ktot
         k  = kfirst + mod(kq-1,ktot)
         call tp2c(a2, va(1,jfirst,k), q(1,jfirst-ng,k,iq),      &
                   cx(1,jfirst-ng,k) , cy(1,jfirst,k),           &
                   im, jm, iord, jord, ng,                       &
                   fx, fy, ffsl(1,k), rcap, acosp,               &
                   flx(1,jfirst,k), mfy(1,jfirst,k),             &
                   cosp, 1, jfirst, jlast )

         do j=jfirst,jlast
            do i=1,im
               q(i,j,k,iq) = q(i,j,k,iq)*dp1(i,j,k) + a2(i,j)
            enddo
         enddo

         if (fill) call fillxy (q(1,jfirst,k,iq),  im, jm, jfirst, jlast,   &
                                acap, cosp, acosp)

         do j=jfirst,jlast
            do i=1,im
               q(i,j,k,iq) = q(i,j,k,iq) / dp2(i,j,k)
            enddo
         enddo
5000  continue
      if( it /= nsplt ) then
!$omp parallel do private(i, j, k) 
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  dp1(i,j,k) = dp2(i,j,k)
               enddo
            enddo
         enddo
      endif
6000  continue

      return
!EOC
 end subroutine trac2d
!-----------------------------------------------------------------------
