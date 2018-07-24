#include <misc.h>
#include <params.h>

module cloudsimulator
!-----------------------------------------------------------------------
!Purpose: to simulate isccp type of clouds based on instantaneous model
!         profiles, including model level clouds
!
!Author:  W.Lin, with isccp simulator codes from Steve Klein
!-----------------------------------------------------------------------
!
!Special Notes:   
!
!        There are three methods to determine the cloud top. Please
!        refer to the comments associated with variable 'top_height'
!        in subroutine 'isccp_cloud_types' for detail. The cloud top
!        may be adjusted using visible or infrared channels in order
!        to compare with ISCCP daytime or nighttime IR clouds.
!        In current version of code, only 'top_height=1' is used,
!        even during night time. The choices of 'top_height' may be
!        explored by defining 'top_height' in subroutine 'ccm_isccp' 
!        before calling 'isccp_cloud_types'. For example, =1 for
!        local daytime and = 3 for local night time. But if the history
!        tape saves time-averaged fields, additional manipulation will
!        be necessary in order to retrieve the daytime and nighttime mean
!        clouds separately. On the other hand, if history tape saving
!        instantaneous fields, the daytime and nighttime clouds can be
!        directly retrieved based on the availability of local insolation.
!        These methods are not included in standard distribution.
!        If interested, please contact wlin@msrc.sunysb.edu. 
!          
!
!public functions/subroutines:
!
!      ccm_isccp:  to be called from tphysbc
!      isccptab:   to be called from cam.F90 to initialize tau-count tables
!
!--------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use cloudsimulatorparms
   use history,       only: outfld, fillvalue

   implicit none

PRIVATE

   real(r8) tautab(0:255) 	        !  ISCCP table for converting count value to 
                                        !  optical thickness
   integer invtau(-20:45000)            !  ISCCP table for converting optical thickness 
                                        !  to count value

!Public functions/subroutines

   public :: isccptab
   public :: ccm_isccp
#ifdef CRM
   public :: crm_isccp
#endif

CONTAINS

   subroutine ccm_isccp (lchnk, ncol, pmid, pint, q, t, ts, concld,   &
           cld, cliqwp, cicewp, rel, rei, emis, cltot, coszrs)
!
   real(r8)    emsfc_lw   !emsfc_lw    longwave emissivity of surface at 10.5 microns
   integer nsubcol!nsubcol  # of columns each grid box is subdivided into
   PARAMETER(nsubcol=50,emsfc_lw=0.99)
   integer, intent(in) :: lchnk                           ! chunk identifier
   integer, intent(in) :: ncol                            ! number of atmospheric columns
   real(r8), intent(in) :: pmid(pcols,pver)  ! Pressure at middle level
   real(r8), intent(in) :: pint(pcols,pverp)  !
   real(r8), intent(in) :: q   (pcols,pver)  !
   real(r8), intent(in) :: t   (pcols,pver)  !
   real(r8), intent(in) :: ts  (pcols)  !
   real(r8), intent(in) :: concld  (pcols,pver)  !convective cloud cover
   real(r8), intent(in) :: cld  (pcols,pver)  !cloud cover
   real(r8), intent(in) :: cliqwp(pcols,pver) ! in-cloud liquid water path
   real(r8), intent(in) :: cicewp(pcols,pver) ! in-cloud ice water path
   real(r8), intent(in) :: rel(pcols,pver)   ! Liquid cloud particle effective radius
   real(r8), intent(in) :: rei(pcols,pver)   ! Ice effective drop size (microns)
   real(r8), intent(in) :: emis(pcols,pver)  ! Cloud longwave emissivity
   real(r8), intent(in) :: cltot(pcols)  ! Total cloud amount
   real(r8), intent(in) :: coszrs(pcols)  ! cosine solar zenith angle (to tell if day or night)

   real(r8)  fq_isccp_s1(pcols,ntau*npres)  !accumulated fq_isccp

   real(r8) tau(pcols,pver)
   real(r8) dtau_s(pver),dtau_c(pver), dem_s(pver),dem_c(pver)
   real(r8) totalcldarea(pcols), meantaucld(pcols), meanptop(pcols), meanttop(pcols), cloudy(pcols)

   integer itau,top_height,overlap,ibox,seed
   integer i,k
   REAL(r8) fq_isccp(ntau,npres)
   real(r8) boxtau(nsubcol), boxptop(nsubcol)
   real abarl         ! A coefficient for extinction optical depth
   real bbarl         ! b coefficient for extinction optical depth
   real abari         ! A coefficient for extinction optical depth
   real bbari         ! b coefficient for extinction optical depth
   real(r8) cldmin    ! the value taken from radcswmx.F90, note: cldmin much less than cldmin from what on cldnrh
   real(r8) cldeps    ! the value taken from radcswmx.F90
   real(r8) pmid1(pver), pint1(pverp), t1(pver), q1(pver), ts1,  &
            cld1(pver), concld1(pver)
   integer itaunpres,it,ip
   character*1 cc
!
!JR Fake fields added for code checking
!
   real(r8) :: fakefld1(pcols)
   real(r8) :: fakefld2(pcols)

   call t_startf ('ccm_isccp')

   cldmin = 1.0e-80_r8
   cldeps = 0.0_r8

   abarl = 2.817e-02
   bbarl = 1.305
   abari = 3.448e-03
   bbari = 2.431

   !compute the optical depth

    do k=1,pver
       do i=1,ncol
          if(cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
             tau(i,k) = (abarl + bbarl/rel(i,k)) * cliqwp(i,k) + &
                        (abari + bbari/rei(i,k)) * cicewp(i,k)
          else
             tau(i,k) = 0.0
          endif
       end do
    end do

   top_height = 1
   overlap = 3

   do i=1,ncol
!
!JR Add coszrs logic to only call the simulator during daytime
!
      if(coszrs(i) > 0.) then        !if cloudy
         if(cltot(i) >= cldmin) then        !if cloudy
            cloudy(i) = 1.0    !cloudy flag, which may be used to derived mean values under cloudy
        		       !conditions. The cloud flag itself is freuency of cloudy condition
	 		       !when average over an accumulation period.
            seed=(pmid(i,pver)-int(pmid(i,pver)))*100+1
            dtau_s(1:pver) = tau(i,1:pver)
            dtau_c(1:pver) = tau(i,1:pver)
            dem_s (1:pver) = emis(i,1:pver)
            dem_c (1:pver) = emis(i,1:pver)
            pmid1 (1:pver) = pmid(i,1:pver)
            pint1 (1:pverp) = pint(i,1:pverp)
            q1    (1:pver) = q   (i,1:pver)
            t1    (1:pver) = t   (i,1:pver)
            ts1            = ts(i)
            cld1  (1:pver) = cld(i,1:pver)
            concld1(1:pver) = concld(i,1:pver)

            call ISCCP_CLOUD_TYPES(pver,nsubcol,seed,pmid1,pint1,q1,cld1,concld1,  &
                                   dtau_s,dtau_c,top_height,overlap, ts1,emsfc_lw, &
                                   t1,dem_s, dem_c,fq_isccp,totalcldarea(i),       &
                                   meanptop(i), meantaucld(i),boxtau,boxptop,meanttop(i))

!save standard ISCCP type of 7x7 clouds
            do ip=1,npres
               do it=1,ntau
                  itaunpres = (ip-1)*ntau+it
                  fq_isccp_s1(i,itaunpres) = fq_isccp(it,ip)
               end do
            end do
! extra index:
!     23-> totalcldarea
!     24-> meanptop
!     25-> meantaucld
!     26-> meanttop
!     27-> cloudy flag
!     print *,lchnk,itaunpres-pverp+1   != 23
!JR got rid of this and wrote them out as separate fields

            fakefld1(i)      = 1.
            fakefld2(i)      = 1.
         else                             !cloud free in the (daytime) grid box
            fq_isccp_s1(i,:) = 0.
            totalcldarea(i)  = 0.
            meanptop(i)      = fillvalue
            meantaucld(i)    = fillvalue
            meanttop(i)      = fillvalue
            cloudy(i)        = 0.

            fakefld1(i)      = 0.
            fakefld2(i)      = fillvalue
         endif
      else                                ! nighttime
         fq_isccp_s1(i,:) = fillvalue
         totalcldarea(i)  = fillvalue
         meanptop(i)      = fillvalue
         meantaucld(i)    = fillvalue
         meanttop(i)      = fillvalue
         cloudy(i)        = fillvalue

         fakefld1(i)      = fillvalue
         fakefld2(i)      = fillvalue
      end if
   end do
!
!JR dont need to call outfld if all points are nighttime
!
   if (any(coszrs(:ncol) > 0.)) then
#ifdef CRM
      call outfld('_FISCCP1',fq_isccp_s1, pcols,lchnk)
      call outfld('_TCLDAR ',totalcldarea,pcols,lchnk)
      call outfld('_MNPTOP ',meanptop    ,pcols,lchnk)
      call outfld('_MNTAU  ',meantaucld  ,pcols,lchnk)
      call outfld('_MNTTOP ',meanttop    ,pcols,lchnk)
      call outfld('_CLOUDY ',cloudy      ,pcols,lchnk)
#else
      call outfld('FISCCP1 ',fq_isccp_s1, pcols,lchnk)
      call outfld('TCLDAREA',totalcldarea,pcols,lchnk)
      call outfld('MEANPTOP',meanptop    ,pcols,lchnk)
      call outfld('MEANTAU ',meantaucld  ,pcols,lchnk)
      call outfld('MEANTTOP',meanttop    ,pcols,lchnk)
      call outfld('CLOUDY  ',cloudy      ,pcols,lchnk)
#ifdef FILLDEBUG
      call outfld('FAKEFLD1',fakefld1    ,pcols,lchnk)
      call outfld('FAKEFLD2',fakefld2    ,pcols,lchnk)
#endif
#endif
   end if

   call t_stopf ('ccm_isccp')
   return
   end subroutine ccm_isccp

!#######################################################################

subroutine isccp_cloud_types(nlev,ncol,seed,pfull,phalf,qv,       &
           cc,conv,dtau_s,dtau_c,top_height,overlap,              &
           skt,emsfc_lw,at,dem_s,dem_c,fq_isccp,                  &
           totalcldarea,meanptop,meantaucld,boxtau,boxptop,meanttop)

!$Id: cloudsimulator.F90,v 1.1.4.2 2003/06/13 15:54:27 hender Exp $

! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.


!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolmax,nlevmax,ncolprint

      parameter(ncolmax=100,nlevmax=pver)

!     -----
!     Input 
!     -----

      INTEGER nlev                      !  number of model levels in column
      INTEGER ncol                      !  number of subcolumns

      INTEGER seed                      !  seed value for random number generator
                                        !  ( see Numerical Recipes Chapter 7)
                                        !  It is recommended that the seed is set
                                        !  to a different value for each model
                                        !  gridbox it is called on, as it is
					!  possible that the choice of the same
					!  seed value every time may introduce some
					!  statistical bias in the results, particularly
					!  for low values of NCOL.

      REAL(r8) pfull(nlev)	      	!  pressure of full model levels (Pascals)
                                        !  pfull(1)    is    top level of model
                                        !  pfull(nlev) is bottom level of model

      REAL(r8) phalf(nlev+1)            !  pressure of half model levels (Pascals)
                                        !  phalf(1)    is    top       of model
                                        !  phalf(nlev+1) is the surface pressure

      REAL(r8) qv(nlev)                 !  water vapor specific humidity (kg vapor/ kg air)
                                        !         on full model levels

      REAL(r8) cc(nlev)                 !  input cloud cover in each model level (fraction) 
                                        !  NOTE:  This is the HORIZONTAL area of each
                                        !         grid box covered by clouds

      REAL(r8) conv(nlev)               !  input convective cloud cover in each model level (fraction) 
                                        !  NOTE:  This is the HORIZONTAL area of each
                                        !         grid box covered by convective clouds

      REAL(r8) dtau_s(nlev)             !  mean 0.67 micron optical depth of stratiform
					!  clouds in each model level
                                        !  NOTE:  this the cloud optical depth of only the
                                        !         cloudy part of the grid box, it is not weighted
                                        !         with the 0 cloud optical depth of the clear
                                        !         part of the grid box

      REAL(r8) dtau_c(nlev)             !  mean 0.67 micron optical depth of convective
					!  clouds in each
                                        !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
					!  1=max
					!  2=rand
					!  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
					!  optical depth to adjust cloud top pressure. Note
					!  that this calculation is most appropriate to compare
					!  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
					!  3 = adjust top height using only the computed
					!  infrared brightness temperature. Note that this
					!  calculation is most appropriate to compare to ISCCP
					!  IR only algortihm (i.e. you can compare to nighttime
					!  ISCCP data with this option)

!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL(r8) skt                      !  skin Temperature (K)
      REAL(r8) emsfc_lw                 !  10.5 micron emissivity of surface (fraction)                                            
      REAL(r8) at(nlev)                 !  temperature in each model level (K)
      REAL(r8) dem_s(nlev)              !  10.5 micron longwave emissivity of stratiform
					!  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
      REAL(r8) dem_c(nlev)              !  10.5 micron longwave emissivity of convective
					!  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
!     ------
!     Output
!     ------

      REAL(r8) fq_isccp(ntau,npres)     !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL(r8) totalcldarea             !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  This should
					!  equal the sum over all entries of fq_isccp
	
	
      ! The following three means are averages over the cloudy areas only.  If no
      ! clouds are in grid box all three quantities should equal zero.	
					
      REAL(r8) meanptop	                !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
      REAL(r8) meanttop                 !  mean cloud top temp (k) - linear averaging
					
      REAL(r8) meantaucld               !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
      
      
      REAL(r8) boxtau(ncol)             !  optical thickness in each column
      
      REAL(r8) boxptop(ncol)            !  cloud top pressure (mb) in each column
					
															
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL frac_out(ncolmax,nlevmax)    ! boxes gridbox divided up into
					! Equivalent of BOX in original version, but
					! indexed by column then row, rather than
					! by row then column

      REAL tca(ncolmax,0:nlevmax)       ! total cloud cover in each model level (fraction)
                                        ! with extra layer of zeroes on top
                                        ! in this version this just contains the values input
                                        ! from cc but replicated across ncol
      REAL cca(ncolmax,nlevmax)         ! convective cloud cover in each model level (fraction)
                                        ! from conv but replicated across ncol

      REAL threshold(ncolmax)           ! pointer to position in gridbox
      REAL maxocc(ncolmax)              ! Flag for max overlapped conv cld
      REAL maxosc(ncolmax)              ! Flag for max overlapped strat cld

      REAL boxpos(ncolmax)              ! ordered pointer to position in gridbox

      REAL threshold_min(ncolmax)       ! minimum value to define range in with new threshold
                                        ! is chosen

      REAL dem(ncolmax),bb              !  working variables for 10.5 micron longwave 
					!  emissivity in part of
					!  gridbox under consideration

!     REAL    ran0 			! type for random number function
      REAL    dtautmp(ncolmax)          ! temporary variable for dtau of layer
      REAL ptrop,attrop,atmax,atmin,btcmin,transmax
      INTEGER ilev,ibox,itrop,ipres,itau,ilev2
      INTEGER acc(nlevmax,ncolmax)
      INTEGER match(nlevmax-1),nmatch,levmatch(ncolmax)
      
      !variables needed for water vapor continuum absorption
      real fluxtop_clrsky,trans_layers_above_clrsky,taumin
      real dem_wv(nlevmax), wtmair, wtmh20, Navo, grav, pstd, t0
      real press, dpress, atmden, rvh20, wk, rhoave, rh20s, rfrgn
      real tmpexp,tauwv
      
      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL tau(ncolmax),tb(ncolmax),ptop(ncolmax),emcld(ncolmax)
      real ttop(ncolmax)
      real fluxtop(ncolmax),trans_layers_above(ncolmax)
      real              fluxtopinit,tauir
      real meanalbedocld 
      REAL albedocld(ncolmax)
      real boxarea
      
!JR   DATA isccp_taumin / 0.3 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!JR Test for validity of array dimensioning

      if (ncol > ncolmax .or. nlev > nlevmax) then
         write(6,*)'ISCCP_CLOUD_TYPES: nlevmax and/or ncolmax too small'
         write(6,*)'nlev, nlevmax=', nlev, nlevmax
         write(6,*)'ncol, ncolmax=', ncol, ncolmax
         call endrun ()
      end if

      call t_startf ('isccp_cloud_types')
      
      ncolprint=0


!     assign 2d tca array using 1d input array cc

      do ilev=0,nlev
        do ibox=1,ncol
	  if (ilev.eq.0) then
	    tca(ibox,ilev)=0
	  else
	    tca(ibox,ilev)=cc(ilev)
	  endif
        enddo
      enddo

!     assign 2d cca array using 1d input array conv

      do ilev=1,nlev
        do ibox=1,ncol
	  cca(ibox,ilev)=conv(ilev)
        enddo
      enddo

      if (ncolprint.ne.0) then
        write (6,'(a)') 'seed:'
        write (6,'(I3.2)') seed

        write (6,'(a)') 'tca_pp_rev:'
        write (6,'(8f5.2)')  &
         ((tca(ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

        write (6,'(a)') 'cca_pp_rev:'
        write (6,'(8f5.2)')  &
         ((cca(ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)
      endif

      if (top_height .eq. 1 .or. top_height .eq. 3) then 

      ptrop=5000.
      atmin = 400.
      atmax = 0.
      do 12 ilev=1,nlev-1
           if ((pfull(ilev)/phalf(nlev+1)) .lt. 0.4 .and.  &
                at(ilev) .gt. at(ilev+1)) then
                ptrop = pfull(ilev+1)
                attrop = at(ilev+1)
                itrop=ilev+1
           end if
           if (at(ilev) .gt. atmax) atmax=at(ilev)
           if (at(ilev) .lt. atmin) atmin=at(ilev)
12    continue

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!
!     find unpermittable data.....
!
      do 13 ilev=1,nlev
            if (cc(ilev) .lt. 0.) then
                write(6,*)  ' error = cloud fraction less than zero'
                CALL ENDRUN ()
            else if (cc(ilev) .gt. 1.) then
                write(6,*)  ' error = cloud fraction greater than 1'
                CALL ENDRUN ()
            end if 
            if (conv(ilev) .lt. 0.) then
                write(6,*) &
                     ' error = convective cloud fraction less than zero'
                CALL ENDRUN ()
            else if (conv(ilev) .gt. 1.) then
                write(6,*) &
                     ' error = convective cloud fraction greater than 1'
                CALL ENDRUN ()
            end if 

            if (dtau_s(ilev) .lt. 0.) then
                write(6,*)   &
                 ' error = stratiform cloud opt. depth less than zero'
                CALL ENDRUN ()
            end if
            if (dem_s(ilev) .lt. 0.) then
                write(6,*)  &
                 ' error = stratiform cloud emissivity less than zero'
                CALL ENDRUN ()
            else if (dem_s(ilev) .gt. 1.) then
                write(6,*)  &
                 ' error = stratiform cloud emissivity greater than 1'
                CALL ENDRUN ()
            end if 

            if (dtau_c(ilev) .lt. 0.) then
                write(6,*)  &
                 ' error = convective cloud opt. depth less than zero'
                CALL ENDRUN ()
            end if
            if (dem_c(ilev) .lt. 0.) then
                write(6,*)  &
                 ' error = convective cloud emissivity less than zero'
                CALL ENDRUN ()
            else if (dem_c(ilev) .gt. 1.) then
                write(6,*)  &
                 ' error = convective cloud emissivity greater than 1'
                CALL ENDRUN ()
            end if
13    continue

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      do ibox=1,ncol
        do ilev=1,nlev
	  frac_out(ibox,ilev)=0.0
        enddo
      enddo

      if (ncolprint.ne.0) then
        write (6,'(a)') 'frac_out_pp_rev:'
        write (6,'(8f5.2)')  &
        ((frac_out(ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

        write (6,'(a)') 'ncol:'
        write (6,'(I3)') ncol
      endif

      do ibox=1,ncol
	boxpos(ibox)=(ibox-.5)/ncol
      enddo

      if (ncolprint.ne.0) then
        write (6,'(a)') 'last_frac_pp:'
        write (6,'(8f5.2)') (tca(ibox,1),ibox=1,ncolprint)
      endif

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
                                  
!     Initialise threshold

        IF (ilev.eq.1) then
          DO ibox=1,ncol
	    ! If max overlap 
	    IF (overlap.eq.1) then
	      ! select pixels spread evenly
	      ! across the gridbox
              threshold(ibox)=boxpos(ibox)
	    ELSE
	      ! select random pixels from the non-convective
	      ! part the gridbox ( some will be converted into
	      ! convective pixels below )
              threshold(ibox)=   &
                 cca(ibox,ilev)+(1-cca(ibox,ilev))*ran0(seed)
            ENDIF
          ENDDO
          if(ncolprint .ne. 0) then
             write (6,'(a)') 'threshold_nsf2:'
             write (6,'(8f5.2)') (threshold(ibox),ibox=1,ncolprint)
          endif
        ENDIF

        IF (ncolprint.ne.0) then
            write (6,'(a)') 'ilev:'
            write (6,'(I2)') ilev
        ENDIF

        DO ibox=1,ncol

          ! All versions
          if (boxpos(ibox).le.cca(ibox,ilev)) then
              maxocc(ibox) = 1
          else
              maxocc(ibox) = 0
          end if

          ! Max overlap
          if (overlap.eq.1) then 
            threshold_min(ibox)=cca(ibox,ilev)
            maxosc(ibox)=1
          endif

          ! Random overlap
          if (overlap.eq.2) then 
            threshold_min(ibox)=cca(ibox,ilev)
            maxosc(ibox)=0
          endif

          ! Max/Random overlap
          if (overlap.eq.3) then 
            threshold_min(ibox)=max(cca(ibox,ilev),   &
               min(tca(ibox,ilev-1),tca(ibox,ilev)))
            if (threshold(ibox).lt.min(tca(ibox,ilev-1),tca(ibox,ilev))  &
                  .and.(threshold(ibox).gt.cca(ibox,ilev))) then
                   maxosc(ibox)= 1
            else
                   maxosc(ibox)= 0
            end if

          endif
    
          ! Reset threshold 

          threshold(ibox)=                                        &
              !if max overlapped conv cloud
              maxocc(ibox) * (                                    &  
                  boxpos(ibox)                                    &          
              ) +                                                 &    
              !else
              (1-maxocc(ibox)) * (                                &  
                  !if max overlapped strat cloud
                  (maxosc(ibox)) * (                              &   
                      !threshold=boxpos
                      threshold(ibox)                             &          
                  ) +                                             &    
                  !else
                  (1-maxosc(ibox)) * (                            &  
                      !threshold_min=random[thrmin,1]
                      threshold_min(ibox)+                        &
                        (1-threshold_min(ibox))*ran0(seed)        &
                 )                                                &
              )

           ENDDO

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
               if (tca(ibox,ilev).gt.threshold(ibox)) then
               frac_out(ibox,ilev)=1
               else
               frac_out(ibox,ilev)=0
               end if               
           ENDDO

!	   Code to partition boxes into startiform and convective parts
!	   goes here

           DO ibox=1,ncol
                if (threshold(ibox).le.cca(ibox,ilev)) then
                    ! = 2 IF threshold le cca(ibox)
                    frac_out(ibox,ilev) = 2 
                else
                    ! = the same IF NOT threshold le cca(ibox) 
                    frac_out(ibox,ilev) = frac_out(ibox,ilev)
                end if
           ENDDO

!         Set last_frac to tca at this level, so as to be tca 
!         from last level next time round

          if (ncolprint.ne.0) then

            write (6,'(a)') 'last_frac:'
            write (6,'(8f5.2)') (tca(ibox,ilev-1),ibox=1,ncolprint)
    
            write (6,'(a)') 'cca:'
            write (6,'(8f5.2)') (cca(ibox,ilev),ibox=1,ncolprint)
    
            write (6,'(a)') 'max_overlap_cc:'
            write (6,'(8f5.2)') (maxocc(ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'max_overlap_sc:'
            write (6,'(8f5.2)') (maxosc(ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_min_nsf2:'
            write (6,'(8f5.2)') (threshold_min(ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_nsf2:'
            write (6,'(8f5.2)') (threshold(ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'frac_out_pp_rev:'
            write (6,'(8f5.2)')    &
             ((frac_out(ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
          endif

200   CONTINUE    !loop over nlev

!
!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
            tau(ibox)=0.
	    albedocld(ibox)=0.
15    continue

      !compute total cloud optical depth for each column     
      do 26 ilev=1,nlev

          
            !increment tau for each of the boxes
            do 16 ibox=1,ncol

                 if (frac_out(ibox,ilev).eq.1) then
                        dtautmp(ibox)= dtau_s(ilev)
                 else if (frac_out(ibox,ilev).eq.2) then
                        dtautmp(ibox)= dtau_c(ilev)
                 else
                        dtautmp(ibox)= 0.
                 end if

                 tau(ibox)=tau(ibox)+dtautmp(ibox)

16          continue 

            if (ncolprint.ne.0) then
                     write(6,'(i2,1X,8(f7.2,1X))')   &
                     ilev,(dtautmp(ibox),ibox=1,ncolprint)
            endif 

26    continue

          if (ncolprint.ne.0) then

              write(6,'(i2,1X,8(f7.2,1X))')          &
                ilev, (tau(ibox),ibox=1,ncolprint)
          endif 
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky and trans_layers_above_clrsky are the clear
!             sky versions of these quantities.


      if (top_height .eq. 1 .or. top_height .eq. 3) then

        
        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644
        wtmh20 = 18.01534
        Navo = 6.023E+23
        grav = 9.806650E+02
        pstd = 1.013250E+06
        t0 = 296.
        if (ncolprint .ne. 0)  & 
               write(6,*)  'ilev   pw (kg/m2)   tauwv      dem_wv'
        do 125 ilev=1,nlev
               !press and dpress are dyne/cm2 = Pascals *10
               press = pfull(ilev)*10.
               dpress = (phalf(ilev+1)-phalf(ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden = dpress/grav
               rvh20 = qv(ilev)*wtmair/wtmh20
               wk = rvh20*Navo*atmden/wtmair
               rhoave = (press/pstd)*(t0/at(ilev))
               rh20s = rvh20*rhoave
               rfrgn = rhoave-rh20s
               tmpexp = exp(-0.02*(at(ilev)-t0))
               tauwv = wk*1.e-20*( (0.0224697*rh20s*tmpexp) +      &
                      (3.41817e-7*rfrgn)         )*0.98
               dem_wv(ilev) = 1. - exp( -1. * tauwv)
               if (ncolprint .ne. 0)                               &
               write(6,'(i2,1X,3(f8.3,3X))') ilev,                 &
                 qv(ilev)*(phalf(ilev+1)-phalf(ilev))/(grav/100.), &
                 tauwv,dem_wv(ilev)
125     continue


        !initialize variables
        fluxtop_clrsky = 0.
        trans_layers_above_clrsky=1.


        do ilev=1,nlev
 
            ! Black body emission at temperature of the layer

	        bb=1 / ( exp(1307.27/at(ilev)) - 1. )
	        !bb= 5.67e-8*at(ilev)**4

	        ! increase TOA flux by flux emitted from layer
	        ! times total transmittance in layers above

                fluxtop_clrsky = fluxtop_clrsky                      &
                  + dem_wv(ilev) * bb * trans_layers_above_clrsky 
            
                ! update trans_layers_above with transmissivity
	        ! from this layer for next time around loop

                trans_layers_above_clrsky=                           &
                  trans_layers_above_clrsky*(1.-dem_wv(ilev))
                   

            if (ncolprint.ne.0) then
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write (6,'(a)') 'emiss_layer,100.*bb,100.*f,total_trans:'
              write (6,'(4(f7.2,1X))') dem_wv(ilev),100.*bb,         &
                   100.*fluxtop_clrsky,trans_layers_above_clrsky
            endif

        enddo   !loop over level
        
        !add in surface emission
        bb=1/( exp(1307.27/skt) - 1. )
        !bb=5.67e-8*skt**4

        fluxtop_clrsky = fluxtop_clrsky + emsfc_lw * bb              &
           * trans_layers_above_clrsky
            
        if (ncolprint.ne.0) then
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emsfc,100.*bb,100.*f,total_trans:'
          write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb,                 &
            100.*fluxtop_clrsky,trans_layers_above_clrsky
	endif
    

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



        if (ncolprint.ne.0) then

            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt,ibox=1,ncolprint)
    
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)')                                      &
             ((at(ilev2),ibox=1,ncolprint),ilev2=1,nlev)

        endif 
        !loop over columns 
        do ibox=1,ncol
      
            fluxtop(ibox)=0.
            trans_layers_above(ibox)=1.
       
        enddo

        do ilev=1,nlev
            
            do ibox=1,ncol

                ! Black body emission at temperature of the layer

	        bb=1 / ( exp(1307.27/at(ilev)) - 1. )
	        !bb= 5.67e-8*at(ilev)**4

	        ! emissivity for point in this layer
                if (frac_out(ibox,ilev).eq.1) then
                dem(ibox)= 1. -                                     &
                           ( (1. - dem_wv(ilev)) * (1. -  dem_s(ilev)) )
                else if (frac_out(ibox,ilev).eq.2) then
                dem(ibox)= 1. -                                     &
                           ( (1. - dem_wv(ilev)) * (1. -  dem_c(ilev)) )
                else
                dem(ibox)=  dem_wv(ilev)
                end if
                

                ! increase TOA flux by flux emitted from layer
	        ! times total transmittance in layers above

                fluxtop(ibox) = fluxtop(ibox)                       &
                  + dem(ibox) * bb                                  &
                  * trans_layers_above(ibox)                         
            
                ! update trans_layers_above with transmissivity
	        ! from this layer for next time around loop

                trans_layers_above(ibox)=                           &
                  trans_layers_above(ibox)*(1.-dem(ibox))

            enddo ! ibox

            if (ncolprint.ne.0) then

              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(ibox),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*bb:'
              write (6,'(8f7.2)') (100.*bb,ibox=1,ncolprint)
        
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
        
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)')                                   &
                (trans_layers_above(ibox),ibox=1,ncolprint)
          endif

        enddo ! ilev

        do ibox=1,ncol

            !add in surface emission

            bb=1/( exp(1307.27/skt) - 1. )
            !bb=5.67e-8*skt**4

            fluxtop(ibox) = fluxtop(ibox)                           &
               + emsfc_lw * bb                                      &
               * trans_layers_above(ibox) 
            
        end do

        if (ncolprint.ne.0) then

          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*bb:'
          write (6,'(8f7.2)') (100.*bb,ibox=1,ncolprint)
    
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
	endif
    
        do ibox=1,ncol

            !now that you have the top of atmosphere radiance account
            !for ISCCP procedures to determine cloud top temperature

            !account for partially transmitting cloud recompute flux 
            !ISCCP would see assuming a single layer cloud
            !note choice here of 2.13, as it is primarily ice
            !clouds which have partial emissivity and need the 
            !adjustment performed in this section
            !
	    !If it turns out that the cloud brightness temperature
	    !is greater than 260K, then the liquid cloud conversion
            !factor of 2.56 is used.
	    !
            !Note that this is discussed on pages 85-87 of 
            !the ISCCP D level documentation (Rossow et al. 1996)
           

            !compute minimum brightness temperature and optical depth
            btcmin = 1. /  ( exp(1307.27/(attrop-5.)) - 1. ) 
            transmax = (fluxtop(ibox)-btcmin)/(fluxtop_clrsky-btcmin)
            taumin = -1. * log(max(min(transmax,0.9999999),0.001))
            
	    !note that the initial setting of tauir is needed so that
	    !tauir has a realistic value should the next if block be
	    !bypassed
            tauir = tau(ibox) / 2.13

            if (top_height .eq. 1 .and. transmax .gt. 0.001 .and.   &
                transmax .le. 0.9999999) then
                    icycle = 1
                    fluxtopinit = fluxtop(ibox)
		    tauir = tau(ibox) / 2.13
10                  emcld(ibox) = 1. - exp(-1. * tauir  )
                    fluxtop(ibox) = fluxtopinit -                   &
                               ((1.-emcld(ibox))*fluxtop_clrsky)
                    fluxtop(ibox)=max(1.E-06,                       &
                               (fluxtop(ibox)/emcld(ibox)))
                    tb(ibox)= 1307.27/ (log(1. + (1./fluxtop(ibox))))
                    if (icycle .eq. 1 .and. tb(ibox) .gt. 260.) then
		         tauir = tau(ibox) / 2.56
			 icycle = 2
			 go to 10
                    end if			 
            end if

            if (tau(ibox) .gt.  (-1.*log(0.9999999))) then 
                
                !cloudy box
                tb(ibox)= 1307.27/ (log(1. + (1./fluxtop(ibox))))
                
                if (top_height .eq. 1 .and. tauir .lt. taumin) then
                         tb(ibox) = attrop - 5. 
			 tau(ibox) = 2.13*taumin
                end if

            else

                !clear sky brightness temperature
                tb(ibox) = 1307.27/(log(1.+(1./fluxtop_clrsky)))

            end if
            
        enddo ! ibox

        if (ncolprint.ne.0) then

          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)')                                      &
      	  (trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)')                                      &
      	  (1.0-trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)')                                      &
      	  (trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(ibox),ibox=1,ncolprint)
	endif
    
      end if
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!



      !compute cloud top pressure
      do 30 ibox=1,ncol
      
               !segregate according to optical thickness
               if (tau(ibox) .le. (-1.*log(0.9999999))) then

                         ptop(ibox)=0.
			 ttop(ibox)=0.
                         levmatch(ibox)=0      

               else 

                     if (top_height .eq. 1 .or. top_height .eq. 3) then  
                                               
                        !find level whose temperature
                        !most closely matches brightness temperature
                        nmatch=0
                        do 29 ilev=1,nlev-1
                        
                            if ((at(ilev)   .ge. tb(ibox) .and.      &
                                 at(ilev+1) .lt. tb(ibox)) .or.      &
                                (at(ilev) .le. tb(ibox) .and.        &
                                 at(ilev+1) .gt. tb(ibox))) then 
     
                                  nmatch=nmatch+1
                                  if(abs(at(ilev)-tb(ibox)) .lt.     &
                                     abs(at(ilev+1)-tb(ibox))) then
                                         match(nmatch)=ilev
                                  else
                                         match(nmatch)=ilev+1
                                  end if
                            end if                        
29                      continue

                        if (nmatch .ge. 1) then
                                 
                            ptop(ibox)=pfull(match(nmatch))
			    ttop(ibox)=at(match(nmatch))
                            levmatch(ibox)=match(nmatch)   
                        else
                                                        
                            if (tb(ibox) .lt. atmin) then
                                 ptop(ibox)=ptrop
 				 ttop(ibox)=atmin
                                 levmatch(ibox)=itrop
                            end if
                            if (tb(ibox) .gt. atmax) then
                                 ptop(ibox)=pfull(nlev)
				 ttop(ibox)=atmax
                                 levmatch(ibox)=nlev
                            end if
                                                                
                        end if
                                                               
                     else
                          ptop(ibox)=0.
			  ttop(ibox)=0.
                          ilev=1
                          do while(ptop(ibox) .eq. 0.                 &
                                    .and. ilev .lt. nlev+1)
                                if (frac_out(ibox,ilev) .ne. 0) then
                                   ptop(ibox)=pfull(ilev)
				   levmatch(ibox)=ilev
                                end if
                                ilev=ilev+1
                          end do
                     end if
               end if

30    continue


!
!     ---------------------------------------------------!


!
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined,
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.
!

      !compute isccp frequencies

      !reset frequencies
      do 38 itau=1,ntau
      do 38 ipres=1,npres
             fq_isccp(itau,ipres)=0.
38    continue

      !reset variables need for averaging cloud properties
      totalcldarea = 0.
      meanalbedocld = 0.
      meanptop = 0.
      meanttop = 0.
      meantaucld = 0.
      boxarea = 1./real(ncol)

      do 39 ibox=1,ncol

            !convert ptop to millibars
            ptop(ibox)=ptop(ibox) / 100.

	    !save for output cloud top pressure and optical thickness
	    boxtau(ibox) = tau(ibox)
	    boxptop(ibox) = ptop(ibox)

      if (tau(ibox) .gt. (-1.*log(0.9999999))                 &
           .and. ptop(ibox) .gt. 0.) then

            !convert optical thickness to albedo
	    albedocld(ibox)=real(invtau(min(nint(100.*tau(ibox)),45000)))
	    
            !contribute to averaging
            totalcldarea = totalcldarea + boxarea
	    meanalbedocld = meanalbedocld + albedocld(ibox)*boxarea
	    meanptop = meanptop + ptop(ibox)*boxarea
	    meanttop = meanttop + ttop(ibox)*boxarea

            !reset itau, ipres
            itau = 0
            ipres = 0

            !determine optical depth category
            if (tau(ibox) .lt. taulim(2)) then
                itau=1
            else if (tau(ibox) .ge. taulim(2)                  &
                                        .and. tau(ibox) .lt. taulim(3)) then
                itau=2
            else if (tau(ibox) .ge. taulim(3) .and. tau(ibox) .lt. taulim(4)) then
                itau=3
            else if (tau(ibox) .ge. taulim(4) .and. tau(ibox) .lt. taulim(5)) then
                itau=4
            else if (tau(ibox) .ge. taulim(5) .and. tau(ibox) .lt. taulim(6)) then
                itau=5
            else if (tau(ibox) .ge. taulim(6) .and. tau(ibox) .lt. taulim(7)) then
                itau=6
            else if (tau(ibox) .ge. taulim(7)) then
                itau=7
            end if

            !determine cloud top pressure category
            if (    ptop(ibox) .gt. prlim(1)  .and.ptop(ibox) .lt. prlim(2)) then
                ipres=1
            else if(ptop(ibox) .ge. prlim(2).and.ptop(ibox) .lt. prlim(3)) then
                ipres=2
            else if(ptop(ibox) .ge. prlim(3).and.ptop(ibox) .lt. prlim(4)) then
                ipres=3
            else if(ptop(ibox) .ge. prlim(4).and.ptop(ibox) .lt. prlim(5)) then
                ipres=4
            else if(ptop(ibox) .ge. prlim(5).and.ptop(ibox) .lt. prlim(6)) then
                ipres=5
            else if(ptop(ibox) .ge. prlim(6).and.ptop(ibox) .lt. prlim(7)) then
                ipres=6
            else if(ptop(ibox) .ge. prlim(7)) then
                ipres=7
            end if

            !update frequencies
            if(ipres .gt. 0.and.itau .gt. 0) then
            fq_isccp(itau,ipres)=fq_isccp(itau,ipres)+ boxarea
            end if

      end if

39    continue


      !compute mean cloud properties
      if (totalcldarea .gt. 0.) then
         meanalbedocld = meanalbedocld / totalcldarea
	 meanptop = meanptop / totalcldarea
         meanttop = meanttop / totalcldarea
	 meantaucld = tautab(min(255,max(1,nint(meanalbedocld))))
      else
!JR New code added here to prevent zeros getting accumulated into 2-d output arrays
         meanptop = fillvalue
         meanttop = fillvalue
         meantaucld = fillvalue
      end if
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
!     to see info replace '3 .eq. 4' with '3 .eq. 3'
!
  
      if (3 .eq. 4) then
            
            !produce character output
            do 4031 ilev=1,nlev
            do 4031 ibox=1,ncol
                   acc(ilev,ibox)=frac_out(ibox,ilev)*2
                   if (levmatch(ibox) .eq. ilev)                      &
                       acc(ilev,ibox)=acc(ilev,ibox)+1
4031         continue

             !print test

	     write (6,*) 'Gridbox decomposition written to unit 9'

             write(9,'(a1)') ' '
             write(9,'(10i5)') (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '

             do ibox=1,ncol
                    write(9,'(40(a1),1x,40(a1))')                     &
                        (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlevmax)  &
                       ,(cchar(acc(ilev,ibox)+1),ilev=1,nlevmax) 
             end do

             if (ncolprint.ne.0) then
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)')                   &
                        'ilev',                                       &
                        'pfull','at',                                 &
                        'cc*100','dem_s','dtau_s',                    &
                        'cchar'

               do 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
                   write(6,'(i2,1X,5(f7.2,1X),50(a1))')               &
                        ilev,                                         &
                        pfull(ilev)/100.,at(ilev),                    &
                        cc(ilev)*100.0,dem_s(ilev),dtau_s(ilev),      &
                        (cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
4012           continue
               write (6,'(a)') 'skt:'
               write (6,'(8f7.2)') skt
                                      
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
	      
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(ibox),ibox=1,ncolprint)
             endif 
      end if 
      call t_stopf ('isccp_cloud_types')
      return
end subroutine isccp_cloud_types

#ifdef CRM
   subroutine crm_isccp (pmid, pint, q, t, ts,    &
           cld, cliqwp, cicewp, rel, rei, emis, coszrs,    &
	   fq_isccp_s1, totalcldarea, lowcldarea, midcldarea, hghcldarea, &
           meantaucld, meanptop, meanttop, cloudy)
!
	include 'crmdims.inc'
   
   real(r8)    emsfc_lw   !emsfc_lw    longwave emissivity of surface at 10.5 microns
   PARAMETER(emsfc_lw=0.99)
   real(r8), intent(in) :: pmid(pver)  ! Pressure at middle level
   real(r8), intent(in) :: pint(pverp)  !
   real(r8), intent(in) :: q   (crm_nx, crm_ny, crm_nz)  !
   real(r8), intent(in) :: t   (crm_nx, crm_ny, crm_nz)  !
   real(r8), intent(in) :: ts    !
   real(r8), intent(in) :: cld  (crm_nx, crm_ny, crm_nz)  !cloud cover
   real(r8), intent(in) :: cliqwp(crm_nx, crm_ny, crm_nz) ! in-cloud liquid water path
   real(r8), intent(in) :: cicewp(crm_nx, crm_ny, crm_nz) ! in-cloud ice water path
   real(r8), intent(in) :: rel(crm_nx, crm_ny, crm_nz)   ! Liquid cloud particle effective radius
   real(r8), intent(in) :: rei(crm_nx, crm_ny, crm_nz)   ! Ice effective drop size (microns)
   real(r8), intent(in) :: emis(crm_nx, crm_ny, crm_nz)  ! Cloud longwave emissivity
   real(r8), intent(in) :: coszrs  ! cosine solar zenith angle (to tell if day or night)

   real(r8), intent(out) :: fq_isccp_s1(ntau*npres)  !accumulated fq_isccp

   real(r8), intent(out) :: totalcldarea !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  This should
					!  equal the sum over all entries of fq_isccp
      ! The following three means are averages over the cloudy areas only.  If no
      ! clouds are in grid box all three quantities should equal zero.
   real(r8), intent(out) :: lowcldarea ! low clouds (below 700 mb)
   real(r8), intent(out) :: midcldarea ! middle clouds (above 700 mb and below 400 mb)
   real(r8), intent(out) :: hghcldarea ! high clouds (above 400 mb)

   real(r8), intent(out) :: lowcldarea ! low clouds (below 700 mb)
   real(r8), intent(out) :: meantaucld !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
   real(r8), intent(out) :: meanptop !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
   real(r8), intent(out) :: meanttop !  mean cloud top temp (k) - linear averaging
   real(r8), intent(out) :: cloudy
   
! local
   integer ncolmax
   parameter (ncolmax=crm_nx*crm_ny)
   real(r8) tau(ncolmax, pver)
   real(r8) q_down(ncolmax, pver)  !
   real(r8) t_down(ncolmax, pver)  !
   real(r8) cld_down  (ncolmax, pver)  !cloud cover
   real(r8) emis_down(ncolmax, pver)  ! Cloud longwave emissivity
   real(r8) cltot  ! Total cloud amount in grid cell
   real(r8) cldcol ! Total cloud amount in crm column

   integer itau,top_height,ibox,seed
   integer i,j,k,m,ii
   REAL(r8) fq_isccp(ntau,npres)
   real abarl         ! A coefficient for extinction optical depth
   real bbarl         ! b coefficient for extinction optical depth
   real abari         ! A coefficient for extinction optical depth
   real bbari         ! b coefficient for extinction optical depth
   real(r8) cldmin    ! the value taken from radcswmx.F90, note: cldmin much less than cldmin from what on cldnrh
   real(r8) cldeps    ! the value taken from radcswmx.F90
   integer itaunpres,it,ip
   character*1 cc
!
!JR Fake fields added for code checking
!
   real(r8) :: fakefld1
   real(r8) :: fakefld2

   call t_startf ('crm_isccp')

   cldmin = 1.0e-80_r8
   cldeps = 0.0_r8

   abarl = 2.817e-02
   bbarl = 1.305
   abari = 3.448e-03
   bbari = 2.431

   !compute the optical depth
    do m=crm_nz+1,pver
       k = pver-m+1
       cld_down(:ncolmax,k)=0.
       emis_down(:ncolmax,k)=0.
       tau(:ncolmax,k)=0.
       ii=0
       do i=1,crm_nx
       do j=1,crm_ny
          ii=ii+1
          t_down(ii,k)=t(i,j,crm_nz)
          q_down(ii,k)=q(i,j,crm_nz)
       end do
       end do
    enddo

    do m=1,crm_nz
       k = pver-m+1
       ii=0
       do i=1,crm_nx
       do j=1,crm_ny
       ii=ii+1
          cld_down(ii,k)=cld(i,j,m)
	  q_down(ii,k)=q(i,j,m)
	  t_down(ii,k)=t(i,j,m)
	  emis_down(ii,k)=emis(i,j,m)
          if(cld_down(ii,k) >= cldmin .and. cld_down(ii,k) >= cldeps) then
             tau(ii,k) = (abarl + bbarl/rel(i,j,m)) * cliqwp(i,j,m) + &
                         (abari + bbari/rei(i,j,m)) * cicewp(i,j,m)
          else
             tau(ii,k) = 0.0
          endif
       end do
       end do
    end do
    cltot=0.
    do i=1,crm_nx
    do j=1,crm_ny
       cldcol=0
       do m=1,crm_nz
          cldcol=max(cldcol,cld(i,j,m))
       end do
       cltot=cltot+cldcol
    end do
    end do
    cltot=cltot/(crm_nx*crm_ny)

    top_height = 1

!
!JR Add coszrs logic to only call the simulator during daytime
!
      if(coszrs > 0.) then        !if cloudy
         if(cltot >= cldmin) then        !if cloudy
            cloudy = 1.0    !cloudy flag, which may be used to derived mean values under cloudy
        		       !conditions. The cloud flag itself is freuency of cloudy condition
	 		       !when average over an accumulation period.

            call ISCCP_CLOUD_TYPES_crm(pver, crm_nx*crm_ny, pmid, pint, q_down, cld_down,  &
                                   tau, top_height, ts, emsfc_lw, &
                                   t_down, emis_down, fq_isccp, totalcldarea,       &
                                   lowcldarea, midcldarea, hghcldarea, &
                                   meanptop, meantaucld, meanttop)

!           save standard ISCCP type of 7x7 clouds
            do ip=1,npres
               do it=1,ntau
                  itaunpres = (ip-1)*ntau+it
                  fq_isccp_s1(itaunpres) = fq_isccp(it,ip)
               end do
            end do
         else                             !cloud free in the (daytime) grid box
            fq_isccp_s1(:ntau*npres) = 0.
            totalcldarea  = 0.
            lowcldarea  = 0.
            midcldarea  = 0.
            hghcldarea  = 0.
            meanptop      = fillvalue
            meantaucld    = fillvalue
            meanttop      = fillvalue
            cloudy        = 0.
         endif
      else                                ! nighttime
         fq_isccp_s1(1:ntau*npres) = fillvalue
         totalcldarea  = fillvalue
         lowcldarea  = fillvalue
         midcldarea  = fillvalue
         hghcldarea  = fillvalue
         meanptop      = fillvalue
         meantaucld    = fillvalue
         meanttop      = fillvalue
         cloudy        = fillvalue
      end if
      


   call t_stopf ('crm_isccp')
   return
   end subroutine crm_isccp

!#######################################################################

subroutine isccp_cloud_types_crm(nlev,ncol,pfull,phalf,qv,       &
           cc,dtau_s,top_height,                              &
           skt,emsfc_lw,at,dem_s,fq_isccp,                  &
           totalcldarea,lowcldarea,midcldarea,hghcldarea, &
           meanptop,meantaucld,meanttop)

!$Id: cloudsimulator.F90,v 1.1.4.2 2003/06/13 15:54:27 hender Exp $

! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.

	include 'crmdims.inc'

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolmax,nlevmax,ncolprint

      parameter(ncolmax=crm_nx*crm_ny,nlevmax=pver)

!     -----
!     Input 
!     -----

      INTEGER nlev                      !  number of model levels in column
      INTEGER ncol                      !  number of subcolumns


      REAL(r8) pfull(nlev)	      	!  pressure of full model levels (Pascals)
                                        !  pfull(1)    is    top level of model
                                        !  pfull(nlev) is bottom level of model

      REAL(r8) phalf(nlev+1)            !  pressure of half model levels (Pascals)
                                        !  phalf(1)    is    top       of model
                                        !  phalf(nlev+1) is the surface pressure

      REAL(r8) qv(ncol,nlev)                 !  water vapor specific humidity (kg vapor/ kg air)
                                        !         on full model levels

      REAL(r8) cc(ncol,nlev)                 !  input cloud cover in each model level (fraction)
                                        !  NOTE:  This is the HORIZONTAL area of each
                                        !         grid box covered by clouds

      REAL(r8) dtau_s(ncol,nlev)             !  mean 0.67 micron optical depth of stratiform
					!  clouds in each model level
                                        !  NOTE:  this the cloud optical depth of only the
                                        !         cloudy part of the grid box, it is not weighted
                                        !         with the 0 cloud optical depth of the clear
                                        !         part of the grid box


      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
					!  optical depth to adjust cloud top pressure. Note
					!  that this calculation is most appropriate to compare
					!  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
					!  3 = adjust top height using only the computed
					!  infrared brightness temperature. Note that this
					!  calculation is most appropriate to compare to ISCCP
					!  IR only algortihm (i.e. you can compare to nighttime
					!  ISCCP data with this option)

!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL(r8) skt                      !  skin Temperature (K)
      REAL(r8) emsfc_lw                 !  10.5 micron emissivity of surface (fraction)
      REAL(r8) at(ncol,nlev)            !  temperature in each model level (K)
      REAL(r8) dem_s(ncol,nlev)         !  10.5 micron longwave emissivity of
					!  clouds in each model level.
!     ------
!     Output
!     ------

      REAL(r8) fq_isccp(ntau,npres)     !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL(r8) totalcldarea             !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  This should
					!  equal the sum over all entries of fq_isccp
	
      REAL(r8) lowcldarea, midcldarea, hghcldarea
	
      ! The following three means are averages over the cloudy areas only.  If no
      ! clouds are in grid box all three quantities should equal zero.	
					
      REAL(r8) meanptop	                !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
      REAL(r8) meanttop                 !  mean cloud top temp (k) - linear averaging
					
      REAL(r8) meantaucld               !  mean optical thickness (dimensionless)
                                        !  linear averaging in albedo performed.
      
      
      REAL(r8) boxtau(ncol)             !  optical thickness in each column
      
      REAL(r8) boxptop(ncol)            !  cloud top pressure (mb) in each column

!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------


      REAL threshold(ncolmax)           ! pointer to position in gridbox
      REAL maxocc(ncolmax)              ! Flag for max overlapped conv cld
      REAL maxosc(ncolmax)              ! Flag for max overlapped strat cld

      REAL boxpos(ncolmax)              ! ordered pointer to position in gridbox

      REAL dem(ncolmax),bb              !  working variables for 10.5 micron longwave 
					!  emissivity in part of
					!  gridbox under consideration

!     REAL    ran0 			! type for random number function
      REAL    dtautmp(ncolmax)          ! temporary variable for dtau of layer
      REAL ptrop(ncolmax),attrop(ncolmax),atmax(ncolmax),atmin(ncolmax),btcmin,transmax
      INTEGER ilev,ibox,itrop(ncolmax),ipres,itau,ilev2
      INTEGER match(nlevmax-1),nmatch,levmatch(ncolmax)
      
      !variables needed for water vapor continuum absorption
      real fluxtop_clrsky(ncolmax),trans_layers_above_clrsky(ncolmax),taumin
      real dem_wv(ncolmax,nlevmax), wtmair, wtmh20, Navo, grav, pstd, t0
      real press, dpress, atmden, rvh20, wk, rhoave, rh20s, rfrgn
      real tmpexp,tauwv

      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL tau(ncolmax),tb(ncolmax),ptop(ncolmax),emcld(ncolmax)
      real ttop(ncolmax)
      real fluxtop(ncolmax),trans_layers_above(ncolmax)
      real              fluxtopinit,tauir
      real meanalbedocld 
      REAL albedocld(ncolmax)
      real boxarea
      
!JR   DATA isccp_taumin / 0.3 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!JR Test for validity of array dimensioning

      if (ncol > ncolmax .or. nlev > nlevmax) then
         write(6,*)'ISCCP_CLOUD_TYPES: nlevmax and/or ncolmax too small'
         write(6,*)'nlev, nlevmax=', nlev, nlevmax
         write(6,*)'ncol, ncolmax=', ncol, ncolmax
         call endrun ()
      end if

      call t_startf ('isccp_cloud_types')
      
      ncolprint=0



      if (top_height .eq. 1 .or. top_height .eq. 3) then

      do ibox=1,ncol
      ptrop(ibox)=5000.
      atmin(ibox) = 400.
      atmax(ibox) = 0.
      do 12 ilev=1,nlev-1
           if ((pfull(ilev)/phalf(nlev+1)) .lt. 0.4 .and.  &
                at(ibox,ilev) .gt. at(ibox,ilev+1)) then
                ptrop(ibox) = pfull(ilev+1)
                attrop(ibox) = at(ibox,ilev+1)
                itrop(ibox)=ilev+1
           end if
           if (at(ibox,ilev) .gt. atmax(ibox)) atmax(ibox)=at(ibox,ilev)
           if (at(ibox,ilev) .lt. atmin(ibox)) atmin(ibox)=at(ibox,ilev)
12    continue
      enddo

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!
!     find unpermittable data.....
!
      do 13 ibox=1,ncol
      do 13 ilev=1,nlev
            if (cc(ibox,ilev) .lt. 0.) then
                write(6,*)  ' error = cloud fraction less than zero'
                CALL ENDRUN ()
            else if (cc(ibox,ilev) .gt. 1.) then
                write(6,*)  ' error = cloud fraction greater than 1'
                CALL ENDRUN ()
            end if 
            if (dtau_s(ibox,ilev) .lt. 0.) then
                write(6,*)   &
                 ' error = cloud opt. depth =',dtau_s(ibox,ilev)
		 print *,'ibox,ilev=',ibox,ilev
                CALL ENDRUN ()
            end if
            if (dem_s(ibox,ilev) .lt. 0.) then
                write(6,*)  &
                 ' error = cloud emissivity less than zero'
                CALL ENDRUN ()
            else if (dem_s(ibox,ilev) .gt. 1.) then
                write(6,*)  &
                 ' error = cloud emissivity greater than 1'
                CALL ENDRUN ()
            end if 

13    continue


!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!
      if (ncolprint.ne.0) then
        write (6,'(a)') 'last_frac_pp:'
        write (6,'(8f5.2)') (cc(ibox,1),ibox=1,ncolprint)
      endif

!     ---------------------------------------------------!
!!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
            tau(ibox)=0.
	    albedocld(ibox)=0.
15    continue

      !compute total cloud optical depth for each column     
      do 26 ilev=1,nlev
            !increment tau for each of the boxes
            do 16 ibox=1,ncol
                 tau(ibox)=tau(ibox)+dtau_s(ibox,ilev)
16          continue

26    continue

          if (ncolprint.ne.0) then
              write(6,'(i2,1X,8(f7.2,1X))')          &
                ilev, (tau(ibox),ibox=1,ncolprint)
          endif 
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky and trans_layers_above_clrsky are the clear
!             sky versions of these quantities.


      if (top_height .eq. 1 .or. top_height .eq. 3) then

        
        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644
        wtmh20 = 18.01534
        Navo = 6.023E+23
        grav = 9.806650E+02
        pstd = 1.013250E+06
        t0 = 296.
        if (ncolprint .ne. 0)  & 
               write(6,*)  'ilev   pw (kg/m2)   tauwv      dem_wv'
        do 125 ilev=1,nlev
               !press and dpress are dyne/cm2 = Pascals *10
               press = pfull(ilev)*10.
               dpress = (phalf(ilev+1)-phalf(ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden = dpress/grav
	       do ibox=1,ncol
                 rvh20 = qv(ibox,ilev)*wtmair/wtmh20
                 wk = rvh20*Navo*atmden/wtmair
                 rhoave = (press/pstd)*(t0/at(ibox,ilev))
                 rh20s = rvh20*rhoave
                 rfrgn = rhoave-rh20s
                 tmpexp = exp(-0.02*(at(ibox,ilev)-t0))
                 tauwv = wk*1.e-20*( (0.0224697*rh20s*tmpexp) +      &
                      (3.41817e-7*rfrgn)         )*0.98
                 dem_wv(ibox,ilev) = 1. - exp( -1. * tauwv)
                 if (ncolprint .ne. 0)                               &
                    write(6,'(i2,1X,3(f8.3,3X))') ilev,                 &
                    qv(ibox,ilev)*(phalf(ilev+1)-phalf(ilev))/(grav/100.), &
                    tauwv,dem_wv(ibox,ilev)
               enddo
125     continue


        !initialize variables
        fluxtop_clrsky(:ncol) = 0.
        trans_layers_above_clrsky(:ncol)=1.


        do ibox=1,ncol
        do ilev=1,nlev
 
            ! Black body emission at temperature of the layer
	        bb=1 / ( exp(1307.27/at(ibox,ilev)) - 1. )
	        !bb= 5.67e-8*at(ilev)**4

	        ! increase TOA flux by flux emitted from layer
	        ! times total transmittance in layers above

                fluxtop_clrsky(ibox) = fluxtop_clrsky(ibox)                      &
                  + dem_wv(ibox,ilev) * bb * trans_layers_above_clrsky(ibox)
            
                ! update trans_layers_above with transmissivity
	        ! from this layer for next time around loop

                trans_layers_above_clrsky(ibox)=                           &
                  trans_layers_above_clrsky(ibox)*(1.-dem_wv(ibox,ilev))
                   

               if (ncolprint.ne.0) then
                 write (6,'(a)') 'ilev:'
                 write (6,'(I2)') ilev
    
                 write (6,'(a)') 'emiss_layer,100.*bb,100.*f,total_trans:'
                 write (6,'(4(f7.2,1X))') dem_wv(ibox,ilev),100.*bb,         &
                   100.*fluxtop_clrsky(ibox),trans_layers_above_clrsky(ibox)
               endif

        enddo   !loop over level
        
        !add in surface emission
        bb=1/( exp(1307.27/skt) - 1. )
        !bb=5.67e-8*skt**4

        fluxtop_clrsky(ibox) = fluxtop_clrsky(ibox) + emsfc_lw * bb              &
           * trans_layers_above_clrsky(ibox)
            
        if (ncolprint.ne.0) then
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emsfc,100.*bb,100.*f,total_trans:'
          write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb,                 &
            100.*fluxtop_clrsky(ibox),trans_layers_above_clrsky(ibox)
	endif
     enddo


        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



        if (ncolprint.ne.0) then

            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt,ibox=1,ncolprint)
    
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)')                                      &
             ((at(ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)

        endif 
        !loop over columns 
        do ibox=1,ncol
      
            fluxtop(ibox)=0.
            trans_layers_above(ibox)=1.
       
        enddo

        do ilev=1,nlev
            
            do ibox=1,ncol

                ! Black body emission at temperature of the layer

	        bb=1 / ( exp(1307.27/at(ibox,ilev)) - 1. )
	        !bb= 5.67e-8*at(ibox,ilev)**4

	        ! emissivity for point in this layer
                dem(ibox)= 1. - ( (1. - dem_wv(ibox,ilev)) * (1. -  dem_s(ibox,ilev)) )

                ! increase TOA flux by flux emitted from layer
	        ! times total transmittance in layers above

                fluxtop(ibox) = fluxtop(ibox)                       &
                  + dem(ibox) * bb                                  &
                  * trans_layers_above(ibox)                         
            
                ! update trans_layers_above with transmissivity
	        ! from this layer for next time around loop

                trans_layers_above(ibox)=                           &
                  trans_layers_above(ibox)*(1.-dem(ibox))

            enddo ! ibox

            if (ncolprint.ne.0) then

              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(ibox),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*bb:'
              write (6,'(8f7.2)') (100.*bb,ibox=1,ncolprint)
        
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
        
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)')                                   &
                (trans_layers_above(ibox),ibox=1,ncolprint)
          endif

        enddo ! ilev

        !add in surface emission

        bb=1/( exp(1307.27/skt) - 1. )
        !bb=5.67e-8*skt**4
        do ibox=1,ncol
            fluxtop(ibox) = fluxtop(ibox)                           &
               + emsfc_lw * bb                                      &
               * trans_layers_above(ibox)
            
        end do

        if (ncolprint.ne.0) then

          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*bb:'
          write (6,'(8f7.2)') (100.*bb,ibox=1,ncolprint)
    
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
	endif
    
        do ibox=1,ncol

            !now that you have the top of atmosphere radiance account
            !for ISCCP procedures to determine cloud top temperature

            !account for partially transmitting cloud recompute flux 
            !ISCCP would see assuming a single layer cloud
            !note choice here of 2.13, as it is primarily ice
            !clouds which have partial emissivity and need the 
            !adjustment performed in this section
            !
	    !If it turns out that the cloud brightness temperature
	    !is greater than 260K, then the liquid cloud conversion
            !factor of 2.56 is used.
	    !
            !Note that this is discussed on pages 85-87 of
            !the ISCCP D level documentation (Rossow et al. 1996)
           

            !compute minimum brightness temperature and optical depth
            btcmin = 1. /  ( exp(1307.27/(attrop(ibox)-5.)) - 1. )
            transmax = (fluxtop(ibox)-btcmin)/(fluxtop_clrsky(ibox)-btcmin)
            taumin = -1. * log(max(min(transmax,0.9999999),0.001))

	    !note that the initial setting of tauir is needed so that
	    !tauir has a realistic value should the next if block be
	    !bypassed
            tauir = tau(ibox) / 2.13

            if (top_height .eq. 1 .and. transmax .gt. 0.001 .and.   &
                transmax .le. 0.9999999) then
                    icycle = 1
                    fluxtopinit = fluxtop(ibox)
		    tauir = tau(ibox) / 2.13
10                  emcld(ibox) = 1. - exp(-1. * tauir  )
                    fluxtop(ibox) = fluxtopinit -                   &
                               ((1.-emcld(ibox))*fluxtop_clrsky(ibox))
                    fluxtop(ibox)=max(1.E-06,                       &
                               (fluxtop(ibox)/emcld(ibox)))
                    tb(ibox)= 1307.27/ (log(1. + (1./fluxtop(ibox))))
                    if (icycle .eq. 1 .and. tb(ibox) .gt. 260.) then
		         tauir = tau(ibox) / 2.56
			 icycle = 2
			 go to 10
                    end if			 
            end if

            if (tau(ibox) .gt.  (-1.*log(0.9999999))) then 
                
                !cloudy box
                tb(ibox)= 1307.27/ (log(1. + (1./fluxtop(ibox))))
                
                if (top_height .eq. 1 .and. tauir .lt. taumin) then
                         tb(ibox) = attrop(ibox) - 5.
			 tau(ibox) = 2.13*taumin
                end if

            else

                !clear sky brightness temperature
                tb(ibox) = 1307.27/(log(1.+(1./fluxtop_clrsky(ibox))))

            end if
            
        enddo ! ibox

        if (ncolprint.ne.0) then

          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100.*fluxtop(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)')                                      &
      	  (trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)')                                      &
      	  (1.0-trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)')                                      &
      	  (trans_layers_above(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(ibox),ibox=1,ncolprint)
	endif
    
      end if
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!



      !compute cloud top pressure
      do 30 ibox=1,ncol
      
               !segregate according to optical thickness
               if (tau(ibox) .le. (-1.*log(0.9999999))) then

                         ptop(ibox)=0.
			 ttop(ibox)=0.
                         levmatch(ibox)=0      

               else 

                     if (top_height .eq. 1 .or. top_height .eq. 3) then  
                                               
                        !find level whose temperature
                        !most closely matches brightness temperature
                        nmatch=0
                        do 29 ilev=1,nlev-1
                        
                            if ((at(ibox,ilev)   .ge. tb(ibox) .and.      &
                                 at(ibox,ilev+1) .lt. tb(ibox)) .or.      &
                                (at(ibox,ilev) .le. tb(ibox) .and.        &
                                 at(ibox,ilev+1) .gt. tb(ibox))) then
     
                                  nmatch=nmatch+1
                                  if(abs(at(ibox,ilev)-tb(ibox)) .lt.     &
                                     abs(at(ibox,ilev+1)-tb(ibox))) then
                                         match(nmatch)=ilev
                                  else
                                         match(nmatch)=ilev+1
                                  end if
                            end if                        
29                      continue

                        if (nmatch .ge. 1) then
                                 
                            ptop(ibox)=pfull(match(nmatch))
			    ttop(ibox)=at(ibox,match(nmatch))
                            levmatch(ibox)=match(nmatch)   
                        else
                                                        
                            if (tb(ibox) .lt. atmin(ibox)) then
                                 ptop(ibox)=ptrop(ibox)
 				 ttop(ibox)=atmin(ibox)
                                 levmatch(ibox)=itrop(ibox)
                            end if
                            if (tb(ibox) .gt. atmax(ibox)) then
                                 ptop(ibox)=pfull(nlev)
				 ttop(ibox)=atmax(ibox)
                                 levmatch(ibox)=nlev
                            end if
                                                                
                        end if
                                                               
                     else
                          ptop(ibox)=0.
			  ttop(ibox)=0.
                          ilev=1
                          do while(ptop(ibox) .eq. 0.                 &
                                    .and. ilev .lt. nlev+1)
                                if (cc(ibox,ilev) .gt. 0.01) then
                                   ptop(ibox)=pfull(ilev)
				   levmatch(ibox)=ilev
                                end if
                                ilev=ilev+1
                          end do
                     end if
               end if

30    continue


!
!     ---------------------------------------------------!


!
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined,
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.
!

      !compute isccp frequencies

      !reset frequencies
      do 38 itau=1,ntau
      do 38 ipres=1,npres
             fq_isccp(itau,ipres)=0.
38    continue

      !reset variables need for averaging cloud properties
      totalcldarea = 0.
      lowcldarea = 0.
      midcldarea = 0.
      hghcldarea = 0.
      meanalbedocld = 0.
      meanptop = 0.
      meanttop = 0.
      meantaucld = 0.
      boxarea = 1./real(ncol)

      do 39 ibox=1,ncol

            !convert ptop to millibars
            ptop(ibox)=ptop(ibox) / 100.

	    !save for output cloud top pressure and optical thickness
	    boxtau(ibox) = tau(ibox)
	    boxptop(ibox) = ptop(ibox)

      if (tau(ibox) .gt. (-1.*log(0.9999999))                 &
           .and. ptop(ibox) .gt. 0.) then

            !convert optical thickness to albedo
	    albedocld(ibox)=real(invtau(min(nint(100.*tau(ibox)),45000)))
	    
            !contribute to averaging
            totalcldarea = totalcldarea + boxarea
            if (ptop(ibox) .ge. 700. ) then
               lowcldarea = lowcldarea + boxarea
            elseif (ptop(ibox) .le. 400.) then
               hghcldarea = hghcldarea + boxarea
            else
               midcldarea = midcldarea + boxarea
            end if
 

	    meanalbedocld = meanalbedocld + albedocld(ibox)*boxarea
	    meanptop = meanptop + ptop(ibox)*boxarea
	    meanttop = meanttop + ttop(ibox)*boxarea

            !reset itau, ipres
            itau = 0
            ipres = 0

            !determine optical depth category
            if (tau(ibox) .lt. taulim(2)) then
                itau=1
            else if (tau(ibox) .ge. taulim(2)                  &
                                        .and. tau(ibox) .lt. taulim(3)) then
                itau=2
            else if (tau(ibox) .ge. taulim(3) .and. tau(ibox) .lt. taulim(4)) then
                itau=3
            else if (tau(ibox) .ge. taulim(4) .and. tau(ibox) .lt. taulim(5)) then
                itau=4
            else if (tau(ibox) .ge. taulim(5) .and. tau(ibox) .lt. taulim(6)) then
                itau=5
            else if (tau(ibox) .ge. taulim(6) .and. tau(ibox) .lt. taulim(7)) then
                itau=6
            else if (tau(ibox) .ge. taulim(7)) then
                itau=7
            end if

            !determine cloud top pressure category
            if (    ptop(ibox) .gt. prlim(1)  .and.ptop(ibox) .lt. prlim(2)) then
                ipres=1
            else if(ptop(ibox) .ge. prlim(2).and.ptop(ibox) .lt. prlim(3)) then
                ipres=2
            else if(ptop(ibox) .ge. prlim(3).and.ptop(ibox) .lt. prlim(4)) then
                ipres=3
            else if(ptop(ibox) .ge. prlim(4).and.ptop(ibox) .lt. prlim(5)) then
                ipres=4
            else if(ptop(ibox) .ge. prlim(5).and.ptop(ibox) .lt. prlim(6)) then
                ipres=5
            else if(ptop(ibox) .ge. prlim(6).and.ptop(ibox) .lt. prlim(7)) then
                ipres=6
            else if(ptop(ibox) .ge. prlim(7)) then
                ipres=7
            end if

            !update frequencies
            if(ipres .gt. 0.and.itau .gt. 0) then
            fq_isccp(itau,ipres)=fq_isccp(itau,ipres)+ boxarea
            end if

      end if

39    continue


      !compute mean cloud properties
      if (totalcldarea .gt. 0.) then
         meanalbedocld = meanalbedocld / totalcldarea
	 meanptop = meanptop / totalcldarea
         meanttop = meanttop / totalcldarea
	 meantaucld = tautab(min(255,max(1,nint(meanalbedocld))))
      else
!JR New code added here to prevent zeros getting accumulated into 2-d output arrays
         meanptop = fillvalue
         meanttop = fillvalue
         meantaucld = fillvalue
      end if
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
!     to see info replace '3 .eq. 4' with '3 .eq. 3'
!
  
      if (3 .eq. 4) then
            

             !print test

	     write (6,*) 'Gridbox decomposition written to unit 9'

             write(9,'(a1)') ' '
             write(9,'(10i5)') (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '


             if (ncolprint.ne.0) then
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)')                   &
                        'ilev',                                       &
                        'pfull','at',                                 &
                        'cc*100','dem_s','dtau_s',                    &
                        'cchar'

               write (6,'(a)') 'skt:'
               write (6,'(8f7.2)') skt
                                      
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
	      
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(ibox),ibox=1,ncolprint)
             endif
      end if 
      call t_stopf ('isccp_cloud_types_crm')
      return
end subroutine isccp_cloud_types_crm
    
#endif

subroutine isccptab
   use pmgrid, only: masterproc
   use mpishorthand
   use filenames, only: isccpdata
   use ioFileMod, only: getfil

   include 'netcdf.inc'

   integer ncid,tautabid,invtauid
   character(len=256) locfn ! local filename
!
   if (masterproc) then
      call getfil (isccpdata, locfn)
      call wrap_open (locfn,NF_NOWRITE,ncid)
      call wrap_inq_varid (ncid,'tautab',tautabid)
      call wrap_inq_varid (ncid,'invtau',invtauid)
      call wrap_get_var_realx (ncid,tautabid,tautab)
      call wrap_get_var_int   (ncid,invtauid,invtau)
      call wrap_close (ncid)
   end if
#if ( defined SPMD )
   call mpibcast (tautab, 256,   mpir8, 0, mpicom)
   call mpibcast (invtau, 45021, mpiint, 0, mpicom)
#endif
   return
end subroutine isccptab

real(r8) function ran0(idum)

!     $Id: cloudsimulator.F90,v 1.1.4.2 2003/06/13 15:54:27 hender Exp $
!     Platform independent random number generator from
!     Numerical Recipies
!     Mark Webb July 1999
      
      integer idum,IA,IM,IQ,IR,k
      real(r8) AM

      parameter (IA=16807, IM=2147483647, AM=1.0/IM, IQ=127773, IR=2836)
      
      if (idum.eq.0) then
	write(6,*) 'idum=',idum
	write(6,*) 'ZERO seed not allowed'
	call endrun ()
      endif

      k=idum/IQ
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran0=am*idum
      return
end function ran0
end module cloudsimulator
