#include <misc.h>
#include <params.h>

module spmd_dyn

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM spectral Eulerian dynamics.
! 
! Author: CCM Core Group
! Modified: P. Worley, September 2002, November 2003
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, masterproc, iam, beglatex, endlatex, numbnd, numlats, numlatsex, &
                           beglat, endlat, begirow, endirow, plev
   use constituents, only: pcnst
   use mpishorthand, only: mpir8, mpicom
   use infnan, only: inf
   implicit none

   private
   public spmdinit_dyn, compute_gsfactors
   save
   integer, public :: npes                 ! Total number of MPI tasks
   integer, public :: nsmps                ! Total number of SMP nodes
   integer, public :: cut(2,0:plat-1)      ! partition for MPI tasks
   integer, public :: cutex(2,0:plat-1)    ! extended partition 
   integer, public :: proc(plat)           ! MPI task id associated with a given lat.
   integer, public :: neighs               ! number of south neighbors to comm guardcells
   integer, public, allocatable :: neighs_proc(:)    ! sorted south process neighbors
   integer, public :: neighn               ! number of north neighbors to comm guardcells
   integer, public, allocatable :: neighn_proc(:)    ! sorted north process neighbors
   integer, public :: npessp               ! number of MPI tasks in spectral space
   integer, public :: bsiz                 ! buffer size
   integer, public :: maxlats              ! max number of lats on any MPI task
   integer, public, allocatable :: nlat_p(:)    ! number of latitudes per MPI task
   integer, public, allocatable :: proc_smp_map(:) ! map of process/SMP node assignments
!
#if (defined MIRROR_DECOMP)
   logical, private :: mirror = .true.     ! flag indicating whether latitudes and their
                                           ! reflections across the equator should assigned 
                                           ! to consecutive processes
#else
   logical, private :: mirror = .false.    ! flag indicating whether latitudes and their
                                           ! reflections across the equator should assigned 
                                           ! to consecutive processes
#endif
!
   real(r8), public, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs

CONTAINS

!========================================================================

   subroutine spmdinit_dyn ()
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processors
! 
! Method: Distribution is S->N for processors 0->npes
! 
! Author: CCM Core Group
! Modified: P. Worley, November 2003 to improve SMP load balance, and to
!           change distribution to 
!             S->E for processors 0,2,..,npes-2
!           and 
!             N->E for processors 1,3,..,npes-1
!           when mirror flag is set (currently a compile-time option)
! 
!-----------------------------------------------------------------------
!     use pspect, only: maxlats
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer procid    ! process id
      integer procids   ! process id SH
      integer procidn   ! process id NH
      integer smpid     ! SMP id
      integer smpids    ! SMP id for SH process
      integer smpidn    ! SMP id for NH process
      integer nlat_base ! minimum number of latitudes per proc
      integer workleft  ! amount of work still to be parcelled out
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer nlat_smp(0:npes-1)  ! number of latitudes per SMP
      integer nproc_smp(0:npes-1) ! number of MPI processes per SMP
      real(r8) avgnlat_proc(0:npes-1) ! average number of latitudes per 
                                      ! MPI process in a given SMP node
      real(r8) minavgnlat_proc        ! minimum average number of latitudes per 
                                      ! MPI process over SMP nodes
      integer neighn_minlat(plat)    ! minimum latitude in north neighbor
      integer neighs_maxlat(plat)    ! maximum latitude in south neighbor
!
!-----------------------------------------------------------------------
!
! Allocate memory for number of lats per proc
!
      allocate (nlat_p (0:npes-1))
      nlat_p(0:npes-1) = 0
!
! Make sure number of PEs and number of latitudes are kosher
!
      call factor (plat, m2, m3, m5)

      if (m2 < 1) then
         write(6,*) 'FACTOR: Problem size is not divisible by 2'
         call endrun
      end if

      if (masterproc) then
         write (6,*) 'Problem factors: 2**',m2,' * 3**',m3,' * 5**',m5
      end if
      call factor (npes, m2, m3, m5)
      
      if (mod(npes,2) /= 0) then
         write(6,*)'SPMDINIT_DYN: nprocs(',npes,') must be a multiple of 2'
         call endrun
      end if
!
! Determine minimum number of latitudes for each process
!
      nlat_base = plat/npes
      do procids=0,npes/2-1
         procidn = npes - procids - 1
         nlat_p(procids) = nlat_base
         nlat_p(procidn) = nlat_p(procids)
      enddo
      maxlats = nlat_base
!
! Calculate initial distribution of latitudes and 
! distribution of processes by SMP
!
      nlat_smp(0:npes-1) = 0
      nproc_smp(0:npes-1) = 0
      do procid=0,npes-1
         smpid = proc_smp_map(procid)
         nproc_smp(smpid) = nproc_smp(smpid) + 1
      enddo
!
      do smpid=0,nsmps-1
         nlat_smp(smpid)     = nlat_base*nlat_smp(smpid)
         avgnlat_proc(smpid) = float(nlat_base)
      enddo
!
! Equi-distribute remaining latitudes across SMPs
! without increasing per process imbalance beyond minimum
!
      workleft = plat - npes*nlat_base
      if (workleft > 0) maxlats = maxlats + 1
      do while (workleft > 0)
!
! (a) Find minimun number of latitudes assigned to an SMP
!
         minavgnlat_proc = avgnlat_proc(0)
         do smpid=1,nsmps-1
            if (minavgnlat_proc > avgnlat_proc(smpid)) then
               minavgnlat_proc = avgnlat_proc(smpid)
            endif
         enddo
!
! (b) Assign an additional latitude to processes with nlat_base
!     latitudes in SMPs with the minimum average number of 
!     latitudes
!
         do procid=npes/2-1,0,-1
            if (mirror) then
               procids = 2*procid
               procidn = procids + 1
            else
               procids = procid
               procidn = npes - procids - 1
            endif
!
            smpids = proc_smp_map(procids)
            smpidn = proc_smp_map(procidn)
            if ((nlat_p(procids) .eq. nlat_base)  .and. &
                ((avgnlat_proc(smpids) .eq. minavgnlat_proc) .or. &
                 (avgnlat_proc(smpidn) .eq. minavgnlat_proc)) .and. &
                (workleft > 0)) then
!
               nlat_p(procids) = nlat_p(procids) + 1
               nlat_smp(smpids) = nlat_smp(smpids) + 1
               avgnlat_proc(smpids) = &
                  float(nlat_smp(smpids))/float(nproc_smp(smpids))
!
               nlat_p(procidn) = nlat_p(procids)
               nlat_smp(smpidn) = nlat_smp(smpidn) + 1
               avgnlat_proc(smpidn) = &
                  float(nlat_smp(smpidn))/float(nproc_smp(smpidn))
!
               workleft = workleft - 2
            endif
         enddo
      end do
!
! Determine latitude assignments
!
      iend = 0
      do procid=0,npes/2-1
         if (mirror) then
            procids = 2*procid
            procidn = procids + 1
         else
            procids = procid
            procidn = npes - procids - 1
         endif
!
         cut(1,procids) = iend + 1
         cut(2,procids) = iend + nlat_p(procids)
         iend = iend + nlat_p(procids)
!
! Assign mirror latitudes
!
         cut(1,procidn) = plat - cut(2,procids) + 1
         cut(2,procidn) = plat - cut(1,procids) + 1
!
! Save local information
!
         if (iam == procids .or. iam == procidn) then
            beglat = cut(1,iam)
            endlat = cut(2,iam)
            numlats = nlat_p(iam)
            begirow = cut(1,procids)
            endirow = cut(2,procids)
         end if
!
      enddo
!
      do procid=0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', &
                      cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                      cut(1,procid),' through ',cut(2,procid)
         end if
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
         cutex(1,procid) = cut(1,procid) - numbnd
         cutex(2,procid) = cut(2,procid) + numbnd
         if (iam == procid) then
            beglatex = cutex(1,procid) + numbnd
            endlatex = cutex(2,procid) + numbnd
            numlatsex = endlatex - beglatex + 1
         end if
      end do
!
! Determine neighbor processors needed for boundary communication.  
! North first.
!
      neighn = 0
      neighn_minlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(1,procid) > cut(2,iam)) .and. &
                (cut(1,procid) <= cut(2,iam)+numbnd)) then
               neighn_minlat(cut(1,procid)) = procid
               neighn = neighn + 1
            endif
         endif
      enddo
!
! Sort north processes by increasing latitude
!
      allocate (neighn_proc (neighn))
      neighn = 0
      do lat=1,plat
         if (neighn_minlat(lat) /= -1) then
            neighn = neighn + 1
            neighn_proc(neighn) = neighn_minlat(lat)
         endif
      enddo
!
! South next.
!
      neighs = 0
      neighs_maxlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(2,procid) < cut(1,iam)) .and. &
                (cut(2,procid) >= cut(1,iam)-numbnd)) then
               neighs_maxlat(cut(2,procid)) = procid
               neighs = neighs + 1
            endif
         endif
      enddo
!
! Sort south processes by decreasing latitude
!
      allocate (neighs_proc (neighs))
      neighs = 0
      do lat=plat,1,-1
         if (neighs_maxlat(lat) /= -1) then
            neighs = neighs + 1
            neighs_proc(neighs) = neighs_maxlat(lat)
         endif
      enddo
!
      if (masterproc) then
         write(6,*)'-----------------------------------------'
         write(6,*)'Number of lats passed north & south = ',numbnd
         write(6,*)'Node  Partition  Extended Partition'
         write(6,*)'-----------------------------------------'
         do procid=0,npes-1
            write(6,200) procid,cut(1,procid),cut(2,procid) ,cutex(1,procid), cutex(2,procid)
200         format(i3,4x,i3,'-',i3,7x,i3,'-',i3)
         end do
      end if
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      call decomp_wavenumbers ()
      call spmdbuf ()

      return
   end subroutine spmdinit_dyn

!========================================================================

   subroutine factor (nitems, m2, m3, m5)
!----------------------------------------------------------------------- 
! 
! Purpose: Factor a given number into powers of 2,3,5
! 
! Method: Brute force application of "mod" function
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nitems      ! Number to be factored into powers of 2,3,5
      integer, intent(out) :: m2,m3,m5   ! Powers of 2, 3, and 5 respectively
!
! Local workspace
!
      integer num                        ! current number to be factored
!
!-----------------------------------------------------------------------
!
      num = nitems
      m2 = 0
      m3 = 0
      m5 = 0
      
2     if (mod(num,2) == 0) then
         m2 = m2 + 1
         num = num/2
         goto 2
      end if
      
3     if (mod(num,3) == 0) then
         m3 = m3 + 1
         num = num/3
         goto 3
      end if
      
5     if (mod(num,5) == 0) then
         m5 = m5 + 1
         num = num/5
         goto 5
      end if
      
      if (num /= 1) then
         write(6,*) 'FACTOR: ',nitems,' has a prime factor other than 2, 3, or 5.  Aborting...'
         call endrun
      end if
      
      return
   end subroutine factor

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Approximately equidistribute both the number of spectral 
!         coefficients and the number of wavenumbers assigned to each 
!         MPI task using a modified version of the mapping due to
!         Barros and Kauranne. 
! 
! Author: P. Worley, September 2002
! 
!-----------------------------------------------------------------------
      use pspect, only: pmmax
      use comspe, only: numm, maxm, locm, nlen, lpspt, lnstart
      use infnan, only: bigint
!
! Local workspace
!
      integer procid      ! processor id
      integer m, lm       ! global and local fourier wavenumber indices
      integer mstride     ! Stride over wavenumbers used in decomposition
      integer begm1       ! Starting Fourier wavenumbers owned by an MPI task
      integer begm2       !  when using Barros & Kauranne decomposition
      integer speccount(0:npes-1)
                          ! number of spectral coefficients assigned to
                          ! each MPI task
!
!-----------------------------------------------------------------------
!
! determine upper bound on number of wavenumbers to be assigned to each 
! process
      if (mod(pmmax,npes) .eq. 0) then
         maxm = pmmax/npes
      else
         maxm = (pmmax/npes) + 1
      endif
      allocate ( locm(1:maxm, 0:npes-1) )
!
! assign wavenumbers to approximately equidistribute the number 
! of spectral coefficients assigned to each processor
      mstride = 2*npes
      npessp = 0
      do procid = 0,npes-1
         numm(procid) = 0
         speccount(procid) = 0
         begm1 = procid + 1
         begm2 = mstride - procid
         do m=begm1,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
         do m=begm2,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
!
         if (numm(procid) .gt. 0) then
            npessp = npessp + 1
         endif
!
      enddo
!
      do procid = 0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', speccount(procid), &
                      ' spectral coefficients and ', numm(procid), &
                      ' m values: ', (locm(lm,procid),lm=1,numm(procid))
         end if
         do lm=numm(procid)+1,maxm
            locm(lm,procid) = bigint
         enddo
      enddo
!
! Calculate number of local spectral coefficients
      lpspt = 0
      do lm=1,numm(iam)
         lpspt = lpspt + nlen(locm(lm,iam))
      enddo
!
! Evaluate displacement info based on truncation params and
! wavenumber assignment
      allocate ( lnstart(1:maxm) )
      lnstart(1) = 0
      do lm=2,numm(iam)
         lnstart(lm) = lnstart(lm-1) + nlen(locm(lm-1,iam))
      enddo
!   
      return
   end subroutine decomp_wavenumbers

!========================================================================

  subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: allocate spmd pack buffers used in pairwise all-all exchanges
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
     use comspe, only: nlen, maxm, locm, numm

     integer maxcount(4),m
     integer length,i,lm
!
! realloc4a max: 8  2 plev*numm*numlats (e.g. tdyn)
!                1  2     *numm*numlats (bpstr)
!                1             *numlats (lat)
!                1                      (numlats)
!
     maxcount(1) = maxlats*(2*maxm*(plev*8 + 1) + 1) + 1
!
! realloc4b max: 8  2 plev*numm*numlats (e.g. vort)
!                4  2     *numm*numlats (e.g. dps)
!
     maxcount(2) = maxlats*(2*maxm*(plev*8 + 4))
!
! realloc5 max: 3                 (len,beglat,numlats)
!               1 numlats         (tmass)
!               5 numlats  *pcnst (e.g. hw1lat)
!               2 4*numlats*pcnst (e.g.hw2al)
!
     maxcount(3) = 3 + maxlats*(6 + (5 + 2*4)*pcnst)
!
! realloc7 max: 2                  (beglat, numlats)
!               3 plev *numlats    (e.g. vmax2d)
!               5      *numlats    (e.g. psurf)
!
     maxcount(4) = maxlats*(3*plev + 5) + 2
     m = maxval(maxcount)
     call mpipack_size (m, mpir8, mpicom, bsiz)
     write(6,*) 'SPMDBUF: Allocating SPMD buffers of size ',bsiz
     allocate(buf1(bsiz/8+1))
     allocate(buf2(bsiz/8+1))
     return
  end subroutine spmdbuf

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = numperlat*(cut(1,p)-1)
     end do
     
  end subroutine compute_gsfactors

#endif

end module spmd_dyn
