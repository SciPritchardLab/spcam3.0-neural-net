#include <misc.h>
#include <preproc.h>

module restDGVMMod

#if (defined DGVM)
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: restDGVMMod
! 
! !DESCRIPTION: 
! Read/Write to/from DGVM info to CLM restart file. 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod, only : masterproc
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restart_dgvm
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_dgvm
!
! !INTERFACE:
  subroutine restart_dgvm (nio, flag)
!
! !DESCRIPTION: 
! Read/write DGVM restart data
!
! !USES:
    use clmtype
    use clmpoint
    use iobinary
    use DGVMMod, only : resetWeightsDGVM
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ci,pi                             !indices 
    type(gridcell_type), pointer :: g            !local pointer to derived subtype
    type(landunit_type), pointer :: l            !local pointer to derived subtype
    type(column_type)  , pointer :: c            !local pointer to derived subtype
    type(pft_type)     , pointer :: p            !local pointer to derived subtype
    integer , pointer, dimension(:) :: ibuf1dp   !pointer to memory to be allocated
    real(r8), pointer, dimension(:) :: rbuf1dp   !pointer to memory to be allocated
    real(r8), pointer, dimension(:) :: rbuf1dc   !pointer to memory to be allocated
    integer begc,endc,numc                       !column 1d indices
    integer begp,endp,nump                       !pft 1d indices
    character(len=16) namec                      !column 1d name 
    character(len=16) namep                      !pft 1d name 
!-----------------------------------------------------------------------

    ! Set up shorthand for derived type names
    begc = cols1d%beg
    endc = cols1d%end
    numc = cols1d%num
    namec = cols1d%name

    begp = pfts1d%beg
    endp = pfts1d%end
    nump = pfts1d%num
    namep = pfts1d%name

   ! Allocate necessary 1d buffers

    allocate (ibuf1dp(nump))
    allocate (rbuf1dp(nump))

    !pft type dgvm physical state variable - wf
    allocate (rbuf1dc(numc))
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)  
    do ci = begc,endc
       c => cpoint(ci)%c
       if (flag == 'read' ) c%cps%wf = rbuf1dc(ci)
       if (flag == 'write') rbuf1dc(ci) = c%cps%wf
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)  
    deallocate (rbuf1dc)

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%t_mo_min  = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%t_mo_min  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - annpsn
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%annpsn = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%annpsn 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - annpsnpot
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%annpsnpot = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%annpsnpot 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft cflux tye - fmicr
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pcf%fmicr = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pcf%fmicr  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - bm_inc
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%bm_inc = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%bm_inc  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - afmicr
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%afmicr = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%afmicr
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - t10min
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%t10min = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%t10min  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - tmomin20
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%tmomin20 = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%tmomin20  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - agdd20
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%agdd20 = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%agdd20  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - itypveg
    if (flag == 'read') call readin (nio, ibuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pps%itype = ibuf1dp(pi)
       if (flag == 'write') ibuf1dp(pi) = p%pps%itype
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - fpcgrid
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%fpcgrid = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%fpcgrid  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - lai_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%lai_ind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%lai_ind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - crownarea
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%crownarea = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%crownarea  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - dphen
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%dphen = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%dphen 
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - leafon
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%leafon = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%leafon  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - leafof
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%leafof = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%leafof  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - firelength
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%firelength = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%firelength
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - litterag
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%litterag = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%litterag  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - litterbg
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%litterbg = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%litterbg  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - cpool_fast
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%cpool_fast = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%cpool_fast  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - cpool_slow
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%cpool_slow = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%cpool_slow  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - k_fast_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%k_fast_ave = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%k_fast_ave  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - k_slow_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read') p%pdgvs%k_slow_ave = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%k_slow_ave  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - litter_decom_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%litter_decom_ave = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%litter_decom_ave  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - nind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%nind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%nind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - lm_ind 
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%lm_ind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%lm_ind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - sm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%sm_ind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%sm_ind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - hm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%hm_ind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%hm_ind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - rm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read' ) p%pdgvs%rm_ind = rbuf1dp(pi)
       if (flag == 'write') rbuf1dp(pi) = p%pdgvs%rm_ind  
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! pft type dgvm physical state - present  
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)  
    do pi = begp,endp
       p => ppoint(pi)%p
       if (flag == 'read') then
          p%pdgvs%present = .false.
          if (rbuf1dp(pi) == 1.0) p%pdgvs%present = .true.
       endif
       if (flag == 'write') then
          rbuf1dp(pi) = 0
          if (p%pdgvs%present) rbuf1dp(pi) = 1.0
       endif
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)  

    ! Determine new pft, column and land properties that result from
    ! restart data input

    if (flag == 'read') call resetWeightsDGVM()

    deallocate (ibuf1dp)
    deallocate (rbuf1dp)

  end subroutine restart_dgvm

#endif

end module restDGVMMod

