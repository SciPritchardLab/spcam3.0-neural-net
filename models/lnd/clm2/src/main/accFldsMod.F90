#include <misc.h>
#include <preproc.h>

module accFldsMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: accFldsMod
!
! !DESCRIPTION: 
! This module contains subroutines that initialize, update and extract
! the user-specified fields over user-defined intervals. Each interval 
! and accumulation type is unique to each field processed. 
! Subroutine [initAccumFlds] defines the fields to be processed
! and the type of accumulation. Subroutine [updateAccumFlds] does 
! the actual accumulation for a given field. Fields are accumulated 
! by calls to subroutine [update_accum_field]. To accumulate a field, 
! it must first be defined in subroutine [initAccumFlds] and then 
! accumulated by calls to [updateAccumFlds]. 
! Four types of accumulations are possible:
!   o average over time interval
!   o running mean over time interval
!   o running accumulation over time interval
! Time average fields are only valid at the end of the averaging interval.
! Running means are valid once the length of the simulation exceeds the
! averaging interval. Accumulated fields are continuously accumulated.
! The trigger value "-99999." resets the accumulation to zero.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initAccFlds      ! Initialization accumulator fields
  public :: initAccClmtype   ! Initialize clmtype variables obtained from accum fields
  public :: updateAccFlds    ! Update accumulator fields
!                              
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
!EOP

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccFlds()
!
! !INTERFACE: 
  subroutine initAccFlds()
!
! !DESCRIPTION: 
! Initializes accumulator and sets up array of accumulated fields
!
! !USES:
    use clmtype
    use accumulMod, only : init_accum_field, print_accum_fields
    use time_manager, only : get_step_size, get_nstep
    use shr_const_mod, only: SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use nanMod, only : bigint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY::
! Created by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES
!
    integer :: nf                        !accumulated field index 
    integer :: i,k                       !do loop indices
    integer :: dtime                     !time step size  
    integer :: pi                        !indices
    integer :: begp, endp                !per-process beginning and ending pft indices
    type(gridcell_type), pointer :: g 	 !local pointer to derived subtype
    type(landunit_type), pointer :: l 	 !local pointer to derived subtype
    type(column_type)  , pointer :: c 	 !local pointer to derived subtype
    type(pft_type)     , pointer :: p 	 !local pointer to derived subtype
    integer, parameter :: not_used = bigint
!------------------------------------------------------------------------

    ! Hourly average of 2m temperature.
	
    dtime = get_step_size()
    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600./dtime), &
         subgrid_type='pft', numlev=1, init_value=0.)

#if (defined DGVM)
    ! 30-day average of 2m temperature.

    call init_accum_field (name='TDA', units='K', &
         desc='30-day average of 2-m temperature', &
         accum_type='timeavg', accum_period=-30, &
         subgrid_type='pft', numlev=1, init_value=0.)

    ! The following are running means. 
    ! The accumulation period is set to 10 days for a 10-day running mean. 

    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20.)

    call init_accum_field (name='FNPSN10', units='UMOL/M2S', &
         desc='10-day running mean net cpy photosynth', & 
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0.)

    call init_accum_field (name='PREC365', units='MM H2O/S', &
         desc='365-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-365, &
         subgrid_type='pft', numlev=1, init_value=0.)

    ! The following are accumulated fields. 
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field). 
    ! Hence, [accper] is not valid.

    call init_accum_field (name='AGDD0', units='K', &
         desc='growing degree-days base 0C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0.)

    call init_accum_field (name='AGDD5', units='K', &
         desc='growing degree-days base -5C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0.)

    call init_accum_field (name='AGDDTW', units='K', &
         desc='growing degree-days base twmax', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0.)

    call init_accum_field (name='AGDD', units='K', &
         desc='growing degree-days base 5C', &
         accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0.)
#endif

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateAccFlds
!
! !INTERFACE:
  subroutine updateAccFlds()
!
! !DESCRIPTION: 
! Update and/or extract accumulated fields
!
! !USES:
    use clmtype
    use clmpoint
    use pftvarcon, only : pftpar
    use shr_const_mod, only: SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use time_manager, only : get_step_size, get_nstep, is_end_curr_day
    use clm_varcon, only : spval
    use accumulMod, only : update_accum_field, extract_accum_field
    use globals, only : ctl
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer :: gi,li,ci,pi               !indices
    integer :: begp, endp                !per-proc beginning and ending pft indices
    real(r8), pointer :: rbufslp(:)      !temporary single level - pft level
    type(gridcell_type), pointer :: g 	 !local pointer to derived subtype
    type(landunit_type), pointer :: l 	 !local pointer to derived subtype
    type(column_type)  , pointer :: c 	 !local pointer to derived subtype
    type(pft_type)     , pointer :: p 	 !local pointer to derived subtype
    integer :: itypveg                   !vegetation type
    integer :: dtime                     !timestep size [seconds]
    integer :: nstep                     !timestep number
    integer :: ier                       !error status    
!------------------------------------------------------------------------

    ! Determine time step size and index

    dtime = get_step_size()
    nstep = get_nstep()

    ! Don't do any accumulation if nstep is zero 
    ! (only applies to coupled or cam mode)

    if (nstep == 0) return

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Determine begnning and ending 1d indices

    begp = pfts1d%beg
    endp = pfts1d%end

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun
    endif

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to 
    ! accext if the time step does not correspond to the end of an 
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = p%pes%t_ref2m
    end do
    call update_accum_field  ('TREFAV', rbufslp, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       if (rbufslp(pi) /= spval) then
          p%pes%t_ref2m_max_inst = max(rbufslp(pi), p%pes%t_ref2m_max_inst)
          p%pes%t_ref2m_min_inst = min(rbufslp(pi), p%pes%t_ref2m_min_inst)
       endif
       if (is_end_curr_day()) then
          p%pes%t_ref2m_max = p%pes%t_ref2m_max_inst
          p%pes%t_ref2m_max_inst = -spval
          p%pes%t_ref2m_min = p%pes%t_ref2m_min_inst
          p%pes%t_ref2m_min_inst = spval
       else if (ctl%secs == int(dtime)) then
          p%pes%t_ref2m_max = spval
          p%pes%t_ref2m_min = spval
       endif
    end do

#if (defined DGVM)
    ! Accumulate and extract TDA 
    ! (accumulates TBOT as 30-day average) 
    ! Also determine t_mo_min 

!$OMP PARALLEL DO PRIVATE (pi,gi,g)
    do pi = begp,endp
       gi = pfts1d%gindex(pi)
       g => gpoint(gi)%g
       rbufslp(pi) = g%a2ls%forc_t
    end do
    call update_accum_field  ('TDA', rbufslp, nstep)
    call extract_accum_field ('TDA', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%t_mo = rbufslp(pi)
       p%pdgvs%t_mo_min = min(p%pdgvs%t_mo_min, rbufslp(pi))
    end do

    ! Accumulate and extract T10 
    !(acumulates TSA as 10-day running mean)

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = p%pes%t_ref2m
    end do
    call update_accum_field  ('T10', rbufslp, nstep)
    call extract_accum_field ('T10', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%t10 = rbufslp(pi)
    end do

    ! Accumulate and extract FNPSN10 
    !(accumulates fpsn-frmf as 10-day running mean)

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = p%pcf%fpsn - p%pcf%frmf
    end do
    call update_accum_field  ('FNPSN10', rbufslp, nstep)
    call extract_accum_field ('FNPSN10', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%fnpsn10 = rbufslp(pi)
    end do

    ! Accumulate and extract PREC365
    ! (accumulates total precipitation as 365-day running mean)

!$OMP PARALLEL DO PRIVATE (pi,gi,g)
    do pi = begp,endp
       gi = pfts1d%gindex(pi)
       g => gpoint(gi)%g
       rbufslp(pi) = g%a2lf%forc_rain + g%a2lf%forc_snow
    end do
    call update_accum_field  ('PREC365', rbufslp, nstep)
    call extract_accum_field ('PREC365', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%prec365 = rbufslp(pi)
    end do

    ! Accumulate growing degree days based on 10-day running mean temperature. 
    ! Accumulate GDD above 0C and -5C using extracted t10 from accumulated variable. 
    ! The trigger to reset the accumulated values to zero is -99999.
    ! agddtw is currently reset at the end of each year in subr. lpj

    ! Accumulate and extract AGDDO

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = (p%pdgvs%t10 - SHR_CONST_TKFRZ) * dtime / SHR_CONST_CDAY
       if (rbufslp(pi) < 0._r8) rbufslp(pi) = -99999.
    end do
    call update_accum_field  ('AGDD0', rbufslp, nstep)
    call extract_accum_field ('AGDD0', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agdd0 = rbufslp(pi)
    end do

    ! Accumulate and extract AGDD5

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = (p%pdgvs%t10 - (SHR_CONST_TKFRZ - 5.0))*dtime / SHR_CONST_CDAY
       if (rbufslp(pi) < 0._r8) rbufslp(pi) = -99999.
    end do
    call update_accum_field  ('AGDD5', rbufslp, nstep)
    call extract_accum_field ('AGDD5', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agdd5 = rbufslp(pi)
    end do

    ! Accumulate and extract AGDDTW

!$OMP PARALLEL DO PRIVATE (pi,p,itypveg)
    do pi = begp,endp
       p => ppoint(pi)%p
       itypveg = p%pps%itype
       rbufslp(pi) = max(0.0, (p%pdgvs%t10 - (SHR_CONST_TKFRZ+pftpar(itypveg,31))) &
            * dtime/SHR_CONST_CDAY) 
    end do
    call update_accum_field  ('AGDDTW', rbufslp, nstep)
    call extract_accum_field ('AGDDTW', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agddtw = rbufslp(pi) 
    end do

    ! Accumulate and extract AGDD

!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp
       p => ppoint(pi)%p
       rbufslp(pi) = max(0.0, (p%pes%t_ref2m - (SHR_CONST_TKFRZ + 5.0)) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDD', rbufslp, nstep)
    call extract_accum_field ('AGDD', rbufslp, nstep)
!$OMP PARALLEL DO PRIVATE (pi,p)
    do pi = begp,endp  
       p => ppoint(pi)%p
       p%pdgvs%agdd = rbufslp(pi) 
    end do
#endif

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccClmtype
!
! !INTERFACE:
  subroutine initAccClmtype
!
! !DESCRIPTION: 
! Initialize clmtype variables that are associated with 
! time accumulated fields. This routine is called in an initial run
! at nstep=0 for cam and csm mode and at nstep=1 for offline mode.  
! This routine is also always called for a restart run and 
! therefore must be called after the restart file is read in
! and the accumulated fields are obtained. 
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clmpoint
    use accumulMod, only : extract_accum_field
    use time_manager, only: get_nstep
    use clm_varctl, only : nsrest
    use clm_varcon, only : spval
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: pi                    !indices
    integer :: begp,endp             !beginning and ending subgrid indices
    integer :: nstep                 !time step
    integer :: ier                   !error status
    real(r8), pointer :: rbufslp(:)  !temporary
    type(pft_type), pointer :: p     !local pointer to derived subtype
!-----------------------------------------------------------------------

    ! Determine begnning and ending 1d indices

    begp = pfts1d%beg
    endp = pfts1d%end

    ! Determine time step

    nstep = get_nstep()

    ! Initialize 2m ref temperature max and min values

    if (nsrest == 0) then
       do pi = begp,endp
          p => ppoint(pi)%p
          p%pes%t_ref2m_max = spval
          p%pes%t_ref2m_min = spval
          p%pes%t_ref2m_max_inst = -spval
          p%pes%t_ref2m_min_inst =  spval
       end do
    end if

#if (defined DGVM)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbufslp'
       call endrun
    endif

    ! Initialize clmtype variables that are to be time accumulated 

    call extract_accum_field ('T10', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%t10 = rbufslp(pi)
    end do

    call extract_accum_field ('TDA', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%t_mo = rbufslp(pi)
    end do

    call extract_accum_field ('AGDD0', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agdd0 = rbufslp(pi)
    end do

    call extract_accum_field ('AGDD5', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agdd5 = rbufslp(pi)
    end do

    call extract_accum_field ('FNPSN10', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%fnpsn10 = rbufslp(pi)
    end do

    call extract_accum_field ('PREC365', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%prec365 = rbufslp(pi)
    end do

    call extract_accum_field ('AGDDTW', rbufslp, nstep)
    do pi = begp,endp
       p => ppoint(pi)%p
       p%pdgvs%agddtw = rbufslp(pi) 
    end do

    call extract_accum_field ('AGDD', rbufslp, nstep)
    do pi = begp,endp  
       p => ppoint(pi)%p
       p%pdgvs%agdd = rbufslp(pi) 
    end do

    deallocate(rbufslp)
#endif

  end subroutine initAccClmtype

end module accFldsMod
