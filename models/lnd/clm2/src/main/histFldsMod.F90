#include <misc.h>
#include <preproc.h>

module histFldsMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: histFldsMod
! 
! !DESCRIPTION: 
! Module containing initialization of clm2 history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history 
! fields.
! 
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private
  include "netcdf.inc"
!
! !PUBLIC TYPES
  public  :: initHistFlds ! Build master field list of all possible history file fields 
!                              
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
!
!EOP
!------------------------------------------------------------------------


contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initHistFlds
!
! !INTERFACE:
  subroutine initHistFlds ()
!
! !DESCRIPTION: 
! Build master field list of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute. A subroutine is called to add each field to the masterlist.
! The array hpindices has the form (gridcell_index, landunit_index, column_index,
! pft_index). The value of type1d massed to masterlist_add_fld determines which of
! the hpindices array values is used to determine the 1d type (gridcell, landunit,
! column or pft) of the output history buffer (i.e. it determines beg1d and end1d
! of the history buffer field). Set default history contents for given field on all 
! tapes by calling [masterlist_make_active] for the appropriatae tape.
! After the masterlist is build, routine [htapes_build] is called for an initial or
! branch run to initialize the actual history tapes.  
!
! !USES
    use clmpoint
    use clmtype
    use clm_varctl, only : nsrest
    use clm_varpar, only : nlevsoi
    use clm_varcon, only : spval
    use histFileMod, only : masterlist_addfld, masterlist_make_active, &
         masterlist_printflds, htapes_build, not_valid
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY
! Created by Mariana Vertenstein 
!
!EOP
!
! LOCAL VARIABLES
    integer :: pi,ci,li,gi   !indices 
    integer :: nf            !masterlist field counter
    integer :: hpindices(4)  !pointer indices into clmtype derived types
    integer :: begg, endg                !per-process beginning and ending gridcell indices
    integer :: begl, endl                !per-process beginning and ending land unit indices
    integer :: begc, endc                !per-process beginning and ending column indices 
    integer :: begp, endp                !per-process beginning and ending pft indices
    type(gridcell_type), pointer :: g 	 !local pointer to derived subtype
    type(landunit_type), pointer :: l 	 !local pointer to derived subtype
    type(column_type)  , pointer :: c 	 !local pointer to derived subtype
    type(pft_type)     , pointer :: p 	 !local pointer to derived subtype
!-----------------------------------------------------------------------

    ! Determine per-process beginning and ending indices

    begp = pfts1d%beg
    endp = pfts1d%end
    begc = cols1d%beg
    endc = cols1d%end
    begl = land1d%beg
    endl = land1d%end
    begg = grid1d%beg
    endg = grid1d%end

    ! Snow properties (will be vertically averaged over the snow profile)
    
    !---snow depth
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_snowdp == not_set) then
       ic_cps_snowdp = pointer_index()
       allocate (clmptr_rs(ic_cps_snowdp)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_snowdp)%val(ci)%p => c%cps%snowdp
       end do
    endif
    hpindices = (/-1, -1, ic_cps_snowdp, not_valid/)
    call masterlist_addfld (fname='SNOWDP', type1d='column', units='m', numlev=1, &
         avgflag='A', long_name='snow height', hpindices=hpindices)
    call masterlist_make_active (name='SNOWDP  ', tape_index=1)

    !---snow age
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_snowage == not_set) then
       ic_cps_snowage = pointer_index()
       allocate (clmptr_rs(ic_cps_snowage)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_snowage)%val(ci)%p => c%cps%snowage
       end do
    endif
    hpindices = (/-1, -1, ic_cps_snowage, not_valid/)
    call masterlist_addfld (fname='SNOWAGE', type1d='column', units='unitless', numlev=1, &
         avgflag='A', long_name='snow age', hpindices=hpindices)
    call masterlist_make_active (name='SNOWAGE ', tape_index=1)

    !---snow fraction
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_frac_sno == not_set) then
       ic_cps_frac_sno = pointer_index()
       allocate (clmptr_rs(ic_cps_frac_sno)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_frac_sno)%val(ci)%p => c%cps%frac_sno
       end do
    endif
    hpindices = (/-1, -1, ic_cps_frac_sno, not_valid/)
    call masterlist_addfld (fname='FSNO', type1d='column', units='unitless', numlev=1, &
         avgflag='A', long_name='fraction of ground covered by snow', hpindices=hpindices)
    call masterlist_make_active (name='FSNO ', tape_index=1)

    ! Temperatures
    
    !---2m air temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_pes_a_t_ref2m == not_set) then
       ic_ces_pes_a_t_ref2m = pointer_index()
       allocate (clmptr_rs(ic_ces_pes_a_t_ref2m)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_ces_pes_a_t_ref2m)%val(ci)%p => c%ces%pes_a%t_ref2m
       end do
    endif
    if (ip_pes_t_ref2m == not_set) then
       ip_pes_t_ref2m = pointer_index()
       allocate (clmptr_rs(ip_pes_t_ref2m)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pes_t_ref2m)%val(pi)%p => p%pes%t_ref2m
       end do
    endif
    hpindices = (/-1, -1, ic_ces_pes_a_t_ref2m, ip_pes_t_ref2m/)
    call masterlist_addfld (fname='TSA', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='2m air temperature', hpindices=hpindices)        
    call masterlist_make_active (name='TSA     ', tape_index=1) 

    !---daily minimum of average 2-m temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pes_t_ref2m_min == not_set) then
       ip_pes_t_ref2m_min = pointer_index()
       allocate (clmptr_rs(ip_pes_t_ref2m_min)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pes_t_ref2m_min)%val(pi)%p => p%pes%t_ref2m_min
       end do
    endif
    hpindices = (/-1, -1, -1, ip_pes_t_ref2m_min/)
    call masterlist_addfld (fname='TREFMNAV', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='daily minimum of average 2-m temperature', hpindices=hpindices)        
    call masterlist_make_active (name='TREFMNAV', tape_index=1) 

    !---daily maximum of average 2-m temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pes_t_ref2m_max == not_set) then
       ip_pes_t_ref2m_max = pointer_index()
       allocate (clmptr_rs(ip_pes_t_ref2m_max)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pes_t_ref2m_max)%val(pi)%p => p%pes%t_ref2m_max
       end do
    endif
    hpindices = (/-1, -1, -1, ip_pes_t_ref2m_max/)
    call masterlist_addfld (fname='TREFMXAV', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='daily maximum of average 2-m temperature', hpindices=hpindices)        
    call masterlist_make_active (name='TREFMXAV', tape_index=1) 

    !---vegetation temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_pes_a_t_veg == not_set) then
       ic_ces_pes_a_t_veg = pointer_index()
       allocate (clmptr_rs(ic_ces_pes_a_t_veg)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_ces_pes_a_t_veg)%val(ci)%p => c%ces%pes_a%t_veg
       end do
    endif
    if (ip_pes_t_veg == not_set) then
       ip_pes_t_veg = pointer_index()
       allocate (clmptr_rs(ip_pes_t_veg)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pes_t_veg)%val(pi)%p => p%pes%t_veg
       end do
    endif
    hpindices = (/-1, -1, ic_ces_pes_a_t_veg, ip_pes_t_veg/)
    call masterlist_addfld (fname='TV', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='vegetation temperature', hpindices=hpindices)    
    call masterlist_make_active (name='TV      ', tape_index=1) 
                                                                       
    !---ground temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_t_grnd == not_set) then
       ic_ces_t_grnd = pointer_index()
       allocate (clmptr_rs(ic_ces_t_grnd)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_ces_t_grnd)%val(ci)%p => c%ces%t_grnd
       end do
    end if
    hpindices = (/-1, -1, ic_ces_t_grnd, not_valid/)
    call masterlist_addfld (fname='TG', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='ground temperature', hpindices=hpindices)        
    call masterlist_make_active (name='TG      ', tape_index=1) 
                                                                       
    !---snow temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_t_snow == not_set) then
       ic_ces_t_snow = pointer_index()
       allocate (clmptr_rs(ic_ces_t_snow)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_ces_t_snow)%val(ci)%p => c%ces%t_snow
          if (c%l%lps%lakpoi) c%ces%t_snow = spval
       end do
    end if
    hpindices = (/-1, -1, ic_ces_t_snow, not_valid/)
    call masterlist_addfld (fname='TSNOW', type1d='column', units='K', numlev=1, &
         avgflag='A', long_name='snow temperature', hpindices=hpindices)          
    call masterlist_make_active (name='TSNOW   ', tape_index=1) 
                                                                       
    !---soil temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_t_soisno == not_set) then
       ic_ces_t_soisno = pointer_index()
       allocate (clmptr_ra(ic_ces_t_soisno)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_ra(ic_ces_t_soisno)%val(ci)%p => c%ces%t_soisno
       end do
    end if
    hpindices = (/-1, -1, ic_ces_t_soisno, not_valid/)
    call masterlist_addfld (fname='TSOI', type1d='column', units='K', numlev=nlevsoi, &
         avgflag='A', long_name='soil temperature', hpindices=hpindices)    
    call masterlist_make_active (name='TSOI    ', tape_index=1) 
                                                                       
    !---lake temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_t_lake == not_set) then
       ic_ces_t_lake = pointer_index()
       allocate (clmptr_ra(ic_ces_t_lake)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_ra(ic_ces_t_lake)%val(ci)%p => c%ces%t_lake
       end do
    end if
    hpindices = (/-1, -1, ic_ces_t_lake, not_valid/)
    call masterlist_addfld (fname='TLAKE', type1d='column', units='K', numlev=nlevsoi, &
         avgflag='A', long_name='lake temperature', hpindices=hpindices)
    call masterlist_make_active (name='TLAKE   ', tape_index=1) 

    ! Specific humidity
    
    !---2m specific humidity
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_ces_pes_a_q_ref2m == not_set) then
       ic_ces_pes_a_q_ref2m = pointer_index()
       allocate (clmptr_rs(ic_ces_pes_a_q_ref2m)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_ces_pes_a_q_ref2m)%val(ci)%p => c%ces%pes_a%q_ref2m
       end do
    endif
    if (ip_pes_q_ref2m == not_set) then
       ip_pes_q_ref2m = pointer_index()
       allocate (clmptr_rs(ip_pes_q_ref2m)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pes_q_ref2m)%val(pi)%p => p%pes%q_ref2m
       end do
    endif
    hpindices = (/-1, -1, ic_ces_pes_a_q_ref2m, ip_pes_q_ref2m/)
    call masterlist_addfld (fname='Q2M', type1d='column', units='kg/kg', numlev=1, &
         avgflag='A', long_name='2m specific humidity', hpindices=hpindices)        
    call masterlist_make_active (name='Q2M     ', tape_index=1) 

    ! Surface radiation                                          

    !---solar rad absorbed by veg
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_sabv == not_set) then
       ic_cef_pef_a_sabv = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_sabv)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_sabv)%val(ci)%p => c%cef%pef_a%sabv
          if (c%l%lps%lakpoi) c%cef%pef_a%sabv = 0.
       end do
    end if
    if (ip_pef_sabv == not_set) then
       ip_pef_sabv = pointer_index()
       allocate (clmptr_rs(ip_pef_sabv)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_sabv)%val(pi)%p => p%pef%sabv
          if (p%c%l%lps%lakpoi) p%pef%sabv = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_sabv, ip_pef_sabv/)
    call masterlist_addfld (fname='SABV', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='solar rad absorbed by veg', hpindices=hpindices)
    call masterlist_make_active (name='SABV    ', tape_index=1)

    !---solar rad absorbed by ground
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_sabg == not_set) then
       ic_cef_pef_a_sabg = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_sabg)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_sabg)%val(ci)%p => c%cef%pef_a%sabg
       end do
    end if
    if (ip_pef_sabg == not_set) then
       ip_pef_sabg = pointer_index()
       allocate (clmptr_rs(ip_pef_sabg)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_sabg)%val(pi)%p => p%pef%sabg
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_sabg, ip_pef_sabg/)
    call masterlist_addfld (fname='SABG', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='solar rad absorbed by ground', hpindices=hpindices)
    call masterlist_make_active (name='SABG    ', tape_index=1)

    !---direct vis incident solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsds_vis_d == not_set) then
       ic_cef_pef_a_fsds_vis_d = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsds_vis_d)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsds_vis_d)%val(ci)%p => c%cef%pef_a%fsds_vis_d
       end do
    end if
    if (ip_pef_fsds_vis_d == not_set) then
       ip_pef_fsds_vis_d = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_vis_d)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_vis_d)%val(pi)%p => p%pef%fsds_vis_d
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsds_vis_d, ip_pef_fsds_vis_d/)
    call masterlist_addfld (fname='FSDSVD', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct vis incident solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSDSVD  ', tape_index=1)

    !---direct nir incident solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsds_nir_d == not_set) then
       ic_cef_pef_a_fsds_nir_d = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsds_nir_d)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsds_nir_d)%val(ci)%p => c%cef%pef_a%fsds_nir_d
       end do
    end if
    if (ip_pef_fsds_nir_d == not_set) then
       ip_pef_fsds_nir_d = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_nir_d)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_nir_d)%val(pi)%p => p%pef%fsds_nir_d
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsds_nir_d, ip_pef_fsds_nir_d/)
    call masterlist_addfld (fname='FSDSND', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct nir incident solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSDSND  ', tape_index=1)

    !---diffuse vis incident solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsds_vis_i == not_set) then
       ic_cef_pef_a_fsds_vis_i = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsds_vis_i)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsds_vis_i)%val(ci)%p => c%cef%pef_a%fsds_vis_i
       end do
    end if
    if (ip_pef_fsds_vis_i == not_set) then
       ip_pef_fsds_vis_i = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_vis_i)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_vis_i)%val(pi)%p => p%pef%fsds_vis_i
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsds_vis_i, ip_pef_fsds_vis_i/)
    call masterlist_addfld (fname='FSDSVI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='diffuse vis incident solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSDSVI  ', tape_index=1)

    !---diffuse nir incident solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsds_nir_i == not_set) then
       ic_cef_pef_a_fsds_nir_i = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsds_nir_i)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsds_nir_i)%val(ci)%p => c%cef%pef_a%fsds_nir_i
       end do
    end if
    if (ip_pef_fsds_nir_i == not_set) then
       ip_pef_fsds_nir_i = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_nir_i)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_nir_i)%val(pi)%p => p%pef%fsds_nir_i
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsds_nir_i, ip_pef_fsds_nir_i/)
    call masterlist_addfld (fname='FSDSNI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='diffuse nir incident solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSDSNI  ', tape_index=1)

    !---direct vis reflected solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsr_vis_d == not_set) then
       ic_cef_pef_a_fsr_vis_d = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsr_vis_d)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsr_vis_d)%val(ci)%p => c%cef%pef_a%fsr_vis_d
       end do
    end if
    if (ip_pef_fsr_vis_d == not_set) then
       ip_pef_fsr_vis_d = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_vis_d)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_vis_d)%val(pi)%p => p%pef%fsr_vis_d
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsr_vis_d, ip_pef_fsr_vis_d/)
    call masterlist_addfld (fname='FSRVD', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct vis reflected solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSRVD   ', tape_index=1)

    !---direct nir reflected solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsr_nir_d == not_set) then
       ic_cef_pef_a_fsr_nir_d = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsr_nir_d)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsr_nir_d)%val(ci)%p => c%cef%pef_a%fsr_nir_d
       end do
    end if
    if (ip_pef_fsr_nir_d == not_set) then
       ip_pef_fsr_nir_d = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_nir_d)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_nir_d)%val(pi)%p => p%pef%fsr_nir_d
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsr_nir_d, ip_pef_fsr_nir_d/)
    call masterlist_addfld (fname='FSRND', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct nir reflected solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSRND   ', tape_index=1)

    !---diffuse vis reflected solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsr_vis_i == not_set) then
       ic_cef_pef_a_fsr_vis_i = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsr_vis_i)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsr_vis_i)%val(ci)%p => c%cef%pef_a%fsr_vis_i
       end do
    end if
    if (ip_pef_fsr_vis_i == not_set) then
       ip_pef_fsr_vis_i = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_vis_i)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_vis_i)%val(pi)%p => p%pef%fsr_vis_i
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsr_vis_i, ip_pef_fsr_vis_i/)
    call masterlist_addfld (fname='FSRVI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='diffuse vis reflected solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSRVI   ', tape_index=1)

    !---diffuse nir reflected solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsr_nir_i == not_set) then
       ic_cef_pef_a_fsr_nir_i = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsr_nir_i)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsr_nir_i)%val(ci)%p => c%cef%pef_a%fsr_nir_i
       end do
    end if
    if (ip_pef_fsr_nir_i == not_set) then
       ip_pef_fsr_nir_i = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_nir_i)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_nir_i)%val(pi)%p => p%pef%fsr_nir_i
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsr_nir_i, ip_pef_fsr_nir_i/)
    call masterlist_addfld (fname='FSRNI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='diffuse nir reflected solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSRNI   ', tape_index=1)

    !---direct vis incident solar radiation at local noon
    !---when defining pointers, set constant values (e.g. over lake) if needed 
    if (ip_pef_fsds_vis_d_ln == not_set) then
       ip_pef_fsds_vis_d_ln = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_vis_d_ln)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_vis_d_ln)%val(pi)%p => p%pef%fsds_vis_d_ln
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pef_fsds_vis_d_ln/)
    call masterlist_addfld (fname='FSDSVDLN', type1d='pft', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', hpindices=hpindices)
    call masterlist_make_active (name='FSDSVDLN', tape_index=1)

    !---direct nir incident solar radiation at local noon
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pef_fsds_nir_d_ln == not_set) then
       ip_pef_fsds_nir_d_ln = pointer_index()
       allocate (clmptr_rs(ip_pef_fsds_nir_d_ln)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsds_nir_d_ln)%val(pi)%p => p%pef%fsds_nir_d_ln
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pef_fsds_nir_d_ln/)
    call masterlist_addfld (fname='FSDSNDLN', type1d='pft', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', hpindices=hpindices)
    call masterlist_make_active (name='FSDSNDLN', tape_index=1)

    !---direct vis reflected solar radiation at local noon
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pef_fsr_vis_d_ln == not_set) then
       ip_pef_fsr_vis_d_ln = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_vis_d_ln)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_vis_d_ln)%val(pi)%p => p%pef%fsr_vis_d_ln
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pef_fsr_vis_d_ln/)
    call masterlist_addfld (fname='FSRVDLN', type1d='pft', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', hpindices=hpindices)
    call masterlist_make_active (name='FSRVDLN ', tape_index=1)

    !---direct nir reflected solar radiation at local noon
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pef_fsr_nir_d_ln == not_set) then
       ip_pef_fsr_nir_d_ln = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr_nir_d_ln)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr_nir_d_ln)%val(pi)%p => p%pef%fsr_nir_d_ln
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pef_fsr_nir_d_ln/)
    call masterlist_addfld (fname='FSRNDLN', type1d='pft', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', hpindices=hpindices)
    call masterlist_make_active (name='FSRNDLN ', tape_index=1)

    !---absobed solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsa == not_set) then
       ic_cef_pef_a_fsa = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsa)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsa)%val(ci)%p => c%cef%pef_a%fsa
       end do
    end if
    if (ip_pef_fsa == not_set) then
       ip_pef_fsa = pointer_index()
       allocate (clmptr_rs(ip_pef_fsa)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsa)%val(pi)%p => p%pef%fsa
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsa, ip_pef_fsa/)
    call masterlist_addfld (fname='FSA', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='absobed solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSA     ', tape_index=1) 

    !---reflected solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_fsr == not_set) then
       ic_cef_pef_a_fsr = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_fsr)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_fsr)%val(ci)%p => c%cef%pef_a%fsr
       end do
    end if
    if (ip_pef_fsr == not_set) then
       ip_pef_fsr = pointer_index()
       allocate (clmptr_rs(ip_pef_fsr)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_fsr)%val(pi)%p => p%pef%fsr
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_fsr, ip_pef_fsr/)
    call masterlist_addfld (fname='FSR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='reflected solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSR     ', tape_index=1) 

    !---surface ndvi
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_pps_a_ndvi == not_set) then
       ic_cps_pps_a_ndvi = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_ndvi)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_ndvi)%val(ci)%p => c%cps%pps_a%ndvi
       end do
    end if
    if (ip_pps_ndvi == not_set) then
       ip_pps_ndvi = pointer_index()
       allocate (clmptr_rs(ip_pps_ndvi)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_ndvi)%val(pi)%p => p%pps%ndvi
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_ndvi, ip_pps_ndvi/)
    call masterlist_addfld (fname='NDVI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='surface ndvi', hpindices=hpindices)
    call masterlist_make_active (name='NDVI    ', tape_index=1) 

    !---net infrared (longwave) radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_lwrad_net == not_set) then
       ic_cef_pef_a_eflx_lwrad_net = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_lwrad_net)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_lwrad_net)%val(ci)%p => c%cef%pef_a%eflx_lwrad_net
       end do
    end if
    if (ip_pef_eflx_lwrad_net == not_set) then
       ip_pef_eflx_lwrad_net = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_lwrad_net)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_lwrad_net)%val(pi)%p => p%pef%eflx_lwrad_net
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lwrad_net, ip_pef_eflx_lwrad_net/)
    call masterlist_addfld (fname='FIRA', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='net infrared (longwave) radiation', hpindices=hpindices)
    call masterlist_make_active (name='FIRA    ', tape_index=1) 

    !---emitted infrared (longwave) radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_lwrad_out == not_set) then
       ic_cef_pef_a_eflx_lwrad_out = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_lwrad_out)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_lwrad_out)%val(ci)%p => c%cef%pef_a%eflx_lwrad_out
       end do
    end if
    if (ip_pef_eflx_lwrad_out == not_set) then
       ip_pef_eflx_lwrad_out = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_lwrad_out)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_lwrad_out)%val(pi)%p => p%pef%eflx_lwrad_out
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lwrad_out, ip_pef_eflx_lwrad_out/)
    call masterlist_addfld (fname='FIRE', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='emitted infrared (longwave) radiation', hpindices=hpindices)
    call masterlist_make_active (name='FIRE    ', tape_index=1) 

    ! Surface energy fluxes                                      

    !---canopy transpiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_lh_vegt == not_set) then
       ic_cef_pef_a_eflx_lh_vegt = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_lh_vegt)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_lh_vegt)%val(ci)%p => c%cef%pef_a%eflx_lh_vegt
          if (c%l%lps%lakpoi) c%cef%pef_a%eflx_lh_vegt = 0.
       end do
    end if
    if (ip_pef_eflx_lh_vegt == not_set) then
       ip_pef_eflx_lh_vegt = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_lh_vegt)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_lh_vegt)%val(pi)%p => p%pef%eflx_lh_vegt
          if (p%c%l%lps%lakpoi) p%pef%eflx_lh_vegt = 0. 
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_vegt, ip_pef_eflx_lh_vegt/)
    call masterlist_addfld (fname='FCTR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='canopy transpiration', hpindices=hpindices)
    call masterlist_make_active (name='FCTR    ', tape_index=1) 

    !---canopy evaporation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_lh_vege == not_set) then
       ic_cef_pef_a_eflx_lh_vege = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_lh_vege)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_lh_vege)%val(ci)%p => c%cef%pef_a%eflx_lh_vege
          if (c%l%lps%lakpoi) c%cef%pef_a%eflx_lh_vege = 0.
       end do
    end if
    if (ip_pef_eflx_lh_vege == not_set) then
       ip_pef_eflx_lh_vege = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_lh_vege)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_lh_vege)%val(pi)%p => p%pef%eflx_lh_vege
          if (p%c%l%lps%lakpoi) p%pef%eflx_lh_vege = 0. 
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_vege, ip_pef_eflx_lh_vege/)
    call masterlist_addfld (fname='FCEV', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='canopy evaporation', hpindices=hpindices)
    call masterlist_make_active (name='FCEV    ', tape_index=1) 

    !---ground evaporation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_lh_grnd == not_set) then
       ic_cef_pef_a_eflx_lh_grnd = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_lh_grnd)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_lh_grnd)%val(ci)%p => c%cef%pef_a%eflx_lh_grnd
       end do
    end if
    if (ip_pef_eflx_lh_grnd == not_set) then
       ip_pef_eflx_lh_grnd = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_lh_grnd)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_lh_grnd)%val(pi)%p => p%pef%eflx_lh_grnd
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_lh_grnd, ip_pef_eflx_lh_grnd/)
    call masterlist_addfld (fname='FGEV', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='ground evaporation', hpindices=hpindices)
    call masterlist_make_active (name='FGEV    ', tape_index=1) 

    !---sensible heat
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_sh_tot == not_set) then
       ic_cef_pef_a_eflx_sh_tot = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_sh_tot)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_sh_tot)%val(ci)%p => c%cef%pef_a%eflx_sh_tot
       end do
    end if
    if (ip_pef_eflx_sh_tot == not_set) then
       ip_pef_eflx_sh_tot = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_sh_tot)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_sh_tot)%val(pi)%p => p%pef%eflx_sh_tot
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_sh_tot, ip_pef_eflx_sh_tot/)
    call masterlist_addfld (fname='FSH', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='sensible heat', hpindices=hpindices)
    call masterlist_make_active (name='FSH     ', tape_index=1) 

    !---sensible heat from veg
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_sh_veg == not_set) then
       ic_cef_pef_a_eflx_sh_veg = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_sh_veg)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_sh_veg)%val(ci)%p => c%cef%pef_a%eflx_sh_veg
          if (c%l%lps%lakpoi) c%cef%pef_a%eflx_sh_veg = 0.
       end do
    end if
    if (ip_pef_eflx_sh_veg == not_set) then
       ip_pef_eflx_sh_veg = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_sh_veg)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_sh_veg)%val(pi)%p => p%pef%eflx_sh_veg
          if (p%c%l%lps%lakpoi) p%pef%eflx_sh_veg = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_sh_veg, ip_pef_eflx_sh_veg/)
    call masterlist_addfld (fname='FSH_V', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='sensible heat from veg', hpindices=hpindices)
    call masterlist_make_active (name='FSH_V   ', tape_index=1)

    !---sensible heat from ground
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_sh_grnd == not_set) then
       ic_cef_pef_a_eflx_sh_grnd = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_sh_grnd)%val(begc:endc))
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_sh_grnd)%val(ci)%p => c%cef%pef_a%eflx_sh_grnd
       end do
    end if
    if (ip_pef_eflx_sh_grnd == not_set) then
       ip_pef_eflx_sh_grnd = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_sh_grnd)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_sh_grnd)%val(pi)%p => p%pef%eflx_sh_grnd
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_sh_grnd, ip_pef_eflx_sh_grnd/)
    call masterlist_addfld (fname='FSH_G', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='sensible heat from ground', hpindices=hpindices)
    call masterlist_make_active (name='FSH_G   ', tape_index=1)

    !---heat flux into soil
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_pef_a_eflx_soil_grnd == not_set) then
       ic_cef_pef_a_eflx_soil_grnd = pointer_index()
       allocate (clmptr_rs(ic_cef_pef_a_eflx_soil_grnd)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_pef_a_eflx_soil_grnd)%val(ci)%p => c%cef%pef_a%eflx_soil_grnd
       end do
    end if
    if (ip_pef_eflx_soil_grnd == not_set) then
       ip_pef_eflx_soil_grnd = pointer_index()
       allocate (clmptr_rs(ip_pef_eflx_soil_grnd)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pef_eflx_soil_grnd)%val(pi)%p => p%pef%eflx_soil_grnd
       end do
    end if
    hpindices = (/-1, -1, ic_cef_pef_a_eflx_soil_grnd, ip_pef_eflx_soil_grnd/)
    call masterlist_addfld (fname='FGR', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='heat flux into soil', hpindices=hpindices)
    call masterlist_make_active (name='FGR     ', tape_index=1) 

    !---snow melt heat flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cef_eflx_snomelt == not_set) then
       ic_cef_eflx_snomelt = pointer_index()
       allocate (clmptr_rs(ic_cef_eflx_snomelt)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cef_eflx_snomelt)%val(ci)%p => c%cef%eflx_snomelt
       end do
    end if
    hpindices = (/-1, -1, ic_cef_eflx_snomelt, not_valid/)
    call masterlist_addfld (fname='FSM', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='snow melt heat flux', hpindices=hpindices)
    call masterlist_make_active (name='FSM     ', tape_index=1) 

    !---zonal surface stress
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cmf_pmf_a_taux == not_set) then
       ic_cmf_pmf_a_taux = pointer_index()
       allocate (clmptr_rs(ic_cmf_pmf_a_taux)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cmf_pmf_a_taux)%val(ci)%p => c%cmf%pmf_a%taux
       end do
    end if
    if (ip_pmf_taux == not_set) then
       ip_pmf_taux = pointer_index()
       allocate (clmptr_rs(ip_pmf_taux)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pmf_taux)%val(pi)%p => p%pmf%taux
       end do
    end if
    hpindices = (/-1, -1, ic_cmf_pmf_a_taux, ip_pmf_taux/)
    call masterlist_addfld (fname='TAUX', type1d='column', units='kg/m/s^2', numlev=1, &
         avgflag='A', long_name='zonal surface stress', hpindices=hpindices)
    call masterlist_make_active (name='TAUX    ', tape_index=1) 

    !---meridional surface stress
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cmf_pmf_a_tauy == not_set) then
       ic_cmf_pmf_a_tauy = pointer_index()
       allocate (clmptr_rs(ic_cmf_pmf_a_tauy)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cmf_pmf_a_tauy)%val(ci)%p => c%cmf%pmf_a%tauy
       end do
    end if
    if (ip_pmf_tauy == not_set) then
       ip_pmf_tauy = pointer_index()
       allocate (clmptr_rs(ip_pmf_tauy)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pmf_tauy)%val(pi)%p => p%pmf%tauy
       end do
    end if
    hpindices = (/-1, -1, ic_cmf_pmf_a_tauy, ip_pmf_tauy/)
    call masterlist_addfld (fname='TAUY', type1d='column', units='kg/m/s^2', numlev=1, &
         avgflag='A', long_name='meridional surface stress', hpindices=hpindices)
    call masterlist_make_active (name='TAUY    ', tape_index=1) 

    ! Vegetation phenology

    !---exposed one-sided leaf area index
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_pps_a_elai == not_set) then
       ic_cps_pps_a_elai = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_elai)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_elai)%val(ci)%p => c%cps%pps_a%elai
       end do
    end if
    if (ip_pps_elai == not_set) then
       ip_pps_elai = pointer_index()
       allocate (clmptr_rs(ip_pps_elai)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_elai)%val(pi)%p => p%pps%elai
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_elai, ip_pps_elai/)
    call masterlist_addfld (fname='ELAI', type1d='column', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='exposed one-sided leaf area index', hpindices=hpindices)
    call masterlist_make_active (name='ELAI    ', tape_index=1) 

    !---exposed one-sided stem area index
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_pps_a_esai == not_set) then
       ic_cps_pps_a_esai = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_esai)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_esai)%val(ci)%p => c%cps%pps_a%esai
       end do
    end if
    if (ip_pps_esai == not_set) then
       ip_pps_esai = pointer_index()
       allocate (clmptr_rs(ip_pps_esai)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_esai)%val(pi)%p => p%pps%esai
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_esai, ip_pps_esai/)
    call masterlist_addfld (fname='ESAI', type1d='column', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='exposed one-sided stem area index', hpindices=hpindices)
    call masterlist_make_active (name='ESAI    ', tape_index=1) 

    ! Canopy physiology                                          

    !---sunlit leaf stomatal resistance
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_pps_a_rssun == not_set) then
       ic_cps_pps_a_rssun = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_rssun)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_rssun)%val(ci)%p => c%cps%pps_a%rssun
          if (c%l%lps%lakpoi) c%cps%pps_a%rssun = spval
       end do
    end if
    if (ip_pps_rssun == not_set) then
       ip_pps_rssun = pointer_index()
       allocate (clmptr_rs(ip_pps_rssun)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_rssun)%val(pi)%p => p%pps%rssun
          if (p%c%l%lps%lakpoi) p%pps%rssun = spval
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_rssun, ip_pps_rssun/)
    call masterlist_addfld (fname='RSSUN', type1d='column', units='s/m', numlev=1, &
         avgflag='M', long_name='sunlit leaf stomatal resistance', hpindices=hpindices)
    call masterlist_make_active (name='RSSUN   ', tape_index=1) 

    !---shaded leaf stomatal resistance
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cps_pps_a_rssha == not_set) then
       ic_cps_pps_a_rssha = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_rssha)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_rssha)%val(ci)%p => c%cps%pps_a%rssha
          if (c%l%lps%lakpoi) c%cps%pps_a%rssha = spval
       end do
    end if
    if (ip_pps_rssha == not_set) then
       ip_pps_rssha = pointer_index()
       allocate (clmptr_rs(ip_pps_rssha)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_rssha)%val(pi)%p => p%pps%rssha
          if (p%c%l%lps%lakpoi) p%pps%rssha = spval
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_rssha, ip_pps_rssha/)
    call masterlist_addfld (fname='RSSHA', type1d='column', units='s/m', numlev=1, &
         avgflag='M', long_name='shaded leaf stomatal resistance', hpindices=hpindices)
    call masterlist_make_active (name='RSSHA   ', tape_index=1) 

    !---transpiration beta factor
    !---when defining pointers, set constant values (e.g. over lake) if neede
    if (ic_cps_pps_a_btran == not_set) then
       ic_cps_pps_a_btran = pointer_index()
       allocate (clmptr_rs(ic_cps_pps_a_btran)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cps_pps_a_btran)%val(ci)%p => c%cps%pps_a%btran
          if (c%l%lps%lakpoi) c%cps%pps_a%btran = spval
       end do
    end if
    if (ip_pps_btran == not_set) then
       ip_pps_btran = pointer_index()
       allocate (clmptr_rs(ip_pps_btran)%val(begp:endp))
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_btran)%val(pi)%p => p%pps%btran
          if (p%c%l%lps%lakpoi) p%pps%btran = spval
       end do
    end if
    hpindices = (/-1, -1, ic_cps_pps_a_btran, ip_pps_btran/)
    call masterlist_addfld (fname='BTRAN', type1d='column', units='s/m', numlev=1, &
         avgflag='A', long_name='transpiration beta factor', hpindices=hpindices)
    call masterlist_make_active (name='BTRAN   ', tape_index=1) 

    !---photosynthesis
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_fpsn == not_set) then
       ip_pcf_fpsn = pointer_index()
       allocate (clmptr_rs(ip_pcf_fpsn)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_fpsn)%val(pi)%p => p%pcf%fpsn
          if (p%c%l%lps%lakpoi) p%pcf%fpsn = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_fpsn/)
    call masterlist_addfld (fname='FPSN', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='photosynthesis', hpindices=hpindices)
    call masterlist_make_active (name='FPSN    ', tape_index=1) 

#if (defined BGC)
    ! Biogeochemical fluxes                                      

    !---microbial respiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_fmicr == not_set) then
       ip_pcf_fmicr = pointer_index()
       allocate (clmptr_rs(ip_pcf_fmicr)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_fmicr)%val(pi)%p => p%pcf%fmicr
          if (p%c%l%lps%lakpoi) p%pcf%fmicr = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_fmicr/)
    call masterlist_addfld (fname='FMICR', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='microbial respiration', hpindices=hpindices)
    call masterlist_make_active (name='FMICR   ', tape_index=1) 

    !---stem maintenance respiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_frms == not_set) then
       ip_pcf_frms = pointer_index()
       allocate (clmptr_rs(ip_pcf_frms)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_frms)%val(pi)%p => p%pcf%frms
          if (p%c%l%lps%lakpoi) p%pcf%frms = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_frms/)
    call masterlist_addfld (fname='FRMS', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='stem maintenance respiration', hpindices=hpindices)
    call masterlist_make_active (name='FRMS    ', tape_index=1) 

    !---root maintenance respiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_frmr == not_set) then
       ip_pcf_frmr = pointer_index()
       allocate (clmptr_rs(ip_pcf_frmr)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_frmr)%val(pi)%p => p%pcf%frmr
          if (p%c%l%lps%lakpoi) p%pcf%frmr = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_frmr/)
    call masterlist_addfld (fname='FRMR', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='root maintenance respiration', hpindices=hpindices)
    call masterlist_make_active (name='FRMR    ', tape_index=1) 

    !---foliage maintenance respiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_frmf == not_set) then
       ip_pcf_frmf = pointer_index()
       allocate (clmptr_rs(ip_pcf_frmf)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_frmf)%val(pi)%p => p%pcf%frmf
          if (p%c%l%lps%lakpoi) p%pcf%frmf = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_frmf/)
    call masterlist_addfld (fname='FRMF', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='foliage maintenance respiration', hpindices=hpindices)
    call masterlist_make_active (name='FRMF    ', tape_index=1) 

    !---growth respiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_frg == not_set) then
       ip_pcf_frg = pointer_index()
       allocate (clmptr_rs(ip_pcf_frg)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_frg)%val(pi)%p => p%pcf%frg
          if (p%c%l%lps%lakpoi) p%pcf%frg = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_frg/)
    call masterlist_addfld (fname='FRG', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='growth respiration', hpindices=hpindices)
    call masterlist_make_active (name='FRG     ', tape_index=1) 

    !---net CO2 flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_fco2 == not_set) then
       ip_pcf_fco2 = pointer_index()
       allocate (clmptr_rs(ip_pcf_fco2)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_fco2)%val(pi)%p => p%pcf%fco2
          if (p%c%l%lps%lakpoi) p%pcf%fco2 = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_fco2/)
    call masterlist_addfld (fname='FCO2', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='net CO2 flux', hpindices=hpindices)
    call masterlist_make_active (name='FCO2    ', tape_index=1) 

    !---net primary production
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pcf_dmi == not_set) then
       ip_pcf_dmi = pointer_index()
       allocate (clmptr_rs(ip_pcf_dmi)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pcf_dmi)%val(pi)%p => p%pcf%dmi
          if (p%c%l%lps%lakpoi) p%pcf%dmi = 0. 
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pcf_dmi/)
    call masterlist_addfld (fname='DMI', type1d='pft', units='umol/m2s', numlev=1, &
         avgflag='A', long_name='net primary production', hpindices=hpindices)
    call masterlist_make_active (name='DMI     ', tape_index=1) 

    !---net primary production
    !---when defining pointers, set constant values (e.g. over lake) if needed
    !---for this case also need to set constant values over lake leading to the
    !---computation of flx_mss_vrt_dst_tot
    if (ip_pdf_flx_mss_vrt_dst_tot == not_set) then
       ip_pdf_flx_mss_vrt_dst_tot = pointer_index()
       allocate (clmptr_rs(ip_pdf_flx_mss_vrt_dst_tot)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdf_flx_mss_vrt_dst_tot)%val(pi)%p => p%pdf%flx_mss_vrt_dst_tot
          if (p%c%l%lps%lakpoi) then
             p%pdf%flx_mss_vrt_dst_tot = 0.
             p%pdf%flx_mss_vrt_dst(1) = 0.
             p%pdf%flx_mss_vrt_dst(2) = 0.
             p%pdf%flx_mss_vrt_dst(3) = 0.
             p%pdf%flx_mss_vrt_dst(4) = 0.
          endif
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdf_flx_mss_vrt_dst_tot/)
    call masterlist_addfld (fname='DSTFLXT', type1d='pft', units='kg/m2/s', numlev=1, &
         avgflag='A', long_name='net primary production', hpindices=hpindices)
    call masterlist_make_active (name='DSTFLXT ', tape_index=1) 

    !---turbulent deposition velocity 1
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdf_vlc_trb_1 == not_set) then
       ip_pdf_vlc_trb_1 = pointer_index()
       allocate (clmptr_rs(ip_pdf_vlc_trb_1)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdf_vlc_trb_1)%val(pi)%p => p%pdf%vlc_trb_1
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdf_vlc_trb_1/)
    call masterlist_addfld (fname='DPVLTRB1', type1d='pft', units='m/s', numlev=1, &
         avgflag='A', long_name='turbulent deposition velocity 1', hpindices=hpindices)
    call masterlist_make_active (name='DPVLTRB1', tape_index=1) 

    !---turbulent deposition velocity 2
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdf_vlc_trb_2 == not_set) then
       ip_pdf_vlc_trb_2 = pointer_index()
       allocate (clmptr_rs(ip_pdf_vlc_trb_2)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdf_vlc_trb_2)%val(pi)%p => p%pdf%vlc_trb_2
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdf_vlc_trb_2/)
    call masterlist_addfld (fname='DPVLTRB2', type1d='pft', units='m/s', numlev=1, &
         avgflag='A', long_name='turbulent deposition velocity 2', hpindices=hpindices)
    call masterlist_make_active (name='DPVLTRB2', tape_index=1) 

    !---turbulent deposition velocity 3
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdf_vlc_trb_3 == not_set) then
       ip_pdf_vlc_trb_3 = pointer_index()
       allocate (clmptr_rs(ip_pdf_vlc_trb_3)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdf_vlc_trb_3)%val(pi)%p => p%pdf%vlc_trb_3
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdf_vlc_trb_3/)
    call masterlist_addfld (fname='DPVLTRB3', type1d='pft', units='m/s', numlev=1, &
         avgflag='A', long_name='turbulent deposition velocity 3', hpindices=hpindices)
    call masterlist_make_active (name='DPVLTRB3', tape_index=1) 

    !---turbulent deposition velocity 4
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdf_vlc_trb_4 == not_set) then
       ip_pdf_vlc_trb_4 = pointer_index()
       allocate (clmptr_rs(ip_pdf_vlc_trb_4)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdf_vlc_trb_4)%val(pi)%p => p%pdf%vlc_trb_4
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdf_vlc_trb_4/)
    call masterlist_addfld (fname='DPVLTRB4', type1d='pft', units='m/s', numlev=1, &
         avgflag='A', long_name='turbulent deposition velocity 4', hpindices=hpindices)
    call masterlist_make_active (name='DPVLTRB4', tape_index=1) 

    !---total VOC flux into atmosphere
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_tot == not_set) then
       ip_pvf_vocflx_tot = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_tot)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_tot)%val(pi)%p => p%pvf%vocflx_tot
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_tot = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_tot/)
    call masterlist_addfld (fname='VOCFLXT', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='total VOC flux into atmosphere', hpindices=hpindices)
    call masterlist_make_active (name='VOCFLXT ', tape_index=1) 

    !---isoprene flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_1 == not_set) then
       ip_pvf_vocflx_1 = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_1)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_1)%val(pi)%p => p%pvf%vocflx_1
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_1 = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_1/)
    call masterlist_addfld (fname='ISOPRENE', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='isoprene flux', hpindices=hpindices)

    !---monoterpene flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_2 == not_set) then
       ip_pvf_vocflx_2 = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_2)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_2)%val(pi)%p => p%pvf%vocflx_2
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_2 = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_2/)
    call masterlist_addfld (fname='MONOTERP', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='monoterpene flux', hpindices=hpindices)
    call masterlist_make_active (name='MONOTERP', tape_index=1) 

    !---other VOC flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_3 == not_set) then
       ip_pvf_vocflx_3 = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_3)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_3)%val(pi)%p => p%pvf%vocflx_3
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_3 = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_3/)
    call masterlist_addfld (fname='OVOC', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='other VOC flux', hpindices=hpindices)
    call masterlist_make_active (name='OVOC    ', tape_index=1) 

    !---other reactive VOC flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_4 == not_set) then
       ip_pvf_vocflx_4 = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_4)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_4)%val(pi)%p => p%pvf%vocflx_4
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_4 = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_4/)
    call masterlist_addfld (fname='ORVOC', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='other reactive VOC flux', hpindices=hpindices)
    call masterlist_make_active (name='ORVOC   ', tape_index=1) 

    !---biogenic CO flux
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pvf_vocflx_5 == not_set) then
       ip_pvf_vocflx_5 = pointer_index()
       allocate (clmptr_rs(ip_pvf_vocflx_5)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pvf_vocflx_5)%val(pi)%p => p%pvf%vocflx_5
          if (p%c%l%lps%lakpoi) p%pvf%vocflx_5 = 0.
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pvf_vocflx_5/)
    call masterlist_addfld (fname='BIOGENCO', type1d='pft', units='uG/M2/H', numlev=1, &
         avgflag='A', long_name='biogenic CO flux', hpindices=hpindices)
    call masterlist_make_active (name='BIOGENCO', tape_index=1) 

#endif

    ! Hydrology

    !---snow depth (liquid water)
    !---when defining pointers, set constant values (e.g. over lake) if neede
    if (ic_cws_h2osno == not_set) then
       ic_cws_h2osno = pointer_index()
       allocate (clmptr_rs(ic_cws_h2osno)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cws_h2osno)%val(ci)%p => c%cws%h2osno
       end do
    end if
    hpindices = (/-1, -1, ic_cws_h2osno, not_valid/)
    call masterlist_addfld (fname='H2OSNO', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='snow depth (liquid water)', hpindices=hpindices)
    call masterlist_make_active (name='H2OSNO  ', tape_index=1) 

    !---intercepted water
    !---when defining pointers, set constant values (e.g. over lake) if neede
    if (ic_cws_pws_a_h2ocan == not_set) then
       ic_cws_pws_a_h2ocan = pointer_index()
       allocate (clmptr_rs(ic_cws_pws_a_h2ocan)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cws_pws_a_h2ocan)%val(ci)%p => c%cws%pws_a%h2ocan
          if (c%l%lps%lakpoi) c%cws%pws_a%h2ocan = 0.
       end do
    end if
    if (ip_pws_h2ocan == not_set) then
       ip_pws_h2ocan = pointer_index()
       allocate (clmptr_rs(ip_pws_h2ocan)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pws_h2ocan)%val(pi)%p =>  p%pws%h2ocan
          if (p%c%l%lps%lakpoi) p%pws%h2ocan = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cws_pws_a_h2ocan, ip_pws_h2ocan/)
    call masterlist_addfld (fname='H2OCAN', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='intercepted water', hpindices=hpindices)
    call masterlist_make_active (name='H2OCAN  ', tape_index=1) 

    !---volumetric soil water
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cws_h2osoi_vol == not_set) then
       ic_cws_h2osoi_vol = pointer_index()
       allocate (clmptr_ra(ic_cws_h2osoi_vol)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_ra(ic_cws_h2osoi_vol)%val(ci)%p => c%cws%h2osoi_vol
       end do
    end if
    hpindices = (/-1, -1, ic_cws_h2osoi_vol, not_valid/)
    call masterlist_addfld (fname='H2OSOI', type1d='column', units='mm3/mm3', numlev=nlevsoi, &
         avgflag='A', long_name='volumetric soil water', hpindices=hpindices)
    call masterlist_make_active (name='H2OSOI  ', tape_index=1) 

    !---soil liquid water
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cws_h2osoi_liq == not_set) then
       ic_cws_h2osoi_liq = pointer_index()
       allocate (clmptr_ra(ic_cws_h2osoi_liq)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_ra(ic_cws_h2osoi_liq)%val(ci)%p => c%cws%h2osoi_liq
       end do
    end if
    hpindices = (/-1, -1, ic_cws_h2osoi_liq, not_valid/)
    call masterlist_addfld (fname='SOILLIQ', type1d='column', units='kg/m2', numlev=nlevsoi, &
         avgflag='A', long_name='soil liquid water', hpindices=hpindices)
    call masterlist_make_active (name='SOILLIQ ', tape_index=1) 

    !---soil liquid ice
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cws_h2osoi_ice == not_set) then
       ic_cws_h2osoi_ice = pointer_index()
       allocate (clmptr_ra(ic_cws_h2osoi_ice)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_ra(ic_cws_h2osoi_ice)%val(ci)%p => c%cws%h2osoi_ice
       end do
    end if
    hpindices = (/-1, -1, ic_cws_h2osoi_ice, not_valid/)
    call masterlist_addfld (fname='SOILICE', type1d='column', units='kg/m2', numlev=nlevsoi, &
         avgflag='A', long_name='soil liquid ice', hpindices=hpindices)
    call masterlist_make_active (name='SOILICE ', tape_index=1) 

    !---snow liquid water
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cws_snowliq == not_set) then
       ic_cws_snowliq = pointer_index()
       allocate (clmptr_rs(ic_cws_snowliq)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cws_snowliq)%val(ci)%p => c%cws%snowliq
       end do
    end if
    hpindices = (/-1, -1, ic_cws_snowliq, not_valid/)
    call masterlist_addfld (fname='SNOWLIQ', type1d='column', units='kg/m2', numlev=1, &
         avgflag='A', long_name='snow liquid water', hpindices=hpindices)
    call masterlist_make_active (name='SNOWLIQ ', tape_index=1) 

    !---snow liquid ice
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cws_snowice == not_set) then
       ic_cws_snowice = pointer_index()
       allocate (clmptr_rs(ic_cws_snowice)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cws_snowice)%val(ci)%p => c%cws%snowice
       end do
    end if
    hpindices = (/-1, -1, ic_cws_snowice, not_valid/)
    call masterlist_addfld (fname='SNOWICE', type1d='column', units='kg/m2', numlev=1, &
         avgflag='A', long_name='snow liquid ice', hpindices=hpindices)
    call masterlist_make_active (name='SNOWICE ', tape_index=1) 

    !---infiltration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_qflx_infl == not_set) then
       ic_cwf_qflx_infl = pointer_index()
       allocate (clmptr_rs(ic_cwf_qflx_infl)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_qflx_infl)%val(ci)%p => c%cwf%qflx_infl
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_qflx_infl, not_valid/)
    call masterlist_addfld (fname='QINFL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='infiltration', hpindices=hpindices)
    call masterlist_make_active (name='QINFL   ', tape_index=1) 

    !---surface runoff
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_qflx_surf == not_set) then
       ic_cwf_qflx_surf = pointer_index()
       allocate (clmptr_rs(ic_cwf_qflx_surf)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_qflx_surf)%val(ci)%p => c%cwf%qflx_surf
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_qflx_surf, not_valid/)
    call masterlist_addfld (fname='QOVER', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='surface runoff', hpindices=hpindices)
    call masterlist_make_active (name='QOVER   ', tape_index=1) 

    !---surface runoff at glaciers, wetlands, lakes
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_qflx_qrgwl == not_set) then
       ic_cwf_qflx_qrgwl = pointer_index()
       allocate (clmptr_rs(ic_cwf_qflx_qrgwl)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_qflx_qrgwl)%val(ci)%p => c%cwf%qflx_qrgwl
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_qflx_qrgwl, not_valid/)
    call masterlist_addfld (fname='QRGWL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='surface runoff at glaciers, wetlands, lakes', hpindices=hpindices)
    call masterlist_make_active (name='QRGWL   ', tape_index=1) 

    !---sub-surface drainage
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_qflx_drain == not_set) then
       ic_cwf_qflx_drain = pointer_index()
       allocate (clmptr_rs(ic_cwf_qflx_drain)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_qflx_drain)%val(ci)%p => c%cwf%qflx_drain
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_qflx_drain, not_valid/)
    call masterlist_addfld (fname='QDRAI', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='sub-surface drainage', hpindices=hpindices)
    call masterlist_make_active (name='QDRAI   ', tape_index=1) 

    !---interception
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_pwf_a_qflx_prec_intr == not_set) then
       ic_cwf_pwf_a_qflx_prec_intr = pointer_index()
       allocate (clmptr_rs(ic_cwf_pwf_a_qflx_prec_intr)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_pwf_a_qflx_prec_intr)%val(ci)%p => c%cwf%pwf_a%qflx_prec_intr
          if (c%l%lps%lakpoi) c%cwf%pwf_a%qflx_prec_intr = 0.
       end do
    end if
    if (ip_pwf_qflx_prec_intr == not_set) then
       ip_pwf_qflx_prec_intr = pointer_index()
       allocate (clmptr_rs(ip_pwf_qflx_prec_intr)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pwf_qflx_prec_intr)%val(pi)%p => p%pwf%qflx_prec_intr
          if (p%c%l%lps%lakpoi) p%pwf%qflx_prec_intr = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_prec_intr, ip_pwf_qflx_prec_intr/)
    call masterlist_addfld (fname='QINTR', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='interception', hpindices=hpindices)
    call masterlist_make_active (name='QINTR   ', tape_index=1) 

    !---throughfall
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_pwf_a_qflx_prec_grnd == not_set) then
       ic_cwf_pwf_a_qflx_prec_grnd = pointer_index()
       allocate (clmptr_rs(ic_cwf_pwf_a_qflx_prec_grnd)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_pwf_a_qflx_prec_grnd)%val(ci)%p => c%cwf%pwf_a%qflx_prec_grnd
       end do
    end if
    if (ip_pwf_qflx_prec_grnd == not_set) then
       ip_pwf_qflx_prec_grnd = pointer_index()
       allocate (clmptr_rs(ip_pwf_qflx_prec_grnd)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pwf_qflx_prec_grnd)%val(pi)%p => p%pwf%qflx_prec_grnd
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_prec_grnd, ip_pwf_qflx_prec_grnd/)
    call masterlist_addfld (fname='QDRIP', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='throughfall', hpindices=hpindices)
    call masterlist_make_active (name='QDRIP   ', tape_index=1) 

    !---snow melt
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_qflx_snomelt == not_set) then
       ic_cwf_qflx_snomelt = pointer_index()
       allocate (clmptr_rs(ic_cwf_qflx_snomelt)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_qflx_snomelt)%val(ci)%p => c%cwf%qflx_snomelt
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_qflx_snomelt, not_valid/)
    call masterlist_addfld (fname='QMELT', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='snow melt', hpindices=hpindices)
    call masterlist_make_active (name='QMELT   ', tape_index=1) 

    !---ground evaporation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_pwf_a_qflx_evap_soi == not_set) then
       ic_cwf_pwf_a_qflx_evap_soi = pointer_index()
       allocate (clmptr_rs(ic_cwf_pwf_a_qflx_evap_soi)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_pwf_a_qflx_evap_soi)%val(ci)%p => c%cwf%pwf_a%qflx_evap_soi
       end do
    end if
    if (ip_pwf_qflx_evap_soi == not_set) then
       ip_pwf_qflx_evap_soi = pointer_index()
       allocate (clmptr_rs(ip_pwf_qflx_evap_soi)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pwf_qflx_evap_soi)%val(pi)%p => p%pwf%qflx_evap_soi
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_evap_soi, ip_pwf_qflx_evap_soi/)
    call masterlist_addfld (fname='QSOIL', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='ground evaporation', hpindices=hpindices)
    call masterlist_make_active (name='QSOIL   ', tape_index=1) 

    !---canopy evaporation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_pwf_a_qflx_evap_can == not_set) then
       ic_cwf_pwf_a_qflx_evap_can = pointer_index()
       allocate (clmptr_rs(ic_cwf_pwf_a_qflx_evap_can)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_pwf_a_qflx_evap_can)%val(ci)%p => c%cwf%pwf_a%qflx_evap_can
          if (c%l%lps%lakpoi) c%cwf%pwf_a%qflx_evap_can = 0.
       end do
    end if
    if (ip_pwf_qflx_evap_can == not_set) then
       ip_pwf_qflx_evap_can = pointer_index()
       allocate (clmptr_rs(ip_pwf_qflx_evap_can)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pwf_qflx_evap_can)%val(pi)%p => p%pwf%qflx_evap_can
          if (p%c%l%lps%lakpoi) p%pwf%qflx_evap_can = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_evap_can, ip_pwf_qflx_evap_can/)
    call masterlist_addfld (fname='QVEGE', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='canopy evaporation', hpindices=hpindices)
    call masterlist_make_active (name='QVEGE   ', tape_index=1) 

    !---canopy transpiration
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwf_pwf_a_qflx_tran_veg == not_set) then
       ic_cwf_pwf_a_qflx_tran_veg = pointer_index()
       allocate (clmptr_rs(ic_cwf_pwf_a_qflx_tran_veg)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwf_pwf_a_qflx_tran_veg)%val(ci)%p => c%cwf%pwf_a%qflx_tran_veg
          if (c%l%lps%lakpoi) c%cwf%pwf_a%qflx_tran_veg = 0.
       end do
    end if
    if (ip_pwf_qflx_tran_veg == not_set) then
       ip_pwf_qflx_tran_veg = pointer_index()
       allocate (clmptr_rs(ip_pwf_qflx_tran_veg)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pwf_qflx_tran_veg)%val(pi)%p => p%pwf%qflx_tran_veg
          if (p%c%l%lps%lakpoi) p%pwf%qflx_tran_veg = 0.
       end do
    end if
    hpindices = (/-1, -1, ic_cwf_pwf_a_qflx_tran_veg, ip_pwf_qflx_tran_veg/)
    call masterlist_addfld (fname='QVEGT', type1d='column', units='mm/s', numlev=1, &
         avgflag='A', long_name='canopy transpiration', hpindices=hpindices)
    call masterlist_make_active (name='QVEGT   ', tape_index=1) 

#if (defined RTM)
    ! RTM River Routing

    !---RTM river flow (maximum subgrid flow)
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_gwf_qchan2 == not_set) then
       ig_gwf_qchan2 = pointer_index()
       allocate (clmptr_rs(ig_gwf_qchan2)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_gwf_qchan2)%val(gi)%p => g%gwf%qchan2
       end do
    endif
    hpindices = (/ig_gwf_qchan2, -1, -1, -1/)
    call masterlist_addfld (fname='QCHANR', type1d='gridcell', units='m3/s', numlev=1, &
         avgflag='A', long_name='RTM river flow (maximum subgrid flow)', hpindices=hpindices)
    call masterlist_make_active (name='QCHANR  ', tape_index=1) 

    !---RTM river discharge into ocean
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_gwf_qchocn2 == not_set) then
       ig_gwf_qchocn2 = pointer_index()
       allocate (clmptr_rs(ig_gwf_qchocn2)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_gwf_qchocn2)%val(gi)%p=> g%gwf%qchocn2
       end do
    endif
    hpindices = (/ig_gwf_qchocn2, -1, -1,-1/)
    call masterlist_addfld (fname='QCHOCNR', type1d='gridcell', units='m3/s', numlev=1, &
         avgflag='A', long_name='RTM river discharge into ocean', hpindices=hpindices)
    call masterlist_make_active (name='QCHOCNR ', tape_index=1) 
#endif

    ! Water and energy balance checks

    !---soil/lake energy conservation error
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cebal_errsoi == not_set) then
       ic_cebal_errsoi = pointer_index()
       allocate (clmptr_rs(ic_cebal_errsoi)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cebal_errsoi)%val(ci)%p => c%cebal%errsoi
       end do
    end if
    hpindices = (/-1, -1, ic_cebal_errsoi, not_valid/)
    call masterlist_addfld (fname='ERRSOI', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='soil/lake energy conservation error', hpindices=hpindices)
    call masterlist_make_active (name='ERRSOI  ', tape_index=1) 

    if (ic_cebal_errseb == not_set) then
       ic_cebal_errseb = pointer_index()
       allocate (clmptr_rs(ic_cebal_errseb)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cebal_errseb)%val(ci)%p => c%cebal%errseb
       end do
    end if
    hpindices = (/-1, -1, ic_cebal_errseb, not_valid/)
    call masterlist_addfld (fname='ERRSEB', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='surface energy conservation error', hpindices=hpindices)
    call masterlist_make_active (name='ERRSEB  ', tape_index=1) 

    !---solar radiation conservation error
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cebal_errsol == not_set) then
       ic_cebal_errsol = pointer_index()
       allocate (clmptr_rs(ic_cebal_errsol)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cebal_errsol)%val(ci)%p => c%cebal%errsol
       end do
    end if
    hpindices = (/-1, -1, ic_cebal_errsol, not_valid/)
    call masterlist_addfld (fname='ERRSOL', type1d='column', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='solar radiation conservation error', hpindices=hpindices)
    call masterlist_make_active (name='ERRSOL  ', tape_index=1) 

    !---total water conservation error
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ic_cwbal_errh2o == not_set) then
       ic_cwbal_errh2o = pointer_index()
       allocate (clmptr_rs(ic_cwbal_errh2o)%val(begc:endc)) 
       do ci = begc,endc
          c => cpoint(ci)%c
          clmptr_rs(ic_cwbal_errh2o)%val(ci)%p => c%cwbal%errh2o
       end do
    end if
    hpindices = (/-1, -1, ic_cwbal_errh2o, not_valid/)
    call masterlist_addfld (fname='ERRH2O', type1d='column', units='mm', numlev=1, &
         avgflag='A', long_name='total water conservation error', hpindices=hpindices)
    call masterlist_make_active (name='ERRH2O  ', tape_index=1) 

    ! Atmospheric forcing                                

    !---atmospheric rain
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2lf_forc_rain == not_set) then
       ig_a2lf_forc_rain = pointer_index()
       allocate (clmptr_rs(ig_a2lf_forc_rain)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2lf_forc_rain)%val(gi)%p => g%a2lf%forc_rain
       end do
    end if
    hpindices = (/ig_a2lf_forc_rain, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='RAIN', type1d='gridcell', units='mm/s', numlev=1, &
         avgflag='A', long_name='atmospheric rain', hpindices=hpindices)
    call masterlist_make_active (name='RAIN    ', tape_index=1) 

    !---atmospheric snow
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2lf_forc_snow == not_set) then
       ig_a2lf_forc_snow = pointer_index()
       allocate (clmptr_rs(ig_a2lf_forc_snow)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2lf_forc_snow)%val(gi)%p => g%a2lf%forc_snow
       end do
    end if
    hpindices = (/ig_a2lf_forc_snow, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='SNOW', type1d='gridcell', units='mm/s', numlev=1, &
         avgflag='A', long_name='atmospheric snow', hpindices=hpindices)
    call masterlist_make_active (name='SNOW    ', tape_index=1) 

    !---atmospheric air temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2ls_forc_t == not_set) then
       ig_a2ls_forc_t = pointer_index()
       allocate (clmptr_rs(ig_a2ls_forc_t)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2ls_forc_t)%val(gi)%p => g%a2ls%forc_t
       end do
    end if
    hpindices = (/ig_a2ls_forc_t, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='TBOT', type1d='gridcell', units='K', numlev=1, &
         avgflag='A', long_name='atmospheric air temperature', hpindices=hpindices)
    call masterlist_make_active (name='TBOT    ', tape_index=1) 

    !---atmospheric air potential temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2ls_forc_th == not_set) then
       ig_a2ls_forc_th = pointer_index()
       allocate (clmptr_rs(ig_a2ls_forc_th)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2ls_forc_th)%val(gi)%p => g%a2ls%forc_th
       end do
    end if
    hpindices = (/ig_a2ls_forc_th, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='THBOT', type1d='gridcell', units='K', numlev=1, &
         avgflag='A', long_name='atmospheric air potential temperature', hpindices=hpindices)
    call masterlist_make_active (name='THBOT   ', tape_index=1) 

    !---atmospheric wind velocity magnitude
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2ls_forc_wind == not_set) then
       ig_a2ls_forc_wind = pointer_index()
       allocate (clmptr_rs(ig_a2ls_forc_wind)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2ls_forc_wind)%val(gi)%p => g%a2ls%forc_wind
       end do
    end if
    hpindices = (/ig_a2ls_forc_wind, not_valid, not_valid, not_valid/)
    call masterlist_addfld (fname='WIND', type1d='gridcell', units='m/s', numlev=1, &
         avgflag='A', long_name='atmospheric wind velocity magnitude', hpindices=hpindices)
    call masterlist_make_active (name='WIND    ', tape_index=1) 

    !---atmospheric specific humidity
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2ls_forc_q == not_set) then
       ig_a2ls_forc_q = pointer_index()
       allocate (clmptr_rs(ig_a2ls_forc_q)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2ls_forc_q)%val(gi)%p => g%a2ls%forc_q
       end do
    end if
    hpindices = (/ig_a2ls_forc_q, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='QBOT', type1d='gridcell', units='kg/kg', numlev=1, &
         avgflag='A', long_name='atmospheric specific humidity', hpindices=hpindices)
    call masterlist_make_active (name='QBOT    ', tape_index=1) 

    !---atmospheric reference height
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2ls_forc_hgt == not_set) then
       ig_a2ls_forc_hgt = pointer_index()
       allocate (clmptr_rs(ig_a2ls_forc_hgt)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2ls_forc_hgt)%val(gi)%p => g%a2ls%forc_hgt
       end do
    end if
    hpindices = (/ig_a2ls_forc_hgt, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='ZBOT', type1d='gridcell', units='m', numlev=1, &
         avgflag='A', long_name='atmospheric reference height', hpindices=hpindices)
    call masterlist_make_active (name='ZBOT    ', tape_index=1) 

    !---atmospheric longwave radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2lf_forc_lwrad == not_set) then
       ig_a2lf_forc_lwrad = pointer_index()
       allocate (clmptr_rs(ig_a2lf_forc_lwrad)%val(begg:endg)) 
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2lf_forc_lwrad)%val(gi)%p => g%a2lf%forc_lwrad
       end do
    end if
    hpindices = (/ig_a2lf_forc_lwrad, not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='FLDS', type1d='gridcell', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='atmospheric longwave radiation', hpindices=hpindices)
    call masterlist_make_active (name='FLDS    ', tape_index=1) 

    !---atmospheric incident solar radiation
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ig_a2lf_forc_solar == not_set) then
       ig_a2lf_forc_solar = pointer_index()
       allocate (clmptr_rs(ig_a2lf_forc_solar)%val(begg:endg))
       do gi = begg,endg
          g => gpoint(gi)%g
          clmptr_rs(ig_a2lf_forc_solar)%val(gi)%p => g%a2lf%forc_solar
       end do
    end if
    hpindices = (/ig_a2lf_forc_solar,  not_valid, not_valid, not_valid/) 
    call masterlist_addfld (fname='FSDS', type1d='gridcell', units='watt/m^2', numlev=1, &
         avgflag='A', long_name='atmospheric incident solar radiation', hpindices=hpindices)
    call masterlist_make_active (name='FSDS    ', tape_index=1) 

#if (defined OFFLINE) && (defined DGVM)
    ! History output of accumulation variables                             

    !---height of top of canopy
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pps_htop == not_set) then
       ip_pps_htop = pointer_index()
       allocate (clmptr_rs(ip_pps_htop)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_htop)%val(pi)%p => p%pps%htop
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pps_htop/)
    call masterlist_addfld (fname='HTOP', type1d='pft', units='m', numlev=1, &
         avgflag='A', long_name='height of top of canopy', hpindices=hpindices)
    call masterlist_make_active (name='HTOP    ', tape_index=1) 

    !---height of bottom of canopy
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pps_hbot == not_set) then
       ip_pps_hbot = pointer_index()
       allocate (clmptr_rs(ip_pps_hbot)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_hbot)%val(pi)%p => p%pps%hbot
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pps_hbot/)
    call masterlist_addfld (fname='HBOT', type1d='pft', units='m', numlev=1, &
         avgflag='A', long_name='height of bottom of canopy', hpindices=hpindices)
    call masterlist_make_active (name='HBOT    ', tape_index=1) 

    !---total one-sided leaf area index
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pps_tlai == not_set) then
       ip_pps_tlai = pointer_index()
       allocate (clmptr_rs(ip_pps_tlai)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_tlai)%val(pi)%p => p%pps%tlai
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pps_tlai/)
    call masterlist_addfld (fname='TLAI', type1d='pft', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='total one-sided leaf area index', hpindices=hpindices)
    call masterlist_make_active (name='TLAI    ', tape_index=1) 

    !---total one-sided stem area index
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pps_tsai == not_set) then
       ip_pps_tsai = pointer_index()
       allocate (clmptr_rs(ip_pps_tsai)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pps_tsai)%val(pi)%p => p%pps%tsai
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pps_tsai/)
    call masterlist_addfld (fname='TSAI', type1d='pft', units='m^2/m^2', &
         numlev=1, avgflag='A', long_name='total one-sided stem area index', hpindices=hpindices)
    call masterlist_make_active (name='TSAI    ', tape_index=1) 

    !---daily average 2-m temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdgvs_t_mo == not_set) then
       ip_pdgvs_t_mo = pointer_index()
       allocate (clmptr_rs(ip_pdgvs_t_mo)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdgvs_t_mo)%val(pi)%p => p%pdgvs%t_mo
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdgvs_t_mo/)
    call masterlist_addfld (fname='TDA', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='daily average 2-m temperature', hpindices=hpindices)
    call masterlist_make_active (name='TDA     ', tape_index=1) 

    !---10-day running mean of 2-m temperature
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdgvs_t10 == not_set) then
       ip_pdgvs_t10 = pointer_index()
       allocate (clmptr_rs(ip_pdgvs_t10)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdgvs_t10)%val(pi)%p => p%pdgvs%t10
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdgvs_t10/)
    call masterlist_addfld (fname='T10', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='10-day running mean of 2-m temperature', hpindices=hpindices)
    call masterlist_make_active (name='T10     ', tape_index=1) 

    !---growing degree-days base 0C
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdgvs_agdd0 == not_set) then
       ip_pdgvs_agdd0 = pointer_index()
       allocate (clmptr_rs(ip_pdgvs_agdd0)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdgvs_agdd0)%val(pi)%p => p%pdgvs%agdd0
       end do
    endif
    hpindices = (/-1, -1, -1, ip_pdgvs_agdd0/)
    call masterlist_addfld (fname='AGDD0', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='growing degree-days base 0C', hpindices=hpindices)
    call masterlist_make_active (name='AGDD0   ', tape_index=1) 

    !---growing degree-days base 5C
    !---when defining pointers, set constant values (e.g. over lake) if needed
    if (ip_pdgvs_agdd5 == not_set) then
       ip_pdgvs_agdd5 = pointer_index()
       allocate (clmptr_rs(ip_pdgvs_agdd5)%val(begp:endp)) 
       do pi = begp,endp
          p => ppoint(pi)%p
          clmptr_rs(ip_pdgvs_agdd5)%val(pi)%p => p%pdgvs%agdd5
       end do
    end if
    hpindices = (/-1, -1, -1, ip_pdgvs_agdd5/)
    call masterlist_addfld (fname='AGDD5', type1d='pft', units='K', numlev=1, &
         avgflag='A', long_name='growing degree-days base 5C', hpindices=hpindices)
    call masterlist_make_active (name='AGDD5   ', tape_index=1) 
#endif

    ! Print masterlist of history fields

    call masterlist_printflds ()

    ! Initialize active history fields if not a restart run
    ! If a restart run, then this information is obtained from the restart history file

    if (nsrest == 0 .or. nsrest == 3) call htapes_build ()

  end subroutine initHistFlds

end module histFldsMod
