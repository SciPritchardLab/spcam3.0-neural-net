#include <params.h>
#ifdef CRM 
subroutine crmic_debug_dump (idstr)
      use buffer
      use phys_grid
      use pmgrid, only: iam

      implicit none
      character(*), intent(in) :: idstr
      integer icol,icrmx,icrmy,icrmz,iz
      character(2) numb
      icol = 2
      icrmx = 12
      icrmy = 1
      icrmz = 13
      iz = 20
!      write (6,*) trim(idstr),'RANK ',rank,': q_crm(icol=',icol,',icrmx=',icrmx,',icrmy=',icrmy,',icrmz=',icrmz,',lchnk=begchunk=',begchunk,')=',q_crm(icol,icrmx,icrmy,icrmz,begchunk)
      if (iam .eq. 3 ) then
        write (numb,'(i2)') iam
        write (6,*) trim (idstr), 'rank',numb,': q_crm=' , q_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': u_crm=' , u_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': v_crm=' , v_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': w_crm=' , w_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': t_crm=' , t_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qp_crm=' , qp_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qn_crm=' , qn_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qrs_crm=' , qrs_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qrl_crm=' , qrl_crm(icol,icrmx,icrmy,icrmz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': fsds_crm=', fsds_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': fsns_crm=', fsns_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': fsdsc_crm=', fsdsc_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': fsntoa_crm=', fsntoa_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': fsntoac_crm=', fsntoac_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': flwds_crm=', flwds_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': flns_crm=', flns_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': flnsc_crm=', flnsc_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': flut_crm=', flut_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': flutc_crm=', flutc_crm(icol,icrmx,icrmy,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qrs1=', qrs1(icol,iz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': qrl1=', qrl1(icol,iz,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': rad_buffer(13)=', rad_buffer (icol,13,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': rad_buffer(1)=', rad_buffer (icol,1,begchunk)
        write (6,*) trim (idstr), 'rank',numb,': rad_buffer(7)=', rad_buffer (icol,7,begchunk)
      
      end if
       end subroutine crmic_debug_dump
#endif
