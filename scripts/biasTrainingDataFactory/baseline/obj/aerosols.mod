  -J     k820309    [          18.0        ta                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/atm/cam/src/physics/cam1/aerosols.F90 AEROSOLS              AEROSOL_INITIALIZE GET_AEROSOL AEROSOL_DIAGNOSTICS AEROSOL_INDIRECT GET_RF_SCALES GET_INT_SCALES AERINT NAER_ALL IDXSUL IDXSSLT IDXOCPHO IDXBCPHO IDXOCPHI IDXBCPHI IDXBG IDXDUSTFIRST NUMDUST IDXCARBONFIRST NUMCARBON AEROSOL_NAME RADFORCE SULSCL_RF CARSCL_RF SSLTSCL_RF DUSTSCL_RF BGSCL_RF TAUBACK SULSCL CARSCL SSLTSCL DUSTSCL          @                                         
       PCOLS PVER PVERP BEGCHUNK ENDCHUNK                                                     
       PLON PLAT PLEV PLEVP MASTERPROC                      @                              
       GET_NCOLS_P SCATTER_FIELD_TO_CHUNK                      @                              
       GET_CURR_CALDAY                                                     
       INF BIGINT                   @                              
       R8 SHR_KIND_R8 %         @                                                   
       #OFFSET              
                                                                                         	                                                      11                                             
                                                      1                                                                                                   2                                                                                                   7                                                                                                   8                                                                                    	               9                                                                                    
               10                                                                                                   11                                                                                                   3                                                                                                   4                                                                                                   7                                                                                                   4                                                                                                                   TWp          n                  	                       CMSUL_V    n                     	                       CMSSLT_V   n                     	                       CMDUST1_V  n                     	                       CMDUST2_V  n                     	                       CMDUST3_V  n                     	                       CMDUST4_V  n                     	                       CMOCPHO_V  n                     	                       CMBCPHO_V  n                     	                       CMOCPHI_V  n                     	                       CMBCPHI_V  n                     	                       CBackgrnd  h  p          p          p            p                                                                                                                                                                                                                                                                                                                                              
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                       
       #         @                                   !                    #MPISHORTHAND!AEROSOL_INITIALIZE%MPIFCMB5 "   #MPISHORTHAND!AEROSOL_INITIALIZE%MPIFCMB9 $   #MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIV1 &   #MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIV2 *   #MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIVC -                                           "                          #AEROSOL_INITIALIZE%MPIFCMB5%MPI_UNWEIGHTED #                                          #                                                            $                          #AEROSOL_INITIALIZE%MPIFCMB9%MPI_WEIGHTS_EMPTY %                                           %                                                            &                          #AEROSOL_INITIALIZE%MPIPRIV1%MPI_BOTTOM '   #AEROSOL_INITIALIZE%MPIPRIV1%MPI_IN_PLACE (   #AEROSOL_INITIALIZE%MPIPRIV1%MPI_STATUS_IGNORE )                                          '                                                           (                                                           )                                p          p            p                                                                          *                          #AEROSOL_INITIALIZE%MPIPRIV2%MPI_STATUSES_IGNORE +   #AEROSOL_INITIALIZE%MPIPRIV2%MPI_ERRCODES_IGNORE ,                                           +                                 p          p          p            p          p                                                                          ,                                p          p            p                                                                          -                          #AEROSOL_INITIALIZE%MPIPRIVC%MPI_ARGVS_NULL .   #AEROSOL_INITIALIZE%MPIPRIVC%MPI_ARGV_NULL /   -                                        .                                 p          p          p            p          p                                  -                                        /                                p          p            p                                  #         @                                   0                    #C 1   #NCOL1 2   #NCOL 3   #PINT 4   #AEROSOLT 5   #SCALE 6             
  @                               1                     
  @                               2                     
  @                               3                     
  @                              4     ø              
    p 	         p          p            p          p                                    D @                              5     P
             
     p ù         p          p          p            p          p          p                                    
  @                              6                   
    p          p            p                          #         @                                   7                   #AEROSOL_DIAGNOSTICS%PHYSICS_STATE 8   #STATE Q   #AEROSOL_MMR R                                                  @                           8     'Ð                   #LCHNK 9   #NCOL :   #PS ;   #PHIS <   #T =   #U >   #V ?   #S @   #OMEGA A   #PMID B   #PDEL C   #RPDEL D   #LNPMID E   #EXNER F   #ZM G   #Q H   #PINT I   #LNPINT J   #ZI K   #TE_INI L   #TE_CUR M   #TW_INI N   #TW_CUR O   #COUNT P                 $                              9                                 $                              :                                $                             ;                             
  p          p            p                                        $                             <            H                 
  p          p            p                                        $                             =     ð                        
  p 	         p          p            p          p                                        $                             >     ð                       
  p 	         p          p            p          p                                        $                             ?     ð                       
  p 	         p          p            p          p                                        $                             @     ð                       
  p 	         p          p            p          p                                        $                             A     ð                    	   
  p 	         p          p            p          p                                        $                             B     ð       &             
   
  p 	         p          p            p          p                                        $                             C     ð       -                
  p 	         p          p            p          p                                        $                             D     ð       5                
  p 	         p          p            p          p                                        $                             E     ð       <                
  p 	         p          p            p          p                                        $                             F     ð       D                
  p 	         p          p            p          p                                        $                             G     ð       K                
  p 	         p          p            p          p                                        $                             H     Ð      S                
  p ù         p          p          p            p          p          p                                        $                             I     ø       i                
  p 	         p          p            p          p                                        $                             J     ø       Hq                
  p 	         p          p            p          p                                        $                             K     ø       y                
  p 	         p          p            p          p                                        $                             L            È                
  p          p            p                                        $                             M                            
  p          p            p                                        $                             N            H                
  p          p            p                                        $                             O                            
  p          p            p                                        $                              P     È                      
                                  Q     Ð             #AEROSOL_DIAGNOSTICS%PHYSICS_STATE 8             
                                 R     P
             
     p ù         p          p          p            p          p          p                          #         @                                   S                
   #AEROSOLS!AEROSOL_INDIRECT!COMCTL T   #AEROSOLS!AEROSOL_INDIRECT!COMCTL_R8 w   #NCOL1 {   #NCOL |   #LCHNK }   #LANDFRAC ~   #PMID    #T    #QM1    #CLD    #ZM    #REL                                                @                          T            "              #AEROSOL_INDIRECT%COMCTL%ITSST U   #AEROSOL_INDIRECT%COMCTL%NSREST V   #AEROSOL_INDIRECT%COMCTL%IRADSW W   #AEROSOL_INDIRECT%COMCTL%IRADLW X   #AEROSOL_INDIRECT%COMCTL%IRADAE Y   #AEROSOL_INDIRECT%COMCTL%NREFRQ Z   #AEROSOL_INDIRECT%COMCTL%ANNCYC [   #AEROSOL_INDIRECT%COMCTL%NLEND \   #AEROSOL_INDIRECT%COMCTL%NLRES ]   #AEROSOL_INDIRECT%COMCTL%NLHST ^   #AEROSOL_INDIRECT%COMCTL%LBRNCH _   #AEROSOL_INDIRECT%COMCTL%AERES `   #AEROSOL_INDIRECT%COMCTL%OZNCYC a   #AEROSOL_INDIRECT%COMCTL%SSTCYC b   #AEROSOL_INDIRECT%COMCTL%ICECYC c   #AEROSOL_INDIRECT%COMCTL%ADIABATIC d   #AEROSOL_INDIRECT%COMCTL%FLXAVE e   #AEROSOL_INDIRECT%COMCTL%IDEAL_PHYS f   #AEROSOL_INDIRECT%COMCTL%NSPLIT g   #AEROSOL_INDIRECT%COMCTL%IORD h   #AEROSOL_INDIRECT%COMCTL%JORD i   #AEROSOL_INDIRECT%COMCTL%KORD j   #AEROSOL_INDIRECT%COMCTL%USE_ETA k   #AEROSOL_INDIRECT%COMCTL%AQUA_PLANET l   #AEROSOL_INDIRECT%COMCTL%DORAMP_SO4 m   #AEROSOL_INDIRECT%COMCTL%DORAMP_SCON n   #AEROSOL_INDIRECT%COMCTL%FULLGRID o   #AEROSOL_INDIRECT%COMCTL%PRINT_STEP_COST p   #AEROSOL_INDIRECT%COMCTL%DOABSEMS q   #AEROSOL_INDIRECT%COMCTL%DOSW r   #AEROSOL_INDIRECT%COMCTL%DOLW s   #AEROSOL_INDIRECT%COMCTL%INDIRECT t   #AEROSOL_INDIRECT%COMCTL%SOM_CONSCHK_FRQ u   #AEROSOL_INDIRECT%COMCTL%ICE_CONSCHK_FRQ v                @                          U                                 @                          V                                @                          W                                @                          X                                @                          Y                                @                          Z                                @                          [                                @                          \                                @                          ]                                 @                          ^     $                           @                          _     (                           @                          `     ,                           @                          a     0                           @                          b     4                           @                          c     8                           @                          d     <                           @                          e     @                           @                          f     D                           @                          g     H                           @                          h     L                           @                          i     P                           @                          j     T                           @                          k     X                           @                          l     \                           @                          m     `                           @                          n     d                           @                          o     h                           @                          p     l                           @                          q     p                           @                          r     t                           @                          s     x                           @                          t     |                           @                          u                                @                          v                                  @                          w                          #AEROSOL_INDIRECT%COMCTL_R8%DIVDAMPN x   #AEROSOL_INDIRECT%COMCTL_R8%PRECC_THRESH y   #AEROSOL_INDIRECT%COMCTL_R8%PRECL_THRESH z                @                         x             
                    @                         y            
                    @                         z            
                 
                                  {                     
                                  |                     
  @                               }                     
                                 ~                   
 &   p          p            p                                    
                                      ð              
 '   p 	         p          p            p          p                                    
                                      ð              
 (   p 	         p          p            p          p                                    
                                      Ð             
 )   p ù         p          p          p            p          p          p                                    
                                      ð              
 *   p 	         p          p            p          p                                    
                                      ð              
 +   p 	         p          p            p          p                                    
                                      ð              
 ,   p 	         p          p            p          p                          #         @                                                       #SCALES              D                                                   
 D    p          p            p                          #         @                                                       #SCALES              D                                                   
 E    p          p            p                          #         @                                                               n      fn#fn      W  b   uapp(AEROSOLS    e  c   J  PPGRID    È  `   J  PMGRID    (  c   J  PHYS_GRID      P   J  TIME_MANAGER    Û  K   J  INFNAN    &  O   J  SHR_KIND_MOD -   u  \       GET_CURR_CALDAY+TIME_MANAGER 4   Ñ  @   a   GET_CURR_CALDAY%OFFSET+TIME_MANAGER      r       NAER_ALL      q       IDXSUL    ô  q       IDXSSLT    e  q       IDXOCPHO    Ö  q       IDXBCPHO    G  q       IDXOCPHI    ¸  r       IDXBCPHI    *  r       IDXBG      q       IDXDUSTFIRST    	  q       NUMDUST    ~	  q       IDXCARBONFIRST    ï	  q       NUMCARBON    `
        AEROSOL_NAME    ß  @       RADFORCE      @       SULSCL_RF    _  @       CARSCL_RF      @       SSLTSCL_RF    ß  @       DUSTSCL_RF      @       BGSCL_RF    _  @       TAUBACK      @       SULSCL    ß  @       CARSCL      @       SSLTSCL    _  @       DUSTSCL #     .      AEROSOL_INITIALIZE O   Í       MPISHORTHAND!AEROSOL_INITIALIZE%MPIFCMB5+MPISHORTHAND=MPIFCMB5 H   M  H     AEROSOL_INITIALIZE%MPIFCMB5%MPI_UNWEIGHTED+MPISHORTHAND O          MPISHORTHAND!AEROSOL_INITIALIZE%MPIFCMB9+MPISHORTHAND=MPIFCMB9 K     H     AEROSOL_INITIALIZE%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPISHORTHAND O   `  Ý     MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIV1+MPISHORTHAND=MPIPRIV1 D   =  H     AEROSOL_INITIALIZE%MPIPRIV1%MPI_BOTTOM+MPISHORTHAND F     H     AEROSOL_INITIALIZE%MPIPRIV1%MPI_IN_PLACE+MPISHORTHAND K   Í  ¤     AEROSOL_INITIALIZE%MPIPRIV1%MPI_STATUS_IGNORE+MPISHORTHAND O   q  º     MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIV2+MPISHORTHAND=MPIPRIV2 M   +  Ä     AEROSOL_INITIALIZE%MPIPRIV2%MPI_STATUSES_IGNORE+MPISHORTHAND M   ï  ¤     AEROSOL_INITIALIZE%MPIPRIV2%MPI_ERRCODES_IGNORE+MPISHORTHAND O     ¯     MPISHORTHAND!AEROSOL_INITIALIZE%MPIPRIVC+MPISHORTHAND=MPIPRIVC H   B  Ä     AEROSOL_INITIALIZE%MPIPRIVC%MPI_ARGVS_NULL+MPISHORTHAND G     ¤     AEROSOL_INITIALIZE%MPIPRIVC%MPI_ARGV_NULL+MPISHORTHAND    ª         GET_AEROSOL    1  @   a   GET_AEROSOL%C "   q  @   a   GET_AEROSOL%NCOL1 !   ±  @   a   GET_AEROSOL%NCOL !   ñ  ´   a   GET_AEROSOL%PINT %   ¥  Ô   a   GET_AEROSOL%AEROSOLT "   y     a   GET_AEROSOL%SCALE $     ¨       AEROSOL_DIAGNOSTICS @   µ  <     AEROSOL_DIAGNOSTICS%PHYSICS_STATE+PHYSICS_TYPES F   ñ  H   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%LCHNK+PHYSICS_TYPES E   9   H   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%NCOL+PHYSICS_TYPES C         a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%PS+PHYSICS_TYPES E   !     a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%PHIS+PHYSICS_TYPES B   ¹!  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%T+PHYSICS_TYPES B   u"  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%U+PHYSICS_TYPES B   1#  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%V+PHYSICS_TYPES B   í#  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%S+PHYSICS_TYPES F   ©$  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%OMEGA+PHYSICS_TYPES E   e%  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%PMID+PHYSICS_TYPES E   !&  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%PDEL+PHYSICS_TYPES F   Ý&  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%RPDEL+PHYSICS_TYPES G   '  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%LNPMID+PHYSICS_TYPES F   U(  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%EXNER+PHYSICS_TYPES C   )  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%ZM+PHYSICS_TYPES B   Í)  Ü   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%Q+PHYSICS_TYPES E   ©*  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%PINT+PHYSICS_TYPES G   e+  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%LNPINT+PHYSICS_TYPES C   !,  ¼   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%ZI+PHYSICS_TYPES G   Ý,     a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%TE_INI+PHYSICS_TYPES G   y-     a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%TE_CUR+PHYSICS_TYPES G   .     a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%TW_INI+PHYSICS_TYPES G   ±.     a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%TW_CUR+PHYSICS_TYPES F   M/  H   a   AEROSOL_DIAGNOSTICS%PHYSICS_STATE%COUNT+PHYSICS_TYPES *   /  o   a   AEROSOL_DIAGNOSTICS%STATE 0   0  Ô   a   AEROSOL_DIAGNOSTICS%AEROSOL_MMR !   Ø0        AEROSOL_INDIRECT 1   î1  @     AEROSOLS!AEROSOL_INDIRECT!COMCTL .   .7  H      AEROSOL_INDIRECT%COMCTL%ITSST /   v7  H      AEROSOL_INDIRECT%COMCTL%NSREST /   ¾7  H      AEROSOL_INDIRECT%COMCTL%IRADSW /   8  H      AEROSOL_INDIRECT%COMCTL%IRADLW /   N8  H      AEROSOL_INDIRECT%COMCTL%IRADAE /   8  H      AEROSOL_INDIRECT%COMCTL%NREFRQ /   Þ8  H      AEROSOL_INDIRECT%COMCTL%ANNCYC .   &9  H      AEROSOL_INDIRECT%COMCTL%NLEND .   n9  H      AEROSOL_INDIRECT%COMCTL%NLRES .   ¶9  H      AEROSOL_INDIRECT%COMCTL%NLHST /   þ9  H      AEROSOL_INDIRECT%COMCTL%LBRNCH .   F:  H      AEROSOL_INDIRECT%COMCTL%AERES /   :  H      AEROSOL_INDIRECT%COMCTL%OZNCYC /   Ö:  H      AEROSOL_INDIRECT%COMCTL%SSTCYC /   ;  H      AEROSOL_INDIRECT%COMCTL%ICECYC 2   f;  H      AEROSOL_INDIRECT%COMCTL%ADIABATIC /   ®;  H      AEROSOL_INDIRECT%COMCTL%FLXAVE 3   ö;  H      AEROSOL_INDIRECT%COMCTL%IDEAL_PHYS /   ><  H      AEROSOL_INDIRECT%COMCTL%NSPLIT -   <  H      AEROSOL_INDIRECT%COMCTL%IORD -   Î<  H      AEROSOL_INDIRECT%COMCTL%JORD -   =  H      AEROSOL_INDIRECT%COMCTL%KORD 0   ^=  H      AEROSOL_INDIRECT%COMCTL%USE_ETA 4   ¦=  H      AEROSOL_INDIRECT%COMCTL%AQUA_PLANET 3   î=  H      AEROSOL_INDIRECT%COMCTL%DORAMP_SO4 4   6>  H      AEROSOL_INDIRECT%COMCTL%DORAMP_SCON 1   ~>  H      AEROSOL_INDIRECT%COMCTL%FULLGRID 8   Æ>  H      AEROSOL_INDIRECT%COMCTL%PRINT_STEP_COST 1   ?  H      AEROSOL_INDIRECT%COMCTL%DOABSEMS -   V?  H      AEROSOL_INDIRECT%COMCTL%DOSW -   ?  H      AEROSOL_INDIRECT%COMCTL%DOLW 1   æ?  H      AEROSOL_INDIRECT%COMCTL%INDIRECT 8   .@  H      AEROSOL_INDIRECT%COMCTL%SOM_CONSCHK_FRQ 8   v@  H      AEROSOL_INDIRECT%COMCTL%ICE_CONSCHK_FRQ 4   ¾@  Ó      AEROSOLS!AEROSOL_INDIRECT!COMCTL_R8 4   A  H      AEROSOL_INDIRECT%COMCTL_R8%DIVDAMPN 8   ÙA  H      AEROSOL_INDIRECT%COMCTL_R8%PRECC_THRESH 8   !B  H      AEROSOL_INDIRECT%COMCTL_R8%PRECL_THRESH '   iB  @   a   AEROSOL_INDIRECT%NCOL1 &   ©B  @   a   AEROSOL_INDIRECT%NCOL '   éB  @   a   AEROSOL_INDIRECT%LCHNK *   )C     a   AEROSOL_INDIRECT%LANDFRAC &   ½C  ´   a   AEROSOL_INDIRECT%PMID #   qD  ´   a   AEROSOL_INDIRECT%T %   %E  Ô   a   AEROSOL_INDIRECT%QM1 %   ùE  ´   a   AEROSOL_INDIRECT%CLD $   ­F  ´   a   AEROSOL_INDIRECT%ZM %   aG  ´   a   AEROSOL_INDIRECT%REL    H  T       GET_RF_SCALES %   iH     a   GET_RF_SCALES%SCALES    ýH  T       GET_INT_SCALES &   QI     a   GET_INT_SCALES%SCALES    åI  H       AERINT 