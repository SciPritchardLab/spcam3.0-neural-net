  �  G   k820309    [          18.0        �ta                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/lnd/clm2/src/main/controlMod.F90 CONTROLMOD                                                     
                                                           
                            @                              
       MAX_TAPES MAX_NAMLEN HIST_EMPTY_HTAPES HIST_DOV2XY HIST_AVGFLAG_PERTAPE HIST_TYPE1D_PERTAPE HIST_NHTFRQ HIST_NDENS HIST_MFILT HIST_FINCL1 HIST_FINCL2 HIST_FINCL3 HIST_FINCL4 HIST_FINCL5 HIST_FINCL6 HIST_FEXCL1 HIST_FEXCL2 HIST_FEXCL3 HIST_FEXCL4 HIST_FEXCL5 HIST_FEXCL6                                                     
       SHR_CONST_CDAY                  � @                              
       R8 SHR_KIND_R8                                                                                                                                                                                                          6                                                                                                    32          D@@                               	                      D@@                               
                         p          p            p                          +          D@@                                                       p          p            p                                  +          D@@                                                        p          p            p                                            D@@                                                        p          p            p                                    D@@                                                        p          p            p                                    D@@                                                        p          p            p                                                                           
                
                      �@        86400.0          @@@                                                    @@@                                                     @@@                                                     @@@                                                    @@@                                                    @@@                                                    @@@                                                     @@@                                                    @@                                                     @@@                                                    @@@                                                    @@@                                                    @@                                                     @@@                                                    @@@                                                    @@@                                                     @@@                              !                      @@@                              "                      @@@                              #                      @@@                              $                      @@                               %                      @@                               &     
                 @@                               '     
                 @@                               (     
                 @@                               )     
                 @@                               *                      @@@                               +                      @@@                               ,                      @@@                               -                      @@@                               .                      @@@                               /                      @@@                              0     
                                                    1                       @@                              2                       @@                              3                                                         4                                
        L            1275068698                                             5                                
        L            1275069467                                             6                                
        L            1275069469                                             7                                
       ) L            1275070505#         @                                   8                    #CAM_CASEID 9   #CAM_CTITLE :   #CAM_IRAD ;   #CAM_NSREST <   #CAM_CRTINIC =   #CAM_NHTFRQ >   #CAM_MFILT ?   #CAM_IRT @             
                               9                    1           
                               :                    1           
                                 ;                     
                                 <                     
                               =                    1           
                                 >                     
                                 ?                     
                                 @           #         @                                   A                                @                                B                      @@@                              C            +           @                               D                         p          p            p                                  #         @                                  E                    #CONTROL_SPMD%NPES F                                            F               �   k      fn#fn      @   J   CLM_VARCTL    K  @   J   SPMDMOD    �  N  J  HISTFILEMOD    �  O   J  SHR_CONST_MOD    (  O   J   SHR_KIND_MOD ,   w  p       R8+SHR_KIND_MOD=SHR_KIND_R8 &   �  q       MAX_TAPES+HISTFILEMOD '   X  r       MAX_NAMLEN+HISTFILEMOD .   �  @       HIST_EMPTY_HTAPES+HISTFILEMOD (   
  �       HIST_DOV2XY+HISTFILEMOD 1   �  �       HIST_AVGFLAG_PERTAPE+HISTFILEMOD 0   :  �       HIST_TYPE1D_PERTAPE+HISTFILEMOD (   �  �       HIST_NHTFRQ+HISTFILEMOD '   j  �       HIST_NDENS+HISTFILEMOD '   �  �       HIST_MFILT+HISTFILEMOD -   �  w       SHR_CONST_CDAY+SHR_CONST_MOD "   		  @       CTITLE+CLM_VARCTL "   I	  @       CASEID+CLM_VARCTL "   �	  @       NSREST+CLM_VARCTL (   �	  @       HIST_CRTINIC+CLM_VARCTL '   	
  @       ARCHIVE_DIR+CLM_VARCTL %   I
  @       MSS_WPASS+CLM_VARCTL #   �
  @       MSS_IRT+CLM_VARCTL "   �
  @       NREVSN+CLM_VARCTL *   	  @       OFFLINE_ATMDIR+CLM_VARCTL #   I  @       FINIDAT+CLM_VARCTL #   �  @       FSURDAT+CLM_VARCTL #   �  @       FPFTCON+CLM_VARCTL '   	  @       FRIVINP_RTM+CLM_VARCTL )   I  @       MKSRF_FVEGTYP+CLM_VARCTL )   �  @       MKSRF_FSOITEX+CLM_VARCTL )   �  @       MKSRF_FSOICOL+CLM_VARCTL )   	  @       MKSRF_FLANWAT+CLM_VARCTL *   I  @       MKSRF_FGLACIER+CLM_VARCTL (   �  @       MKSRF_FURBAN+CLM_VARCTL &   �  @       MKSRF_FLAI+CLM_VARCTL /   	  @       MKSRF_OFFLINE_FGRID+CLM_VARCTL /   I  @       MKSRF_OFFLINE_EDGEN+CLM_VARCTL /   �  @       MKSRF_OFFLINE_EDGEE+CLM_VARCTL /   �  @       MKSRF_OFFLINE_EDGES+CLM_VARCTL /   	  @       MKSRF_OFFLINE_EDGEW+CLM_VARCTL 2   I  @       MKSRF_OFFLINE_FNAVYORO+CLM_VARCTL "   �  @       CONCHK+CLM_VARCTL     �  @       IRAD+CLM_VARCTL "   	  @       WRTDIA+CLM_VARCTL (   I  @       CSM_DOFLXAVE+CLM_VARCTL &   �  @       RTM_NSTEPS+CLM_VARCTL #   �  @       PERTLIM+CLM_VARCTL "   	  @       MASTERPROC+PMGRID #   I  @       RPNTDIR+CLM_VARCTL #   �  @       RPNTFIL+CLM_VARCTL +   �  z       MPI_CHARACTER+MPISHORTHAND )   C  z       MPI_INTEGER+MPISHORTHAND )   �  z       MPI_LOGICAL+MPISHORTHAND '   7  z       MPI_REAL8+MPISHORTHAND    �  �       CONTROL_INIT (   t  L   a   CONTROL_INIT%CAM_CASEID (   �  L   a   CONTROL_INIT%CAM_CTITLE &     @   a   CONTROL_INIT%CAM_IRAD (   L  @   a   CONTROL_INIT%CAM_NSREST )   �  L   a   CONTROL_INIT%CAM_CRTINIC (   �  @   a   CONTROL_INIT%CAM_NHTFRQ '     @   a   CONTROL_INIT%CAM_MFILT %   X  @   a   CONTROL_INIT%CAM_IRT    �  H       CONTROL_PRINT    �  @       MKFSURDAT       @       RPNTPATH    `  �       RUNTYP    �  _       CONTROL_SPMD 0   [  @     CONTROL_SPMD%NPES+SPMD_DYN=NPES 