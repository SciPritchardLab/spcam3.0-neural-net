  �5  h   k820309    [          18.0        ��#a                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/lnd/clm2/src/main/lp_coupling.F90 LP_COUPLING              LP_COUPLING_INIT LP_COUPLING_FINALIZE ALLTOALL_CLUMP_TO_CHUNK_INIT ALLTOALL_CLUMP_TO_CHUNK ALLTOALL_CHUNK_TO_CLUMP                                                     
       MPIR8 MPICOM          @       �   @                              
       NPES                                                     
       IAM                      @                              
       GET_NCLUMPS GET_CLUMP_OWNER_ID GET_CLUMP_NCELLS_ID GET_CLUMP_COORD_ID GET_CLUMP_GCELL_INFO          @       �   @                              
       GET_CHUNK_COORD_OWNER_P                  � @                              
       R8 SHR_KIND_R8            @@                                                      @@                                                                                      	            %         @                               
                            %         @                                                          #CID              
                                             %         @                                                          #CID              
                                             #         @                                                      #CID    #NCELLS    #LONS    #LATS              
                                                       
                                                                                                                p          5 O p            5 O p                                                                                         	    p          5 O p            5 O p                          #         @                                                      #CID    #CELL    #GI              
                                                       
                                                                                                    #         @                                                        #         @                                                       #         @                                                     #ALLTOALL_CLUMP_TO_CHUNK_INIT%NPES    #ALLTOALL_CLUMP_TO_CHUNK_INIT%CHUNK_BUF_NRECS    #ALLTOALL_CLUMP_TO_CHUNK_INIT%BLOCK_BUF_NRECS    #ALLTOALL_CLUMP_TO_CHUNK_INIT%NLCOLS    #ALLTOALL_CLUMP_TO_CHUNK_INIT%NGCOLS    #ALLTOALL_CLUMP_TO_CHUNK_INIT%NCHUNKS     #ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK !   #ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK "   #ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM #   #SRFFLX2D 1                                                          @                           #     '�                   #NAME $   #ASDIR %   #ASDIF &   #ALDIR '   #ALDIF (   #LWUP )   #LHF *   #SHF +   #WSX ,   #WSY -   #TREF .   #TS /   #CFLX 0                �                              $                                        �                              %                             
  p          p            p                                       �                              &            X                 
  p          p            p                                       �                              '            �                 
  p          p            p                                       �                              (            �                 
  p          p            p                                       �                              )                            
  p          p            p                                       �                              *            X                
  p          p            p                                       �                              +            �                
  p          p            p                                       �                              ,            �             	   
  p          p            p                                       �                              -                         
   
  p          p            p                                       �                              .            X                
  p          p            p                                       �                              /            �                
  p          p            p                                       �                              0            �                
  p 	         p          p            p          p                                                                                                                                                                                                                                                                                                                                                                                                    !                                                       "                     
D     �                           1             �            5 r "     & 5 r "   5 r !         5 r !   5 r "   p                          #ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM #   #         @                                   2                  #ALLTOALL_CLUMP_TO_CHUNK%NPES 3   #ALLTOALL_CLUMP_TO_CHUNK%CHUNK_BUF_NRECS 4   #ALLTOALL_CLUMP_TO_CHUNK%BLOCK_BUF_NRECS 5   #ALLTOALL_CLUMP_TO_CHUNK%NLCOLS 6   #ALLTOALL_CLUMP_TO_CHUNK%NGCOLS 7   #ALLTOALL_CLUMP_TO_CHUNK%NCHUNKS 8   #ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK 9   #ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK :   #ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM ;   #SRFFLX2D I                                                                                                                           @                           ;     '�                   #NAME <   #ASDIR =   #ASDIF >   #ALDIR ?   #ALDIF @   #LWUP A   #LHF B   #SHF C   #WSX D   #WSY E   #TREF F   #TS G   #CFLX H                �                              <                                        �                              =                             
  p          p            p                                       �                              >            X                 
  p          p            p                                       �                              ?            �                 
  p          p            p                                       �                              @            �                 
  p          p            p                                       �                              A                            
  p          p            p                                       �                              B            X                
  p          p            p                                       �                              C            �                
  p          p            p                                       �                              D            �             	   
  p          p            p                                       �                              E                         
   
  p          p            p                                       �                              F            X                
  p          p            p                                       �                              G            �                
  p          p            p                                       �                              H            �                
  p 	         p          p            p          p                                                                   3                                                     4                                                     5                                                     6                                                     7                                                     8                                                       9                                                       :                     
D     �                           I             �            5 r :     & 5 r :   5 r 9         5 r 9   5 r :   p                          #ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM ;   #         @                                   J                  #ALLTOALL_CHUNK_TO_CLUMP%NPES K   #ALLTOALL_CHUNK_TO_CLUMP%CHUNK_BUF_NRECS L   #ALLTOALL_CHUNK_TO_CLUMP%BLOCK_BUF_NRECS M   #ALLTOALL_CHUNK_TO_CLUMP%NLCOLS N   #ALLTOALL_CHUNK_TO_CLUMP%NGCOLS O   #ALLTOALL_CHUNK_TO_CLUMP%NCHUNKS P   #ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK Q   #ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK R   #ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE S   #SRF_STATE f                                                     @                           S     '@                   #TBOT T   #ZBOT U   #UBOT V   #VBOT W   #QBOT X   #PBOT Y   #FLWDS Z   #PRECSC [   #PRECSL \   #PRECC ]   #PRECL ^   #SOLL _   #SOLS `   #SOLLD a   #SOLSD b   #SRFRAD c   #THBOT d   #TSSUB e                �                              T                              
  p          p            p                                       �                              U            @                 
  p          p            p                                       �                              V            �                 
  p          p            p                                       �                              W            �                 
  p          p            p                                       �                              X                             
  p          p            p                                       �                              Y            @                
  p          p            p                                       �                              Z            �                
  p          p            p                                       �                              [            �                
  p          p            p                                       �                              \                          	   
  p          p            p                                       �                              ]            @             
   
  p          p            p                                       �                              ^            �                
  p          p            p                                       �                              _            �                
  p          p            p                                       �                              `                             
  p          p            p                                       �                              a            @                
  p          p            p                                       �                              b            �                
  p          p            p                                       �                              c            �                
  p          p            p                                       �                              d                             
  p          p            p                                       �                              e             @                
  p 	         p          p            p          p                                                                   K                                                     L                                                     M                                                     N                                                     O                                                     P                                                       Q                                                       R                     
      �                           f             @           5 r R     & 5 r R   5 r Q         5 r Q   5 r R   p                          #ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE S      �   m      fn#fn !     �   b   uapp(LP_COUPLING    �  M   J  MPISHORTHAND    �  E   J  SPMD_DYN    "  D   J  PMGRID    f  �   J  LND_GRID      X   J  PHYS_GRID    Y  O   J  SHR_KIND_MOD #   �  @       MPIR8+MPISHORTHAND $   �  @       MPICOM+MPISHORTHAND    (  @       NPES+SPMD_DYN %   h  P       GET_NCLUMPS+LND_GRID ,   �  Y       GET_CLUMP_OWNER_ID+LND_GRID 0     @   a   GET_CLUMP_OWNER_ID%CID+LND_GRID -   Q  Y       GET_CLUMP_NCELLS_ID+LND_GRID 1   �  @   a   GET_CLUMP_NCELLS_ID%CID+LND_GRID ,   �  q       GET_CLUMP_COORD_ID+LND_GRID 0   [  @   a   GET_CLUMP_COORD_ID%CID+LND_GRID 3   �  @   a   GET_CLUMP_COORD_ID%NCELLS+LND_GRID 1   �  �   a   GET_CLUMP_COORD_ID%LONS+LND_GRID 1     �   a   GET_CLUMP_COORD_ID%LATS+LND_GRID .   #  c       GET_CLUMP_GCELL_INFO+LND_GRID 2   �  @   a   GET_CLUMP_GCELL_INFO%CID+LND_GRID 3   �  @   a   GET_CLUMP_GCELL_INFO%CELL+LND_GRID 1   	  @   a   GET_CLUMP_GCELL_INFO%GI+LND_GRID !   F	  H       LP_COUPLING_INIT %   �	  H       LP_COUPLING_FINALIZE -   �	        ALLTOALL_CLUMP_TO_CHUNK_INIT @   �  �      ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM+COMSRF E   �  P   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%NAME+COMSRF F   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ASDIR+COMSRF F   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ASDIF+COMSRF F   4  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ALDIR+COMSRF F   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ALDIF+COMSRF E   l  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%LWUP+COMSRF D     �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%LHF+COMSRF D   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%SHF+COMSRF D   @  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%WSX+COMSRF D   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%WSY+COMSRF E   x  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%TREF+COMSRF C     �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%TS+COMSRF E   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%CFLX+COMSRF @   l  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%NPES+SPMD_DYN=NPES W   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%CHUNK_BUF_NRECS+PHYS_GRID=CHUNK_BUF_NRECS W   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%BLOCK_BUF_NRECS+PHYS_GRID=BLOCK_BUF_NRECS E   ,  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%NLCOLS+PHYS_GRID=NLCOLS E   l  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%NGCOLS+PHYS_GRID=NGCOLS G   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%NCHUNKS+PHYS_GRID=NCHUNKS =   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK+PPGRID =   ,  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK+PPGRID 6   l  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX2D (   ^        ALLTOALL_CLUMP_TO_CHUNK ;   x  �      ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM+COMSRF @   H  P   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%NAME+COMSRF A   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ASDIR+COMSRF A   4  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ASDIF+COMSRF A   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ALDIR+COMSRF A   l  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ALDIF+COMSRF @     �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%LWUP+COMSRF ?   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%LHF+COMSRF ?   @  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%SHF+COMSRF ?   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%WSX+COMSRF ?   x  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%WSY+COMSRF @      �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%TREF+COMSRF >   �   �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%TS+COMSRF @   L!  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%CFLX+COMSRF ;   "  @     ALLTOALL_CLUMP_TO_CHUNK%NPES+SPMD_DYN=NPES R   H"  @     ALLTOALL_CLUMP_TO_CHUNK%CHUNK_BUF_NRECS+PHYS_GRID=CHUNK_BUF_NRECS R   �"  @     ALLTOALL_CLUMP_TO_CHUNK%BLOCK_BUF_NRECS+PHYS_GRID=BLOCK_BUF_NRECS @   �"  @     ALLTOALL_CLUMP_TO_CHUNK%NLCOLS+PHYS_GRID=NLCOLS @   #  @     ALLTOALL_CLUMP_TO_CHUNK%NGCOLS+PHYS_GRID=NGCOLS B   H#  @     ALLTOALL_CLUMP_TO_CHUNK%NCHUNKS+PHYS_GRID=NCHUNKS 8   �#  @     ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK+PPGRID 8   �#  @     ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK+PPGRID 1   $  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX2D (   �$  �      ALLTOALL_CHUNK_TO_CLUMP =   �&       ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE+COMSRF B   �'  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%TBOT+COMSRF B   y(  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%ZBOT+COMSRF B   )  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%UBOT+COMSRF B   �)  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%VBOT+COMSRF B   M*  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%QBOT+COMSRF B   �*  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PBOT+COMSRF C   �+  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%FLWDS+COMSRF D   !,  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECSC+COMSRF D   �,  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECSL+COMSRF C   Y-  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECC+COMSRF C   �-  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECL+COMSRF B   �.  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLL+COMSRF B   -/  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLS+COMSRF C   �/  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLLD+COMSRF C   e0  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLSD+COMSRF D   1  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SRFRAD+COMSRF C   �1  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%THBOT+COMSRF C   92  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%TSSUB+COMSRF ;   �2  @     ALLTOALL_CHUNK_TO_CLUMP%NPES+SPMD_DYN=NPES R   53  @     ALLTOALL_CHUNK_TO_CLUMP%CHUNK_BUF_NRECS+PHYS_GRID=CHUNK_BUF_NRECS R   u3  @     ALLTOALL_CHUNK_TO_CLUMP%BLOCK_BUF_NRECS+PHYS_GRID=BLOCK_BUF_NRECS @   �3  @     ALLTOALL_CHUNK_TO_CLUMP%NLCOLS+PHYS_GRID=NLCOLS @   �3  @     ALLTOALL_CHUNK_TO_CLUMP%NGCOLS+PHYS_GRID=NGCOLS B   54  @     ALLTOALL_CHUNK_TO_CLUMP%NCHUNKS+PHYS_GRID=NCHUNKS 8   u4  @     ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK+PPGRID 8   �4  @     ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK+PPGRID 2   �4  �   a   ALLTOALL_CHUNK_TO_CLUMP%SRF_STATE 