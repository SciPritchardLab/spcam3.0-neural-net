  �   >   k820309    [          18.0        �ta                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/atm/cam/src/physics/cam1/moistconvection.F90 MOISTCONVECTION              MFINTI CMFMCA CP GRAV RGRAV RGAS LIMCNV                                                     
       R8 SHR_KIND_R8 #         @                                                      #MOISTCONVECTION!MFINTI!COMHYB    #RAIR    #CPAIR    #GRAVIT    #LATVAP    #RHOWTR                                                      @                               �                    #MFINTI%COMHYB%HYAI    #MFINTI%COMHYB%HYAM    #MFINTI%COMHYB%HYBI    #MFINTI%COMHYB%HYBM    #MFINTI%COMHYB%HYBD    #MFINTI%COMHYB%HYPI 	   #MFINTI%COMHYB%HYPM 
   #MFINTI%COMHYB%HYPD    #MFINTI%COMHYB%PS0    #MFINTI%COMHYB%PSR    #MFINTI%COMHYB%PRSFAC    #MFINTI%COMHYB%NPRLEV              �   @       �                                             
      p          p            p                                            �   @       �                                     �       
      p          p            p                                            �   @       �                                     �      
      p          p            p                                            �   @       �                                     �      
      p          p            p                                            �   @       �                                     �      
      p          p            p                                            �   @       �                  	                   �      
      p          p            p                                            �   @       �                  
                   �      
      p          p            p                                            �   @       �                                     �      
      p          p            p                                            �   @       �                       �      
                 �   @       �                       �      
                 �   @       �                       �      
                 �   @        �                       �                       
                                      
                
                                      
                
                                      
                
                                      
                
                                      
      #         @                                                      #CMFMCA%CHUNK_BUF_NRECS    #CMFMCA%BLOCK_BUF_NRECS    #CMFMCA%NLCOLS    #CMFMCA%NGCOLS    #CMFMCA%NCHUNKS    #CMFMCA%NPES    #CMFMCA%ENDCHUNK    #CMFMCA%BEGCHUNK    #LCHNK    #NCOL    #NSTEP     #ZTODT !   #PMID "   #PDEL #   #RPDEL $   #ZM %   #TPERT &   #QPERT '   #PHIS (   #PBLHT )   #T *   #Q +   #CMFDT ,   #DQ -   #CMFMC .   #CMFDQR /   #CMFSL 0   #CMFLQ 1   #PRECC 2   #QC 3   #CNT 4   #CNB 5   #ICWMR 6   #RLIQ 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                       
  @                                                    
                                                        
  @                              !     
                
  @                              "     �              
 	   p 	         p          p            p          p                                    
                                 #     �              
 
   p 	         p          p            p          p                                    
                                 $     �              
    p 	         p          p            p          p                                    
                                 %     �              
    p 	         p          p            p          p                                    
                                 &                   
    p          p            p                                    
                                 '                   
    p 	         p          p            p          p                                    
                                 (                   
    p          p            p                                    
                                 )                   
    p          p            p                                    
                                 *     �              
    p 	         p          p            p          p                                    
                                 +     �             
    p �         p          p          p            p          p          p                                    D                                ,     �              
     p 	         p          p            p          p                                    D                                -     �             
     p �         p          p          p            p          p          p                                    D                                .     �              
     p 	         p          p            p          p                                    D                                /     �              
     p 	         p          p            p          p                                    D                                0     �              
     p 	         p          p            p          p                                    D                                1     �              
     p 	         p          p            p          p                                    D                                2                   
     p          p            p                                    D                                3     �              
     p 	         p          p            p          p                                    D                                4                   
     p          p            p                                    D                                5                   
     p          p            p                                    D                                6     �              
     p 	         p          p            p          p                                    D                                7                   
     p          p            p                                     @                               8     
                  @                               9     
                  @                               :     
                  @                               ;     
                  @@                               <               �   |      fn#fn %     8   b   uapp(MOISTCONVECTION    T  O   J  SHR_KIND_MOD    �  �       MFINTI .   j  r  �   MOISTCONVECTION!MFINTI!COMHYB #   �  �      MFINTI%COMHYB%HYAI #   �  �      MFINTI%COMHYB%HYAM #   $  �      MFINTI%COMHYB%HYBI #   �  �      MFINTI%COMHYB%HYBM #   l  �      MFINTI%COMHYB%HYBD #     �      MFINTI%COMHYB%HYPI #   �  �      MFINTI%COMHYB%HYPM #   X  �      MFINTI%COMHYB%HYPD "   �  H      MFINTI%COMHYB%PS0 "   D	  H      MFINTI%COMHYB%PSR %   �	  H      MFINTI%COMHYB%PRSFAC %   �	  H      MFINTI%COMHYB%NPRLEV    
  @   a   MFINTI%RAIR    \
  @   a   MFINTI%CPAIR    �
  @   a   MFINTI%GRAVIT    �
  @   a   MFINTI%LATVAP      @   a   MFINTI%RHOWTR    \  S      CMFMCA A   �  @     CMFMCA%CHUNK_BUF_NRECS+PHYS_GRID=CHUNK_BUF_NRECS A   �  @     CMFMCA%BLOCK_BUF_NRECS+PHYS_GRID=BLOCK_BUF_NRECS /   /  @     CMFMCA%NLCOLS+PHYS_GRID=NLCOLS /   o  @     CMFMCA%NGCOLS+PHYS_GRID=NGCOLS 1   �  @     CMFMCA%NCHUNKS+PHYS_GRID=NCHUNKS *   �  @     CMFMCA%NPES+SPMD_DYN=NPES 0   /  @     CMFMCA%ENDCHUNK+PPGRID=ENDCHUNK 0   o  @     CMFMCA%BEGCHUNK+PPGRID=BEGCHUNK    �  @   a   CMFMCA%LCHNK    �  @   a   CMFMCA%NCOL    /  @   a   CMFMCA%NSTEP    o  @   a   CMFMCA%ZTODT    �  �   a   CMFMCA%PMID    c  �   a   CMFMCA%PDEL      �   a   CMFMCA%RPDEL    �  �   a   CMFMCA%ZM      �   a   CMFMCA%TPERT      �   a   CMFMCA%QPERT    �  �   a   CMFMCA%PHIS    [  �   a   CMFMCA%PBLHT    �  �   a   CMFMCA%T    �  �   a   CMFMCA%Q    w  �   a   CMFMCA%CMFDT    +  �   a   CMFMCA%DQ    �  �   a   CMFMCA%CMFMC    �  �   a   CMFMCA%CMFDQR    g  �   a   CMFMCA%CMFSL      �   a   CMFMCA%CMFLQ    �  �   a   CMFMCA%PRECC    c  �   a   CMFMCA%QC      �   a   CMFMCA%CNT    �  �   a   CMFMCA%CNB    ?  �   a   CMFMCA%ICWMR    �  �   a   CMFMCA%RLIQ    �  @       CP    �  @       GRAV       @       RGRAV    G   @       RGAS    �   @       LIMCNV 