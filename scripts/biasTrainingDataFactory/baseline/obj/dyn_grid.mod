  H  $   k820309    [          18.0        ��#a                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/atm/cam/src/dynamics/sld/dyn_grid.F90 DYN_GRID                                                     
       PLEV                                                     
       R8 SHR_KIND_R8                                                                                                                                                                                                          30%         @                                                           #GLON    #GLAT              
                                                       
                                             #         @                                                       #GLON 	   #GLAT 
   #CNT    #BLOCKID    #BCID              
                                  	                     
                                  
                     
                                                      D                                                          p          5 � p        r        5 � p        r                               D                                                          p          5 � p        r        5 � p        r                      #         @                                                       #BLOCK_FIRST    #BLOCK_LAST                                           D                                                       D                                             %         @                                                           #BLOCKID              
                                             %         @                                                           #BLOCKID    #BCID              
                                                       
                                             #         @                                                       #BLOCKID    #BCID    #LVLSIZ    #LEVELS              
                                                       
                                                       
                                                      D                                                          p          5 � p        r        5 � p        r                      %         @                                                           #BLOCKID    #BCID              
                                                       
                                             %         @                                                           #BLOCKID    #BCID               
                                                       
                                              %         @                                !                          #GET_BLOCK_OWNER_D%NPES "   #BLOCKID #                                            "                      
                                  #              �   n      fn#fn      E   J  PMGRID    S  O   J   SHR_KIND_MOD ,   �  p       R8+SHR_KIND_MOD=SHR_KIND_R8      r       PLEV+PMGRID &   �  d       GET_BLOCK_COORD_CNT_D +   �  @   a   GET_BLOCK_COORD_CNT_D%GLON +   (  @   a   GET_BLOCK_COORD_CNT_D%GLAT "   h  |       GET_BLOCK_COORD_D '   �  @   a   GET_BLOCK_COORD_D%GLON '   $  @   a   GET_BLOCK_COORD_D%GLAT &   d  @   a   GET_BLOCK_COORD_D%CNT *   �  �   a   GET_BLOCK_COORD_D%BLOCKID '   X  �   a   GET_BLOCK_COORD_D%BCID #     �       GET_BLOCK_BOUNDS_D /   �  @   a   GET_BLOCK_BOUNDS_D%BLOCK_FIRST .   �  @   a   GET_BLOCK_BOUNDS_D%BLOCK_LAST $     ]       GET_BLOCK_COL_CNT_D ,   o  @   a   GET_BLOCK_COL_CNT_D%BLOCKID $   �  g       GET_BLOCK_LVL_CNT_D ,     @   a   GET_BLOCK_LVL_CNT_D%BLOCKID )   V  @   a   GET_BLOCK_LVL_CNT_D%BCID #   �  w       GET_BLOCK_LEVELS_D +   	  @   a   GET_BLOCK_LEVELS_D%BLOCKID (   M	  @   a   GET_BLOCK_LEVELS_D%BCID *   �	  @   a   GET_BLOCK_LEVELS_D%LVLSIZ *   �	  �   a   GET_BLOCK_LEVELS_D%LEVELS    �
  g       GET_LON_D "   �
  @   a   GET_LON_D%BLOCKID    (  @   a   GET_LON_D%BCID    h  g       GET_LAT_D "   �  @   a   GET_LAT_D%BLOCKID      @   a   GET_LAT_D%BCID "   O  y       GET_BLOCK_OWNER_D 5   �  @     GET_BLOCK_OWNER_D%NPES+SPMD_DYN=NPES *     @   a   GET_BLOCK_OWNER_D%BLOCKID 