  �  =   k820309    [          18.0        yta                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/TimeMgmt/ESMF_TimeMod.F ESMF_TIMEMOD                                                     
                                                           
                                                              u #ESMF_TIMEINITIS    #ESMF_TIMEINITUNDEFINED    #ESMF_TIMECOPYINIT                                                           u #ESMF_TIMESETIS                                                           u #ESMF_TIMEGETIS                                                           u #ESMF_TIMEINCREMENTIS                                                           u #ESMF_TIMEDECREMENTIS 	                                                
                                                       0               � @                                '                     #DAY    #TOD                 � D                                                             � D                                                        #ESMF_TOD                   � @                               '                    #TYPE    #SEC    #MSEC                 � D                                                             � D                                                            � D                                               &         @@   X                                                         #DAYS    #SECONDS    #RC    #ESMF_TIME              
@ @                                                    
@ @                                                    F @                                           &         @@   X                                                         #RC    #ESMF_TIME              F @                                           &         @@   X                                                         #ORIG    #RC    #ESMF_TIME              
@ @                                                   #ESMF_TIME              F @                                           #         @      X                                                 #TIME    #DAYS    #SECONDS    #RC              D @                                                    #ESMF_TIME              
@ @                                                    
@ @                                                    F @                                           #         @      X                                                 #TIME    #DAYS    #SECONDS    #RC              
@ @                                                   #ESMF_TIME              D @                                                     D @                                                     F @                                           &         @@   X                                                         #TIME     #DAYS !   #SECONDS "   #RC #   #ESMF_TIME              
@ @                                                    #ESMF_TIME              
@ @                               !                     
@ @                               "                     F @                               #            &         @@   X                             	                            #TIME $   #DAYS %   #SECONDS &   #RC '   #ESMF_TIME              
@ @                               $                    #ESMF_TIME              
@ @                               %                     
@ @                               &                     F @                               '            %         @@                               (                    
       #TIME )   #RC *             
@ @                               )                    #ESMF_TIME              F @                               *            #         @                                   +                    #TIME ,   #ORIG -   #RC .             D @                               ,                     #ESMF_TIME              
@ @                               -                    #ESMF_TIME              F @                               .            #         @                                   /                    #EARLYTIME 0   #LATETIME 1   #DIFF 2   #ISLATER 3   #RC 4             
@ @                               0                    #ESMF_TIME              
@ @                               1                    #ESMF_TIME              D @                               2                     #ESMF_TIME              D @                               3                      F @                               4            #         @                                   5                    #TIME 6   #RC 7             
@ @                               6                    #ESMF_TIME              F @                               7               �   �      fn#fn "   "  @   J   ESMF_BASICUTILMOD    b  @   J   ESMF_TODMOD "   �  �       gen@ESMF_TIMEINIT !   *  T       gen@ESMF_TIMESET !   ~  T       gen@ESMF_TIMEGET '   �  Z       gen@ESMF_TIMEINCREMENT '   ,  Z       gen@ESMF_TIMEDECREMENT /   �  q       ESMF_SUCCESS+ESMF_BASICUTILMOD    �  b       ESMF_TIME    Y  H   !   ESMF_TIME%DAY    �  ^   !   ESMF_TIME%TOD %   �  m       ESMF_TOD+ESMF_TODMOD /   l  H   %   ESMF_TOD%TYPE+ESMF_TODMOD=TYPE -   �  H   %   ESMF_TOD%SEC+ESMF_TODMOD=SEC /   �  H   %   ESMF_TOD%MSEC+ESMF_TODMOD=MSEC     D  ~       ESMF_TIMEINITIS %   �  @   a   ESMF_TIMEINITIS%DAYS (     @   a   ESMF_TIMEINITIS%SECONDS #   B  @   a   ESMF_TIMEINITIS%RC '   �  g       ESMF_TIMEINITUNDEFINED *   �  @   a   ESMF_TIMEINITUNDEFINED%RC "   )  q       ESMF_TIMECOPYINIT '   �  W   a   ESMF_TIMECOPYINIT%ORIG %   �  @   a   ESMF_TIMECOPYINIT%RC    1	  q       ESMF_TIMESETIS $   �	  W   a   ESMF_TIMESETIS%TIME $   �	  @   a   ESMF_TIMESETIS%DAYS '   9
  @   a   ESMF_TIMESETIS%SECONDS "   y
  @   a   ESMF_TIMESETIS%RC    �
  q       ESMF_TIMEGETIS $   *  W   a   ESMF_TIMEGETIS%TIME $   �  @   a   ESMF_TIMEGETIS%DAYS '   �  @   a   ESMF_TIMEGETIS%SECONDS "     @   a   ESMF_TIMEGETIS%RC %   A  �       ESMF_TIMEINCREMENTIS *   �  W   a   ESMF_TIMEINCREMENTIS%TIME *      @   a   ESMF_TIMEINCREMENTIS%DAYS -   `  @   a   ESMF_TIMEINCREMENTIS%SECONDS (   �  @   a   ESMF_TIMEINCREMENTIS%RC %   �  �       ESMF_TIMEDECREMENTIS *   h  W   a   ESMF_TIMEDECREMENTIS%TIME *   �  @   a   ESMF_TIMEDECREMENTIS%DAYS -   �  @   a   ESMF_TIMEDECREMENTIS%SECONDS (   ?  @   a   ESMF_TIMEDECREMENTIS%RC !     b       ESMF_TIMEGETDAYS &   �  W   a   ESMF_TIMEGETDAYS%TIME $   8  @   a   ESMF_TIMEGETDAYS%RC    x  d       ESMF_TIMECOPY #   �  W   a   ESMF_TIMECOPY%TIME #   3  W   a   ESMF_TIMECOPY%ORIG !   �  @   a   ESMF_TIMECOPY%RC    �  �       ESMF_TIMEDIFF (   N  W   a   ESMF_TIMEDIFF%EARLYTIME '   �  W   a   ESMF_TIMEDIFF%LATETIME #   �  W   a   ESMF_TIMEDIFF%DIFF &   S  @   a   ESMF_TIMEDIFF%ISLATER !   �  @   a   ESMF_TIMEDIFF%RC    �  Z       ESMF_TIMEPRINT $   -  W   a   ESMF_TIMEPRINT%TIME "   �  @   a   ESMF_TIMEPRINT%RC 