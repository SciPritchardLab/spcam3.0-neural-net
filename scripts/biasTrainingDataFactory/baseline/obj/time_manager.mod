  L  à   k820309    [          18.0        ö#a                                                                                                          
       /home1/08110/tg874091/repos/spcam3.0-neural-net/models/atm/cam/src/control/time_manager.F90 TIME_MANAGER       $       TIMEMGR_PRESET TIMEMGR_INIT ADVANCE_TIMESTEP GET_STEP_SIZE GET_NSTEP GET_CURR_DATE GET_PREV_DATE GET_START_DATE GET_REF_DATE GET_PERP_DATE GET_CURR_TIME GET_CURR_CALDAY IS_FIRST_STEP IS_FIRST_RESTART_STEP IS_END_CURR_DAY IS_END_CURR_MONTH IS_LAST_STEP IS_PERPETUAL TIMEMGR_WRITE_RESTART TIMEMGR_READ_RESTART TIMEMGR_RESTART CALENDAR DTIME NESTEP NELAPSE START_YMD START_TOD STOP_YMD STOP_TOD REF_YMD REF_TOD PERPETUAL_YMD PERPETUAL_RUN IC_YMD IC_TOD TM_AQUA_PLANET                                                     
       MASTERPROC                                                     
  !     ESMF_ERRHANDLERSETTYPE ESMF_ERR_RETURN ESMF_ERRPRINT ESMF_SUCCESS ESMF_TIME ESMF_TIMEINIT ESMF_TIMEGET ESMF_TIMEGETDAYS ESMF_TIMEINCREMENT ESMF_TIMEDECREMENT ESMF_DATE ESMF_DATEINIT ESMF_GREGORIAN ESMF_NO_LEAP ESMF_DATEGET ESMF_DATEINCREMENTSEC ESMF_DATEINCREMENTDAY ESMF_DATEDECREMENT ESMF_DATEDIFF ESMF_DATEGETFLTDAYOFYEAR ESMF_TIMEMGR ESMF_TIMEMGRINIT ESMF_TIMEMGRADVANCE ESMF_TIMEMGRGETNSTEP ESMF_TIMEMGRGETSTEPSIZE ESMF_TIMEMGRGETSTARTDATE ESMF_TIMEMGRGETBASEDATE ESMF_TIMEMGRLASTSTEP ESMF_TIMEMGRGETCURRDATE ESMF_TIMEMGRGETPREVDATE ESMF_DATEISLATER ESMF_TIMEMGRRESTARTWRITE ESMF_TIMEMGRRESTARTREAD                      @                              
       TO_UPPER                      @                              
       DYCORE_IS                                                     
       MPICOM MPIINT MPILOG                                                     
       R8 SHR_KIND_R8                                                       u #ESMF_TIMEINITIS    #ESMF_TIMEINITUNDEFINED    #ESMF_TIMECOPYINIT    &         @   @                                                       #DAYS    #SECONDS 	   #RC 
   #ESMF_TIME              
                                                       
                                  	                                                      
            &         @   @                                                       #RC    #ESMF_TIME                                                           &         @   @                                                        #ORIG    #RC    #ESMF_TIME              
                                                      #ESMF_TIME                                                                                                                u #ESMF_TIMEGETIS    #         @     @                                               #TIME    #DAYS    #SECONDS    #RC              
                                                      #ESMF_TIME                                                                                                                                                                                                                               u #ESMF_TIMEINCREMENTIS    &         @   @                                                        #TIME    #DAYS    #SECONDS    #RC    #ESMF_TIME              
                                                      #ESMF_TIME              
                                                       
                                                                                                                                                          u #ESMF_TIMEDECREMENTIS    &         @   @                                                        #TIME    #DAYS    #SECONDS    #RC    #ESMF_TIME              
                                                      #ESMF_TIME              
                                                       
                                                                                                                                                          u #ESMF_DATEINITIS     #ESMF_DATEINITUNDEFINED &   #ESMF_DATECOPYINIT (   &         @   @                                                       #TYPE !   #YEARMMDD "   #TOD #   #RC $   #ESMF_DATE %             
                                  !                     
                                  "                     
                                  #                                                      $            &         @   @                            &                           #RC '   #ESMF_DATE %                                              '            &         @   @                           (                           #ORIG )   #RC *   #ESMF_DATE %             
                                  )                   #ESMF_DATE %                                              *                                                                 u #ESMF_DATEGETIS +   #         @     @                           +                    #DATE ,   #YEARMMDD -   #TOD .   #RC /             
                                  ,                   #ESMF_DATE %                                              -                                                       .                                                       /                                                                  u #ESMF_TIMEMGRINITSTD 0   #ESMF_TIMEMGRINITNOBASESTD 7   #ESMF_TIMEMGRINITIS <   #ESMF_TIMEMGRINITNOBASEIS G   &         @   @                           0     È                     #STEPSIZE 1   #STARTDATE 2   #STOPDATE 3   #BASEDATE 4   #RC 5   #ESMF_TIMEMGR 6             
                                  1                    #ESMF_TIME              
                                  2                   #ESMF_DATE %             
                                  3                   #ESMF_DATE %             
                                  4                   #ESMF_DATE %                                              5            &         @   @                            7     È                     #STEPSIZE 8   #STARTDATE 9   #STOPDATE :   #RC ;   #ESMF_TIMEMGR 6             
                                  8                    #ESMF_TIME              
                                  9                   #ESMF_DATE %             
                                  :                   #ESMF_DATE %                                              ;            &         @   @                            <     È                  
   #STEPDAYS =   #STEPSECS >   #STARTCALENDARDATE ?   #STARTTOD @   #STOPCALENDARDATE A   #STOPTOD B   #BASECALENDARDATE C   #BASETOD D   #TYPE E   #RC F   #ESMF_TIMEMGR 6             
                                  =                     
                                  >                     
                                  ?                     
                                  @                     
                                  A                     
                                  B                     
                                  C                     
                                  D                     
                                  E                                                      F            &         @   @                            G     È                     #STEPDAYS H   #STEPSECS I   #STARTCALENDARDATE J   #STARTTOD K   #STOPCALENDARDATE L   #STOPTOD M   #TYPE N   #RC O   #ESMF_TIMEMGR 6             
                                  H                     
                                  I                     
                                  J                     
                                  K                     
                                  L                     
                                  M                     
                                  N                                                      O                                                                 u #ESMF_TIMEMGRGETSTEPSIZESTD P   #ESMF_TIMEMGRGETSTEPSIZEIS T   #         @     @                            P                    #TIMEMGR Q   #STEPSIZE R   #RC S             
                                  Q     È             #ESMF_TIMEMGR 6                                              R                     #ESMF_TIME                                               S            #         @     @                           T                    #TIMEMGR U   #DAYS V   #SECONDS W   #RC X             
                                  U     È             #ESMF_TIMEMGR 6                                              V                                                       W                                                       X                                                                 u #ESMF_TIMEMGRRESTARTWRITEIS Y   #         @     @                           Y                    #TIMEMGR Z   #TYPE [   #NSTEP \   #STEPDAYS ]   #STEPSEC ^   #STARTYYMMDD _   #STARTSEC `   #STOPYYMMDD a   #STOPSEC b   #BASEYYMMDD c   #BASESEC d   #CURRYYMMDD e   #CURRSEC f   #RC g             
                                  Z     È             #ESMF_TIMEMGR 6                                              [                                                       \                                                       ]                                                       ^                                                       _                                                       `                                                       a                                                       b                                                       c                                                       d                                                       e                                                       f                                                       g                                                                  u #ESMF_TIMEMGRRESTARTREADIS h   &         @   @                           h     È                     #TYPE i   #NSTEP j   #STEPDAYS k   #STEPSEC l   #STARTYYMMDD m   #STARTSEC n   #STOPYYMMDD o   #STOPSEC p   #BASEYYMMDD q   #BASESEC r   #CURRYYMMDD s   #CURRSEC t   #RC u   #ESMF_TIMEMGR 6             
                                  i                     
                                  j                     
                                  k                     
                                  l                     
                                  m                     
                                  n                     
                                  o                     
                                  p                     
                                  q                     
                                  r                     
                                  s                     
                                  t                                                      u                           @                               '                     #DAY v   #TOD w                 D                             v                                 D                              w                          #ESMF_TOD x                  À @                         x     '                    #TYPE y   #SEC z   #MSEC {                 D                            y                                 D                            z                                D                            {                                 @                          %     '                    #CALENDAR |   #YEAR    #MONTH    #DAY    #TOD    #JULIANDAY    #DAYOFYEAR                  D                              |     à                      #ESMF_CALENDAR }                  À @                         }     'à                    #TYPE ~   #DIM    #DIMRUNNINGSUM    #DIY                  D                            ~                                 D                                                           p          p            p                                        D                                        p                   p          p            p                                        D                                 Ø                           D                                  à                           D                                  è                           D                                  ð                           D                                          ø              #ESMF_TOD x                 D                                                            D                                                             @                           6     'È                   #NSTEP    #STEPSIZE    #STARTDATE    #STOPDATE    #BASEDATE    #CURRDATE    #PREVDATE                  D                                                              D                                                         #ESMF_TIME                  D                                          (              #ESMF_DATE %                 D                                          H             #ESMF_DATE %                 D                                          h             #ESMF_DATE %                 D                                                       #ESMF_DATE %                 D                                          ¨             #ESMF_DATE %   $        @                                                           #STR    H r      5 O p                              
                                                    1 %         @                                                          #NAME                                                                  1            @@                                                      @@                                                      @@                                           #         @                                                        #         @                                                        #         @                                                       %         @                                                             %         @                                                             #         @                                                      #YR    #MON    #DAY    #TOD     #OFFSET ¡             D                                                       D                                                       D                                                       D @                                                      
 @                               ¡           #         @                                   ¢                    #YR £   #MON ¤   #DAY ¥   #TOD ¦             D                                 £                      D                                 ¤                      D                                 ¥                      D @                               ¦            #         @                                   §                    #YR ¨   #MON ©   #DAY ª   #TOD «             D                                 ¨                      D                                 ©                      D                                 ª                      D @                               «            #         @                                   ¬                    #YR ­   #MON ®   #DAY ¯   #TOD °             D                                 ­                      D                                 ®                      D                                 ¯                      D @                               °            #         @                                   ±                    #YR ²   #MON ³   #DAY ´   #TOD µ             D                                 ²                      D                                 ³                      D                                 ´                      D @                               µ            #         @                                   ¶                    #DAYS ·   #SECONDS ¸             D @                               ·                      D @                               ¸            %         @                                ¹                    
       #OFFSET º             
 @                               º           %         @                                 »                            %         @                                 ¼                            %         @                                 ½                            %         @                                 ¾                            %         @                                 ¿                            %         @                                 À                            #         @                                   Á                    #FTN_UNIT Â             
                                  Â           #         @                                   Ã                    #FTN_UNIT Ä             
                                  Ä           #         @                                   Å                               @@                              Æ                       @@                               Ç                      @                                È                                                       É                      @@                               Ê                      @@                               Ë                       @                               Ì                       @                               Í                       @                               Î                       @                               Ï                       @                               Ð                                                       Ñ                       @                               Ò                                                       Ó                                                       Ô                           @                                LEN        q      fn#fn "     á  b   uapp(TIME_MANAGER    ò  K   J  PMGRID !   =    J  ESMF_TIMEMGMTMOD    Ù  I   J  STRING_UTILS    "  J   J  DYCORE    l  U   J  MPISHORTHAND    Á  O   J  SHR_KIND_MOD /            gen@ESMF_TIMEINIT+ESMF_TIMEMOD -     ~      ESMF_TIMEINITIS+ESMF_TIMEMOD 2     @   a   ESMF_TIMEINITIS%DAYS+ESMF_TIMEMOD 5   V  @   a   ESMF_TIMEINITIS%SECONDS+ESMF_TIMEMOD 0     @   a   ESMF_TIMEINITIS%RC+ESMF_TIMEMOD 4   Ö  g      ESMF_TIMEINITUNDEFINED+ESMF_TIMEMOD 7   =	  @   a   ESMF_TIMEINITUNDEFINED%RC+ESMF_TIMEMOD /   }	  q      ESMF_TIMECOPYINIT+ESMF_TIMEMOD 4   î	  W   a   ESMF_TIMECOPYINIT%ORIG+ESMF_TIMEMOD 2   E
  @   a   ESMF_TIMECOPYINIT%RC+ESMF_TIMEMOD .   
  T       gen@ESMF_TIMEGET+ESMF_TIMEMOD ,   Ù
  q      ESMF_TIMEGETIS+ESMF_TIMEMOD 1   J  W   a   ESMF_TIMEGETIS%TIME+ESMF_TIMEMOD 1   ¡  @   a   ESMF_TIMEGETIS%DAYS+ESMF_TIMEMOD 4   á  @   a   ESMF_TIMEGETIS%SECONDS+ESMF_TIMEMOD /   !  @   a   ESMF_TIMEGETIS%RC+ESMF_TIMEMOD 4   a  Z       gen@ESMF_TIMEINCREMENT+ESMF_TIMEMOD 2   »        ESMF_TIMEINCREMENTIS+ESMF_TIMEMOD 7   C  W   a   ESMF_TIMEINCREMENTIS%TIME+ESMF_TIMEMOD 7     @   a   ESMF_TIMEINCREMENTIS%DAYS+ESMF_TIMEMOD :   Ú  @   a   ESMF_TIMEINCREMENTIS%SECONDS+ESMF_TIMEMOD 5     @   a   ESMF_TIMEINCREMENTIS%RC+ESMF_TIMEMOD 4   Z  Z       gen@ESMF_TIMEDECREMENT+ESMF_TIMEMOD 2   ´        ESMF_TIMEDECREMENTIS+ESMF_TIMEMOD 7   <  W   a   ESMF_TIMEDECREMENTIS%TIME+ESMF_TIMEMOD 7     @   a   ESMF_TIMEDECREMENTIS%DAYS+ESMF_TIMEMOD :   Ó  @   a   ESMF_TIMEDECREMENTIS%SECONDS+ESMF_TIMEMOD 5     @   a   ESMF_TIMEDECREMENTIS%RC+ESMF_TIMEMOD /   S         gen@ESMF_DATEINIT+ESMF_DATEMOD -   Û        ESMF_DATEINITIS+ESMF_DATEMOD 2   c  @   a   ESMF_DATEINITIS%TYPE+ESMF_DATEMOD 6   £  @   a   ESMF_DATEINITIS%YEARMMDD+ESMF_DATEMOD 1   ã  @   a   ESMF_DATEINITIS%TOD+ESMF_DATEMOD 0   #  @   a   ESMF_DATEINITIS%RC+ESMF_DATEMOD 4   c  g      ESMF_DATEINITUNDEFINED+ESMF_DATEMOD 7   Ê  @   a   ESMF_DATEINITUNDEFINED%RC+ESMF_DATEMOD /   
  q      ESMF_DATECOPYINIT+ESMF_DATEMOD 4   {  W   a   ESMF_DATECOPYINIT%ORIG+ESMF_DATEMOD 2   Ò  @   a   ESMF_DATECOPYINIT%RC+ESMF_DATEMOD .     T       gen@ESMF_DATEGET+ESMF_DATEMOD ,   f  q      ESMF_DATEGETIS+ESMF_DATEMOD 1   ×  W   a   ESMF_DATEGETIS%DATE+ESMF_DATEMOD 5   .  @   a   ESMF_DATEGETIS%YEARMMDD+ESMF_DATEMOD 0   n  @   a   ESMF_DATEGETIS%TOD+ESMF_DATEMOD /   ®  @   a   ESMF_DATEGETIS%RC+ESMF_DATEMOD 5   î  ®       gen@ESMF_TIMEMGRINIT+ESMF_TIMEMGRMOD 4     £      ESMF_TIMEMGRINITSTD+ESMF_TIMEMGRMOD =   ?  W   a   ESMF_TIMEMGRINITSTD%STEPSIZE+ESMF_TIMEMGRMOD >     W   a   ESMF_TIMEMGRINITSTD%STARTDATE+ESMF_TIMEMGRMOD =   í  W   a   ESMF_TIMEMGRINITSTD%STOPDATE+ESMF_TIMEMGRMOD =   D  W   a   ESMF_TIMEMGRINITSTD%BASEDATE+ESMF_TIMEMGRMOD 7     @   a   ESMF_TIMEMGRINITSTD%RC+ESMF_TIMEMGRMOD :   Û        ESMF_TIMEMGRINITNOBASESTD+ESMF_TIMEMGRMOD C   p  W   a   ESMF_TIMEMGRINITNOBASESTD%STEPSIZE+ESMF_TIMEMGRMOD D   Ç  W   a   ESMF_TIMEMGRINITNOBASESTD%STARTDATE+ESMF_TIMEMGRMOD C     W   a   ESMF_TIMEMGRINITNOBASESTD%STOPDATE+ESMF_TIMEMGRMOD =   u  @   a   ESMF_TIMEMGRINITNOBASESTD%RC+ESMF_TIMEMGRMOD 3   µ  û      ESMF_TIMEMGRINITIS+ESMF_TIMEMGRMOD <   °  @   a   ESMF_TIMEMGRINITIS%STEPDAYS+ESMF_TIMEMGRMOD <   ð  @   a   ESMF_TIMEMGRINITIS%STEPSECS+ESMF_TIMEMGRMOD E   0  @   a   ESMF_TIMEMGRINITIS%STARTCALENDARDATE+ESMF_TIMEMGRMOD <   p  @   a   ESMF_TIMEMGRINITIS%STARTTOD+ESMF_TIMEMGRMOD D   °  @   a   ESMF_TIMEMGRINITIS%STOPCALENDARDATE+ESMF_TIMEMGRMOD ;   ð  @   a   ESMF_TIMEMGRINITIS%STOPTOD+ESMF_TIMEMGRMOD D   0  @   a   ESMF_TIMEMGRINITIS%BASECALENDARDATE+ESMF_TIMEMGRMOD ;   p  @   a   ESMF_TIMEMGRINITIS%BASETOD+ESMF_TIMEMGRMOD 8   °  @   a   ESMF_TIMEMGRINITIS%TYPE+ESMF_TIMEMGRMOD 6   ð  @   a   ESMF_TIMEMGRINITIS%RC+ESMF_TIMEMGRMOD 9   0  Ø      ESMF_TIMEMGRINITNOBASEIS+ESMF_TIMEMGRMOD B     @   a   ESMF_TIMEMGRINITNOBASEIS%STEPDAYS+ESMF_TIMEMGRMOD B   H  @   a   ESMF_TIMEMGRINITNOBASEIS%STEPSECS+ESMF_TIMEMGRMOD K     @   a   ESMF_TIMEMGRINITNOBASEIS%STARTCALENDARDATE+ESMF_TIMEMGRMOD B   È  @   a   ESMF_TIMEMGRINITNOBASEIS%STARTTOD+ESMF_TIMEMGRMOD J      @   a   ESMF_TIMEMGRINITNOBASEIS%STOPCALENDARDATE+ESMF_TIMEMGRMOD A   H   @   a   ESMF_TIMEMGRINITNOBASEIS%STOPTOD+ESMF_TIMEMGRMOD >      @   a   ESMF_TIMEMGRINITNOBASEIS%TYPE+ESMF_TIMEMGRMOD <   È   @   a   ESMF_TIMEMGRINITNOBASEIS%RC+ESMF_TIMEMGRMOD <   !         gen@ESMF_TIMEMGRGETSTEPSIZE+ESMF_TIMEMGRMOD ;   !  k      ESMF_TIMEMGRGETSTEPSIZESTD+ESMF_TIMEMGRMOD C   ò!  Z   a   ESMF_TIMEMGRGETSTEPSIZESTD%TIMEMGR+ESMF_TIMEMGRMOD D   L"  W   a   ESMF_TIMEMGRGETSTEPSIZESTD%STEPSIZE+ESMF_TIMEMGRMOD >   £"  @   a   ESMF_TIMEMGRGETSTEPSIZESTD%RC+ESMF_TIMEMGRMOD :   ã"  t      ESMF_TIMEMGRGETSTEPSIZEIS+ESMF_TIMEMGRMOD B   W#  Z   a   ESMF_TIMEMGRGETSTEPSIZEIS%TIMEMGR+ESMF_TIMEMGRMOD ?   ±#  @   a   ESMF_TIMEMGRGETSTEPSIZEIS%DAYS+ESMF_TIMEMGRMOD B   ñ#  @   a   ESMF_TIMEMGRGETSTEPSIZEIS%SECONDS+ESMF_TIMEMGRMOD =   1$  @   a   ESMF_TIMEMGRGETSTEPSIZEIS%RC+ESMF_TIMEMGRMOD =   q$  `       gen@ESMF_TIMEMGRRESTARTWRITE+ESMF_TIMEMGRMOD ;   Ñ$       ESMF_TIMEMGRRESTARTWRITEIS+ESMF_TIMEMGRMOD C   Ô%  Z   a   ESMF_TIMEMGRRESTARTWRITEIS%TIMEMGR+ESMF_TIMEMGRMOD @   .&  @   a   ESMF_TIMEMGRRESTARTWRITEIS%TYPE+ESMF_TIMEMGRMOD A   n&  @   a   ESMF_TIMEMGRRESTARTWRITEIS%NSTEP+ESMF_TIMEMGRMOD D   ®&  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STEPDAYS+ESMF_TIMEMGRMOD C   î&  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STEPSEC+ESMF_TIMEMGRMOD G   .'  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STARTYYMMDD+ESMF_TIMEMGRMOD D   n'  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STARTSEC+ESMF_TIMEMGRMOD F   ®'  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STOPYYMMDD+ESMF_TIMEMGRMOD C   î'  @   a   ESMF_TIMEMGRRESTARTWRITEIS%STOPSEC+ESMF_TIMEMGRMOD F   .(  @   a   ESMF_TIMEMGRRESTARTWRITEIS%BASEYYMMDD+ESMF_TIMEMGRMOD C   n(  @   a   ESMF_TIMEMGRRESTARTWRITEIS%BASESEC+ESMF_TIMEMGRMOD F   ®(  @   a   ESMF_TIMEMGRRESTARTWRITEIS%CURRYYMMDD+ESMF_TIMEMGRMOD C   î(  @   a   ESMF_TIMEMGRRESTARTWRITEIS%CURRSEC+ESMF_TIMEMGRMOD >   .)  @   a   ESMF_TIMEMGRRESTARTWRITEIS%RC+ESMF_TIMEMGRMOD <   n)  _       gen@ESMF_TIMEMGRRESTARTREAD+ESMF_TIMEMGRMOD :   Í)       ESMF_TIMEMGRRESTARTREADIS+ESMF_TIMEMGRMOD ?   Ý*  @   a   ESMF_TIMEMGRRESTARTREADIS%TYPE+ESMF_TIMEMGRMOD @   +  @   a   ESMF_TIMEMGRRESTARTREADIS%NSTEP+ESMF_TIMEMGRMOD C   ]+  @   a   ESMF_TIMEMGRRESTARTREADIS%STEPDAYS+ESMF_TIMEMGRMOD B   +  @   a   ESMF_TIMEMGRRESTARTREADIS%STEPSEC+ESMF_TIMEMGRMOD F   Ý+  @   a   ESMF_TIMEMGRRESTARTREADIS%STARTYYMMDD+ESMF_TIMEMGRMOD C   ,  @   a   ESMF_TIMEMGRRESTARTREADIS%STARTSEC+ESMF_TIMEMGRMOD E   ],  @   a   ESMF_TIMEMGRRESTARTREADIS%STOPYYMMDD+ESMF_TIMEMGRMOD B   ,  @   a   ESMF_TIMEMGRRESTARTREADIS%STOPSEC+ESMF_TIMEMGRMOD E   Ý,  @   a   ESMF_TIMEMGRRESTARTREADIS%BASEYYMMDD+ESMF_TIMEMGRMOD B   -  @   a   ESMF_TIMEMGRRESTARTREADIS%BASESEC+ESMF_TIMEMGRMOD E   ]-  @   a   ESMF_TIMEMGRRESTARTREADIS%CURRYYMMDD+ESMF_TIMEMGRMOD B   -  @   a   ESMF_TIMEMGRRESTARTREADIS%CURRSEC+ESMF_TIMEMGRMOD =   Ý-  @   a   ESMF_TIMEMGRRESTARTREADIS%RC+ESMF_TIMEMGRMOD '   .  b       ESMF_TIME+ESMF_TIMEMOD /   .  H   %   ESMF_TIME%DAY+ESMF_TIMEMOD=DAY /   Ç.  ^   %   ESMF_TIME%TOD+ESMF_TIMEMOD=TOD %   %/  m      ESMF_TOD+ESMF_TODMOD /   /  H   %   ESMF_TOD%TYPE+ESMF_TODMOD=TYPE -   Ú/  H   %   ESMF_TOD%SEC+ESMF_TODMOD=SEC /   "0  H   %   ESMF_TOD%MSEC+ESMF_TODMOD=MSEC '   j0  £       ESMF_DATE+ESMF_DATEMOD 9   1  c   %   ESMF_DATE%CALENDAR+ESMF_DATEMOD=CALENDAR /   p1        ESMF_CALENDAR+ESMF_CALENDARMOD 9   ï1  H   %   ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE 7   72     %   ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM K   Ó2     %   ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM 7   o3  H   %   ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY 1   ·3  H   %   ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR 3   ÿ3  H   %   ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH /   G4  H   %   ESMF_DATE%DAY+ESMF_DATEMOD=DAY /   4  ^   %   ESMF_DATE%TOD+ESMF_DATEMOD=TOD ;   í4  H   %   ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY ;   55  H   %   ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR -   }5  °       ESMF_TIMEMGR+ESMF_TIMEMGRMOD 9   -6  H   %   ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP ?   u6  _   %   ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD=STEPSIZE A   Ô6  _   %   ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD=STARTDATE ?   37  _   %   ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD=STOPDATE ?   7  _   %   ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD=BASEDATE ?   ñ7  _   %   ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD=CURRDATE ?   P8  _   %   ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD=PREVDATE &   ¯8         TO_UPPER+STRING_UTILS *   <9  L   a   TO_UPPER%STR+STRING_UTILS !   9  Z       DYCORE_IS+DYCORE &   â9  L   a   DYCORE_IS%NAME+DYCORE $   .:  @       MPICOM+MPISHORTHAND $   n:  @       MPIINT+MPISHORTHAND $   ®:  @       MPILOG+MPISHORTHAND    î:  H       TIMEMGR_PRESET    6;  H       TIMEMGR_INIT !   ~;  H       ADVANCE_TIMESTEP    Æ;  P       GET_STEP_SIZE    <  P       GET_NSTEP    f<  w       GET_CURR_DATE !   Ý<  @   a   GET_CURR_DATE%YR "   =  @   a   GET_CURR_DATE%MON "   ]=  @   a   GET_CURR_DATE%DAY "   =  @   a   GET_CURR_DATE%TOD %   Ý=  @   a   GET_CURR_DATE%OFFSET    >  k       GET_PREV_DATE !   >  @   a   GET_PREV_DATE%YR "   È>  @   a   GET_PREV_DATE%MON "   ?  @   a   GET_PREV_DATE%DAY "   H?  @   a   GET_PREV_DATE%TOD    ?  k       GET_START_DATE "   ó?  @   a   GET_START_DATE%YR #   3@  @   a   GET_START_DATE%MON #   s@  @   a   GET_START_DATE%DAY #   ³@  @   a   GET_START_DATE%TOD    ó@  k       GET_REF_DATE     ^A  @   a   GET_REF_DATE%YR !   A  @   a   GET_REF_DATE%MON !   ÞA  @   a   GET_REF_DATE%DAY !   B  @   a   GET_REF_DATE%TOD    ^B  k       GET_PERP_DATE !   ÉB  @   a   GET_PERP_DATE%YR "   	C  @   a   GET_PERP_DATE%MON "   IC  @   a   GET_PERP_DATE%DAY "   C  @   a   GET_PERP_DATE%TOD    ÉC  _       GET_CURR_TIME #   (D  @   a   GET_CURR_TIME%DAYS &   hD  @   a   GET_CURR_TIME%SECONDS     ¨D  \       GET_CURR_CALDAY '   E  @   a   GET_CURR_CALDAY%OFFSET    DE  P       IS_FIRST_STEP &   E  P       IS_FIRST_RESTART_STEP     äE  P       IS_END_CURR_DAY "   4F  P       IS_END_CURR_MONTH    F  P       IS_LAST_STEP    ÔF  P       IS_PERPETUAL &   $G  V       TIMEMGR_WRITE_RESTART /   zG  @   a   TIMEMGR_WRITE_RESTART%FTN_UNIT %   ºG  V       TIMEMGR_READ_RESTART .   H  @   a   TIMEMGR_READ_RESTART%FTN_UNIT     PH  H       TIMEMGR_RESTART    H  @       CALENDAR    ØH  @       DTIME    I  @       NESTEP    XI  @       NELAPSE    I  @       START_YMD    ØI  @       START_TOD    J  @       STOP_YMD    XJ  @       STOP_TOD    J  @       REF_YMD    ØJ  @       REF_TOD    K  @       PERPETUAL_YMD    XK  @       PERPETUAL_RUN    K  @       IC_YMD    ØK  @       IC_TOD    L  @       TM_AQUA_PLANET .   XL  <      TO_UPPER%LEN+STRING_UTILS=LEN 