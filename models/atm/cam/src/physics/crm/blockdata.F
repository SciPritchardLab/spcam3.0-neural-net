
	block data

	include 'grid.inc'
	include 'params.inc'
    
      	data ggr	/9.80616/
      	data cp 	/1004.64/
      	data diffelq	/2.21e-05/
      	data lcond	/2.5104e+06/
      	data lfus	/0.3336e+06/
      	data lsub	/2.8440e+06/ 
      	data rv		/461./
      	data rgas	/287./     
      	data therco	/2.40e-02/
      	data muelq	/1.717e-05/
      	data tbgmin	/261.16/  
!old  	data tbgmin	/253.16/  
      	data tbgmax	/270.16/   
!old   	data tbgmax	/273.16/   
      	data tprmin	/268.16/  
      	data tprmax	/283.16/  
      	data tgrmin	/223.16/  
      	data tgrmax	/283.16/  
      	data a_rain	/842./ 
      	data b_rain	/0.8/
      	data a_snow	/4.84/ 
      	data b_snow	/0.25/ 
      	data a_grau	/94.5/ 
!      	data a_grau	/127.94/ 
      	data b_grau	/0.5/ 
      	data qcw0	/1.e-3/
      	data qci0l	/1.e-5/   ! over land
      	data qci0o	/1.e-5/   ! over ocean
      	data alphaelq	/0.001/
      	data betaelq	/0.001/
      	data nzeror	/8.e+06/
      	data nzeros	/3.e+06/
      	data nzerog	/4.e+06/
!      	data nzerog	/4.e+04/
      	data rhor	/1000./   
      	data rhos	/100./
      	data rhog	/400./
!      	data rhog	/917./
      	data erccoef	/1.0/
      	data esccoef	/1.0/
      	data esicoef	/0.1/
      	data egccoef	/1.0/
      	data egicoef	/0.1/
      	data qp_threshold /1.e-8/

	end
