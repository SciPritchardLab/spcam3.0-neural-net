clear; close all;
% note Jerry you need to update the path below to wherever you put the
% ps2dpres.m and makeprintfig.m utilities on your own computer:
% note also the code assumes the zonalmean_*.nc pre-processed output files
% exist locally to this same .m file.  So put them in the same place. 

addpath ('~/Library/matlab_utils'); % The ps2dpres.m util is the main dependency.. needed for correct pressure-weighting in the RMSE.

%fid = fopen ('../mapping.txt'); % Jordan's mapping file.
% Full list of trials.
%A = textscan (fid,'%d %d %f %f %d %f %f %f %d %d %s','headerlines',1); 
% Note these lists included many tuning tests that were not tested
% prognostically:
%val_loss = A{4}; 
%all_runids = A{1}; 
%runids = [ 00001  00018  00035  00052  00069  00086  00103 00002  00019  00036  00053  00070  00087  00104 00003  00020  00037  00054  00071  00088  00105 00004  00021  00038  00055  00072  00089  00106 00005  00022  00039  00056  00073  00090  00107 00006  00023  00040  00057  00074  00091  00108 00007  00024  00041  00058  00075  00092  00008  00025  00042  00059  00076  00093  00009  00026  00043  00060  00077  00094  00010  00027  00044  00061  00078  00095  00011  00028  00045  00062  00079  00096  00012  00029  00046  00063  00080  00097 00013  00030  00047  00064  00081  00098 00014  00031  00048  00065  00082  00099 00015  00032  00049  00066  00083  00100 00016  00033  00050  00067  00084  00101 00017  00034  00051  00068  00085  00102];
%if (length(runids)~= length(all_runids))
%    error ('unexpected.'); 
%end
runids = 1:108; 

% Trim non-existent runs (online results were only available from a subset of full
% SHERPA offline results)
isvalid = ones(length(runids),1); 
for k=1:length(runids); 
    fid = fopen (sprintf ('./zonalmean_%05d_Q_rcat.nc',runids(k)),'r'); 
    if (fid == -1)
        isvalid(k) = 0;
    else
        fclose (fid); 
    end
end
runids = runids(isvalid==1); 
%val_loss = val_loss(runids); 


%% Purpose: Calculation+ visualization of the area- and pressure-weighted tropospheric RMSE metric.
% Depends on  files of this format to have been preprocessed, e.g. sprintf ('./results/zonalmean_monthly_rcat_%s_%s.nc',varname,ctrlstr)
% where varnames = PS,T,Q
% See template_preprocess_script.csh which could be run on stampede2 with
% some minor updates.



n=length(runids); 
varnames = {'T','Q'}; % Master variables to analyze; depends on preprocessed results.
nvar = length(varnames); 
ctrlstr = 'control';

CCC = jet;
doposter = 1; % Activates printer-friendly figure formatting for my CSSI poster.
if (doposter)
    hf= makeprintfig (20,7);
    hfnum = get (hf,'Number'); 
    fontsize = 18; 
end
for kvar = 1:nvar % Cycle over 
    figure (hfnum); 
    ha = subplot (1,nvar,kvar); 
    if (doposter)
        set (ha,'Fontsize',fontsize); 
    end
    varname = varnames{kvar}; 
    runid = 1;
    % Note files of this format are assumed to pre-exist.
    % suffix "0" means the CONTROL simulation we use as baseline to assess
    % error.
    
    f0 = sprintf ('./zonalmean_%s_%s_rcat.nc',ctrlstr,varname)
    % Read in the zonal mean climatology:
    data0 = ncread(f0,varname); 
    lat = ncread(f0,'lat'); 
    lev =ncread(f0,'lev');
    iz = 13:30; % Indices of the vertical levels of the troposphere where we will
    % restrict the analysis (important to ignore upper atm, lev indicies 1:12
    
    % It is important to read in the surface pressure separately...
    fp0 = sprintf ('./zonalmean_%s_PS_rcat.nc',ctrlstr)
    ncid = netcdf.open (fp0,'NC_NOWRITE'); 
    %... in order to calculate the pressure thicknesses of each model layer
    %(needed for mass weighting the RMSE so different layers are weighted
    %differenting) 
    dpres = ps2dpres(ncid);  % Uses my utility ps2dpres... could be easily replicated in python.
    % The puchline is that the PS(x,y,t) field plus two auxilliary static
    % vectors (hybi, hyai) determine the 3D pressure field -- p (x,y,z,t).
    dp = squeeze(nanmean(permute(dpres(iz,:,:),[3 2 1]))); % Compute the time mean of dp.
    % We will use the control sim's dp for analyzing the test sims too
    % (rough but fine to 1st order).
    
    [x,y] = meshgrid (lev(iz),lat); 
    coslat = cos(y*pi/180); % for equal-area-weighting of the RMSE.
    
    A0 = squeeze(nanmean(permute(data0(:,iz,:),[3 1 2]))); % Time mean of the control.
    rmsemin = 1e10; 
    
    for k=1:n % loop through the runs we want to analze
        runid = runids(k);
        % Extract the preprocessed zonal mean monthly data for the current
        % run:
        f = sprintf ('./zonalmean_%05d_%s_rcat.nc',runid,varname); 
        data = ncread(f,varname); 
        nt = size (data,3); 
        for it = 1:nt
            % Pull out an individual month.
            A = squeeze(data(:,iz,it)); 
            % Calculate the bias compared to the control simulation's
            % time-mean climatology (assumes CTRL doesn't drift, verified).            
            sqbias = (A-A0).^2; 
            %Mass- and area- weighted mean squared error below ~ 150 hPa
            % Store the result in a vector that records monthly RMSE.
            rmse(it) = sqrt(nansum(nansum(sqbias.*dp.*coslat))/nansum(nansum(dp.*coslat))); 
        end
        % Just for saving stats of the most skillful model
        % (best rmse + duration exceeds 10 months)
        %if (nt > 10 & rmse (nt) < rmsemin)
        if (rmse(nt) < rmsemin)
            rmsemin = rmse(nt);
            bestbias(:,:) = A-A0; 
            ibest = k; ntbest = nt; 
        end
        
%         % Some random shit to colorize lines based on their offline
%         % validation loss:
%         colindex = ceil(lossindex(k)*64);
%         if (colindex < 1) colindex = 1; end
%         if (colindex > 64) colindex = 64; end
% 
%         thecolor = CCC(colindex,:);
%         
%         % Plot the rmse for the current run as a time series:
%         h = plot (1:nt,rmse,'Color',thecolor);  hold on; 
         h = plot (1:nt,rmse);  hold on; 
         if (doposter)
            set (h,'Linewidth',2.0);
          end
        
%         if (lossindex(k) < prctile (lossindex,10))
%             if (doposter)
%                 set (h,'Linewidth',5.0); 
%             else
%                 set (h,'Linewidth',2.0); 
%             end
%         end
        clear rmse; 
        grid on;
        
    end
    % By the en of this loop
        
%    colormap (CCC); caxis ([min(val_loss) max(val_loss)]); 
hc = colorbar ('vert'); 

    hx = xlabel ('Sim-month'); 
    hy = ylabel (sprintf('RMSE of %s (p >= %.1f hPa)',varname,lev(iz(1)))); 
    xlim ([1 27]); 
    if (doposter)
        set (hc,'Fontsize',fontsize);
        set (hx,'Fontsize',fontsize);
        set (hy,'Fontsize',fontsize);
        print (gcf, '-dpdf','Posterfig.pdf'); 
    end
    hf = makeprintfig(5,5); 
    subplot (2,1,1); 
    imagesc(lat,lev(iz),A0'); axis ('xy'); shading flat; 
    set (gca, 'ydir','reverse'); 
    title (sprintf ('Baseline %s',varname)); 
    colorbar ('vert'); 
    xlabel ('Latitude (deg N)'); ylabel ('Pressure (hPa)'); 
  %  clims = caxis; 
%     subplot (3,1,2,'Fontsize',fontsize); 
%     imagesc(lat,lev(iz),A0'+bestbias'); axis ('xy'); shading flat; 
%     xlabel ('Latitude (deg N)'); ylabel ('Pressure (hPa)'); 
%     set (gca, 'ydir','reverse'); 
%     title (sprintf ('Best model=(runid=%05d) at final month=%d, RMSE=%.5e',runids(ibest),ntbest,rmsemin)); 
%     colorbar ('vert');
%     caxis (clims); 
%     
    subplot (2,1,2,'Fontsize',fontsize); 
    xlabel ('Latitude (deg N)'); ylabel ('Pressure (hPa)'); 
            [x,y] = meshgrid (lat,lev(iz)); 

    if (strcmp(varname,'T'))        
        imagesc(lat,lev(iz),bestbias'); axis ('xy'); shading flat; 
        caxis ([-5 5]); hold on; contour (x,y,bestbias',[-2 2],'k'); 
        title (sprintf ('Temperature bias of best model (2 deg contoured)')); 
    elseif (strcmp(varname,'Q'))
        imagesc (lat,lev(iz),bestbias'./A0'*100); axis ('xy'); shading flat; 
        caxis ([-100 100]); hold on; contour (x,y,bestbias',[0 0],'k'); 
        title (sprintf ('Humidity bias of best model (%% diff)')); 
        contour(x,y,bestbias'./A0'*100,[-20 20],'k'); axis ('xy'); shading flat; 

    else
        imagesc(lat,lev,bestbias'); axis ('xy'); shading flat; 
        title (sprintf ('Bias of best model')); 
    end
    set (gca, 'ydir','reverse'); 
        xlabel ('Latitude (deg N)'); ylabel ('Pressure (hPa)'); 

    colorbar ('vert');
    print (gcf, '-dpdf',sprintf('best_climatology_%s.pdf',varname)); 

    figure (hf); 

    
end


    