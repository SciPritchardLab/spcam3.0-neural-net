clear; close all;
dothis = 0;
if (dothis)
varnames = {'SPDT','SPDQ','QRL','QRS'}; nvar = length(varnames); 
for kvar = 1:nvar
  varname = char(varnames(kvar)); 
  fid = fopen (sprintf('%s_percentile_declarations_f90.txt',varname),'w'); 
  filename = sprintf('./scratch/%s_rcat.nc',varname); 
  aux = ncread(filename,varname); 
  for k=1:30
    aux2 = squeeze(aux(:,:,k,:)); 
    X = reshape (aux2,[1 prod(size(aux2))]); 
    pvals(k,:) = prctile(X,[1 99]);
    disp (sprintf ('Did %s lev %d',varname,k)); 
  end
  fprintf (fid,sprintf ('%s_1st_percentile = (/ %s /)\n',varname,vec2f90str(pvals(:,1)))); 
  fprintf (fid,sprintf ('%s_99th_percentile = (/ %s /)\n',varname,vec2f90str(pvals(:,2)))); 
end
end
varnames = {'PRECT','FLUT'}; nvar = length(varnames); %2D
for kvar = 1:nvar
  varname = char(varnames(kvar)); 
  fid = fopen (sprintf('%s_percentile_declarations_f90.txt',varname),'w');
  filename = sprintf('./scratch/%s_rcat.nc',varname);   
  aux2 = ncread(filename,varname);
  X = reshape (aux2,[1 prod(size(aux2))]);
  ppvals = prctile(X,[1 99]);
  fprintf (fid,sprintf ('%s_1st_percentile = %e\n',varname,ppvals(1))); 
  fprintf (fid,sprintf ('%s_99th_percentile = %e\n',varname,ppvals(2))); 
end
