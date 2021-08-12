function [dpres] = ps2dpres (ncid,suff)
% Purpose: construct pressure from PS, hyam, and hybm
% Second argument optional

if (nargin == 2)
    ps = netcdf.getVar(ncid,netcdf.inqVarID(ncid,['PS' suff])); 
else
    ps = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'PS')); 
end
    
hyai = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hyai')); 
hybi = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hybi')); 
nz = length(hyai)-1;

szdat = size(ps); 
n = prod(szdat); 
ps = reshape (ps,[1 n]); 

p_interface = nan ([nz+1 n]); 
for k = 1:nz+1
    p_interface (k,:) = 1e5*hyai(k) + ps*hybi(k); 
end
dpres = reshape (diff(p_interface),[nz szdat]); 

    

