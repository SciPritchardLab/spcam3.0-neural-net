subroutine chkdims (fileid, name, varid, lonid, latid, verbose)
  implicit none
  include 'netcdf.inc'

  integer fileid, varid, lonid, latid
  logical verbose
  character*(*) name

  integer ret
  integer ndims, dimids(nf_max_dims)

  ret = nf_inq_varid (fileid, name, varid)
  if (ret.eq.NF_NOERR) then
    dimids(:) = -999
    ret = nf_inq_varndims (fileid, varid, ndims)
    ret = nf_inq_vardimid (fileid, varid, dimids)
    if (ret.ne.NF_NOERR) then
      write(6,*)'NF_INQ_VAR failed for ',name
      stop 999
    end if
    if (dimids(1).ne.lonid .or. dimids(2).ne.latid) then
      write(6,*)'1st 2 dims of ', name,' must be lon by lat'
      stop 999
    end if
    if (verbose) write(6,*)'Reducing ',name
  else
    stop 999
  end if
end subroutine chkdims
