!-------------------------------------------------------------------------------
!
! Purpose:
!
! Wrapper routines for the netCDF library for input and output data.
!
! Author: Jim Rosinski
!
! $Id: wrap_nf.f90,v 1.1.2.1 2003/10/22 22:53:11 rosinski Exp $
!
!-------------------------------------------------------------------------------

!===============================================================================

   subroutine wrap_redef (nfid)
   implicit none
   include 'netcdf.inc'
   
   integer, intent(in):: nfid

   integer ret      ! NetCDF return code

   ret = nf_redef (nfid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_redef

!===============================================================================

   subroutine wrap_enddef (nfid)
   implicit none
   include 'netcdf.inc'
   
   integer, intent(in):: nfid

   integer ret      ! NetCDF return code

   ret = nf_enddef (nfid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_enddef

!===============================================================================

   subroutine wrap_inq_dimid (nfid, dimname, dimid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension id
!
!-------------------------------------------------------------------------------
   include 'netcdf.inc'
   
   integer, intent(in):: nfid
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname

   integer ret      ! NetCDF return code

   ret = nf_inq_dimid (nfid, dimname, dimid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_inq_dimid

!===============================================================================

   subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension length for a given dimension
!
!-------------------------------------------------------------------------------
   include 'netcdf.inc'
   
   integer, intent(in)::  nfid
   integer, intent(in)::  dimid 
   integer, intent(out):: dimlen
   
   integer ret      ! NetCDF return code

   ret = nf_inq_dimlen (nfid, dimid, dimlen)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_inq_dimlen

!===============================================================================

   subroutine wrap_inq_varid (nfid, varname, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable ID
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'
   
   integer, intent(in):: nfid
   integer, intent(out):: varid
   character*(*), intent(in):: varname
   
   integer ret      ! NetCDF return code

   ret = nf_inq_varid (nfid, varname, varid)
   if (ret/=NF_NOERR) then
     write(6,*)'wrap_inq_varid: id for ',trim(varname),' not found'
     call handle_error (ret)
   end if
   end subroutine wrap_inq_varid

!===============================================================================

   subroutine wrap_get_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the attribute text from the given variable ID and attribute name
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
   character*(*), intent(out):: atttext

   integer ret      ! NetCDF return code

   ret = nf_get_att_text (nfid, varid, attname, atttext)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_get_att_text

!===============================================================================

   subroutine wrap_put_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given attribute text to variable ID.
!
! This routine violates the convetion that the wrapper codes take an identical
! set of arguments as the netcdf library code.  The length of the character
! argument is computed inside the wrapper.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
   character*(*), intent(in):: atttext

   integer ret      ! NetCDF return code
   integer siz

   siz = len_trim(atttext)
   ret = nf_put_att_text (nfid, varid, attname, siz, atttext)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_att_text

!===============================================================================

   subroutine wrap_put_att_double (nfid, varid, attname, xtype, len, &
                                   attval)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given real attribute to the variable id
!
!-------------------------------------------------------------------------------
   use prec
   implicit none
   include 'netcdf.inc'
   
   integer , intent(in):: nfid
   integer , intent(in):: varid
   integer , intent(in):: xtype
   integer , intent(in):: len
   character*(*) , intent(in):: attname
   real(r8) , intent(in):: attval

   integer ret      ! NetCDF return code

   ret = nf_put_att_double (nfid, varid, attname, xtype, len, attval)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_att_double
!===============================================================================

   subroutine wrap_def_dim (nfid, dimname, len, dimid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the input dimension
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'
   integer, intent(in):: nfid
   integer, intent(in):: len
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname
   
   integer ret      ! NetCDF return code

   ret = nf_def_dim (nfid, dimname, len, dimid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_def_dim

!===============================================================================

   subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the given variable
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in)::xtype
   integer, intent(in)::nvdims
   integer, intent(out)::varid
   integer, intent(in):: vdims(nvdims+1)
   character*(*), intent(in):: name
   
   integer ret      ! NetCDF return code

   ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_def_var

!===============================================================================

   subroutine wrap_get_var_double (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   use prec
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_var_double (nfid, varid, arr)

   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_DOUBLE: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_double

!===============================================================================

   subroutine wrap_get_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given integer variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_var_int (nfid, varid, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_int

!===============================================================================

   subroutine wrap_get_vara_double (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given real variable from a input file
!
!-------------------------------------------------------------------------------
   use prec
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in)::varid
   integer, intent(in)::start(*)
   integer, intent(in)::count(*)
   real(r8), intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_vara_double (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_DOUBLE: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_double

!===============================================================================

   subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given integer variable from a input file.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   integer, intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_vara_int (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_int

!===============================================================================

   subroutine wrap_open (path, omode, ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Open a netCDF file
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   character*(*), intent(in):: path
   integer, intent(in):: omode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code

   ret = nf_open (path, omode, ncid)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_OPEN: nf_open failed for file ', trim(path)
     call handle_error (ret)
   end if
   end subroutine wrap_open

!===============================================================================

   subroutine wrap_close (ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Close netCDF file
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: ncid

   integer ret      ! NetCDF return code

   ret = nf_close (ncid)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
     call handle_error (ret)
   end if
   end subroutine wrap_close

!===============================================================================

   subroutine wrap_put_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a integer variable on output file.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_put_var_int (nfid, varid, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_var_int

!===============================================================================

   subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of a integer variable on a output file.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   integer, intent(in):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_put_vara_int (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_int

!===============================================================================

   subroutine wrap_put_vara_text (nfid, varid, start, count, text)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of the given text variable to output file.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   character*(*), intent(in):: text(*)

   integer ret      ! NetCDF return code

   ret = nf_put_vara_text (nfid, varid, start, count, text)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_text

!===============================================================================

   subroutine wrap_put_vara_double (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Output the given portion of the real array.
!
!-------------------------------------------------------------------------------
   use prec
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   real(r8), intent(in):: arr(*)

   integer ret      ! NetCDF return code
   ret = nf_put_vara_double (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_double

!===============================================================================

   subroutine wrap_put_var_double (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put the given real variable to output file.
!
!-------------------------------------------------------------------------------
   use prec
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(in):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_put_var_double (nfid, varid, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_put_var_double

!===============================================================================

   subroutine handle_error(ret)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Handle netCDF errors.
!
!-------------------------------------------------------------------------------
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: ret
   
   write(6,*)nf_strerror(ret)
   call abort
   end subroutine handle_error

!===============================================================================
