module mmap_m

   use ISO_C_BINDING
   use base_m
   implicit none

   ! Fortran2003 bindings for C functions.

   interface
      function c_mmap(fd,lrdwr,length,offset) bind(C,name='c_mmap')
         use ISO_C_BINDING
         implicit none
         type(C_PTR)             :: c_mmap
         integer(C_INT),value    :: fd
         integer(C_INT),value    :: lrdwr
         integer(C_SIZE_T),value :: length
         integer(C_SIZE_T),value :: offset
      end function c_mmap
   end interface

   contains

   function mmap_dbl(fname,nelem,lrdwr) result(array)
      character(len=*),intent(in) :: fname
      integer*8                   :: nelem
      integer,intent(in)          :: lrdwr
      integer(C_INT)              :: fd,crdwr
      real(dbl)                   :: dum(1)
      real(dbl),pointer           :: array(:)
      type(C_PTR)                 :: cptr
      integer(C_SIZE_T)           :: length,offset
      call c_open(fd,fname,lrdwr)
      crdwr = lrdwr
      length = size(dum)*nelem
      offset = 0
      cptr = c_mmap(fd,crdwr,length,offset)
      call C_F_POINTER(cptr,array,[nelem])
   end function mmap_dbl

end module mmap_m
