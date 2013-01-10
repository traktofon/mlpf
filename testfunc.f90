! vim: set ts=3 sw=3 :
module testfunc

   use base
   implicit none

   contains

   function coulomb(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: f
      real(dbl)            :: r
      r = 0.d0
      do f=1,size(x)
         r = r + x(f)**2
      enddo
      v = 1/sqrt(r)
   end function coulomb

end module testfunc
