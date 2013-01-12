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


   function gauss(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: f
      real(dbl)            :: r
      r = 0.d0
      do f=1,size(x)
         r = r + x(f)**2
      enddo
      v = exp(-r)
   end function gauss


   function coulomb2(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: f
      real(dbl)            :: r
      real(dbl),parameter  :: eps = 1.d-8
      integer              :: np
      np = size(x)/2
      r = 0.d0
      do f=1,np
         r = r + (x(np+f) - x(f))**2
      enddo
      v = 1/sqrt(r+eps)
   end function coulomb2


   function coulomb3(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: p1,p2,k1,k2,i
      real(dbl)            :: r(3)
      real(dbl),parameter  :: eps = 1.d-6
      ! x = (x1,x2,x3,y1,y2,y3,z1,z2,z3)
      r(1) = sqrt((x(1)-x(2))**2 + (x(4)-x(5))**2 + (x(7)-x(8))**2) + eps
      r(2) = sqrt((x(2)-x(3))**2 + (x(5)-x(6))**2 + (x(8)-x(9))**2) + eps
      r(3) = sqrt((x(3)-x(1))**2 + (x(6)-x(4))**2 + (x(9)-x(7))**2) + eps
      ! x = (x1,y1,z1,x2,y2,z2,x3,y3,z3)
      !r(1) = sqrt((x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2) + eps
      !r(2) = sqrt((x(4)-x(7))**2 + (x(5)-x(8))**2 + (x(6)-x(9))**2) + eps
      !r(3) = sqrt((x(7)-x(1))**2 + (x(8)-x(2))**2 + (x(9)-x(3))**2) + eps
      v = 1.d0/r(1) + 1.d0/r(2) + 1.d0/r(3)
   end function coulomb3


end module testfunc
