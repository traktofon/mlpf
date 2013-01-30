! vim: set ts=3 sw=3 :
module hfcomod

   use base
   implicit none

   contains

   function hfco1(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl)            :: rch,rcf,rco,ohco,ofco,phi
      ! rch ohco rcf ofco rco phi
      !   1    2   3    4   5   6
      rch = x(1)
      rcf = x(3)
      rco = x(5)
      ohco = x(2)
      ofco = x(4)
      phi = x(6)
      call hfco(rch,rcf,rco,ohco,ofco,phi,v)
      return
   end function hfco1

end module hfcomod
