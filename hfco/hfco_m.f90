! vim: set ts=3 sw=3 :
module hfco_m

   use base_m
   implicit none

   contains

   function hfco1(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl)            :: rch,rcf,rco,ohco,ofco,phi
      rch  = x(1)
      rcf  = x(2)
      ohco = x(4)
      ofco = x(5)
      rco  = x(6)
      phi  = x(3)
      call hfco(rch,rcf,rco,ohco,ofco,phi,v)
      return
   end function hfco1

end module hfco_m
