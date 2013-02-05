! vim: set ts=3 sw=3 :
module pes3cmod

   use base
   implicit none

   contains

   function pes3c(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl)            :: dr1,dr2,dr4,r3(3),zpp,u1,u2,dihedral,angle1,angle2
      real(dbl),parameter  :: d0 = 1.6d0
      dr1      = x(1)
      dr2      = x(2)
      dr4      = x(3)
      zpp      = x(4)
      r3(1)    = x(5)
      r3(2)    = x(7)
      dihedral = x(8)
      u1       = x(6)
      u2       = x(9)
      r3(3)    = zpp*(dr4-2.d0*d0)
      angle1   = dacos(u1)
      angle2   = dacos(u2)
      call pes3cvpd(dr1,dr2,dr4,r3,angle1,angle2,dihedral,v)
      return
   end function pes3c

end module pes3cmod
