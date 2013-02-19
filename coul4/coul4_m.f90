! vim: set ts=3 sw=3 :
module coul4_m

   use base_m
   implicit none

   contains

   function coul4(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl),parameter  :: eps = 1.d-2
      real(dbl),parameter  :: tfac = 0.5d0
      real(dbl)            :: R,r1,r2,th1,th2,phi
      real(dbl)            :: xA,xB,xC,yC,xD,yD,zD
      real(dbl)            :: rAB,rAC,rAD,rBC,rBD,rCD
      R   = x(1)
      r1  = x(2)
      r2  = x(3)
      th1 = x(4)
      th2 = x(5)
      phi = x(6)
      xA = -tfac*R
      xB = (1.d0-tfac)*R
      xC = r1*cos(th1)
      yC = r1*sin(th1)
      xD = r2*cos(th2)*cos(phi)
      yD = r2*sin(th2)*cos(phi)
      zD = r2*sin(phi)
      rAB = eps + R
      rAC = eps + sqrt((xC-xA)**2 + yC**2)
      rAD = eps + sqrt((xD-xA)**2 + yD**2 + zD**2)
      rBC = eps + sqrt((xC-xB)**2 + yC**2)
      rBD = eps + sqrt((xD-xB)**2 + yD**2 + zD**2)
      rCD = eps + sqrt((xD-xC)**2 + (yD-yC)**2 + zD**2)
      v = 1.d0/rAB - 1.d0/rAC - 1.d0/rAD - 1.d0/rBC - 1.d0/rBD + 1.d0/rCD
   end function coul4

end module coul4_m
