! vim: set ts=3 sw=3 :
module zundelmod

   use base
   implicit none

   contains

   function zundel(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: q(0:14), v
      q( 0) = x( 1) ! x
      q( 1) = x( 4) ! y
      q( 2) = x( 8) ! z
      q( 3) = x(12) ! R
      q( 4) = x( 5) ! alpha
      q( 5) = x(13) ! R1_a
      q( 6) = x(14) ! R2_a
      q( 7) = x(15) ! theta_a (va)
      q( 8) = x( 2) ! beta_a (ua)
      q( 9) = x( 6) ! lambda_a
      q(10) = x( 9) ! R1_b
      q(11) = x(10) ! R2_b
      q(12) = x(11) ! theta_b (vb)
      q(13) = x( 3) ! beta_b (ub)
      q(14) = x( 7) ! lambda_b
      call potpoly(v,q)
      return
   end function zundel

end module zundelmod
