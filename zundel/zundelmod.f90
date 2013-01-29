! vim: set ts=3 sw=3 :
module zundelmod

   use base
   implicit none

   contains

   function zundel(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: q(0:14), v
      ! dofs in x are ordered as follows:
      !   z  R  a  x  y  la  lb  ua  ub  r1a  r2a  va  r1b  r2b  vb
      ! # 1  2  3  4  5   6   7   8   9   10   11  12   13   14  15
      q( 0) = x( 4) ! x
      q( 1) = x( 5) ! y
      q( 2) = x( 1) ! z
      q( 3) = x( 2) ! R
      q( 4) = x( 3) ! alpha
      q( 5) = x(10) ! R1_a
      q( 6) = x(11) ! R2_a
      q( 7) = x(12) ! theta_a
      q( 8) = x( 8) ! beta_a
      q( 9) = x( 6) ! lambda_a
      q(10) = x(13) ! R1_b
      q(11) = x(14) ! R2_b
      q(12) = x(15) ! theta_b
      q(13) = x( 9) ! beta_b
      q(14) = x( 7) ! lambda_b
      call potpoly(v,q)
      return
   end function zundel

end module zundelmod
