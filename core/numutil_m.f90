module numutil_m

   use base_m
   implicit none
   private

   public :: deql, dleq, dgeq

   real(dbl),parameter :: eps = 1.0d-10

   contains

   !--------------------------------------------------------------------
   elemental function deql(x,y)
   ! Checks whether x and y are "equal" (within precision eps).
   !--------------------------------------------------------------------
      real(dbl),intent(in) :: x,y
      logical              :: deql
      deql = (abs(x-y) <= 0.5d0*eps*(abs(x)+abs(y)))
    end function deql


   !--------------------------------------------------------------------
   elemental function dleq(x,y)
   ! Checks whether x is smaller than or "equal" to y.
   !--------------------------------------------------------------------
      real(dbl),intent(in) :: x,y
      logical              :: dleq
      dleq = (x<y) .or. deql(x,y)
   end function dleq


   !--------------------------------------------------------------------
   elemental function dgeq(x,y)
   ! Checks whether x is greather than or "equal" to y.
   !--------------------------------------------------------------------
      real(dbl),intent(in) :: x,y
      logical              :: dgeq
      dgeq = (x>y) .or. deql(x,y)
   end function dgeq

end module numutil_m
