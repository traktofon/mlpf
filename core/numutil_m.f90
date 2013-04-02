module numutil_m

   use base_m
   implicit none
   private

   public :: deql

   real(dbl),parameter :: eps = 1.0d-10

   contains

   !--------------------------------------------------------------------
   elemental function deql(x,y)
   ! Checks whether x and y are "equal" (within precision eps)
   !--------------------------------------------------------------------
      real(dbl),intent(in) :: x,y
      logical              :: deql
      deql = (abs(x-y) <= 0.5d0*eps*(abs(x)+abs(y)))
    end function deql

end module numutil_m
