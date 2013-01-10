! vim: set ts=3 sw=3 :
module dof

   use base
   implicit none

   type :: dof_t
      integer               :: gdim  ! number of grid points
      real(dbl),allocatable :: x(:)  ! coordinates of the grid points
      character(len=16)     :: label ! name of this DOF
   end type dof_t

   type :: dof_tp
      type(dof_t),pointer :: p
   end type dof_tp

end module dof
