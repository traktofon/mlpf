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

   contains

   function new_dof(lbl,gdim,xi,xf) result (p)
      implicit none
      character(len=16),intent(in) :: lbl
      integer,intent(in)           :: gdim
      real(dbl),intent(in)         :: xi,xf
      type(dof_t),pointer          :: p
      integer                      :: g
      allocate(p)
      p%label = lbl
      p%gdim = gdim
      allocate(p%x(gdim))
      do g=1,gdim
         p%x(g) = xi + (g-1)*(xf-xi)/(gdim-1)
      enddo
   end function new_dof

end module dof
