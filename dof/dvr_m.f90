module dvr_m

   use dof_m
   use base
   implicit none

   type,extends(dof_t),abstract :: dvr_t
      real(dbl),allocatable :: trafo(:,:)
      real(dbl),allocatable :: d1mat(:,:)
      real(dbl),allocatable :: d2mat(:,:)
      contains
      procedure :: init_dvr
   end type dvr_t

   contains

   subroutine init_dvr(dof)
      class(dvr_t),intent(inout) :: dof
      allocate(dof%x(dof%gdim))
      allocate(dof%w(dof%gdim))
      allocate(dof%trafo(dof%gdim, dof%gdim))
      allocate(dof%d1mat(dof%gdim, dof%gdim))
      allocate(dof%d2mat(dof%gdim, dof%gdim))
   end subroutine init_dvr

end module dvr_m
