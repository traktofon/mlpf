module meta_dof_m

   use dvr_ho_m
   use dvr_sin_m
   use dvr_exp_m
   use dof_m
   use base_m
   implicit none

   contains

   subroutine init_doftyps
      call init_doftyp_ho
      call init_doftyp_sin
      call init_doftyp_exp
   end subroutine init_doftyps


   ! for compatibility with old test programs
   function new_dof(label,gdim,xi,xf) result (dof)
      class(dof_t),pointer        :: dof
      character(len=*),intent(in) :: label
      integer,intent(in)          :: gdim
      real(dbl),intent(in)        :: xi,xf
      call unpickle_sin(dof,gdim,[0],[xi,xf])
      call dof%set_label(label)
      call dof%init
   end function new_dof

end module meta_dof_m
