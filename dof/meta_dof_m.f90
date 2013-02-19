module meta_dof_m

   use dvr_ho_m
   use dvr_sin_m
   use dvr_exp_m
   implicit none

   contains

   subroutine init_doftyps
      call init_doftyp_ho
      call init_doftyp_sin
      call init_doftyp_exp
   end subroutine init_doftyps

end module meta_dof_m
