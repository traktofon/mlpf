program mlpf2npot

   use vtree_m
   use dof_m
   use meta_dof_m
   use mlpf_io_m
   use hiertuck_m
   use logging_m
   use base_m
   implicit none

   character(len=c5) :: mlpfname,buffer
   integer :: ierr,modc
   type(vtree_t),pointer :: tree
   type(dof_tp),pointer :: dofs(:)
   logical,allocatable :: modeexpand(:)
   real(dbl),pointer :: core(:)

   call init_doftyps

   ! log to stdout
   call set_logger("data", LOGLEVEL_DEBUG, 6)

   call get_command_argument(1, mlpfname, status=ierr)
   if (ierr == -1) &
      call stopnow("filename too long")
   if (ierr /= 0) &
      call stopnow("missing 1st arg = filename")

   call get_command_argument(2, buffer, status=ierr)
   if (ierr /= 0) &
      call stopnow("missing 2nd arg = modc")
   read(buffer,*,iostat=ierr) modc
   if (ierr /= 0 .or. modc <= 0) &
      call stopnow("2nd arg = modc must be a positive integer")

   call mlpf_from_file(mlpfname, tree, dofs)

   if (modc > tree%numleaves) &
      call stopnow("MLPF has less modes than modc")
   allocate(modeexpand(tree%numleaves))
   modeexpand(:) = .false.
   modeexpand(modc) = .true.

   call expand_ht(tree, core, modeexpand)

end program mlpf2npot
