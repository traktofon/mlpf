! vim: set ts=3 sw=3 :
module mlpf_io_m

   use vtree_m
   use dof_m
   use dof_io_m
   use hiertuck_m
   use base_m
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine mlpf_from_file(fname,t,dofs)
   !--------------------------------------------------------------------
      character(len=*),intent(in) :: fname     ! the filename
      type(vtree_t),pointer       :: t         ! the MLPF tree
      type(dof_tp),pointer        :: dofs(:)   ! the DOF definitions
      integer                     :: lun,ierr,m
      type(vnode_t),pointer       :: no
      real(dbl)                   :: versnum

      open(newunit=lun, file=trim(fname), status="old", form="unformatted", iostat=ierr)
      if (ierr /= 0) &
         call stopnow('error opening file "'//trim(fname)//'"')

      read (lun) versnum
      dofs => rddvrdef(lun,versnum)
      t => load_vtree_def(lun)
      do m=1,t%numleaves
         no => t%leaves(m)%p
         call init_vleaf(no,dofs)
      enddo
      call load_vtree_data(t,lun)

   end subroutine mlpf_from_file

end module mlpf_io_m

