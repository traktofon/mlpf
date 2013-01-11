! vim: set ts=3 sw=3 :
module hiertuck

   use dof
   use tree
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine init_leaves(t,dofs)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      type(dof_tp),intent(in) :: dofs(:)
      type(node_t),pointer    :: no
      integer :: l,ndofs,f,idof
      do l=1,t%numleaves
         no => t%leaves(l)%p
         ndofs = no%nmodes
         allocate(no%ndim(ndofs))
         do f=1,ndofs
            idof = no%dofs(f)
            no%ndim(f) = dofs(idof)%p%gdim
         enddo
         no%plen = product(no%ndim)
      enddo
   end subroutine

end module hiertuck
