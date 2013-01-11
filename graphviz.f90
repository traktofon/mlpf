! vim: set ts=3 sw=3 :
module graphviz

   use dof
   use tree
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine mkdot(iout,t,dofs)
   ! Creates an input file for graphviz, for visualising the tree.
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)      :: iout    ! output channel
      type(tree_t),intent(in) :: t       ! the tree
      type(dof_tp),intent(in) :: dofs(:) ! the DOFs
      type(node_t),pointer    :: cnode  
      integer                 :: i,f,idof 
      ! output header
      write(iout,'(5(a,/),a)') &
         'digraph G {', &
         '// customize the following three lines to your liking', &
         'graph [layout=dot rankdir=BT nodesep=0.1 ranksep=0.3]', &
         'node  [fontname=Arial width=0.25 height=0.25 margin=0.05]', &
         'edge  [fontname=Arial]', &
         '// the rest describes the ML tree'
      ! loop over nodes
      do i=1,t%numnodes
         cnode => t%preorder(i)%p
         ! output this node (label by number)
         write (iout,'(a,i0,a,i0,a)') &
            'n', cnode%num, &
            ' [shape=circle label="', cnode%num, &
            '" margin=0]'
         if (cnode%isleaf) then
            ! loop over primitive subnodes
            do f=1,cnode%nmodes
               idof = cnode%dofs(f)
               ! output subnode
               write (iout,'(a,i0,3a)') 'f', idof, &
                  ' [shape=box,label="', trim(dofs(idof)%p%label), '"]'
               ! output edge
               write (iout,'(a,i0,a,i0,a)') &
                  'f', idof, ' -> n', cnode%num, &
                  ' [dir=none]'
            end do
         endif
         ! output edge to parent
         if (associated(cnode%parent)) then
            write (iout,'(a,i0,a,i0,a)') &
               'n', cnode%num, ' -> n', cnode%parent%num, &
               ' [dir=none]'
         endif
      end do
      ! output footer
      write(iout,'(a)') '}'
   end subroutine mkdot

end module graphviz
