! vim: set ts=3 sw=3 :
program test

   use dof
   use tree
   use graphviz
   implicit none

   integer,parameter         :: ndofs = 17
   integer,parameter         :: ncomb = 2
   integer,parameter         :: gdim = 25
   type(dof_tp),allocatable  :: dofs(:)  
   type(node_tp),allocatable :: nodes(:) 
   type(tree_t),pointer      :: t
   integer                   :: f,g,nmodes,nleft,m

   ! Make DOF grids.
   allocate(dofs(ndofs))
   do f=1,ndofs
      allocate(dofs(f)%p)
      dofs(f)%p%gdim = gdim
      write (dofs(f)%p%label, '(a,i0)') 'dof#',f
      allocate(dofs(f)%p%x(gdim))
      do g=1,gdim
         dofs(f)%p%x(g) = dble(g)
      enddo
   enddo

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf(dofs(f:f))
   enddo

   ! Combine modes
   do while (nmodes > 1)
      nleft  = mod(nmodes,ncomb)
      nmodes = nmodes/ncomb
      do m=1,nmodes
         nodes(m)%p => make_node(nodes(ncomb*(m-1)+1:ncomb*m))
      enddo
      if (nleft>1) then
         nodes(nmodes+1)%p => make_node(nodes(ncomb*nmodes+1:ncomb*nmodes+nleft))
         nmodes = nmodes+1
      elseif (nleft==1) then
         nodes(nmodes+1)%p => nodes(ncomb*nmodes+1)%p
         nmodes = nmodes+1
      endif
   enddo

   ! Build tree.
   t => make_tree(nodes(1)%p)
   
   ! Examine everything.
   write (*,*) 'Tree has:'
   write (*,*) t%numnodes, ' nodes'
   write (*,*) t%numdofs, ' dofs'
   write (*,*) t%numlayers, ' layers'
   write (*,*) 'Levels are:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%preorder(m)%p%layer
   enddo
   write (*,*)
   write (*,*) 'Leaves?'
   do m=1,t%numnodes
      write (*,'(x,l)',advance="no") t%preorder(m)%p%isleaf
   enddo
   write (*,*)
   write (*,*) 'Pre-order is:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%preorder(m)%p%num
   enddo
   write (*,*)
   write (*,*) 'Post-order is:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%postorder(m)%p%num
   enddo
   write (*,*)
   call mkdot(42,t)

end program test
