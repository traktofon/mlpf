! vim: set ts=3 sw=3 :
program test

   use dof
   use tree
   use graphviz
   use genpot
   use testfunc
   use tuckerdecomp
   implicit none

   integer,parameter         :: ndofs = 6
   integer,parameter         :: ncomb = 2
   integer,parameter         :: gdim = 15
   real(dbl),parameter       :: accuracy = 1.d-8
   type(dof_tp),allocatable  :: dofs(:)  
   type(node_tp),allocatable :: nodes(:) 
   type(tree_t),pointer      :: t
   integer                   :: f,g,nmodes,nleft,m,vdim
   integer,allocatable       :: ndim(:),mdim(:)
   type(basis_t),allocatable :: basis(:)
   real(dbl),allocatable     :: v(:)
   type(node_t),pointer      :: no

   ! Make DOF grids.
   allocate(dofs(ndofs))
   do f=1,ndofs
      allocate(dofs(f)%p)
      dofs(f)%p%gdim = gdim
      write (dofs(f)%p%label, '(a,i0)') '#',f
      allocate(dofs(f)%p%x(gdim))
      do g=1,gdim
         dofs(f)%p%x(g) = dble(g)
      enddo
   enddo

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
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
   write (*,*) t%numleaves, ' leaves'
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
   write (*,*) 'Leaves are:'
   do m=1,t%numleaves
      write (*,'(x,i0)',advance="no") t%leaves(m)%p%num
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
   call mkdot(42,t,dofs)
   write (*,*) 'Wrote graphviz input file to channel 42.'
   write (*,*)

   ! Generate potential.
   vdim = 1
   do f=1,ndofs
      vdim = vdim*dofs(f)%p%gdim
   enddo
   allocate(v(vdim))
   write (*,*) 'Generating potential, size =',vdim,'...'
   call buildpot(coulomb,dofs,v)

   ! Generate initial Potfit (basis tensors + core tensor)
   nmodes = t%numdofs
   allocate(ndim(nmodes))
   allocate(mdim(nmodes))
   allocate(basis(nmodes))
   ndim(:) = gdim
   mdim(:) = gdim
   write (*,*) 'Computing basis tensors...'
   call compute_basis(v, ndim, accuracy, mdim, basis)
   write (*,*) mdim

   ! Do the hierarchical Tucker decomposition.

end program test
