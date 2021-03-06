! vim: set ts=3 sw=3 :
module vtree_m

   use base_m
   implicit none
   private

   public :: vnode_t, vnode_tp, vtree_t
   public :: make_vleaf, make_vnode, make_vtree, dispose_vtree, &
             leaf_shape, set_maxnbasis, &
             dump_vtree_def, load_vtree_def, &
             dump_vtree_data, load_vtree_data, &
             vleaf_modenum

   type :: vnode_t
      !--- Local node-related data ---
      logical                :: isleaf               ! .true. => children are dofs
      integer                :: nmodes               ! number of submodes/dofs
      type(vnode_tp),pointer :: modes(:) => null()   ! pointers to submodes
      integer,pointer        :: dofnums(:) => null() ! list of dofs
      type(vnode_t),pointer  :: parent   => null()   ! pointer to parent node
      !--- Tree-related data ---
      type(vtree_t),pointer  :: tree => null()       ! pointer to tree that this node belongs to
      integer                :: num                  ! internal number of this node
      integer                :: postnum              ! post-order number of this node
      integer                :: layer                ! layer (level) of this node inside the tree
      !--- MLPF-related data ----
      integer,pointer        :: ndim(:) => null()    ! number of basis tensors/grid points of the children
      integer                :: plen                 ! product(ndim)
      integer                :: nbasis               ! actual number of basis tensors at this node
      integer                :: maxnbasis = 0        ! maximum number of basis tensors at this node
      real(dbl),pointer      :: wghts(:) => null()   ! basis weights, dim=(nbasis)
      real(dbl),pointer      :: basis(:,:) => null() ! basis tensors, dim=(plen,nbasis)
   end type vnode_t


   type :: vnode_tp
      type(vnode_t),pointer  :: p => null()
   end type vnode_tp


   type :: vtree_t
      type(vnode_t),pointer  :: topnode => null()      ! the top node
      integer                :: numnodes               ! total number of nodes in tree
      integer                :: numdofs                ! total number of DOFs in tree
      integer                :: numleaves              ! total number of leaves in tree
      integer                :: numlayers              ! number of layers in tree
      type(vnode_tp),pointer :: preorder(:) => null()  ! nodes in pre-order
      type(vnode_tp),pointer :: postorder(:) => null() ! nodes in post-order
      type(vnode_tp),pointer :: leaves(:) => null()    ! leaf nodes
   end type vtree_t

   contains


   !--------------------------------------------------------------------
   function vleaf_modenum(no) result(num)
   ! Given a vnode _no_, returns _num_ such that _no_ is the
   ! _num_'s leaf in the vtree, or zero.
   !--------------------------------------------------------------------
      type(vnode_t),intent(in) :: no
      integer                  :: num
      type(vtree_t),pointer    :: tree
      integer                  :: m
      tree => no%tree
      num = 0
      do m=1,tree%numleaves
         if (tree%leaves(m)%p%num == no%num) then
            num = m
            return
         endif
      enddo
   end function vleaf_modenum


   !--------------------------------------------------------------------
   function make_vleaf(dofnums) result(leaf)
   !--------------------------------------------------------------------
      type(vnode_t),pointer :: leaf
      integer,intent(in)    :: dofnums(:)
      integer               :: ndofs,f
      ndofs = size(dofnums)
      allocate(leaf)
      ! Mark node as leaf.
      leaf%isleaf = .true.
      ! Link the dofs into the node.
      leaf%nmodes = ndofs
      allocate(leaf%dofnums(ndofs))
      do f=1,ndofs
         leaf%dofnums(f) = dofnums(f)
      enddo
   end function make_vleaf


   !--------------------------------------------------------------------
   function make_vnode(children) result(node)
   !--------------------------------------------------------------------
      type(vnode_t),pointer        :: node       
      type(vnode_tp),intent(inout) :: children(:)
      integer                      :: nmodes,m   
      nmodes = size(children)
      allocate(node)
      ! Mark node as non-leaf.
      node%isleaf = .false.
      ! Link the submodes into the node, and set their parent.
      node%nmodes = nmodes
      allocate(node%modes(nmodes))
      do m=1,nmodes
         node%modes(m)%p => children(m)%p
         children(m)%p%parent => node
      enddo
   end function make_vnode


   !--------------------------------------------------------------------
   function make_vtree(topnode) result(tree)
   !--------------------------------------------------------------------
      type(vtree_t),pointer              :: tree                
      type(vnode_t),target,intent(inout) :: topnode             
      integer                            :: numnodes,numdofs,numleaves,numlayers
      integer                            :: idx,m,f
      type(vnode_t),pointer              :: node
      allocate(tree)
      ! Set the top node.
      tree%topnode => topnode
      ! Count number of nodes and dofs in whole tree.
      call count_nodes(topnode,numnodes,numdofs,numleaves,numlayers)
      tree%numnodes  = numnodes
      tree%numdofs   = numdofs
      tree%numleaves = numleaves
      tree%numlayers = numlayers
      ! Record the pre-/post-order sequence of the nodes.
      allocate(tree%preorder(numnodes))
      idx = 1
      call set_preorder(topnode,tree%preorder,idx)
      allocate(tree%postorder(numnodes))
      idx = 1
      call set_postorder(topnode,tree%postorder,idx)
      ! Record leaves.
      allocate(tree%leaves(numleaves))
      f = 1
      do m=1,numnodes
         node => tree%preorder(m)%p
         if (node%isleaf) then
            tree%leaves(f)%p => node
            f = f+1
         endif
      enddo
      ! Set tree-related data in nodes.
      do m=1,numnodes
         node => tree%preorder(m)%p
         node%tree => tree
      enddo
      call set_layer(topnode,1)
   end function make_vtree


   !--------------------------------------------------------------------
   recursive subroutine count_nodes(node,numnodes,numdofs,numleaves,numlayers)
   !--------------------------------------------------------------------
   ! Counts the total number of nodes, dofs, leaves, and layers below
   ! a certain node, itself included.
   !--------------------------------------------------------------------
      implicit none
      type(vnode_t),intent(in) :: node            
      integer,intent(out)      :: numnodes,numdofs,numleaves,numlayers
      type(vnode_t),pointer    :: child           
      integer                  :: m,chnodes,chdofs,chleaves,chlayers
      if (node%isleaf) then
         numnodes  = 1
         numdofs   = node%nmodes
         numleaves = 1
         numlayers = 1
      else
         numnodes  = 1
         numdofs   = 0
         numleaves = 0
         numlayers = 1
         do m=1,node%nmodes
            child => node%modes(m)%p
            call count_nodes(child,chnodes,chdofs,chleaves,chlayers)
            numnodes  = numnodes  + chnodes
            numdofs   = numdofs   + chdofs
            numleaves = numleaves + chleaves
            numlayers = max(numlayers,chlayers+1)
         enddo
      endif
   end subroutine count_nodes


   !--------------------------------------------------------------------
   recursive subroutine set_preorder(node,seq,idx)
   !--------------------------------------------------------------------
      type(vnode_t),target,intent(inout) :: node
      type(vnode_tp),intent(inout)       :: seq(:)
      integer,intent(inout)              :: idx
      type(vnode_t),pointer              :: child
      integer                            :: m
      ! First record this node.
      seq(idx)%p => node
      node%num = idx
      idx = idx+1
      ! Then record all children, if any.
      if (.not. node%isleaf) then
         do m=1,node%nmodes
            child => node%modes(m)%p
            call set_preorder(child,seq,idx)
         enddo
      endif
   end subroutine set_preorder


   !--------------------------------------------------------------------
   recursive subroutine set_postorder(node,seq,idx)
   !--------------------------------------------------------------------
      type(vnode_t),target,intent(inout) :: node
      type(vnode_tp),intent(inout)       :: seq(:)
      integer,intent(inout)              :: idx
      type(vnode_t),pointer              :: child
      integer                            :: m
      ! First record all children, if any.
      if (.not. node%isleaf) then
         do m=1,node%nmodes
            child => node%modes(m)%p
            call set_postorder(child,seq,idx)
         enddo
      endif
      ! Then record this node.
      seq(idx)%p => node
      node%postnum = idx
      idx = idx+1
   end subroutine set_postorder


   !--------------------------------------------------------------------
   recursive subroutine set_layer(node,layer)
   !--------------------------------------------------------------------
      type(vnode_t),intent(inout) :: node 
      integer,intent(in)          :: layer
      type(vnode_t),pointer       :: child
      integer                     :: m    
      ! Set the layer of this node.
      node%layer = layer
      ! Set the layer of all its children.
      if (.not. node%isleaf) then
         do m=1,node%nmodes
            child => node%modes(m)%p
            call set_layer(child,layer+1)
         enddo
      endif
   end subroutine set_layer


   !-------------------------------------------------------------------
   subroutine leaf_shape(t,vdim,vlen)
   ! Returns the potential's shape according to the leaf-layer.
   !-------------------------------------------------------------------
      type(vtree_t),intent(in)       :: t
      integer,intent(out)            :: vdim(:)
      real(dbl),intent(out),optional :: vlen
      integer                        :: nmodes,m
      nmodes = t%numleaves
      do m=1,nmodes
         vdim(m) = t%leaves(m)%p%plen
      enddo
      if (present(vlen)) then
         vlen = 1.d0
         do m=1,nmodes
            vlen = vlen * vdim(m)
         enddo
      endif
   end subroutine leaf_shape


   !--------------------------------------------------------------------
   subroutine set_maxnbasis(node,maxnb)
   !--------------------------------------------------------------------
      type(vnode_t),intent(inout) :: node
      integer,intent(in)          :: maxnb
      node%maxnbasis = maxnb
   end subroutine set_maxnbasis
      

   !--------------------------------------------------------------------
   subroutine dump_vtree_def(t,lun)
   ! Write a serialized description of the tree structure.
   !--------------------------------------------------------------------
      type(vtree_t),intent(in) :: t
      integer,intent(in)       :: lun
      integer                  :: m,f
      type(vnode_t),pointer    :: no
      integer*4                :: nmodes1
      logical*4                :: isleaf1
      integer*4,allocatable    :: modes(:)
      ! Write number of nodes in tree.
      write(lun) t%numnodes
      ! Write the nodes in post-order. This makes deserializing easier.
      do m=1,t%numnodes
         no => t%postorder(m)%p
         ! Write leaf flag and number of dofs/submodes
         nmodes1 = no%nmodes
         isleaf1 = no%isleaf
         write(lun) isleaf1, nmodes1
         allocate(modes(no%nmodes))
         if (no%isleaf) then
            ! Leaf: write the DOF numbers
            modes = no%dofnums
         else
            ! Non-leaf: write the post-order numbers of the child nodes
            do f=1,no%nmodes
               modes(f) = no%modes(f)%p%postnum
            enddo
         endif
         write(lun) (modes(f), f=1,no%nmodes)
         deallocate(modes)
      enddo
   end subroutine dump_vtree_def


   !--------------------------------------------------------------------
   function load_vtree_def(lun) result(t)
   ! Read a serialized description of the tree structure and
   ! construct a tree from it.
   !--------------------------------------------------------------------
      integer,intent(in)         :: lun
      type(vtree_t),pointer      :: t
      integer                    :: numnodes,nmodes,m,f
      logical                    :: isleaf
      type(vnode_tp),allocatable :: nodes(:),children(:)
      integer,allocatable        :: modes(:)
      type(vnode_t),pointer      :: topnode
      ! Read total number of nodes and allocate pointer array.
      read(lun) numnodes
      allocate(nodes(numnodes))
      ! The nodes are coming in post-order, i.e. all children come
      ! before their parent.
      do m=1,numnodes
         ! Read leaf flag and number of dofs/submodes
         read(lun) isleaf,nmodes
         ! Read the dof/submode numbers
         allocate(modes(nmodes))
         read(lun) (modes(f), f=1,nmodes)
         if (isleaf) then
            ! Leaf: modes contains the DOF numbers
            nodes(m)%p => make_vleaf(modes)
         else
            ! Non-leaf: modes contains the post-order numbers of the child nodes
            ! Gather pointers to the children in another array...
            allocate(children(nmodes))
            do f=1,nmodes
               children(f)%p => nodes(modes(f))%p
            enddo
            ! ... and construct a new parent node for these children.
            nodes(m)%p => make_vnode(children)
            deallocate(children)
         endif
         deallocate(modes)
      enddo
      ! Last node is the topnode.
      topnode => nodes(numnodes)%p
      deallocate(nodes)
      ! Link tree structure.
      t => make_vtree(topnode)
   end function load_vtree_def



   !--------------------------------------------------------------------
   subroutine dump_vtree_data(t,lun)
   !--------------------------------------------------------------------
      type(vtree_t),pointer :: t
      integer,intent(in)    :: lun
      type(vnode_t),pointer :: no
      integer               :: m,g,j
      logical*4             :: lwghts
      integer*4             :: plen1,nbasis1
      do m=1,t%numnodes
         no => t%postorder(m)%p
         lwghts = associated(no%wghts)
         plen1 = no%plen
         nbasis1 = no%nbasis
         write(lun) plen1, nbasis1, lwghts
         if (lwghts) &
            write(lun) (no%wghts(j), j=1,no%nbasis)
         write(lun) ((no%basis(g,j), g=1,no%plen), j=1,no%nbasis)
      enddo
   end subroutine dump_vtree_data
   

   !--------------------------------------------------------------------
   subroutine load_vtree_data(t,lun)
   !--------------------------------------------------------------------
      type(vtree_t),pointer :: t
      integer,intent(in)    :: lun
      type(vnode_t),pointer :: no
      integer               :: m,g,j,f,plen,nbasis
      logical               :: lwghts
      do m=1,t%numnodes
         no => t%postorder(m)%p
         if (.not. no%isleaf) then
            allocate(no%ndim(no%nmodes))
            do f=1,no%nmodes
               no%ndim(f) = no%modes(f)%p%nbasis
            enddo
            no%plen = product(no%ndim)
         endif
         read(lun,err=500) plen,nbasis,lwghts
         if (plen /= no%plen) &
            call stopnow("tree data doesn't match system tree")
         no%nbasis = nbasis
         if (lwghts) then
            allocate(no%wghts(nbasis))
            read(lun,err=500) (no%wghts(j), j=1,nbasis)
         endif
         allocate(no%basis(plen,nbasis))
         read(lun,err=500) ((no%basis(g,j), g=1,plen), j=1,nbasis)
      enddo
      return
 500  call stopnow("error reading tree data")
   end subroutine load_vtree_data


   !--------------------------------------------------------------------
   recursive subroutine dispose_vnode(node)
   !--------------------------------------------------------------------
      type(vnode_t),pointer :: node
      integer               :: m
      if (node%isleaf) then
         deallocate(node%dofnums)
      else
         do m=1,node%nmodes
            call dispose_vnode(node%modes(m)%p)
         enddo
         deallocate(node%modes)
      endif
      if (associated(node%basis)) deallocate(node%basis)
      if (associated(node%wghts)) deallocate(node%wghts)
      if (associated(node%ndim))  deallocate(node%ndim)
      deallocate(node)
   end subroutine dispose_vnode


   !--------------------------------------------------------------------
   subroutine dispose_vtree(tree)
   !--------------------------------------------------------------------
      type(vtree_t),pointer :: tree
      deallocate(tree%preorder)
      deallocate(tree%postorder)
      deallocate(tree%leaves)
      call dispose_vnode(tree%topnode)
      deallocate(tree)
   end subroutine dispose_vtree

end module vtree_m
