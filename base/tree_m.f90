! vim: set ts=3 sw=3 :
module tree_m

   use base_m
   use logging_m
   implicit none
   private

   public :: node_t, node_tp, tree_t
   public :: make_leaf, make_node, make_tree, examine_tree, leaf_shape, set_maxnbasis

   type :: node_t
      !--- Local node-related data ---
      logical               :: isleaf               ! .true. => children are dofs
      integer               :: nmodes               ! number of submodes/dofs
      type(node_tp),pointer :: modes(:) => null()   ! pointers to submodes
      integer,pointer       :: dofs(:)  => null()   ! list of dofs
      type(node_t),pointer  :: parent   => null()   ! pointer to parent node
      !--- Tree-related data ---
      type(tree_t),pointer  :: tree => null()       ! pointer to tree that this node belongs to
      integer               :: num                  ! internal number of this node
      integer               :: layer                ! layer (level) of this node inside the tree
      !--- MLPF-related data ----
      integer,pointer       :: ndim(:) => null()    ! number of basis tensors/grid points of the children
      integer               :: plen                 ! product(ndim)
      integer               :: nbasis               ! actual number of basis tensors at this node
      integer               :: maxnbasis = 0        ! maximum number of basis tensors at this node
      real(dbl),pointer     :: basis(:,:) => null() ! basis tensors, dim=(plen,nbasis)
   end type node_t


   type :: node_tp
      type(node_t),pointer  :: p => null()
   end type node_tp


   type :: tree_t
      type(node_t),pointer  :: topnode => null()      ! the top node
      integer               :: numnodes               ! total number of nodes in tree
      integer               :: numdofs                ! total number of DOFs in tree
      integer               :: numleaves              ! total number of leaves in tree
      integer               :: numlayers              ! number of layers in tree
      type(node_tp),pointer :: preorder(:) => null()  ! nodes in pre-order
      type(node_tp),pointer :: postorder(:) => null() ! nodes in post-order
      type(node_tp),pointer :: leaves(:) => null()    ! leaf nodes
   end type tree_t

   contains


   !--------------------------------------------------------------------
   function make_leaf(dofs) result(leaf)
   !--------------------------------------------------------------------
      type(node_t),pointer :: leaf
      integer,intent(in)   :: dofs(:)
      integer              :: ndofs,f
      ndofs = size(dofs)
      allocate(leaf)
      ! Mark node as leaf.
      leaf%isleaf = .true.
      ! Link the dofs into the node.
      leaf%nmodes = ndofs
      allocate(leaf%dofs(ndofs))
      do f=1,ndofs
         leaf%dofs(f) = dofs(f)
      enddo
   end function make_leaf


   !--------------------------------------------------------------------
   function make_node(children) result(node)
   !--------------------------------------------------------------------
      type(node_t),pointer        :: node       
      type(node_tp),intent(inout) :: children(:)
      integer                     :: nmodes,m   
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
   end function make_node


   !--------------------------------------------------------------------
   function make_tree(topnode) result(tree)
   !--------------------------------------------------------------------
      type(tree_t),pointer              :: tree                
      type(node_t),target,intent(inout) :: topnode             
      integer                           :: numnodes,numdofs,numleaves,numlayers
      integer                           :: idx,m,f
      type(node_t),pointer              :: node
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
         node%num  =  m
      enddo
      call set_layer(topnode,1)
   end function make_tree


   !--------------------------------------------------------------------
   recursive subroutine count_nodes(node,numnodes,numdofs,numleaves,numlayers)
   !--------------------------------------------------------------------
   ! Counts the total number of nodes, dofs, leaves, and layers below
   ! a certain node, itself included.
   !--------------------------------------------------------------------
      implicit none
      type(node_t),intent(in) :: node            
      integer,intent(out)     :: numnodes,numdofs,numleaves,numlayers
      type(node_t),pointer    :: child           
      integer                 :: m,chnodes,chdofs,chleaves,chlayers
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
      type(node_t),target,intent(in) :: node  
      type(node_tp),intent(inout)    :: seq(:)
      integer,intent(inout)          :: idx   
      type(node_t),pointer           :: child 
      integer                        :: m     
      ! First record this node.
      seq(idx)%p => node
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
      type(node_t),target,intent(in) :: node  
      type(node_tp),intent(inout)    :: seq(:)
      integer,intent(inout)          :: idx   
      type(node_t),pointer           :: child 
      integer                        :: m     
      ! First record all children, if any.
      if (.not. node%isleaf) then
         do m=1,node%nmodes
            child => node%modes(m)%p
            call set_postorder(child,seq,idx)
         enddo
      endif
      ! Then record this node.
      seq(idx)%p => node
      idx = idx+1
   end subroutine set_postorder


   !--------------------------------------------------------------------
   recursive subroutine set_layer(node,layer)
   !--------------------------------------------------------------------
      type(node_t),intent(inout) :: node 
      integer,intent(in)         :: layer
      type(node_t),pointer       :: child
      integer                    :: m    
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
      type(tree_t),intent(in) :: t
      integer,allocatable     :: vdim(:)
      integer,intent(out)     :: vlen
      integer                 :: nmodes,m
      nmodes = t%numleaves
      allocate(vdim(nmodes))
      do m=1,nmodes
         vdim(m) = t%leaves(m)%p%plen
      enddo
      vlen = product(vdim)
   end subroutine leaf_shape


   !--------------------------------------------------------------------
   subroutine set_maxnbasis(node,maxnb)
   !--------------------------------------------------------------------
      type(node_t),intent(inout) :: node
      integer,intent(in)         :: maxnb
      node%maxnbasis = maxnb
   end subroutine set_maxnbasis
      

   !-------------------------------------------------------------------
   subroutine examine_tree(t)
   ! Print various information about the tree structure.
   ! Mainly for debugging.
   !--------------------------------------------------------------------
      type(tree_t),intent(in) :: t
      integer                 :: m
      integer,save            :: logid=0
      character(len=400)      :: msg
      character,parameter     :: NL = ACHAR(10)

      call get_logger(logid,"tree")

      write (msg,'(a,4(a,i4,a))') &
         'Tree has:', &
         NL, t%numnodes, ' nodes', &
         NL, t%numdofs, ' dofs', &
         NL, t%numleaves, ' leaves', &
         NL, t%numlayers, ' layers'
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(x,i0))') &
         'Levels are:', NL, &
         (t%preorder(m)%p%layer, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(x,l))') &
         'Leaves?', NL, &
         (t%preorder(m)%p%isleaf, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(x,i0))') &
         'Leaves are:', NL, &
         (t%leaves(m)%p%num, m=1,t%numleaves)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(x,i0))') &
         'Pre-order is:', NL, &
         (t%preorder(m)%p%num, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(x,i0))') &
         'Post-order is:', NL, &
         (t%postorder(m)%p%num, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)
   end subroutine examine_tree

end module tree_m
