! vim: set ts=3 sw=3 :
module tree_m

   use base_m
   use logging_m
   implicit none
   private

   public :: node_t, node_tp, tree_t
   public :: make_leaf, make_node, make_tree, &
             examine_tree, leaf_shape, set_maxnbasis, &
             dump_tree_def, load_tree_def, dump_tree_data, load_tree_data

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
      real(dbl),pointer     :: wghts(:) => null()   ! basis weights, dim=(nbasis)
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
         node => tree%postorder(m)%p
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
      type(tree_t),intent(in)      :: t
      integer,intent(out)          :: vdim(:)
      integer,intent(out),optional :: vlen
      integer                      :: nmodes,m
      nmodes = t%numleaves
      do m=1,nmodes
         vdim(m) = t%leaves(m)%p%plen
      enddo
      if (present(vlen))  vlen = product(vdim)
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


   !--------------------------------------------------------------------
   subroutine dump_tree_def(t,lun)
   ! Write a serialized description of the tree structure.
   !--------------------------------------------------------------------
      type(tree_t),intent(in) :: t
      integer,intent(in)      :: lun
      integer                 :: m,f
      type(node_t),pointer    :: no
      ! Write number of nodes in tree.
      write(lun) t%numnodes
      ! Write the nodes in post-order. This makes deserializing easier.
      do m=1,t%numnodes
         no => t%postorder(m)%p
         ! Write leaf flag and number of dofs/submodes
         write(lun) no%isleaf, no%nmodes
         if (no%isleaf) then
            ! Leaf: write the DOF numbers
            write(lun) (no%dofs(f), f=1,no%nmodes)
         else
            ! Non-leaf: write the post-order numbers of the child nodes
            write(lun) (no%modes(f)%p%num, f=1,no%nmodes)
         endif
      enddo
   end subroutine dump_tree_def


   !--------------------------------------------------------------------
   function load_tree_def(lun) result(t)
   ! Read a serialized description of the tree structure and
   ! construct a tree from it.
   !--------------------------------------------------------------------
      integer,intent(in)        :: lun
      type(tree_t),pointer      :: t
      integer                   :: numnodes,nmodes,m,f
      logical                   :: isleaf
      type(node_tp),allocatable :: nodes(:),children(:)
      integer,allocatable       :: modes(:)
      type(node_t),pointer      :: topnode
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
            nodes(m)%p => make_leaf(modes)
         else
            ! Non-leaf: modes contains the post-order numbers of the child nodes
            ! Gather pointers to the children in another array...
            allocate(children(nmodes))
            do f=1,nmodes
               children(f)%p => nodes(modes(f))%p
            enddo
            ! ... and construct a new parent node for these children.
            nodes(m)%p => make_node(children)
            deallocate(children)
         endif
         deallocate(modes)
      enddo
      ! Last node is the topnode.
      topnode => nodes(numnodes)%p
      deallocate(nodes)
      ! Link tree structure.
      t => make_tree(topnode)
   end function load_tree_def



   !--------------------------------------------------------------------
   subroutine dump_tree_data(t,lun)
   !--------------------------------------------------------------------
      type(tree_t),pointer :: t
      integer,intent(in)   :: lun
      type(node_t),pointer :: no
      integer              :: m,g,j
      do m=1,t%numnodes
         no => t%postorder(m)%p
         write(lun) no%plen, no%nbasis
         write(lun) (no%wghts(j), j=1,no%nbasis)
         write(lun) ((no%basis(g,j), g=1,no%plen), j=1,no%nbasis)
      enddo
   end subroutine dump_tree_data
   

   !--------------------------------------------------------------------
   subroutine load_tree_data(t,lun)
   !--------------------------------------------------------------------
      type(tree_t),pointer :: t
      integer,intent(in)   :: lun
      type(node_t),pointer :: no
      integer              :: m,g,j,f,plen,nbasis
      do m=1,t%numnodes
         no => t%postorder(m)%p
         if (.not. no%isleaf) then
            allocate(no%ndim(no%nmodes))
            do f=1,no%nmodes
               no%ndim(f) = no%modes(f)%p%nbasis
            enddo
            no%plen = product(no%ndim)
         endif
         read(lun,err=500) plen,nbasis
         if (plen /= no%plen) &
            call stopnow("tree data doesn't match system tree")
         no%nbasis = nbasis
         allocate(no%wghts(nbasis))
         allocate(no%basis(plen,nbasis))
         read(lun,err=500) (no%wghts(j), j=1,nbasis)
         read(lun,err=500) ((no%basis(g,j), g=1,plen), j=1,nbasis)
      enddo
      return
 500  call stopnow("error reading tree data")
   end subroutine load_tree_data


end module tree_m
