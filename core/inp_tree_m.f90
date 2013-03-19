!=======================================================================
module inp_tree_m
!=======================================================================

   use tree_m
   use dof_m
   use base_m
   implicit none
   private

   public :: make_inp_leaf, make_inp_node, dispose_inp_node, inp2node

   type,public :: inp_node_t
      logical                   :: isleaf
      integer                   :: nmodes
      type(inp_node_tp),pointer :: modes(:)  => null()
      character(len=c1),pointer :: labels(:) => null()
      type(inp_node_t),pointer  :: parent    => null()
      integer                   :: val = 0
   end type inp_node_t

   type,public :: inp_node_tp
      type(inp_node_t),pointer  :: p => null()
   end type inp_node_tp

   contains


   !--------------------------------------------------------------------
   function make_inp_leaf(labels) result(leaf)
   !--------------------------------------------------------------------
      type(inp_node_t),pointer     :: leaf
      character(len=c1),intent(in) :: labels(:)
      integer                      :: ndofs,f
      ndofs = size(labels)
      allocate(leaf)
      leaf%isleaf = .true.
      leaf%nmodes = ndofs
      allocate(leaf%labels(ndofs))
      do f=1,ndofs
         leaf%labels(f) = labels(f)
      enddo
   end function make_inp_leaf


   !--------------------------------------------------------------------
   function make_inp_node(children) result(node)
   !--------------------------------------------------------------------
      type(inp_node_t),pointer        :: node
      type(inp_node_tp),intent(inout) :: children(:)
      integer                         :: nmodes,m
      nmodes = size(children)
      allocate(node)
      node%isleaf = .false.
      node%nmodes = nmodes
      allocate(node%modes(nmodes))
      do m=1,nmodes
         node%modes(m)%p => children(m)%p
         children(m)%p%parent => node
      enddo
   end function make_inp_node


   !--------------------------------------------------------------------
   recursive subroutine dispose_inp_node(node)
   !--------------------------------------------------------------------
      type(inp_node_t),pointer :: node
      integer                  :: m
      if (node%isleaf) then
         deallocate(node%labels)
      else
         do m=1,node%nmodes
            call dispose_inp_node(node%modes(m)%p)
         enddo
         deallocate(node%modes)
      endif
      deallocate(node)
   end subroutine dispose_inp_node


   !--------------------------------------------------------------------
   recursive function inp2node(inpn,dofs) result(node)
   ! Converts an inp_node_t to a node_t.
   !--------------------------------------------------------------------
      type(node_t),pointer      :: node
      type(inp_node_t),pointer  :: inpn
      type(dof_tp),intent(in)   :: dofs(:)
      type(node_tp),allocatable :: children(:)
      integer,allocatable       :: dofnums(:)
      integer                   :: nmodes,m,f
      character(len=c1)         :: label
      if (inpn%isleaf) then
         ! Leaf: map modelabels to dof-numbers
         nmodes = inpn%nmodes
         allocate(dofnums(nmodes))
         do m=1,nmodes
            label = inpn%labels(m)
            f = find_dofnum_by_label(label, dofs)
            if (f==0) &
               call stopnow("modelabel not found: "//trim(label))
            dofnums(m) = f
         enddo
         node => make_leaf(dofnums)
         deallocate(dofnums)
      else
         ! Non-leaf: convert children first
         nmodes = inpn%nmodes
         allocate(children(nmodes))
         do m=1,nmodes
            children(m)%p => inp2node(inpn%modes(m)%p, dofs)
         enddo
         node => make_node(children)
         deallocate(children)
      endif
      call set_maxnbasis(node, inpn%val)
   end function inp2node

end module inp_tree_m
