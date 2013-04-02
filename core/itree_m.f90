!=======================================================================
module itree_m
!=======================================================================

   use vtree_m
   use dof_m
   use base_m
   implicit none
   private

   public :: make_ileaf, make_inode, dispose_inode, i2vnode

   type,public :: inode_t
      logical                   :: isleaf
      integer                   :: nmodes
      type(inode_tp),pointer    :: modes(:)  => null()
      character(len=c1),pointer :: labels(:) => null()
      type(inode_t),pointer     :: parent    => null()
      integer                   :: val = 0
   end type inode_t

   type,public :: inode_tp
      type(inode_t),pointer  :: p => null()
   end type inode_tp

   contains


   !--------------------------------------------------------------------
   function make_ileaf(labels) result(leaf)
   !--------------------------------------------------------------------
      type(inode_t),pointer        :: leaf
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
   end function make_ileaf


   !--------------------------------------------------------------------
   function make_inode(children) result(node)
   !--------------------------------------------------------------------
      type(inode_t),pointer        :: node
      type(inode_tp),intent(inout) :: children(:)
      integer                      :: nmodes,m
      nmodes = size(children)
      allocate(node)
      node%isleaf = .false.
      node%nmodes = nmodes
      allocate(node%modes(nmodes))
      do m=1,nmodes
         node%modes(m)%p => children(m)%p
         children(m)%p%parent => node
      enddo
   end function make_inode


   !--------------------------------------------------------------------
   recursive subroutine dispose_inode(node)
   !--------------------------------------------------------------------
      type(inode_t),pointer :: node
      integer               :: m
      if (node%isleaf) then
         deallocate(node%labels)
      else
         do m=1,node%nmodes
            call dispose_inode(node%modes(m)%p)
         enddo
         deallocate(node%modes)
      endif
      deallocate(node)
   end subroutine dispose_inode


   !--------------------------------------------------------------------
   recursive function i2vnode(inpn,dofs) result(node)
   ! Converts an inp_node_t to a node_t.
   !--------------------------------------------------------------------
      type(vnode_t),pointer      :: node
      type(inode_t),pointer      :: inpn
      type(dof_tp),intent(in)    :: dofs(:)
      type(vnode_tp),allocatable :: children(:)
      integer,allocatable        :: dofnums(:)
      integer                    :: nmodes,m,f
      character(len=c1)          :: label
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
         node => make_vleaf(dofnums)
         deallocate(dofnums)
      else
         ! Non-leaf: convert children first
         nmodes = inpn%nmodes
         allocate(children(nmodes))
         do m=1,nmodes
            children(m)%p => i2vnode(inpn%modes(m)%p, dofs)
         enddo
         node => make_vnode(children)
         deallocate(children)
      endif
      call set_maxnbasis(node, inpn%val)
   end function i2vnode

end module itree_m
