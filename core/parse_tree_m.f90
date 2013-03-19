!=======================================================================
module parse_tree_m
!=======================================================================
!
! leaf :~ modelabel ( "," modelabel )*
! node :~ ( leaf | "(" node node* ")" ) ( "=" INTEGER )?
!
!=======================================================================

   use tokenize_m
   use strutil_m
   use base_m
   implicit none
   private

   public :: parse_tree, dispose_inp_tree

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

   integer,parameter :: maxndof = 10

   contains


   !--------------------------------------------------------------------
   function make_leaf(labels) result(leaf)
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
   end function make_leaf


   !--------------------------------------------------------------------
   function parse_leaf(t) result(leaf)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: t
      type(inp_node_t),pointer        :: leaf
      character(len=maxtoklen)        :: token
      integer                         :: nf
      character(len=c1)               :: labels(maxndof)
      nf = 0
      do while (.not. t%is_stopped())
         ! read current token
         token = t%get()
         ! store this as modelabel
         nf = nf+1
         if (nf > maxndof) &
            call t%error("too many DOFs in mode")
         if (len_trim(token) > c1) &
            call t%error("modelabel too long")
         labels(nf) = trim(token)
         ! continue parsing
         call t%gofwd
         token = t%get()
         if (token /= ",") exit
         call t%gofwd
      enddo
      ! create leaf 
      leaf => make_leaf(labels(1:nf))
   end function parse_leaf


   !--------------------------------------------------------------------
   function make_node(children) result(node)
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
   end function make_node


   !--------------------------------------------------------------------
   recursive function parse_node(t) result(node)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: t
      type(inp_node_t),pointer        :: node
      character(len=maxtoklen)        :: token
      integer                         :: nm
      type(inp_node_tp)               :: ms(maxndof)
      if (t%is_stopped()) &
         call t%error("expected start or end of mode definition")
      token = t%get()
      if (token == "(") then
         nm = 0
         call t%gofwd
         do
            nm = nm+1
            if (nm > maxndof) &
               call t%error("too many submodes in mode")
            ms(nm)%p => parse_node(t)
            token = t%get()
            if (token == ")") then
               call t%gofwd
               exit
            endif
         enddo
         node => make_node(ms(1:nm))
      else
         node => parse_leaf(t)
      endif
      token = t%get()
      if (token == "=") then
         call t%gofwd
         node%val = parse_int(t)
      endif
   end function parse_node


   !--------------------------------------------------------------------
   function parse_tree(t,topnode) result(flag)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: t
      type(inp_node_t),pointer        :: topnode
      character(len=maxtoklen)        :: token
      logical                         :: flag

      token = t%get()
      call ucase(token)
      if (token /= "TREE-SECTION") then
         flag = .false.
         return
      endif
      call t%clear_stop
      call t%add_stop("END-TREE-SECTION")
      call t%gofwd

      topnode => parse_node(t)
      flag = .true.

      if (t%stopreason() /= STOPREASON_STOPWORD) &
         call t%error("expected end of section")
      call t%clear_stop
      call t%gofwd
   end function parse_tree


   !--------------------------------------------------------------------
   subroutine dispose_inp_tree(topnode)
   !--------------------------------------------------------------------
      type(inp_node_t),pointer :: topnode
      call dispose_node(topnode)
   end subroutine dispose_inp_tree


   !--------------------------------------------------------------------
   recursive subroutine dispose_node(node)
   !--------------------------------------------------------------------
      type(inp_node_t),pointer :: node
      integer                  :: m
      if (node%isleaf) then
         deallocate(node%labels)
      else
         do m=1,node%nmodes
            call dispose_node(node%modes(m)%p)
         enddo
         deallocate(node%modes)
      endif
      deallocate(node)
   end subroutine dispose_node


end module parse_tree_m
