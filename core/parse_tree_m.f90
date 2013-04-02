!=======================================================================
module parse_tree_m
!=======================================================================
!
! leaf :~ modelabel ( "," modelabel )*
! node :~ ( leaf | "(" node node* ")" ) ( "=" INTEGER )?
!
!=======================================================================

   use tokenize_m
   use itree_m
   use base_m
   implicit none
   private

   public :: parse_tree

   integer,parameter :: maxndof = 10

   contains

   !--------------------------------------------------------------------
   function parse_leaf(t) result(leaf)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: t
      type(inode_t),pointer           :: leaf
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
      leaf => make_ileaf(labels(1:nf))
   end function parse_leaf


   !--------------------------------------------------------------------
   recursive function parse_node(t) result(node)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: t
      type(inode_t),pointer           :: node
      character(len=maxtoklen)        :: token
      integer                         :: nm
      type(inode_tp)                  :: ms(maxndof)
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
         node => make_inode(ms(1:nm))
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
      type(inode_t),pointer           :: topnode
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

end module parse_tree_m
