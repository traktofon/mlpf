module parsetree_m

! leaf ::= MODELABEL ( "," MODELABEL )*
! node ::= ( leaf | "(" node node* ")" ) ( "=" INTEGER )?

   use tokenize_m
   use tree_m
   use meta_dof_m
   use dof_m
   use base_m
   implicit none

   type(dof_tp) :: dofs(99)
   integer      :: ndof = 0

   contains


   ! placeholder
   function get_dofnum_by_label(lbl) result (f)
      character(len=*),intent(in) :: lbl
      integer                     :: f
      ndof = ndof+1
      dofs(ndof)%p => new_dof(trim(lbl), 10, 0.d0, 1.d0)
      f = ndof
   end function get_dofnum_by_label


   function parse_leaf(t) result(leaf)
      type(tokenizer_t),intent(inout) :: t
      type(node_t),pointer            :: leaf
      character(len=maxtoklen)        :: token
      integer                         :: f,nf
      integer                         :: fs(16)
      nf = 0
      do while (.not. t%is_stopped())
         ! read current token
         token = t%get()
         ! match with modelabels
         f = get_dofnum_by_label(token)
         if (f==0) call t%error("expected valid modelabel")
         ! store dof number
         nf = nf+1
         fs(nf) = f
         ! continue parsing
         call t%gofwd
         token = t%get()
         if (token /= ",") exit
         call t%gofwd
      enddo
      ! create leaf 
      leaf => make_leaf(fs(1:nf))
   end function parse_leaf


   recursive function parse_node(t) result(node)
      type(tokenizer_t),intent(inout) :: t
      type(node_t),pointer            :: node
      character(len=maxtoklen)        :: token
      integer                         :: nm,maxnb,ierr
      type(node_tp)                   :: ms(16)
      if (t%is_stopped()) &
         call t%error("expected start or end of mode definition")
      token = t%get()
      if (token == "(") then
         nm = 0
         call t%gofwd
         do
            nm = nm+1
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
         token = t%get()
         read (token,*,iostat=ierr) maxnb
         if (ierr/=0) call t%error("expected integer")
         call set_maxnbasis(node,maxnb)
         call t%gofwd
      endif
   end function parse_node


end module parsetree_m
