module parsetree_m

! leaf ::= MODELABEL ( "," MODELABEL )*
! node ::= ( leaf | "(" node node* ")" ) ( "=" INTEGER )?

   use tokenize_m
   use tree_m
   use meta_dof_m
   use dof_m
   use base_m
   implicit none


   contains


   function parse_leaf(t,dofs) result(leaf)
      type(tokenizer_t),intent(inout) :: t
      type(dof_tp),intent(in)         :: dofs(:)
      type(node_t),pointer            :: leaf
      character(len=maxtoklen)        :: token
      integer                         :: f,nf
      integer                         :: fs(16)
      nf = 0
      do while (.not. t%is_stopped())
         ! read current token
         token = t%get()
         ! match with modelabels
         f = find_dofnum_by_label(token,dofs)
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


   recursive function parse_node(t,dofs) result(node)
      type(tokenizer_t),intent(inout) :: t
      type(dof_tp),intent(in)         :: dofs(:)
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
            ms(nm)%p => parse_node(t,dofs)
            token = t%get()
            if (token == ")") then
               call t%gofwd
               exit
            endif
         enddo
         node => make_node(ms(1:nm))
      else
         node => parse_leaf(t,dofs)
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