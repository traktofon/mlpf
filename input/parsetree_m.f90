module parsetree_m

! leaf ::= MODELABEL [ "," MODELABEL ]*
! node ::= { leaf | "(" node [ node ]* ")" } [ "=" INTEGER ]

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


   function parse_leaf() result(leaf)
      type(node_t),pointer     :: leaf
      character(len=maxtoklen) :: token
      integer                  :: f,nf
      integer                  :: fs(16)
      logical                  :: lend
      nf = 0
      lend = .false.
      do while (.not.lend)
         ! read current token
         call get_token(token)
         ! match with modelabels
         f = get_dofnum_by_label(token)
         if (f==0) call parserror("expected valid modelabel")
         ! store dof number
         nf = nf+1
         fs(nf) = f
         ! continue parsing
         call next_token(lend)
         call get_token(token)
         if (token /= ",") exit
         call next_token(lend)
      enddo
      ! create leaf 
      leaf => make_leaf(fs(1:nf))
   end function parse_leaf


   recursive function parse_node() result(node)
      type(node_t),pointer     :: node
      character(len=maxtoklen) :: token
      integer                  :: nm,maxnb,ierr
      type(node_tp)            :: ms(16)
      logical                  :: lend
      call get_token(token)
      if (token == "(") then
         nm = 0
         call next_token(lend)
         do
            nm = nm+1
            ms(nm)%p => parse_node()
            call get_token(token)
            if (token == ")") then
               call next_token(lend)
               exit
            endif
         enddo
         node => make_node(ms(1:nm))
      else
         node => parse_leaf()
      endif
      call get_token(token)
      if (token == "=") then
         call next_token(lend)
         call get_token(token)
         read (token,*,iostat=ierr) maxnb
         if (ierr/=0) call parserror("expected integer")
         call set_maxnbasis(node,maxnb)
         call next_token(lend)
      endif
   end function parse_node


   subroutine parserror(msg)
      character(len=*),intent(in) :: msg
      character(len=maxtoklen)    :: token
      call get_token(token)
      write (*,'(4a)') 'Parse Error at "',trim(token),'": ',trim(msg)
      stop 1
   end subroutine parserror


end module parsetree_m
