! vim: set ts=3 sw=3 :
module hiertuck

   use dof
   use tree
   use tuckerdecomp
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine init_leaves(t,dofs)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      type(dof_tp),intent(in) :: dofs(:)
      type(node_t),pointer    :: no
      integer :: l,ndofs,f,idof
      do l=1,t%numleaves
         no => t%leaves(l)%p
         ndofs = no%nmodes
         allocate(no%ndim(ndofs))
         do f=1,ndofs
            idof = no%dofs(f)
            no%ndim(f) = dofs(idof)%p%gdim
         enddo
         no%plen = product(no%ndim)
      enddo
   end subroutine


   !--------------------------------------------------------------------
   subroutine compute_ht(t,v,vdim,limit)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      real(dbl),intent(inout) :: v(:)
      integer,intent(inout)   :: vdim(:)
      real(dbl),intent(in)    :: limit
      integer                 :: l,m,nc,d1,d2,f,i
      integer                 :: rank,mdim,vlen
      type(node_t),pointer    :: no
      integer                 :: udim(size(vdim))
      integer                 :: xmode(size(vdim))
      type(node_tp)           :: xnode(size(vdim))
      type(basis_t)           :: basis(size(vdim))

      ! On entry, v has rank = t%numleaves.
      rank = size(vdim)

      ! Loop over layers from bottom to top.
      ! The bottom-most layer is skipped as the initial potfit has
      ! already been computed before.
      do l = t%numlayers-1, 1, -1
         write (*,'(/,a,i0,a)') '*** LAYER ',l,' ***'
         write (*,'(a,99(x,i0))') 'vdim =', (vdim(i), i=1,rank)
         vlen = product(vdim(1:rank))
         write (*,'(a,i0)') 'vlen = ', vlen

         d1 = 1 ! counting dimensions of v, without mode-combination
         d2 = 1 ! counting dimensions of v, with mode-combination
         nc = 0 ! counting mode-combinations
         ! vdim(:) is the shape of v without mode-combination
         ! udim(:) is the shape of v with mode-combination

         ! Go through all nodes in the current layer, and all leaves above.
         do m = 1, t%numnodes
            no => t%preorder(m)%p
            if (no%layer > l) cycle
            if (no%layer < l .and. .not.no%isleaf) cycle
            ! Dimensions of v corresponding to leaf nodes must be
            ! left alone, i.e. the tensor will not be projected along
            ! such dimensions.
            if (no%isleaf) then
               udim(d2) = vdim(d1) 
               basis(d2)%btyp = btyp_unit
               d1 = d1+1
            ! Dimensions of v corresponding to internal nodes will
            ! be mode-combined.
            else
               ! Now we can set information in the node about the basis
               ! sizes of its children.  These sizes had been determined
               ! when doing the previous layer.
               allocate(no%ndim(no%nmodes))
               do f=1,no%nmodes
                  no%ndim(f) = no%modes(f)%p%nbasis ! = vdim(d1+f-1)
               enddo
               no%plen = product(no%ndim)
               ! Store the mode-combined dimension.
               udim(d2) = no%plen
               d1 = d1 + no%nmodes
               ! Remember to compute the basis for this mode/node.
               nc = nc+1
               xmode(nc) = d2
               xnode(nc)%p => no
            endif
            d2 = d2+1
         enddo

         ! Now forget about the old shape of v.
         rank = d2-1
         vdim(1:rank) = udim(1:rank)
         write (*,'(a,99(x,i0))') 'vdim =', (vdim(i), i=1,rank)

         ! For the top layer, there can be only one mode left,
         ! i.e. rank=1.  Instead of computing a basis for this
         ! mode, we store the left-over core tensor as the "basis".
         if (l==1) then
            no => t%topnode
            no%nbasis = 1
            allocate(no%basis(no%plen,1))
            no%basis(:,1) = v(1:vlen)
            write (*,'(/,a,g22.15)') '||v|| = ', sqrt(sum(v(1:vlen)**2))

         ! Otherwise, compute the basis tensors for the combined modes.
         else
            do m = 1,nc
               ! Recall the mode number, and the node.
               d2 = xmode(m)
               no => xnode(m)%p
               ! Set maximum allowed basis size.
               ! TODO: in future this might be stored in the tree
               mdim = vdim(d2)
               ! Compute the basis and store it in the node.
               call compute_basis(v(1:vlen), vdim(1:rank), d2, limit, mdim, no%basis)
               no%nbasis = mdim
               ! Add this mode's basis to the list, for later projection.
               basis(d2)%btyp = btyp_rect
               basis(d2)%b => no%basis
            enddo
            ! Project the previous core tensor onto the bases.
            call compute_core(v(1:vlen),vdim(1:rank),basis)
            ! v and vdim have been overwritten.
            write (*,'(a,99(x,i0))') 'vdim =', (vdim(i), i=1,rank)
         endif
      enddo
   end subroutine compute_ht

end module hiertuck
