! vim: set ts=3 sw=3 :
module hiertuck

   use dof
   use tree
   use modeutil
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
   subroutine bla(t,v,vdim,limit)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      real(dbl),intent(inout) :: v(:)
      integer,intent(inout)   :: vdim(:)
      real(dbl),intent(in)    :: limit
      integer                 :: l,m,nc,d1,d2,f
      integer                 :: rank,mdim
      type(node_t),pointer    :: no
      integer                 :: udim(size(vdim))
      integer                 :: todo(size(vdim))
      type(node_tp)           :: nnn(size(vdim))
      type(basis_t)           :: basis(size(vdim))

      ! On entry, v has rank = t%numleaves.
      rank = size(vdim)
      ! Loop over layers from bottom to top.
      do l = t%numlayers-1, 1, -1
         write (*,*) 'LAYER',l
         write (*,*) 'vdim =', vdim(1:rank)
         d1 = 1 ! counting dimensions of v without mode-combination
         d2 = 1 ! counting dimensions of v with mode-combination
         nc = 0 ! counting mode-combinations
         ! Go through all nodes in the current layer, and all leaves above.
         do m = 1, t%numnodes
            no => t%preorder(m)%p
            if (no%layer > l) cycle
            if (no%layer < l .and. .not.no%isleaf) cycle
            if (no%isleaf) then
               udim(d2) = vdim(d1)
               basis(d2)%btyp = btyp_unit
               d1 = d1+1
            else
               allocate(no%ndim(no%nmodes))
               do f=1,no%nmodes
                  no%ndim(f) = no%modes(f)%p%nbasis ! = vdim(c+f-1)
               enddo
               no%plen = product(no%ndim)
               udim(d2) = no%plen
               d1 = d1 + no%nmodes
               nc = nc+1
               todo(nc) = d2
               nnn(nc)%p => no
            endif
            d2 = d2+1
         enddo
         rank = d2-1
         vdim(1:rank) = udim(1:rank)
         write (*,*) 'vdim =', vdim(1:rank)
         write (*,*) 'todo =', todo(1:nc)
         ! TODO: l==1 => just store core tensor
         do m = 1,nc
            d2 = todo(m)
            mdim = vdim(d2)
            call compute_basis(v,vdim(1:rank),d2,limit,mdim,nnn(m)%p%basis)
            nnn(m)%p%nbasis = mdim
            basis(d2)%btyp = btyp_rect
            basis(d2)%b => nnn(m)%p%basis
            udim(m) = mdim
         enddo
         write (*,*) 'udim =', udim(1:rank)
         call compute_core(v,vdim(1:rank),basis)
      enddo
      print *,v(1)
   end subroutine bla

end module hiertuck
