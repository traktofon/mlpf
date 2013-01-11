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
      real(dbl),intent(in),target    :: v(:)
      integer,intent(in)      :: vdim(:)
      real(dbl),intent(in)    :: limit
      integer                 :: l,m,nc,d1,d2,f
      integer                 :: rank,mdim
      type(node_t),pointer    :: no
      real(dbl),allocatable   :: u1(:),u2(:)
      integer                 :: u1dim(size(vdim)),u2dim(size(vdim))
      integer                 :: todo(size(vdim))
      type(node_tp)           :: nnn(size(vdim))
      type(basis_t)           :: basis(size(vdim))

      ! On entry, v has rank = t%numleaves.
      allocate(u1(size(v)))
      u1 = v
      u1dim =  vdim
      rank = size(vdim)
      ! Loop over layers from bottom to top.
      do l = t%numlayers-1, 1, -1
         write (*,*) 'LAYER',l
         write (*,*) 'u1dim =', u1dim(1:rank)
         d1 = 1 ! counting dimensions of u1
         d2 = 1 ! counting dimensions of u2
         nc = 0 ! counting mode-combinations
         ! Go through all nodes in the current layer, and all leaves above.
         do m = 1, t%numnodes
            no => t%preorder(m)%p
            if (no%layer > l) cycle
            if (no%layer < l .and. .not.no%isleaf) cycle
            if (no%isleaf) then
               u2dim(d2) = u1dim(d1)
               basis(d2)%btyp = btyp_unit
               d1 = d1+1
            else
               allocate(no%ndim(no%nmodes))
               do f=1,no%nmodes
                  no%ndim(f) = no%modes(f)%p%nbasis ! = u1dim(c+f-1)
               enddo
               no%plen = product(no%ndim)
               u2dim(d2) = no%plen
               d1 = d1 + no%nmodes
               nc = nc+1
               todo(nc) = d2
               nnn(nc)%p => no
            endif
            d2 = d2+1
         enddo
         rank = d2-1
         write (*,*) 'u2dim =', u2dim(1:rank)
         write (*,*) 'todo  =', todo(1:nc)
         u1dim(1:rank) = u2dim(1:rank)
         ! TODO: l==1 => just store core tensor
         do m = 1,nc
            d2 = todo(m)
            mdim = u2dim(d2)
            call compute_basis(u1,u2dim(1:rank),d2,limit,mdim,nnn(m)%p%basis)
            nnn(m)%p%nbasis = mdim
            basis(d2)%btyp = btyp_rect
            basis(d2)%b => nnn(m)%p%basis
            u1dim(m) = mdim
         enddo
         write (*,*) 'u1dim =', u1dim(1:rank)
         call compute_core(u1,u2dim(1:rank),basis,u2,u1dim(1:rank))
         deallocate(u1)
         allocate(u1(size(u2)))
         u1 = u2
         deallocate(u2)
      enddo
      print *,u1(1)
   end subroutine bla

end module hiertuck
