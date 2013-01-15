! vim: set ts=3 sw=3 :
program test

   use dof
   use tree
   use graphviz
   use genpot
   use testfunc
   use tuckerdecomp
   use hiertuck
   use linear
   implicit none

   integer,parameter         :: ndofs = 6
   integer,parameter         :: ncomb = 2
   integer,parameter         :: gdim = 11
   real(dbl),parameter       :: accuracy = 1.d-6
   type(dof_tp),allocatable  :: dofs(:)  
   type(node_tp),allocatable :: nodes(:) 
   type(tree_t),pointer      :: t
   integer                   :: f,g,nmodes,nleft,ll,m,i,vlen,mdim,nmod1
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v1(:),v0(:)
   type(basis_t),allocatable :: basis(:)
   real(dbl)                 :: vnorm,limit,ee2,error2

   ! Make DOF grids.
   allocate(dofs(ndofs))
   do f=1,ndofs
      allocate(dofs(f)%p)
      dofs(f)%p%gdim = gdim
      write (dofs(f)%p%label, '(a,i0)') '#',f
      allocate(dofs(f)%p%x(gdim))
      do g=1,gdim
         dofs(f)%p%x(g) = dble(g-1)/max(dble(gdim-1),1.d0)
      enddo
   enddo

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
   enddo

   ! Combine modes
   ll = 0
   do while (nmodes > 1)
      nleft  = mod(nmodes,ncomb)
      nmodes = nmodes/ncomb
      if (mod(ll,4)==0) then
         ! extra nodes on right
         do m=1,nmodes
            nodes(m)%p => make_node(nodes(ncomb*(m-1)+1:ncomb*m))
         enddo
         if (nleft>1) then
            nodes(nmodes+1)%p => make_node(nodes(ncomb*nmodes+1:ncomb*nmodes+nleft))
            nmodes = nmodes+1
         elseif (nleft==1) then
            nodes(nmodes+1)%p => nodes(ncomb*nmodes+1)%p
            nmodes = nmodes+1
         endif
      elseif (mod(ll,2)==1) then
         ! extra modes in middle
         m=1
         nmod1 = nmodes/2
         nmod1 = nmod1 + mod(ll/2,2)*mod(nmodes,2)
         do i=1,nmod1
            nodes(m)%p => make_node(nodes( ncomb*(i-1)+1 : ncomb*i ))
            m=m+1
         enddo
         if (nleft>1) then
            nodes(m)%p => make_node(nodes( ncomb*nmod1+1 : ncomb*nmod1+nleft ))
            m=m+1
         elseif (nleft==1) then
            nodes(m)%p => nodes(ncomb*nmod1+1)%p
            m=m+1
         endif
         do i = nmod1+1, nmodes
            nodes(m)%p => make_node(nodes( nleft+ncomb*(i-1)+1 : nleft+ncomb*i ))
            m=m+1
         enddo
         if (nleft/=0) nmodes=nmodes+1
      else
         ! extra nodes on left
         m=1
         if (nleft>1) then
            nodes(m)%p => make_node(nodes(1:nleft))
            m=2
         elseif (nleft==1) then
            ! nodes(1)%p => nodes(1)%p
            m=2
         endif
         do i=1,nmodes
            nodes(m)%p => make_node(nodes(nleft+ncomb*(i-1)+1 : nleft+ncomb*i))
            m=m+1
         enddo
         if (nleft/=0) nmodes=nmodes+1
      endif
      ll = ll+1
   enddo

   ! Build tree.
   t => make_tree(nodes(1)%p)
   
   ! Examine everything.
   write (*,*) 'Tree has:'
   write (*,*) t%numnodes, ' nodes'
   write (*,*) t%numdofs, ' dofs'
   write (*,*) t%numleaves, ' leaves'
   write (*,*) t%numlayers, ' layers'
   write (*,*) 'Levels are:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%preorder(m)%p%layer
   enddo
   write (*,*)
   write (*,*) 'Leaves?'
   do m=1,t%numnodes
      write (*,'(x,l)',advance="no") t%preorder(m)%p%isleaf
   enddo
   write (*,*)
   write (*,*) 'Leaves are:'
   do m=1,t%numleaves
      write (*,'(x,i0)',advance="no") t%leaves(m)%p%num
   enddo
   write (*,*)
   write (*,*) 'Pre-order is:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%preorder(m)%p%num
   enddo
   write (*,*)
   write (*,*) 'Post-order is:'
   do m=1,t%numnodes
      write (*,'(x,i0)',advance="no") t%postorder(m)%p%num
   enddo
   write (*,*)

   ! Generate potential.
   vlen = 1
   do f=1,ndofs
      vlen = vlen*dofs(f)%p%gdim
   enddo
   allocate(v(vlen))
   write (*,*) 'Generating potential, size =',vlen,'...'
   call buildpot(coulomb,dofs,v,vnorm)
   write (*,'(a,g22.15)') '||v|| = ', vnorm
   allocate(v0(vlen))
   v0 = v
   limit = (accuracy*vnorm)**2
   write (*,'(a,es22.15)') 'limit = ', limit

   ! Generate initial Potfit (basis tensors + core tensor)
   call init_leaves(t,dofs)
   nmodes = t%numleaves
   allocate(vdim(nmodes))
   allocate(basis(nmodes))
   do m=1,nmodes
      vdim(m) = t%leaves(m)%p%plen
   enddo
   write (*,*) 'Computing basis tensors...'
   error2 = 0.d0
   do m=1,nmodes
      mdim = vdim(m)
      call compute_basis(v, vdim, m, limit, mdim, t%leaves(m)%p%basis, ee2)
      t%leaves(m)%p%nbasis = mdim
      basis(m)%btyp = btyp_rect
      basis(m)%b => t%leaves(m)%p%basis
      write (*,'(a,i0,a,i0,a,es8.2)') '  mode ',m,' needs ',mdim,' basis tensors, err^2 = ',ee2
      error2 = error2 + ee2
   enddo
   write (*,*) 'Computing core tensor...'
   call contract_core(v,vdim,basis)
   vlen = product(vdim)
   write (*,*)

   ! Do the hierarchical Tucker decomposition.
   write (*,*) 'Generating HT decomposition...'
   call compute_ht(t, v(1:vlen), vdim, limit, ee2)
   error2 = error2 + ee2
   write (*,'(a,es22.15)') 'err^2 = ', error2
   write (*,'(a,es22.15)') 'error = ', sqrt(error2)
   write (*,'(a,es22.15)') 'accu. = ', sqrt(error2)/vnorm

   call mkdot(42,t,dofs)
   call flush(42)
   write (*,*)
   write (*,*) 'Wrote graphviz input file to channel 42.'

   ! Expand it again.
   write (*,*)
   write (*,*) 'Expanding HT representation...'
   call expand_ht(t,v1)
   write (*,*)
   write (*,*) 'Comparing original and expanded tensor:'
   write (*,*)
   call compare(v0,v1)

end program test
