! vim: set ts=3 sw=3 :
program test_pes3c

   use dof
   use tree
   use graphviz
   use genpot
   use testfunc
   use tuckerdecomp
   use hiertuck
   use linear
   implicit none

   integer,parameter         :: ndofs = 9
   real(dbl),parameter       :: accuracy = 1.d-3
   real(dbl),parameter       :: gfac = 1.5
   type(dof_tp)              :: dofs(ndofs)
   type(node_tp),allocatable :: nodes(:)
   type(tree_t),pointer      :: t
   integer                   :: f,nmodes,m,vlen,mdim
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   type(basis_t),allocatable :: basis(:)
   real(dbl)                 :: vnorm,vmax,vmin
   real(dbl)                 :: limit,ee2,error2
   real(dbl)                 :: xi,xf
   integer                   :: gdim
   character(len=16)         :: lbl

   ! Make DOF grids.
   f=0

   f = f+1
   lbl  = "dr1"
   gdim = int(11/gfac)
   xi   = 1.4d0
   xf   = 2.4d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "dr2"
   gdim = int(11/gfac)
   xi   = 1.4d0
   xf   = 2.4d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "dr4"
   gdim = int(11/gfac)
   xi   = 4.2d0
   xf   = 5.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "zpp"
   gdim = int(15/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r3x"
   gdim = int(9/gfac)
   xi   = -0.8d0
   xf   =  0.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r3y"
   gdim = int(9/gfac)
   xi   = -0.8d0
   xf   =  0.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "phi"
   gdim = int(15/gfac)
   xi   = 0.0d0
   xf   = 6.28318530718*(gdim-1)/gdim
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "u1"
   gdim = int(11/gfac)
   xi   = -0.8d0
   xf   =  0.35d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "u2"
   gdim = int(11/gfac)
   xi   = -0.35d0
   xf   =  0.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
   enddo

   ! Combine modes
   m = 0
   ! (dr1,dr2)
   m = m+1
   nodes(m)%p => make_node(nodes(1:2))
   ! (dr4,zpp)
   m = m+1
   nodes(m)%p => make_node(nodes(3:4))
   ! (r3x,r3y,phi)
   m = m+1
   nodes(m)%p => make_node(nodes(5:7))
   ! (u1,u2)
   m = m+1
   nodes(m)%p => make_node(nodes(8:9))

   nmodes = m
   m = 0
   ! ( (dr1,dr2) , (dr4,zpp) )
   m = m+1
   nodes(m)%p => make_node(nodes(1:2))
   ! ( (r3x,r3y,phi) , (u1,u2) )
   m = m+1
   nodes(m)%p => make_node(nodes(3:4))

   nmodes = m
   m = 0
   ! ( ( (dr1,dr2) , (dr4,zpp) ) , ( (r3x,r3y,phi) , (u1,u2) ) )
   m = m+1
   nodes(m)%p => make_node(nodes(1:2))

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
   call buildpot(pes3c,dofs,v,vnorm,vmax,vmin)
   write (*,'(a,g22.15)') '||v|| = ', vnorm
   write (*,'(a,g22.15)') 'v_max = ', vmax
   write (*,'(a,g22.15)') 'v_min = ', vmin
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
      call compute_basis_svd(v, vdim, m, limit, mdim, t%leaves(m)%p%basis, ee2)
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
   ! v was destroyed
   deallocate(v)

   call mkdot(42,t,dofs)
   call flush(42)
   write (*,*)
   write (*,*) 'Wrote graphviz input file to channel 42.'

   ! Expand it again.
   write (*,*)
   write (*,*) 'Expanding HT representation...'
   call expand_ht(t,v)
   write (*,*)
   write (*,*) 'Comparing original and expanded tensor:'
   write (*,*)
   call compare(v0,v)

end program test_pes3c
