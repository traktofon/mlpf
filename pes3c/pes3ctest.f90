! vim: set ts=3 sw=3 :
program test_pes3c

   use dof
   use tree
   use graphviz
   use genpot
   use pes3cmod
   use tuckerdecomp
   use hiertuck
   use linear
   implicit none

   integer,parameter         :: ndofs = 9
   real(dbl),parameter       :: acc = 2.27817e-05 ! 5.0 cm^-1 (target RMSE)
   real(dbl),parameter       :: gfac = 2.0
   type(dof_tp)              :: dofs(ndofs)
   type(node_tp),allocatable :: nodes(:)
   type(node_t),pointer      :: no
   type(tree_t),pointer      :: t
   integer                   :: f,nmodes,m,vlen,vlen0,mdim
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   type(basis_t),allocatable :: basis(:)
   real(dbl)                 :: vnorm,vmax,vmin
   real(dbl)                 :: limit,layerlimit,esq,acesq
   real(dbl)                 :: xi,xf
   integer                   :: gdim
   character(len=16)         :: lbl
   integer                   :: logid_progress = 0
   integer                   :: logid_data = 0
   character(len=160)        :: msg
   integer                   :: idot

   ! Set up logging
   call get_logger(logid_progress, "progress")
   call get_logger(logid_data, "data")

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
      call init_leaf(nodes(f)%p, dofs)
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
   
   ! Check it.
   call examine_tree(t,6)

   ! Generate potential.
   call write_log(logid_progress, LOGLEVEL_INFO, "Generating potential...")
   call leaf_shape(t,vdim,vlen)
   allocate(v(vlen))
   call buildpot(pes3c,dofs,v,vnorm,vmax,vmin)

   ! Keep a copy for later verification.
   allocate(v0(vlen))
   v0 = v
   vlen0 = vlen

   ! Set error limit.
   limit = vlen0 * acc**2
   write (msg,'(a,es22.15)') 'Total err^2 limit = ', limit
   call write_log(logid_data, LOGLEVEL_INFO, msg)

   ! Generate initial Potfit (basis tensors + core tensor).
   call write_log(logid_progress, LOGLEVEL_INFO, 'Generating initial potfit...')
   call potfit_from_v(t, v, vdim, limit, acesq)
   ! v is the new core tensor, size has changed!
   vlen = product(vdim)

   ! Do the hierarchical Tucker decomposition.
   call write_log(logid_progress, LOGLEVEL_INFO, 'Generating HT decomposition...')
   call compute_ht(t, v(1:vlen), vdim, limit-acesq, esq)
   ! v was destroyed
   deallocate(v)

   ! Error information.
   acesq = acesq + esq
   write (msg,'(a,es22.15)') 'Estimated total err^2 = ', acesq
   call write_log(logid_data, LOGLEVEL_INFO, msg)
   write (msg,'(a,es22.15)') 'Estimated RMSE <= ', sqrt(acesq/vlen0)
   call write_log(logid_data, LOGLEVEL_INFO, msg)

   ! Produce graphviz input file.
   call write_log(logid_progress, LOGLEVEL_INFO, 'Writing input file for graphviz...')
   call open_logfile(idot,"tree.dot")
   call mkdot(idot,t,dofs)
   call flush(idot)
   call close_logfile(idot)

   ! Expand the HT decomposition again.
   call write_log(logid_progress, LOGLEVEL_INFO, 'Expanding HT representation...')
   call expand_ht(t,v)

   ! And compare.
   call write_log(logid_progress, LOGLEVEL_INFO, 'Comparing original and expanded tensor...')
   call compare(v0,v)

end program test_pes3c
