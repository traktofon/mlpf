! vim: set ts=3 sw=3 :
program test_pes3c

   use dof
   use tree
   use graphviz
   use genpot
   use zundelmod
   use tuckerdecomp
   use hiertuck
   use linear
   implicit none

   integer,parameter         :: ndofs = 15
   real(dbl),parameter       :: acc = 2.27817e-05 ! 5.0 cm^-1 (target RMSE)
   real(dbl),parameter       :: gfac = 5.0
   type(dof_tp)              :: dofs(ndofs)
   type(node_tp),allocatable :: nodes(:)
   type(tree_t),pointer      :: t
   integer                   :: f,nmodes,m,vlen,vlen0
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   real(dbl)                 :: vnorm,vmax,vmin,dnorm,dmax
   real(dbl)                 :: limit,esq,acesq
   real(dbl)                 :: xi,xf
   integer                   :: gdim
   character(len=16)         :: lbl
   integer                   :: logid_progress = 0
   integer                   :: logid_data = 0
   character(len=160)        :: msg
   integer                   :: idot,ilog

   ! Set up logging
   !call open_logfile(ilog,"output")
   !call set_logger("data", LOGLEVEL_INFO, ilog)
   !call set_logger("tree", LOGLEVEL_INFO, ilog)
   !call set_logger("progress", LOGLEVEL_INFO)
   call get_logger(logid_progress, "progress")
   call get_logger(logid_data, "data")

   ! Make DOF grids.
   f=0

   f = f+1
   lbl  = "z"
   gdim = int(19/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "R"
   gdim = int(20/gfac)
   xi   = 3.9d0
   xf   = 6.0d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "a"
   gdim = int(17/gfac)
   xi   = 0.0d0
   xf   = 6.28318530718*(gdim-1)/gdim
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "x"
   gdim = int(5/gfac)
   xi   = -0.8d0
   xf   =  0.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "y"
   gdim = int(5/gfac)
   xi   = -0.8d0
   xf   =  0.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "la"
   gdim = int(19/gfac)
   xi   = 1.3415926535897d0
   xf   = 4.9415926535897d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "lb"
   gdim = int(19/gfac)
   xi   = -1.8d0
   xf   =  1.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "ua"
   gdim = int(9/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "ub"
   gdim = int(9/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r1a"
   gdim = int(9/gfac)
   xi   =  0.5d0
   xf   =  1.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r2a"
   gdim = int(9/gfac)
   xi   =  2.2d0
   xf   =  3.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "va"
   gdim = int(7/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r1b"
   gdim = int(9/gfac)
   xi   =  0.5d0
   xf   =  1.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "r2b"
   gdim = int(9/gfac)
   xi   =  2.2d0
   xf   =  3.8d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   f = f+1
   lbl  = "vb"
   gdim = int(7/gfac)
   xi   = -0.5d0
   xf   =  0.5d0
   dofs(f)%p => new_dof(lbl,gdim,xi,xf)

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
      call init_leaf(nodes(f)%p, dofs)
   enddo

   ! Combine modes
   ! (z,R)
   nodes(1)%p => make_node(nodes(1:2))
   ! (a,x,y)
   nodes(2)%p => make_node(nodes(3:5))
   ! (la,lb)
   nodes(3)%p => make_node(nodes(6:7))
   ! (ua,ub)
   nodes(4)%p => make_node(nodes(8:9))
   ! (r1a,r2a,va)
   nodes(5)%p => make_node(nodes(10:12))
   ! (r1b,r2b,vb)
   nodes(6)%p => make_node(nodes(13:15))

   ! ( (z,R) , (a,x,y) )
   nodes(1)%p => make_node(nodes(1:2))
   ! ( (la,lb) , (ua,ub) )
   nodes(2)%p => make_node(nodes(3:4))
   ! (r1a,r2a,va)
   nodes(3)%p => nodes(5)%p
   ! (r1b,r2b,vb)
   nodes(4)%p => nodes(6)%p

   ! ( ( (z,R) , (a,x,y) ) , ( (la,lb) , (ua,ub) ) )
   nodes(1)%p => make_node(nodes(1:2))
   ! ( (r1a,r2a,va) , (r1b,r2b,vb) )
   nodes(2)%p => make_node(nodes(3:4))

   ! ( ( ( (z,R) , (a,x,y) ) , ( (la,lb) , (ua,ub) ) ) , ( (r1a,r2a,va) , (r1b,r2b,vb) ) )
   nodes(1)%p => make_node(nodes(1:2))

   ! Build tree.
   t => make_tree(nodes(1)%p)
   
   ! Check it.
   call examine_tree(t)

   ! Generate potential.
   call write_log(logid_progress, LOGLEVEL_INFO, "Generating potential...")
   call h5o2_init(.true., .false., 1.5d0, 1.d0) ! from Oriol's "test"
   call leaf_shape(t,vdim,vlen)
   allocate(v(vlen))
   call buildpot(zundel,dofs,v,vnorm,vmax,vmin)

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
   call compare(v0,v,dnorm,dmax)
   write (msg,'(a,es22.15)') 'Maximum error = ', dmax
   call write_log(logid_data, LOGLEVEL_INFO, msg)
   write (msg,'(a,es22.15)') 'RMSE = ', dnorm/sqrt(1.d0*vlen0)
   call write_log(logid_data, LOGLEVEL_INFO, msg)

   ! Clean up.
   call close_logfile(ilog)

end program test_pes3c
