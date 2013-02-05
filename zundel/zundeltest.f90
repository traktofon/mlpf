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
   real(dbl),parameter       :: acc = 4.55634e-05 ! 10.0 cm^-1 (target RMSE)
   type(dof_tp)              :: dofs(ndofs)
   type(node_tp),allocatable :: nodes(:)
   type(tree_t),pointer      :: t
   integer                   :: f,nmodes,vlen,vlen0
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   real(dbl)                 :: vnorm,vmax,vmin,dnorm,dmax
   real(dbl)                 :: limit,esq,acesq
   integer                   :: logid_progress = 0
   integer                   :: logid_data = 0
   character(len=160)        :: msg
   integer                   :: idot,ilog

   ! Set up logging
   call open_logfile(ilog,"output")
   call set_logger("data", LOGLEVEL_INFO, ilog)
   call set_logger("tree", LOGLEVEL_DEBUG, ilog)
   call set_logger("progress", LOGLEVEL_INFO)
   call get_logger(logid_progress, "progress")
   call get_logger(logid_data, "data")

   ! Make DOF grids.
   ! x,ua,ub - y,a - z,r1b - r2b,vb - R,r1a - r2a,va
   dofs( 1)%p => new_dof("x"  , 3, -0.8d0, 0.8d0)
   dofs( 4)%p => new_dof("y"  , 3, -0.8d0, 0.8d0)
   dofs( 8)%p => new_dof("z"  , 3, -0.5d0, 0.5d0)
   dofs(12)%p => new_dof("R"  , 3,  3.9d0, 6.0d0)
   dofs( 5)%p => new_dof("a"  , 3,  0.0d0, 6.283185d0)
   dofs(13)%p => new_dof("r1a", 3,  0.5d0, 1.8d0)
   dofs(14)%p => new_dof("r2a", 3,  2.2d0, 3.8d0)
   dofs(15)%p => new_dof("va" , 3, -0.5d0, 0.5d0)
   dofs( 2)%p => new_dof("ua" , 3, -0.5d0, 0.5d0)
   dofs( 6)%p => new_dof("la" , 3,  1.3415926535897d0, 4.9415926535897d0)
   dofs( 9)%p => new_dof("r1b", 3,  0.5d0, 1.8d0)
   dofs(10)%p => new_dof("r2b", 3,  2.2d0, 3.8d0)
   dofs(11)%p => new_dof("vb" , 3, -0.5d0, 0.5d0)
   dofs( 3)%p => new_dof("ub" , 3, -0.5d0, 0.5d0)
   dofs( 7)%p => new_dof("lb" , 3, -1.8d0, 1.8d0)

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
      call init_leaf(nodes(f)%p, dofs)
   enddo

   ! Combine modes
   nodes(1)%p => nodes(1)%p
   nodes(2)%p => nodes(2)%p
   nodes(3)%p => nodes(3)%p
   nodes(4)%p => make_node(nodes(4:5))
   nodes(5)%p => make_node(nodes(6:7))
   nodes(6)%p => make_node(nodes(8:9))
   nodes(7)%p => make_node(nodes(10:11))
   nodes(8)%p => make_node(nodes(12:13))
   nodes(9)%p => make_node(nodes(14:15))

   nodes(1)%p => make_node(nodes(1:3))
   nodes(2)%p => make_node(nodes(4:5))
   nodes(3)%p => make_node(nodes(6:7))
   nodes(4)%p => make_node(nodes(8:9))

   nodes(1)%p => make_node(nodes(1:2))
   nodes(2)%p => make_node(nodes(3:4))

   nodes(1)%p => make_node(nodes(1:2))

   ! Build tree.
   t => make_tree(nodes(1)%p)
   
   ! Check it.
   call examine_tree(t)

   ! Generate potential.
   call write_log(logid_progress, LOGLEVEL_INFO, "Generating potential...")
   call h5o2_init(.true., .true., 1.5d0, 1.d0)
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
