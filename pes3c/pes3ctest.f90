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
   call open_logfile(ilog,"output")
   call set_logger("data", LOGLEVEL_INFO, ilog)
   call set_logger("tree", LOGLEVEL_INFO, ilog)
   call set_logger("progress", LOGLEVEL_INFO)
   call get_logger(logid_progress, "progress")
   call get_logger(logid_data, "data")

   ! Make DOF grids.
   ! 11/11/11/15/9/9/15/11/11
   dofs(1)%p => new_dof("dr1", 7, 1.4d0, 2.4d0)
   dofs(2)%p => new_dof("dr2", 7, 1.4d0, 2.4d0)
   dofs(3)%p => new_dof("dr4", 7, 4.2d0, 5.5d0)
   dofs(4)%p => new_dof("zpp",10,-0.5d0, 0.5d0)
   dofs(5)%p => new_dof("r3x", 6,-0.8d0, 0.8d0)
   dofs(7)%p => new_dof("r3y", 6,-0.8d0, 0.8d0)
   dofs(8)%p => new_dof("phi",10, 0.0d0, 6.283185d0)
   dofs(6)%p => new_dof("u1" , 7,-0.8d0, 0.35d0)
   dofs(9)%p => new_dof("u2" , 7,-0.35d0,0.8d0)

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
      call init_leaf(nodes(f)%p, dofs)
   enddo

   ! Combine modes
   nodes(1)%p => make_node(nodes(1:2))
   nodes(2)%p => make_node(nodes(3:4))
   nodes(3)%p => make_node(nodes(5:6))
   nodes(4)%p => make_node(nodes(7:9))

   nodes(1)%p => make_node(nodes(1:2))
   nodes(2)%p => make_node(nodes(3:4))

   nodes(1)%p => make_node(nodes(1:2))

   ! Build tree.
   t => make_tree(nodes(1)%p)
   
   ! Check it.
   call examine_tree(t)

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
   call compare(v0,v,dnorm,dmax)
   write (msg,'(a,es22.15)') 'Maximum error = ', dmax
   call write_log(logid_data, LOGLEVEL_INFO, msg)
   write (msg,'(a,es22.15)') 'RMSE = ', dnorm/sqrt(1.d0*vlen0)
   call write_log(logid_data, LOGLEVEL_INFO, msg)

   ! Clean up.
   call close_logfile(ilog)

end program test_pes3c
