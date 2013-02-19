! vim: set ts=3 sw=3 :
program test_coul4

   use meta_dof_m
   use tree_m
   use graphviz_m
   use genpot_m
   use coul4_m
   use tuckerdecomp_m
   use hiertuck_m
   use linear_m
   implicit none

   integer,parameter         :: ndofs = 6
   real(dbl),parameter       :: acc = 2.27817e-05 ! 5.0 cm^-1 (target RMSE)
   type(dof_tp)              :: dofs(ndofs)
   type(node_tp),allocatable :: nodes(:)
   type(tree_t),pointer      :: t
   integer                   :: f,nmodes,vlen,vlen0
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   real(dbl)                 :: vnorm,vmax,vmin,dnorm,dmax
   real(dbl)                 :: limit,esq,acesq
   real(dbl)                 :: xi,xf
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
   dofs(1)%p => new_dof( "R"  , 24 , 0.1d0, 5.0d0 )
   dofs(2)%p => new_dof( "r1" , 24 , 0.1d0, 5.d0 )
   dofs(3)%p => new_dof( "r2" , 24 , 0.1d0, 5.d0 )
   dofs(4)%p => new_dof( "th1", 24 , 0.0d0, 3.1415926d0 )
   dofs(5)%p => new_dof( "th2", 24 , 0.0d0, 3.1415926d0 )
   dofs(6)%p => new_dof( "phi", 24 , 0.0d0, 6.2831853d0 )

   ! Make leaf nodes.
   nmodes = ndofs
   allocate(nodes(nmodes))
   do f=1,nmodes
      nodes(f)%p => make_leaf( (/ f /) )
      call init_leaf(nodes(f)%p, dofs)
   enddo

   ! Combine modes
   nodes(1)%p => nodes(1)%p
   nodes(2)%p => make_node(nodes(2:3))
   nodes(3)%p => make_node(nodes(4:5))
   nodes(4)%p => nodes(6)%p

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
   call buildpot(coul4,dofs,v,vnorm,vmax,vmin)

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

end program test_coul4
