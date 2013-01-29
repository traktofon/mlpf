! vim: set ts=3 sw=3 :
program test

   use logging
   use dof
   use tree
   use graphviz
   use genpot
   use tuckerdecomp
   use hiertuck
   use linear
   use testfunc
   implicit none

   integer,parameter         :: ndofs = 12
   integer,parameter         :: ncomb = 2
   integer,parameter         :: gdim1 = 3
   integer,parameter         :: gdim2 = 5
   real(dbl),parameter       :: acc = 1.d-3 ! target RMSE
   type(dof_tp),allocatable  :: dofs(:)
   type(node_tp),allocatable :: nodes(:)
   type(node_t),pointer      :: no
   type(tree_t),pointer      :: t
   integer                   :: f,gdim,g,nmodes,nleft,ll,m,i,vlen,vlen0,mdim,nmod1
   integer,allocatable       :: vdim(:)
   real(dbl),allocatable     :: v(:),v0(:)
   type(basis_t),allocatable :: basis(:)
   real(dbl)                 :: vnorm,vmax,vmin
   real(dbl)                 :: limit,layerlimit,esq,acesq
   integer                   :: logid_progress = 0
   integer                   :: logid_data = 0
   character(len=160)        :: msg
   integer                   :: idot

   ! Set up logging
   call get_logger(logid_progress, "progress")
   call get_logger(logid_data, "data")

   ! Make DOF grids.
   allocate(dofs(ndofs))
   do f=1,ndofs
      gdim = gdim1 + mod(f,gdim2-gdim1+1)
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
      call init_leaf(nodes(f)%p, dofs)
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
   
   ! Check it.
   call examine_tree(t,6)

   ! Generate potential.
   call write_log(logid_progress, LOGLEVEL_INFO, "Generating potential...")
   call leaf_shape(t,vdim,vlen)
   allocate(v(vlen))
   call buildpot(coulombn,dofs,v,vnorm,vmax,vmin)

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

end program test
