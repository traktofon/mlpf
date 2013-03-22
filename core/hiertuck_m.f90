! vim: set ts=3 sw=3 :
module hiertuck_m

   use logging_m
   use dof_m
   use tree_m
   use tuckerdecomp_m
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine init_leaf(no,dofs)
   !--------------------------------------------------------------------
      implicit none
      type(node_t),intent(inout)  :: no
      type(dof_tp),intent(in)     :: dofs(:)
      integer                     :: ndofs,f,idof
      ndofs = no%nmodes
      allocate(no%ndim(ndofs))
      do f=1,ndofs
         idof = no%dofs(f)
         no%ndim(f) = dofs(idof)%p%gdim
      enddo
      no%plen = product(no%ndim)
   end subroutine init_leaf


   !--------------------------------------------------------------------
   subroutine potfit_from_v(t,v,vdim,limit,acesq)
   ! Generate initial Potfit (basis tensors + core tensor) from an
   ! in-memory tensor v.
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      real(dbl),intent(inout) :: v(:)
      integer,intent(inout)   :: vdim(:)
      real(dbl),intent(in)    :: limit      ! total allowed err^2
      real(dbl),intent(out)   :: acesq      ! accumulated squared error (estimate)
      type(node_t),pointer    :: no
      real(dbl)               :: layerlimit ! allowed err^2 for each node of this layer
      real(dbl)               :: esq
      integer,save            :: logid_data=0
      integer,save            :: logid_progress=0
      character(len=160)      :: msg
      integer                 :: nmodes,m,mdim
      type(basis_t)           :: basis(size(vdim))

      ! Set up logging.
      call get_logger(logid_data,"data")
      call get_logger(logid_progress,"progress")

      ! Initialize error control.
      ! So far, no nodes have been processed, so we set the allowed err^2 for
      ! each node by dividing the total allowed err^2 by the number of nodes
      ! that will have to be processed. Note that the top node is not processed.
      acesq = 0.d0
      layerlimit = limit/(t%numnodes-1) ! divide limit by number of unprocessed nodes
      write (msg,'(a,es22.15)') 'initial potfit: err^2 limit = ', layerlimit
      call write_log(logid_data, LOGLEVEL_INFO, msg)

      ! Set up the shape array.
      nmodes = t%numleaves
      do m=1,nmodes
         vdim(m) = t%leaves(m)%p%plen
      enddo

      ! First compute the basis tensors.
      call write_log(logid_progress, LOGLEVEL_INFO, '  Computing basis tensors...')
      do m=1,nmodes
         no => t%leaves(m)%p
         ! Set the maximum number of basis tensors we want.
         mdim = vdim(m)
         if (no%maxnbasis > 0)  mdim=min(mdim,no%maxnbasis)
         ! Compute basis for this mode, store the basis tensors directly in the node.
         call compute_basis_svd(v, vdim, m, layerlimit, mdim, no%basis, esq)
         no%nbasis = mdim
         ! Add the basis to the list, for later contraction.
         basis(m)%btyp = btyp_rect
         basis(m)%b => no%basis
         ! Keep track of error.
         write (msg,'(a,i0,a,i0,a,es8.2)') '  mode ',m,' needs ',mdim,' basis tensors, err^2 = ',esq
         call write_log(logid_data, LOGLEVEL_INFO, msg)
         acesq = acesq + esq
      enddo

      ! Now compute the core tensor.
      call write_log(logid_progress, LOGLEVEL_INFO, '  Computing core tensor...')
      call contract_core(v,vdim,basis)
      ! v and vdim have been overwritten!

   end subroutine potfit_from_v



   !--------------------------------------------------------------------
   subroutine compute_ht(t,v,vdim,limit,acesq)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      real(dbl),intent(inout) :: v(:)
      integer,intent(inout)   :: vdim(:)
      real(dbl),intent(in)    :: limit
      real(dbl),intent(out)   :: acesq
      integer                 :: l,m,nc,d1,d2,f,i
      integer                 :: order,mdim,vlen
      type(node_t),pointer    :: no
      integer                 :: udim(size(vdim))
      integer                 :: xmode(size(vdim))
      type(node_tp)           :: xnode(size(vdim))
      type(basis_t)           :: basis(size(vdim))
      real(dbl)               :: esq,limitleft,layerlimit
      integer                 :: nodesleft
      logical                 :: lhosvd
      integer,save            :: logid_data=0
      integer,save            :: logid_progress=0
      character(len=160)      :: msg

      ! Set up logging.
      call get_logger(logid_data,"data")
      call get_logger(logid_progress,"progress")

      ! On entry, v has order = t%numleaves (number of dimensions)
      order = size(vdim)

      ! Initialize error tracking.
      acesq = 0.d0                             ! accumulated squared error (estimate)
      nodesleft = t%numnodes - t%numleaves - 1 ! leaves have been processed already, and topnode won't be processed
      limitleft = limit

      ! Loop over layers from bottom to top.
      ! The bottom-most layer is skipped as the initial potfit has
      ! already been computed before.
      do l = t%numlayers-1, 1, -1
         vlen = product(vdim(1:order))

         write (msg,'(a,i0,a)') '  compute_ht: LAYER ',l
         call write_log(logid_progress, LOGLEVEL_INFO, msg)
         write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)
         write (msg,'(a,i0,a,i0)') '  layer ',l,': vlen = ', vlen
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)

         d1 = 1 ! counting dimensions of v, without mode-combination
         d2 = 1 ! counting dimensions of v, with mode-combination
         nc = 0 ! counting mode-combinations
         ! vdim(:) is the shape of v without mode-combination
         ! udim(:) is the shape of v with mode-combination

         ! Go through all nodes in the current layer, and all leaves above.
         do m = 1, t%numnodes
            no => t%preorder(m)%p
            if (no%layer > l) cycle
            if (no%layer < l .and. .not.no%isleaf) cycle
            ! Dimensions of v corresponding to leaf nodes must be
            ! left alone, i.e. the tensor will not be projected along
            ! such dimensions.
            if (no%isleaf) then
               udim(d2) = vdim(d1) 
               basis(d2)%btyp = btyp_unit
               d1 = d1+1
            ! Dimensions of v corresponding to internal nodes will
            ! be mode-combined.
            else
               ! Now we can set information in the node about the basis
               ! sizes of its children.  These sizes had been determined
               ! when doing the previous layer.
               allocate(no%ndim(no%nmodes))
               do f=1,no%nmodes
                  no%ndim(f) = no%modes(f)%p%nbasis ! = vdim(d1+f-1)
               enddo
               no%plen = product(no%ndim)
               ! Store the mode-combined dimension.
               udim(d2) = no%plen
               d1 = d1 + no%nmodes
               ! Remember to compute the basis for this mode/node.
               nc = nc+1
               xmode(nc) = d2
               xnode(nc)%p => no
            endif
            d2 = d2+1
         enddo

         ! Now forget about the old shape of v.
         order = d2-1
         vdim(1:order) = udim(1:order)
         write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)

         ! For the top layer, there can be only one mode left,
         ! i.e. order=1.  Instead of computing a basis for this
         ! mode, we store the left-over core tensor as the "basis".
         if (l==1) then
            no => t%topnode
            no%nbasis = 1
            allocate(no%basis(no%plen,1))
            no%basis(:,1) = v(1:vlen)
            write (msg,'(a,g22.15)') '||v~|| = ', sqrt(sum(v(1:vlen)**2))
            call write_log(logid_data, LOGLEVEL_DEBUG, msg)

         ! Otherwise, compute the basis tensors for the combined modes.
         else
            ! Set the targeted accuracy for this layer.
            if (l==2 .and. t%topnode%nmodes==2) then
               ! normal SVD
               lhosvd = .false.
               layerlimit = limitleft
            else
               ! HOSVD
               lhosvd = .true.
               layerlimit = limitleft/nodesleft
            endif
            write (msg,'(a,i0,a,es22.15)') 'LAYER ',l,': err^2 limit = ', layerlimit
            call write_log(logid_data, LOGLEVEL_INFO, msg)
            do m = 1,nc
               ! Recall the mode number, and the node.
               d2 = xmode(m)
               no => xnode(m)%p
               ! Set maximum allowed basis size.
               mdim = vdim(d2)
               if (no%maxnbasis > 0)  mdim = min(mdim, no%maxnbasis)
               ! Compute the basis and store it in the node.
               call compute_basis_svd(v(1:vlen), vdim(1:order), d2, layerlimit, mdim, no%basis, esq)
               write (msg,'(a,i0,a,i0,a,es8.2)') '  node ',no%num,' needs ',mdim,' basis tensors, err^2 = ',esq
               call write_log(logid_data, LOGLEVEL_INFO, msg)
               no%nbasis = mdim
               ! Add this mode's basis to the list, for later projection.
               basis(d2)%btyp = btyp_rect
               basis(d2)%b => no%basis
               ! Accumulate estimated error^2.
               if (lhosvd .or. m==1) then
                  acesq = acesq + esq
               endif
               nodesleft = nodesleft-1
            enddo
            ! Project the previous core tensor onto the bases.
            call contract_core(v(1:vlen),vdim(1:order),basis)
            ! v and vdim have been overwritten.
            write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
            call write_log(logid_data, LOGLEVEL_DEBUG, msg)
            ! Update the error tracking.
            limitleft = limitleft - acesq
         endif
      enddo
   end subroutine compute_ht


   !--------------------------------------------------------------------
   subroutine expand_ht(t,v)
   !--------------------------------------------------------------------
      implicit none
      type(tree_t),intent(in) :: t
      real(dbl),pointer       :: v(:)
      type(node_t),pointer    :: no
      integer                 :: vdim(t%numdofs)
      integer                 :: xdim(t%numdofs)
      type(basis_t)           :: basis(t%numdofs)
      integer                 :: order,l,m,d1,d2,i,f
      integer,save            :: logid_data=0
      integer,save            :: logid_progress=0
      character(len=160)      :: msg

      ! Set up logging.
      call get_logger(logid_data,"data")
      call get_logger(logid_progress,"progress")

      ! Start with top-level core tensor.
      no => t%topnode
      order = no%nmodes
      vdim(1:order)  = no%ndim(:)
      allocate(v(no%plen))
      v(:) = no%basis(:,1)

      ! Go through remaining layers from top to bottom.
      do l = 2, t%numlayers
         write (msg,'(a,i0,a)') '  expand_ht: LAYER ',l
         call write_log(logid_progress, LOGLEVEL_INFO, msg)
         d1 = 1
         d2 = 1
         ! Go through all nodes in the current layer, and all leaves above.
         do m = 1, t%numnodes
            no => t%preorder(m)%p
            if (no%layer > l) cycle
            ! Don't expand leaves until the end.
            if (l < t%numlayers .and. no%isleaf) then
               basis(d2)%btyp = btyp_unit
               xdim(d1) = vdim(d2)
               d1 = d1+1
               d2 = d2+1
            ! Expand inner nodes of this layer,
            ! as well as leaves at the bottom layer.
            elseif (no%layer == l .or. no%isleaf) then
               basis(d2)%btyp = btyp_rect
               basis(d2)%b => no%basis
               do f=1,no%nmodes
                  xdim(d1) = no%ndim(f)
                  d1 = d1+1
               enddo
               d2 = d2+1
            endif
         enddo
         write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)
         call expand_core(v, vdim(1:order), basis(1:order))
         write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)
         order = d1-1
         vdim(1:order) = xdim(1:order)
         write (msg,'(a,i0,a,99(x,i0))') '  layer ',l,': vdim =', (vdim(i), i=1,order)
         call write_log(logid_data, LOGLEVEL_DEBUG, msg)
      enddo

   end subroutine expand_ht


end module hiertuck_m
