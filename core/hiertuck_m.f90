! vim: set ts=3 sw=3 :
module hiertuck_m

   use logging_m
   use dof_m
   use dof_io_m
   use vtree_m
   use graphviz_m
   use tuckerdecomp_m
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine init_vleaf(no,dofs,lcheck)
   ! A leaf contains one or more DOFs.
   !--------------------------------------------------------------------
      implicit none
      type(vnode_t),intent(inout) :: no
      type(dof_tp),intent(in)     :: dofs(:)
      logical,optional            :: lcheck
      integer                     :: ndofs,f,idof
      character(len=80)           :: msg
      logical                     :: lc
      lc = .true.
      if (present(lcheck)) lc = lcheck
      ndofs = no%nmodes
      ! Check that the DOFs are actually consecutive.
      if (lc) then
         do f=2,ndofs
            if (no%dofnums(f) - no%dofnums(f-1) /= 1) then
               write(msg,'(a,i0,a)') &
                  'DOFs in node ',no%num,' are not consecutive!'
               call stopnow(msg)
            endif
         enddo
      endif
      ! Copy the DOFs' dimensions into the node.
      allocate(no%ndim(ndofs))
      do f=1,ndofs
         idof = no%dofnums(f)
         no%ndim(f) = dofs(idof)%p%gdim
      enddo
      ! Calculate size of product grid.
      no%plen = product(no%ndim)
   end subroutine init_vleaf


   !--------------------------------------------------------------------
   subroutine potfit_from_v(t,v,vdim,limit,acesq)
   ! Generate initial Potfit (basis tensors + core tensor) from an
   ! in-memory tensor v.
   !--------------------------------------------------------------------
      implicit none
      type(vtree_t),intent(in) :: t
      real(dbl),intent(inout)  :: v(:)
      integer,intent(inout)    :: vdim(:)
      real(dbl),intent(in)     :: limit      ! total allowed err^2
      real(dbl),intent(out)    :: acesq      ! accumulated squared error (estimate)
      type(vnode_t),pointer    :: no
      real(dbl)                :: layerlimit ! allowed err^2 for each node of this layer
      real(dbl)                :: esq
      integer,save             :: logid_data=0
      integer,save             :: logid_progress=0
      character(len=160)       :: msg
      integer                  :: nmodes,m,mdim
      type(basis_t)            :: basis(size(vdim))

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
         ! Compute basis for this mode, store the basis weights/tensors directly in the node.
         call compute_basis_svd(v, vdim, m, layerlimit, mdim, no%wghts, no%basis, esq)
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
   subroutine potfit_from_npot(npfile,dofs,t,v,vdim)
   !--------------------------------------------------------------------
      character(len=c5),intent(in) :: npfile
      type(dof_tp),pointer         :: dofs(:)
      type(vtree_t),intent(in)     :: t
      real(dbl),pointer            :: v(:)
      integer,intent(inout)        :: vdim(:)
      type(dof_tp),pointer         :: npdofs(:)
      type(vtree_t),pointer        :: nptree
      integer                      :: modc,modcdim
      real(dbl),pointer            :: dtens(:)
      integer,pointer              :: dtensdim(:)
      real(dbl)                    :: limit,ee2
      type(vnode_t),pointer        :: cnode
      integer                      :: idot

      ! Read the natpot file:
      !   DVR def. + Mode def. + Contr.mode + D-tensor
      call load_natpot(npfile,npdofs,nptree,modc,dtens,dtensdim)

      ! To get a full Tucker decomposition, we need the basis tensors
      ! for the contracted mode.
      cnode => nptree%leaves(modc)%p
      modcdim = 0
      limit = 0.d0
      call compute_basis_svd(dtens,dtensdim,modc,limit,modcdim,cnode%wghts,cnode%basis,ee2)

      ! TODO: C-tensor

      ! DEBUG: dump natpot tree
      open(newunit=idot,file="natpot.dot",form="formatted",status="unknown")
      call mkdot(idot,nptree,npdofs)
      close(idot)
      stop 1

   end subroutine potfit_from_npot


   !--------------------------------------------------------------------
   subroutine compute_ht(t,v,vdim,limit,acesq)
   !--------------------------------------------------------------------
      implicit none
      type(vtree_t),intent(in) :: t
      real(dbl),intent(inout)  :: v(:)
      integer,intent(inout)    :: vdim(:)
      real(dbl),intent(in)     :: limit
      real(dbl),intent(out)    :: acesq
      integer                  :: l,m,nc,d1,d2,f,i
      integer                  :: order,mdim,vlen
      type(vnode_t),pointer    :: no
      integer                  :: udim(size(vdim))
      integer                  :: xmode(size(vdim))
      type(vnode_tp)           :: xnode(size(vdim))
      type(basis_t)            :: basis(size(vdim))
      real(dbl)                :: esq,limitleft,layerlimit
      integer                  :: nodesleft
      logical                  :: lhosvd
      integer,save             :: logid_data=0
      integer,save             :: logid_progress=0
      character(len=160)       :: msg

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
            allocate(no%wghts(1))
            allocate(no%basis(no%plen,1))
            no%wghts(1) = 1.d0
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
               call compute_basis_svd(v(1:vlen), vdim(1:order), d2, layerlimit, mdim, no%wghts, no%basis, esq)
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
      type(vtree_t),intent(in) :: t
      real(dbl),pointer        :: v(:)
      type(vnode_t),pointer    :: no
      integer                  :: vdim(t%numdofs)
      integer                  :: xdim(t%numdofs)
      type(basis_t)            :: basis(t%numdofs)
      integer                  :: order,l,m,d1,d2,i,f
      integer,save             :: logid_data=0
      integer,save             :: logid_progress=0
      character(len=160)       :: msg

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


   !--------------------------------------------------------------------
   subroutine load_natpot(fname,dofs,nptree,modc,dtens,dtensdim)
   !--------------------------------------------------------------------
      character(len=c5),intent(in) :: fname
      type(dof_tp),pointer         :: dofs(:)
      type(vtree_t),pointer        :: nptree
      integer,intent(out)          :: modc
      real(dbl),pointer            :: dtens(:)
      integer,pointer              :: dtensdim(:)
      integer                      :: lun,lcount,i
      real(dbl)                    :: fver
      integer                      :: ndof,f,nmode,m,g
      integer                      :: dimbef,dimaft,dimmodc
      integer,allocatable          :: onedpot(:),potdim(:)
      type(vnode_t),pointer        :: no

      nullify(nptree)
      open(newunit=lun, file=trim(fname), status="old", form="unformatted", err=510)
      ! File version.
      read(unit=lun,err=500) fver
      ! Skip obsolete information.
      read(unit=lun,err=500)
      ! DVR information.
      dofs => rddvrdef(lun,fver)
      ! Skip textual information.
      read(unit=lun,err=500)
      read(unit=lun,err=500) lcount
      do i=1,lcount
         read(unit=lun,err=500)
      enddo

      ! Number of contracted mode, and subtracted 1D-potentials.
      ndof = size(dofs)
      allocate(onedpot(ndof))
      read(unit=lun,err=500) modc, (onedpot(f), f=1,ndof)
      read(unit=lun,err=500) ! skip lpconm
      ! Leaf structure (DOFs of each primitive mode)
      nptree => rdgrddef(lun)
      nmode = nptree%numleaves
      do m=1,nmode
         no => nptree%leaves(m)%p
         call init_vleaf(no,dofs,lcheck=.false.)
      enddo

      ! Number of natural potentials for each mode.
      allocate(potdim(nmode))
      read(unit=lun,err=500) (potdim(m), m=1,nmode)
      ! Ignore the subtracted 1D-potentials.
      do f=1,ndof
         if (onedpot(f) > 0) read(unit=lun,err=500)
      enddo

      ! Read natural potentials and D-tensor.
      dimbef = product(potdim(1:modc-1))
      dimaft = product(potdim(modc+1:nmode))
      dimmodc = nptree%leaves(modc)%p%plen
      allocate(dtensdim(nmode))
      allocate(dtens(dimbef*dimmodc*dimaft))
      do m=1,nmode
         no => nptree%leaves(m)%p
         if (m==modc) then ! contracted mode -> D-tensor
            dtensdim(m) = dimmodc
            no%nbasis = dimmodc
            call rddtens(dtens,dimbef,dimmodc,dimaft)
         else ! other mode
            dtensdim(m) = potdim(m)
            no%nbasis = potdim(m)
            allocate(no%basis(no%plen,no%nbasis))
            do i=1,no%nbasis
               read(lun,err=500) (no%basis(g,i), g=1,no%plen)
            enddo
         endif
      enddo

      ! Clean up.
      deallocate(potdim)
      deallocate(onedpot)
      close(lun)
      return

 500  call stopnow("error reading file: "//trim(fname))
 510  call stopnow("cannot open file: "//trim(fname))

      contains

      subroutine rddtens(v,vdim,gdim,ndim)
         integer   :: vdim,gdim,ndim
         real(dbl) :: v(vdim,gdim,ndim)
         integer   :: i,j,g
         do i=1,ndim
            do j=1,vdim
               read(unit=lun,err=600) (v(j,g,i), g=1,gdim)
            enddo
         enddo
         return
 600     call stopnow("error reading D-tensor from file: "//trim(fname))
      end subroutine rddtens

   end subroutine load_natpot

end module hiertuck_m
