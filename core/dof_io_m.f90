module dof_io_m

   use dof_m
   use base_m
   use vtree_m
   implicit none

   contains

   !--------------------------------------------------------------------
   function rddvrdef(lun,fver) result (dofs)
   !--------------------------------------------------------------------
      type(dof_tp),pointer            :: dofs(:)
      integer,intent(in)              :: lun
      real(dbl),intent(in)            :: fver
      integer                         :: ndof,f,mbaspar,i
      character(len=c1),allocatable   :: modelabel(:)
      integer,allocatable             :: gdim(:),basis(:), ipbaspar(:,:)
      real(dbl),allocatable           :: rpbaspar(:,:), xend(:,:)
      type(doftyp_t),pointer          :: doftyp
      procedure(unpickle_dof),pointer :: unpickle

      read(unit=lun,err=500) ndof
      allocate(modelabel(ndof))
      allocate(gdim(ndof))
      allocate(basis(ndof))
      read(unit=lun,err=500) (modelabel(f), f=1,ndof)
      read(unit=lun,err=500) (gdim(f), f=1,ndof)
      read(unit=lun,err=500) mbaspar
      read(unit=lun,err=500) (basis(f), f=1,ndof)
      read(unit=lun,err=500) ! skip ldvr, can be derived from basis
      allocate(ipbaspar(mbaspar,ndof))
      allocate(rpbaspar(mbaspar,ndof))
      allocate(xend(2,ndof))
      read(unit=lun,err=500) ((ipbaspar(i,f), i=1,mbaspar), f=1,ndof)
      read(unit=lun,err=500) ((rpbaspar(i,f), i=1,mbaspar), f=1,ndof)
      read(unit=lun,err=500) ((xend(i,f), i=1,2), f=1,ndof)
      ! skip obsolete data
      if (fver >= 8.3016d0 .and. fver < 8.40d0) then
         read(unit=lun,err=500)
         read(unit=lun,err=500)
         read(unit=lun,err=500)
      endif

      allocate(dofs(ndof))
      do f=1,ndof
         doftyp => find_doftyp_by_id(basis(f))
         if (.not.associated(doftyp)) &
            call stopnow("rddvrdef: unknown basis type")
         unpickle => doftyp%unpickle
         call unpickle(dofs(f)%p, gdim(f), ipbaspar(:,f), rpbaspar(:,f))
         call dofs(f)%p%set_label(modelabel(f))
      enddo

      deallocate(xend,rpbaspar,ipbaspar,basis,gdim,modelabel)
      return

  500 call stopnow("rddvrdef: error reading dvr information")

   end function rddvrdef



   !--------------------------------------------------------------------
   subroutine wrdvrdef(lun,dofs)
   !--------------------------------------------------------------------
      integer,intent(in)              :: lun
      type(dof_tp),intent(in)         :: dofs(:)
      class(dof_t),pointer            :: dof
      integer                         :: ndof,f,i
      integer,parameter               :: mbaspar = 6
      integer,allocatable             :: basis(:), ipbaspar(:,:)
      real(dbl),allocatable           :: rpbaspar(:,:), xend(:,:)

      ndof = size(dofs)
      write(unit=lun,err=500) ndof
      write(unit=lun,err=500) (dofs(f)%p%label, f=1,ndof)
      write(unit=lun,err=500) (dofs(f)%p%gdim, f=1,ndof)
      write(unit=lun,err=500) mbaspar
      allocate(basis(ndof))
      allocate(ipbaspar(mbaspar,ndof))
      allocate(rpbaspar(mbaspar,ndof))
      allocate(xend(2,ndof))
      do f=1,ndof
         dof => dofs(f)%p
         call dof%pickle(basis(f), ipbaspar(:,f), rpbaspar(:,f))
         xend(1,f) = dof%x(1)
         xend(2,f) = dof%x(dof%gdim)
      end do
      write(unit=lun,err=500) (basis(f), f=1,ndof)
      write(unit=lun,err=500) (ildvr(basis(f)), f=1,ndof)
      write(unit=lun,err=500) ((ipbaspar(i,f), i=1,mbaspar), f=1,ndof)
      write(unit=lun,err=500) ((rpbaspar(i,f), i=1,mbaspar), f=1,ndof)
      write(unit=lun,err=500) ((xend(i,f), i=1,2), f=1,ndof)
      deallocate(xend,rpbaspar,ipbaspar,basis)
      return

  500 call stopnow("wrdvrdef: error writing dvr information")

   end subroutine wrdvrdef



   !--------------------------------------------------------------------
   pure function ildvr(id)
   !--------------------------------------------------------------------
      integer,intent(in) :: id
      integer            :: ildvr
      select case (id)
         case (-1,0,4,6)
            ildvr = 0
         case default
            ildvr = 1
      end select
   end function ildvr



   !--------------------------------------------------------------------
   function rdgrddef(lun) result(tree)
   !--------------------------------------------------------------------
      integer,intent(in)         :: lun
      type(vtree_t),pointer      :: tree
      type(vnode_tp),allocatable :: leaves(:)
      type(vnode_t),pointer      :: topnode
      integer                    :: ndof,nmode,m,f
      integer,allocatable        :: nmoddof(:),dofnums(:)

      nullify(tree)
      read(lun,err=500) ! skip dentype
      ! Number of DOFs and number of primitive modes (i.e. leaves)
      read(lun,err=500) ndof, nmode
      read(lun,err=500) ! skip nstate,npacket,npackts,feb,meb,fpb,mpb
      ! Number of DOFs inside each mode.
      allocate(nmoddof(nmode))
      read(lun,err=500) (nmoddof(m), m=1,nmode)
      ! The DOF-numbers for each mode.
      allocate(dofnums(maxval(nmoddof)))
      allocate(leaves(nmode))
      do m=1,nmode
         read(lun,err=500) (dofnums(f), f=1,nmoddof(m))
         leaves(m)%p => make_vleaf(dofnums(1:nmoddof(m)))
      enddo
      read(lun,err=500) ! skip dofspf
      read(lun,err=500) ! skip lmult,leb,lmulpack
      ! Clean up.
      deallocate(dofnums)
      deallocate(nmoddof)
      ! Make a tree out of the leaves.
      topnode => make_vnode(leaves)
      tree => make_vtree(topnode)
      deallocate(leaves)
      return

  500 call stopnow("rddvrdef: error reading dvr information")

   end function rdgrddef

end module dof_io_m
