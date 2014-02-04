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
      integer*4                       :: ndof,f,mbaspar,i
      character(len=c1),allocatable   :: modelabel(:)
      integer*4,allocatable           :: gdim(:),basis(:),ipbaspar(:,:)
      integer,allocatable             :: gdim1(:),basis1(:),ipbaspar1(:,:)
      real(dbl),allocatable           :: rpbaspar(:,:), xend(:,:)
      type(doftyp_t),pointer          :: doftyp
      procedure(unpickle_dof),pointer :: unpickle

      read(unit=lun,err=500) ndof
      allocate(modelabel(ndof))
      allocate(gdim(ndof))
      allocate(gdim1(ndof))
      allocate(basis(ndof))
      allocate(basis1(ndof))
      read(unit=lun,err=500) (modelabel(f), f=1,ndof)
      read(unit=lun,err=500) (gdim(f), f=1,ndof)
      read(unit=lun,err=500) mbaspar
      read(unit=lun,err=500) (basis(f), f=1,ndof)
      read(unit=lun,err=500) ! skip ldvr, can be derived from basis
      allocate(ipbaspar(mbaspar,ndof))
      allocate(ipbaspar1(mbaspar,ndof))
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
         basis1(f) = basis(f)
         doftyp => find_doftyp_by_id(basis1(f))
         if (.not.associated(doftyp)) &
            call stopnow("rddvrdef: unknown basis type")
         unpickle => doftyp%unpickle
         gdim1(f) = gdim(f)
         ipbaspar1(:,f) = ipbaspar(:,f)
         call unpickle(dofs(f)%p, gdim1(f), ipbaspar1(:,f), rpbaspar(:,f))
         call dofs(f)%p%set_label(modelabel(f))
      enddo

      deallocate(xend,rpbaspar,ipbaspar1,ipbaspar,basis1,basis,gdim1,gdim,modelabel)
      return

  500 call stopnow("rddvrdef: error reading dvr information")

   end function rddvrdef



   !--------------------------------------------------------------------
   subroutine wrdvrdef(lun,dofs)
   !--------------------------------------------------------------------
      integer,intent(in)              :: lun
      type(dof_tp),intent(in)         :: dofs(:)
      class(dof_t),pointer            :: dof
      integer*4                       :: ndof,f,i
      integer*4,parameter             :: mbaspar = 6
      integer*4,allocatable           :: basis(:), ipbaspar(:,:)
      integer,allocatable             :: basis1(:), ipbaspar1(:,:)
      real(dbl),allocatable           :: rpbaspar(:,:), xend(:,:)
      integer*4,allocatable           :: gdim(:)

      ndof = size(dofs)
      write(unit=lun,err=500) ndof
      write(unit=lun,err=500) (dofs(f)%p%label, f=1,ndof)
      allocate(gdim(ndof))
      do f=1,ndof
         gdim(f) = dofs(f)%p%gdim
      enddo
      write(unit=lun,err=500) (gdim(f), f=1,ndof)
      write(unit=lun,err=500) mbaspar
      allocate(basis(ndof))
      allocate(basis1(ndof))
      allocate(ipbaspar(mbaspar,ndof))
      allocate(ipbaspar1(mbaspar,ndof))
      allocate(rpbaspar(mbaspar,ndof))
      allocate(xend(2,ndof))
      do f=1,ndof
         dof => dofs(f)%p
         call dof%pickle(basis1(f), ipbaspar1(:,f), rpbaspar(:,f))
         xend(1,f) = dof%x(1)
         xend(2,f) = dof%x(dof%gdim)
         basis(f) = basis1(f)
         ipbaspar(:,f) = ipbaspar1(:,f)
      end do
      write(unit=lun,err=500) (basis(f), f=1,ndof)
      write(unit=lun,err=500) (ildvr(basis1(f)), f=1,ndof)
      write(unit=lun,err=500) ((ipbaspar(i,f), i=1,mbaspar), f=1,ndof)
      write(unit=lun,err=500) ((rpbaspar(i,f), i=1,mbaspar), f=1,ndof)
      write(unit=lun,err=500) ((xend(i,f), i=1,2), f=1,ndof)
      deallocate(xend,rpbaspar,ipbaspar1,ipbaspar,basis1,basis)
      return

  500 call stopnow("wrdvrdef: error writing dvr information")

   end subroutine wrdvrdef



   !--------------------------------------------------------------------
   pure function ildvr(id)
   !--------------------------------------------------------------------
      integer,intent(in) :: id
      integer*4          :: ildvr
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
      integer*4                  :: ndof,nmode,m,f
      integer*4,allocatable      :: nmoddof(:),dofnums(:)
      integer,allocatable        :: dofnums1(:)

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
      allocate(dofnums1(maxval(nmoddof)))
      allocate(leaves(nmode))
      do m=1,nmode
         read(lun,err=500) (dofnums(f), f=1,nmoddof(m))
         do f=1,nmoddof(m)
            dofnums1(f) = dofnums(f)
         enddo
         leaves(m)%p => make_vleaf(dofnums1(1:nmoddof(m)))
      enddo
      read(lun,err=500) ! skip dofspf
      read(lun,err=500) ! skip lmult,leb,lmulpack
      ! Clean up.
      deallocate(dofnums1)
      deallocate(dofnums)
      deallocate(nmoddof)
      ! Make a tree out of the leaves.
      topnode => make_vnode(leaves)
      tree => make_vtree(topnode)
      deallocate(leaves)
      return

  500 call stopnow("rdgrddef: error reading grid information")

   end function rdgrddef



   !--------------------------------------------------------------------
   subroutine wrgrddef(lun,tree)
   !--------------------------------------------------------------------
      integer,intent(in)       :: lun
      type(vtree_t),intent(in) :: tree
      integer*4                :: ndof,nmode,ifalse
      integer*4                :: dentype,nstate,npacket,npackts,feb,meb,fpb,mpb
      integer*4,allocatable    :: dofspf(:)
      integer                  :: m,f,i
      type(vnode_t),pointer    :: no

      ! density type = wavefunction
      dentype = 0
      write(lun,err=500) dentype
      ! number of DOFs and primitive modes
      ndof = tree%numdofs
      nmode = tree%numleaves
      write(lun,err=500) ndof,nmode
      ! dummy data about electronic states and packets;
      ! official natpot files seem to have zeroes here
      nstate = 0
      npacket = 0
      npackts = 0
      feb = 0
      meb = 0
      fpb = 0
      mpb = 0
      write(lun,err=500) nstate,npacket,npackts,feb,meb,fpb,mpb
      ! DOF numbers for each mode
      write(lun,err=500) (INT(tree%leaves(m)%p%nmodes,4), m=1,nmode)
      do m=1,nmode
         no => tree%leaves(m)%p
         write(lun,err=500) (INT(no%dofnums(i),4), i=1,no%nmodes)
      enddo
      ! mode numbers for each DOF
      allocate(dofspf(ndof))
      do f=1,ndof
         modeloop: do m=1,nmode
            no => tree%leaves(m)%p
            do i=1,no%nmodes
               if (f == no%dofnums(i)) then
                  dofspf(f) = m
                  cycle modeloop
               endif
            enddo
         enddo modeloop
      enddo
      write(lun,err=500) (dofspf(f), f=1,ndof)
      deallocate(dofspf)
      ! lmult, leb, lmulpack (as integers)
      ifalse = 0
      write(lun,err=500) ifalse,ifalse,ifalse

      return

  500 call stopnow("wrgrddef: error writing grid information")

   end subroutine wrgrddef

end module dof_io_m
