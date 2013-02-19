module dof_io_m

   use dof_m
   use base_m
   implicit none

   contains

   function rddvrdef(lun,fver) result (dofs)
      implicit none
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
            call stopnow("rddvrdef : unknown basis type")
         unpickle => doftyp%unpickle
         call unpickle(dofs(f)%p, gdim(f), modelabel(f), ipbaspar(:,f), rpbaspar(:,f))
      enddo

      deallocate(xend,rpbaspar,ipbaspar,basis,gdim,modelabel)
      return

  500 call stopnow("rddvrdef : error reading dvr information")

   end function rddvrdef

end module dof_io_m
