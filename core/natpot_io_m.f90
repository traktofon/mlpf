! vim: set ts=3 sw=3 :
!=======================================================================
module natpot_io_m
!=======================================================================

   use version_m
   use dof_m
   use dof_io_m
   use vtree_m
   implicit none

   contains !===========================================================

   !--------------------------------------------------------------------
   subroutine write_natpot(fname,dofs,tree,modc,core)
   !--------------------------------------------------------------------
      character(len=*),intent(in) :: fname   ! the output filename
      type(dof_tp),intent(in)     :: dofs(:) ! the DVR definitions
      type(vtree_t),intent(in)    :: tree    ! the MLPF tree
      integer,intent(in)          :: modc    ! number of the contracted mode
      real(dbl),intent(in)        :: core(:) ! the contracted C-tensor
      integer                     :: lun,ierr,ndof,nmode,f,m,dimbef,dimaft
      integer*4                   :: idum,lcount
      integer*4,allocatable       :: onedpot(:),pdim(:)
      logical*4                   :: lpconm
      character(len=80)           :: line
      type(vnode_t),pointer       :: no

      open(newunit=lun, file=trim(fname), form="unformatted", status="unknown", iostat=ierr)
      if (ierr /= 0) &
         call stopnow('error opening output file "'//trim(fname)//'"')

      ! --- HEADER ---

      ! file version
      write(lun,err=500) mctdh_compat_versnum
      ! dummy data
      idum = 0
      write(lun,err=500) idum,idum,idum
      ! DVR definition
      call wrdvrdef(lun,dofs)
      ! textual information
      line = "created by mlpf2npot"
      write(lun,err=500) line
      lcount = 1
      line = "TODO: here should be more information"
      write(lun,err=500) lcount
      write(lun,err=500) line

      ! number of contracted mode, and subtracted 1D-potentials
      ndof = size(dofs)
      allocate(onedpot(ndof))
      onedpot = 0
      write(lun,err=500) INT(modc,4), (onedpot(f), f=1,ndof)
      deallocate(onedpot)

      ! are there combined modes?
      nmode = tree%numleaves
!     lpconm = .false.
!     do m=1,nmode
!        no => tree%leaves(m)%p
!        if (no%nmodes > 1) then
!           lpconm = .true.
!           exit
!        endif
!     enddo
      lpconm = ANY([ (tree%leaves(m)%p%nmodes > 1 , m=1,nmode) ])
      write(lun,err=500) lpconm

      ! grid definition
      call wrgrddef(lun,tree)

      ! potdim
      allocate(pdim(nmode))
      do m=1,nmode
         no => tree%leaves(m)%p
         if (m == modc) then
            pdim(m) = no%plen
         else
            pdim(m) = no%nbasis
         endif
      enddo
      print *, "pdim = (",pdim,")" !!! DEBUG
      write(lun,err=500) (pdim(m), m=1,nmode)

      ! --- DATA ---

      do m=1,nmode
         no => tree%leaves(m)%p
         if (m == modc) then
            ! write contracted C-tensor
            dimbef = product(pdim(1:modc-1))
            dimaft = product(pdim(modc+1:nmode))
            call vecout(core,dimbef,no%plen,dimaft)
         else
            ! write natural potentials
            call vecout(no%basis,1,no%plen,no%nbasis)
         endif
      enddo

      ! --- PARAMETERS ---

      ! TODO

      deallocate(pdim)
      close(lun)
      return

  500 call stopnow("error writing natpot header")

      contains

      subroutine vecout(c,ni,nj,nk)
         integer,intent(in)   :: ni,nj,nk
         real(dbl),intent(in) :: c(ni,nj,nk)
         integer              :: i,j,k
         do k=1,nk
            do i=1,ni
               write(lun,err=520) (c(i,j,k), j=1,nj)
            enddo
         enddo
         return
  520    call stopnow("error writing natpot data")
      end subroutine vecout

   end subroutine write_natpot


end module natpot_io_m

