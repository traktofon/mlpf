! vim: set ts=3 sw=3 :
module tuckerdecomp

   use base
   use modeutil
   use linear
   implicit none

   integer,parameter :: btyp_unit = 1
   integer,parameter :: btyp_rect = 2

   type :: basis_t
      integer           :: btyp  
      real(dbl),pointer :: b(:,:) => null()
   end type basis_t

   contains

   !--------------------------------------------------------------------
   subroutine compute_basis(v, gdim, m, limit, mdim, basis, ee2)
   ! Computes the 1-dim. basis tensors for a Tucker decomposition of the
   ! tensor v along its m-th dimension.
   !--------------------------------------------------------------------
   ! The number of basis tensors to be returned is specified as follows:
   ! at most mdim basis tensors are returned, but if the accuracy
   ! limit is reached with less basis tensors, mdim is reduced accordingly.
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in)  :: v(:)       ! input tensor
      integer,intent(in)    :: gdim(:)    ! shape of v
      integer,intent(in)    :: m          ! mode for which basis should be computed
      real(dbl),intent(in)  :: limit      ! accuracy limit for keeping the basis tensors
      integer,intent(inout) :: mdim       ! in:maximum/out:actual number of basis tensors
      real(dbl),pointer     :: basis(:,:) ! the computed basis tensors (allocated here)
      real(dbl),intent(out) :: ee2        ! error estimate (squared), should be < limit
      integer               :: vd,gd,nd,lwork,info,nw,i
      real(dbl)             :: lworkr
      real(dbl),allocatable :: dmat(:,:),eval(:),work(:)

      ! Allocate and build this mode's density matrix.
      call vgn_shape(m,gdim,vd,gd,nd)
      allocate(dmat(gd,gd))
      call build_dmat(v,vd,gd,nd,dmat)
      ! Diagonalize the density matrix.
      allocate(eval(gd))
      call dsyev('V', 'U', gd, dmat, gd, eval, lworkr, -1, info) ! workspace query
      lwork = int(lworkr)
      allocate(work(lwork))
      call dsyev('V', 'U', gd, dmat, gd, eval, work, lwork, info)
      if (info /= 0) then
         write (*,*) "ERROR: DSYEV returned info = ",info
         stop 1
      endif
      deallocate(work)
      ! eval contains the natural weights in ascending order.
      ! Determine how many basis tensors to keep.
      call get_basis_size(eval, limit, nw, ee2)
      nw = min(mdim,nw)
      ! Copy the important basis tensors.
      allocate(basis(gd,nw))
      do i=1,nw
         basis(:,i) = dmat(:,gd-i+1)
      enddo
      mdim = nw
      ! Free unneeded memory.
      deallocate(eval)
      deallocate(dmat)
   end subroutine compute_basis



   !--------------------------------------------------------------------
   subroutine contract_core(v, vdim, basis)
   ! Contracts the core tensor for a Tucker decomposition of the tensor v
   ! using the given list of 1-dim. basis tensors.
   !--------------------------------------------------------------------
   ! v and vdim will be overwritten!
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(inout)  :: v(:)     ! input tensor/core tensor
      integer,intent(inout)    :: vdim(:)  ! initial/final shape of v
      type(basis_t),intent(in) :: basis(:) ! 1-dim. basis tensors for all modes
      integer                  :: nmodes,m,vloc,vlen,ulen,vd,gd,nd
      real(dbl),allocatable    :: u(:)     ! temporary tensor
      integer                  :: udim(size(vdim))

      nmodes = size(vdim)
      vloc = 0
      vlen = size(v)
      ! TODO: effort can be saved by performing the projections in an
      !       order depending on the basis sizes.
      do m=1,nmodes
         ! Project the input tensor onto the basis of mode m.
         ! But if the basis is unit, then nothing needs to be done.
         if (basis(m)%btyp == btyp_unit) cycle
         ! To avoid unnecessary copying, we keep track of where the
         ! current tensor is stored, and store the transformed tensor
         ! in an unused array.
         if (vloc==0) then
            ! Source is v, target is u.
            udim = vdim
            udim(m) = size(basis(m)%b,2) ! number of basis columns
            ulen = product(udim)
            allocate(u(ulen))
            call vgn_shape(m, vdim, vd, gd, nd)
            call matrix_tensor_tp(basis(m)%b, v(1:vlen), u(1:ulen), vd, gd, nd, udim(m))
            vloc = 1
         elseif (vloc==1) then
            ! Source is u, target is v.
            vdim = udim
            vdim(m) = size(basis(m)%b,2) ! number of basis columns
            vlen = product(vdim)
            call vgn_shape(m, udim, vd, gd, nd)
            call matrix_tensor_tp(basis(m)%b, u(1:ulen), v(1:vlen), vd, gd, nd, vdim(m))
            deallocate(u)
            vloc = 0
         endif
      enddo
      ! Copy result to v, if necessary
      if (vloc==1) then
         ! Source is u.
         vdim = udim
         v(1:ulen) = u(:)
         deallocate(u)
      endif
   end subroutine contract_core



   !--------------------------------------------------------------------
   subroutine expand_core(v, vdim, basis, u1, u1dim)
   ! Expands the core tensor given by v (shape in vdim) by using the
   ! set of 1-dim. basis tensors.  Resulting tensor is returned in
   ! u1 (allocate here, shape in u1dim).
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in)     :: v(:)
      integer,intent(in)       :: vdim(:)
      type(basis_t),intent(in) :: basis(:)
      real(dbl),allocatable    :: u1(:)
      integer,intent(out)      :: u1dim(size(vdim))
      integer                  :: nmodes,m,vloc,plen,vd,gd,nd
      real(dbl),allocatable    :: u2(:)
      integer                  :: u2dim(size(vdim))

      nmodes = size(vdim)
      vloc = 0
      do m=1,nmodes
         if (basis(m)%btyp == btyp_unit) cycle
         if (vloc==0) then
            ! Source is v, target is u1.
            u1dim = vdim
            u1dim(m) = size(basis(m)%b,1)
            plen = product(u1dim)
            allocate(u1(plen))
            call vgn_shape(m, vdim, vd, gd, nd)
            call matrix_tensor_nt(basis(m)%b, v, u1, vd, gd, nd, u1dim(m))
            vloc = 1
         elseif (vloc==1) then
            ! Source is u1, target is u2.
            u2dim = u1dim
            u2dim(m) = size(basis(m)%b,1)
            plen = product(u2dim)
            allocate(u2(plen))
            call vgn_shape(m, u1dim, vd, gd, nd)
            call matrix_tensor_nt(basis(m)%b, u1, u2, vd, gd, nd, u2dim(m))
            deallocate(u1)
            vloc = 2
         else
            ! Source is u2, target is u1.
            u1dim = u2dim
            u1dim(m) = size(basis(m)%b,1)
            plen = product(u1dim)
            allocate(u1(plen))
            call vgn_shape(m, u2dim, vd, gd, nd)
            call matrix_tensor_nt(basis(m)%b, u2, u1, vd, gd, nd, u1dim(m))
            deallocate(u2)
            vloc = 1
         endif
      enddo
      ! Copy result to u1, if necessary
      if (vloc==0) then
         ! Nothing happened. This shouldn't happen.
         u1dim = vdim
         plen = product(u1dim)
         allocate(u1(plen))
         u1 = v
      elseif (vloc==2) then
         ! Move data from u2 to u1.
         u1dim = u2dim
         plen = product(u1dim)
         allocate(u1(plen))
         u1 = u2
         deallocate(u2)
      endif
   end subroutine expand_core



   !--------------------------------------------------------------------
   subroutine build_dmat(v, vd, gd, nd, dm)
   ! Generate the density matrix along the middle dimension of a
   ! 3-dim. tensor.
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in)  :: v(vd,gd,nd)
      integer,intent(in)    :: vd,gd,nd
      real(dbl),intent(out) :: dm(gd,gd)
      integer               :: a,b,vi,ni
      dm = 0.d0
      do ni=1,nd
         do b=1,gd
            do a=b,gd
               do vi=1,vd
                  dm(a,b) = dm(a,b) + v(vi,a,ni)*v(vi,b,ni)
               enddo
            enddo
         enddo
      enddo
      do b=1,gd
         do a=1,b-1
            dm(a,b) = dm(b,a)
         enddo
      enddo
   end subroutine build_dmat



   !--------------------------------------------------------------------
   subroutine get_basis_size(wghts, limit, bsz, ee2)
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in) :: wghts(:) ! list of weights in ascending order
      real(dbl),intent(in) :: limit    ! parameter for determing how many weights to keep
      integer,intent(out)  :: bsz      ! result: number of weights to keep
      real(dbl),intent(out):: ee2      ! error estimate based on neglected weights
      integer   :: nwghts,i
      real(dbl) :: wsum
      ! Sum up all weights until limit is reached.
      ! The remaining weights are those we want to keep.
      nwghts = size(wghts)
      wsum = 0.d0
      i = 0
      do while (wsum < limit .and. i < nwghts)
         i = i+1
         ee2 = wsum
         wsum = wsum + max(0.d0,wghts(i)) ! ignore negative weights
      enddo
      bsz = nwghts - i + 1
   end subroutine get_basis_size

end module tuckerdecomp
