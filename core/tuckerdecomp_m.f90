! vim: set ts=3 sw=3 :
module tuckerdecomp_m

   use base_m
   use modeutil_m
   use linear_m
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
   ! This routine builds a potential density matrix and diagonalizes it.
   ! THIS IS NOT VERY ACCURATE. Eigenvalues that are small compared to
   ! the maximum eigenvalue will be noisy and might become negative.
   ! (I tried to improve the accuracy of the generated density matrix by
   ! using Kahan summation, but it didn't help.)
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
      ! eval contains the natural weights in *ascending* order.
      ! Determine how many basis tensors to keep.
      call get_basis_size(eval,nw)
      nw  = min(mdim,nw)
      ee2 = basis_error_sq(eval,nw)
      ! Copy the important basis tensors.
      allocate(basis(gd,nw))
      do i=1,nw
         basis(:,i) = dmat(:,gd-i+1)
      enddo
      mdim = nw
      ! Free unneeded memory.
      deallocate(eval)
      deallocate(dmat)

   contains

      subroutine get_basis_size(wghts, bsz)
         implicit none
         real(dbl),intent(in) :: wghts(:) ! list of weights in ascending order
         integer,intent(out)  :: bsz      ! result: number of weights to keep
         integer              :: nwghts,i
         real(dbl)            :: wsum
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

      pure function basis_error_sq(wghts, bsz) result (esq)
         implicit none
         real(dbl),intent(in) :: wghts(:)
         integer,intent(in)   :: bsz
         real(dbl)            :: esq
         integer              :: i
         esq = 0.d0
         do i = 1, size(wghts)-bsz
            esq = esq + max(0.d0,wghts(i))
         enddo
      end function basis_error_sq

   end subroutine compute_basis



   !--------------------------------------------------------------------
   subroutine compute_basis_svd(v, gdim, m, limit, mdim, basis, ee2)
   !--------------------------------------------------------------------
   ! Computes the 1-dim. basis tensors for a Tucker decomposition of the
   ! tensor v along its m-th dimension.
   !--------------------------------------------------------------------
   ! This routine builds a matrizication of the tensor (i.e. reshaping
   ! the tensor into [ rest x gdim(m) ]) and computes an SVD to obtain
   ! the basis from the right singular vectors.  While mathematically
   ! equivalent to the density matrix approach, the SVD approach is far
   ! more accurate and can be faster (for large tensors, and when using
   ! an optimized LAPACK).
   ! Alternatively, one might build the transposed matricization and use
   ! the left singular vectors. This works too, but proved to be slower
   ! (~2x) in my tests (using ACML).
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
      real(dbl),allocatable :: vmat(:,:)  ! matricization of v
      real(dbl),allocatable :: sval(:),work(:)
      integer               :: vd,gd,nd,lwork,info,nw,i
      real(dbl)             :: rdum,lworkr
      ! Allocate and build matricization of v.
      call vgn_shape(m,gdim,vd,gd,nd)
      allocate(vmat(vd*nd,gd))
      call build_vmat(v,vmat)
      ! Compute SVD and left singular vectors.
      allocate(sval(min(gd,vd*nd)))
      call dgesvd('N', 'O', vd*nd, gd, vmat, vd*nd, sval, &
                  rdum, 1, rdum, 1, lworkr, -1, info) ! workspace query
      lwork = int(lworkr)
      allocate(work(lwork))
      call dgesvd('N', 'O', vd*nd, gd, vmat, vd*nd, sval, &
                  rdum, 1, rdum, 1, work, lwork, info)
      if (info /= 0) then
         write (*,*) "ERROR: DGESVD returned info = ",info
         stop 1
      endif
      deallocate(work)
      ! sval contains the singular values in *descending* order.
      ! Determine how many basis tensors to keep.
      call get_basis_size(sval,nw)
      nw  = min(mdim,nw)
      ee2 = basis_error_sq(sval,nw)
      ! Copy the important basis tensors.
      allocate(basis(gd,nw))
      do i=1,nw
         basis(:,i) = vmat(i,:)
      enddo
      mdim = nw
      ! Free unneeded memory.
      deallocate(sval)
      deallocate(vmat)

   contains
   
      subroutine build_vmat(x,xmat)
         real(dbl),intent(in)  :: x(vd,gd,nd)
         real(dbl),intent(out) :: xmat(vd,nd,gd)
         integer               :: i,j,k
         if (nd==1) then
            xmat(:,1,:) = x(:,:,1)
         else
            do k=1,nd
               do j=1,gd
                  do i=1,vd
                     xmat(i,k,j) = x(i,j,k)
                  enddo
               enddo
            enddo
         endif
      end subroutine build_vmat

      subroutine get_basis_size(sval, bsz)
         implicit none
         real(dbl),intent(in) :: sval(:) ! list of singular values in descending order
         integer,intent(out)  :: bsz     ! number of weights to keep
         integer              :: i
         real(dbl)            :: wsum
         ! Sum up all weights until limit is reached.
         ! The remaining weights are those we want to keep.
         wsum = 0.d0
         i = size(sval)
         do while (wsum < limit .and. i > 0)
            wsum = wsum + max(0.d0,sval(i))**2 ! ignore negative singular values
            i = i-1
         enddo
         bsz = i+1
      end subroutine get_basis_size

      pure function basis_error_sq(sval, bsz) result (esq)
         implicit none
         real(dbl),intent(in) :: sval(:)
         integer,intent(in)   :: bsz
         real(dbl)            :: esq
         integer              :: i
         esq = 0.d0
         do i = size(sval), bsz+1, -1
            esq = esq + max(0.d0,sval(i))**2
         enddo
      end function basis_error_sq

   end subroutine compute_basis_svd



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
   subroutine expand_core(v, vdim, basis)
   ! Expands the core tensor given by v (shape in vdim) by using the
   ! set of 1-dim. basis tensors.  
   !--------------------------------------------------------------------
   ! v will be reallocated, and vdim will be overwritten!
   !--------------------------------------------------------------------
      implicit none
      real(dbl),pointer        :: v(:)
      integer,intent(inout)    :: vdim(:)
      type(basis_t),intent(in) :: basis(:)
      real(dbl),allocatable    :: u(:)
      integer                  :: udim(size(vdim))
      integer                  :: nmodes,m,vloc,plen,vd,gd,nd

      nmodes = size(vdim)
      vloc = 0
      do m=1,nmodes
         if (basis(m)%btyp == btyp_unit) cycle
         if (vloc==0) then
            ! Source is v, target is u.
            udim = vdim
            udim(m) = size(basis(m)%b,1)
            plen = product(udim)
            allocate(u(plen))
            call vgn_shape(m, vdim, vd, gd, nd)
            call matrix_tensor_nt(basis(m)%b, v, u, vd, gd, nd, udim(m))
            deallocate(v)
            vloc = 1
         elseif (vloc==1) then
            ! Source is u, target is v.
            vdim = udim
            vdim(m) = size(basis(m)%b,1)
            plen = product(vdim)
            allocate(v(plen))
            call vgn_shape(m, udim, vd, gd, nd)
            call matrix_tensor_nt(basis(m)%b, u, v, vd, gd, nd, vdim(m))
            deallocate(u)
            vloc = 0
         endif
      enddo
      ! Copy result to v, if necessary
      if (vloc==1) then
         ! Source is u.
         plen = size(u)
         allocate(v(plen))
         v = u
         vdim = udim
         deallocate(u)
      endif
   end subroutine expand_core



   !--------------------------------------------------------------------
   subroutine build_dmat(v, vd, gd, nd, dm)
   ! Generate the density matrix along the middle dimension of a
   ! 3-dim. tensor.
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)    :: vd,gd,nd
      real(dbl),intent(in)  :: v(vd,gd,nd)
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


end module tuckerdecomp_m
