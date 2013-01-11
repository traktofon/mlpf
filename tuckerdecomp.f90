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
   subroutine compute_basis(v, gdim, m, limit, mdim, basis)
   ! Computes the rank-1 basis tensors for a Tucker decomposition of the
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
      call get_basis_size(eval, limit, nw)
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
   subroutine compute_core(v, vdim, basis, u, udim)
   ! Computes the core tensor for a Tucker decomposition of the tensor v.
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in)     :: v(:)     ! input tensor
      integer,intent(in)       :: vdim(:)  ! shape of v
      type(basis_t),intent(in) :: basis(:) ! rank-1 basis tensors for all modes
      real(dbl),pointer        :: u(:)     ! output core tensor (allocated here)
      integer,intent(in)       :: udim(:)  ! shape of u
      integer                  :: nmodes,m,vloc,plen,vd,gd,nd            
      real(dbl),allocatable    :: tmp1(:),tmp2(:)                        
      integer                  :: tmp1dim(size(vdim)),tmp2dim(size(vdim))

      nmodes = size(vdim)
      vloc = 0
      do m=1,nmodes
         ! Project the input tensor onto the basis of mode m.
         ! But if the basis is unit, then nothing needs to be done.
         if (basis(m)%btyp == btyp_unit) cycle
         ! To avoid unnecessary copying, we keep track of where the
         ! current tensor is stored, and store the transformed tensor
         ! in an unused array.
         if (vloc==0) then
            ! Source is v, target is tmp1.
            tmp1dim = vdim
            tmp1dim(m) = udim(m)
            plen = product(tmp1dim)
            allocate(tmp1(plen))
            call vgn_shape(m, vdim, vd, gd, nd)
            call matrix_tensor(basis(m)%b, v, tmp1, vd, gd, nd, udim(m))
            vloc = 1
         elseif (vloc==1) then
            ! Source is tmp1, target is tmp2.
            tmp2dim = tmp1dim
            tmp2dim(m) = udim(m)
            plen = product(tmp2dim)
            allocate(tmp2(plen))
            call vgn_shape(m, tmp1dim, vd, gd, nd)
            call matrix_tensor(basis(m)%b, tmp1, tmp2, vd, gd, nd, udim(m))
            deallocate(tmp1)
            vloc = 2
         else
            ! Source is tmp2, target is tmp1.
            tmp1dim = tmp2dim
            tmp1dim(m) = udim(m)
            plen = product(tmp1dim)
            allocate(tmp1(plen))
            call vgn_shape(m, tmp2dim, vd, gd, nd)
            call matrix_tensor(basis(m)%b, tmp2, tmp1, vd, gd, nd, udim(m))
            deallocate(tmp2)
            vloc = 1
         endif
      enddo
      ! Copy result to u.
      plen = product(udim)
      allocate(u(plen))
      if (vloc==0) then
         ! Source is v... shouldn't actually happen.
         u = v
      elseif (vloc==1) then
         ! Source is tmp1.
         u = tmp1
         deallocate(tmp1)
      else
         ! Source is tmp2
         u = tmp2
         deallocate(tmp2)
      endif
   end subroutine compute_core



   !--------------------------------------------------------------------
   subroutine build_dmat(v, vd, gd, nd, dm)
   ! Generate the density matrix along the middle dimension of a
   ! rank-3 tensor.
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in)  :: v(vd,gd,nd)
      integer,intent(in)    :: vd,gd,nd
      real(dbl),intent(out) :: dm(gd,gd)
      integer               :: a,b,vi,ni
      do b=1,gd
         do a=b,gd
            dm(a,b) = 0.d0
            do ni=1,nd
               do vi=1,vd
                  dm(a,b) = dm(a,b) + v(vi,a,ni)*v(vi,b,ni)
               enddo
            enddo
            if (a /= b)  dm(b,a) = dm(a,b)
         enddo
      enddo
   end subroutine build_dmat



   !--------------------------------------------------------------------
   subroutine get_basis_size(wghts, limit, bsz)
   !--------------------------------------------------------------------
      implicit none
      real(dbl),intent(in) :: wghts(:) ! list of weights in ascending order
      real(dbl),intent(in) :: limit    ! parameter for determing how many weights to keep
      integer,intent(out)  :: bsz      ! result: number of weights to keep
      integer   :: nwghts,i
      real(dbl) :: wsum
      ! Sum up all weights until limit is reached.
      ! The remaining weights are those we want to keep.
      nwghts = size(wghts)
      wsum = 0.d0
      i = 0
      do while (wsum < limit .and. i < nwghts)
         i = i+1
         wsum = wsum + max(0.d0,wghts(i)) ! ignore negative weights
      enddo
      bsz = nwghts - i + 1
   end subroutine get_basis_size

end module tuckerdecomp
