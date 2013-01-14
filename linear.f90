! vim: set ts=3 sw=3 :
module linear

   use base
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine matrix_tensor_tp(m, v, u, n1, n2, n3, nn)
   ! Multiplies the transpose of a rectangular matrix to a tensor,
   ! u = m^T * v.
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)   :: n1,n2,n3,nn
      real(dbl),intent(in) :: m(n2,nn)
      real(dbl),intent(in) :: v(n1,n2,n3)
      real(dbl),intent(out):: u(n1,nn,n3)
      integer              :: i,j,jj,k
      u = 0.d0
      do k=1,n3
         do jj=1,nn
            do j=1,n2
               do i=1,n1
                  u(i,jj,k) = u(i,jj,k) + m(j,jj)*v(i,j,k)
               enddo
            enddo
         enddo
      enddo
   end subroutine matrix_tensor_tp


   !--------------------------------------------------------------------
   subroutine matrix_tensor_nt(m, v, u, n1, n2, n3, nn)
   ! Multiplies a rectangular matrix to a tensor,
   ! u = m * v.
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)   :: n1,n2,n3,nn
      real(dbl),intent(in) :: m(nn,n2)
      real(dbl),intent(in) :: v(n1,n2,n3)
      real(dbl),intent(out):: u(n1,nn,n3)
      integer              :: i,j,jj,k
      u = 0.d0
      do k=1,n3
         do j=1,n2
            do jj=1,nn
               do i=1,n1
                  u(i,jj,k) = u(i,jj,k) + m(jj,j)*v(i,j,k)
               enddo
            enddo
         enddo
      enddo
   end subroutine matrix_tensor_nt


end module linear
