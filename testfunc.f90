! vim: set ts=3 sw=3 :
module testfunc

   use base
   implicit none

   contains

   function coulomb(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: f
      real(dbl)            :: r
      real(dbl),parameter  :: eps = 1.d-8
      r = 0.d0
      do f=1,size(x)
         r = r + x(f)**2
      enddo
      v = 1/sqrt(r+eps)
   end function coulomb


   function gauss(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      integer              :: f
      real(dbl)            :: r
      r = 0.d0
      do f=1,size(x)
         r = r + x(f)**2
      enddo
      v = exp(-r)
   end function gauss


   function coulombn(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl),allocatable:: r(:),qq(:)
      real(dbl),parameter  :: eps = 1.d-6
      integer              :: np,i,j,k,nk,c
      integer              :: qi,qj
      np = size(x)/3
      ! x = (x_1 ... x_np y_1 ... y_np z_1 ... z_np)
      nk = (np*(np-1))/2
      allocate(r(nk))
      allocate(qq(nk))
      r = eps
      k = 1
      do i=1,np
         do j=i+1,np
            do c=0,2 ! x,y,z
               r(k) = r(k) + (x(np*c+i) - x(np*c+j))**2
            enddo
            qi = (-1)**i
            qj = (-1)**j
            qq(k) = dble(qi*qj)
            k = k+1
         enddo
      enddo
      v = sum(qq/sqrt(r))
      deallocate(r)
   end function coulombn


   function pes3c(x) result(v)
      real(dbl),intent(in) :: x(:)
      real(dbl)            :: v
      real(dbl)            :: dr1,dr2,dr4,r3(3),zpp,u1,u2,dihedral,angle1,angle2
      real(dbl),parameter  :: d0 = 1.6d0
      dr1      = x(1)
      dr2      = x(2)
      dr4      = x(3)
      zpp      = x(4)
      r3(1)    = x(5)
      r3(2)    = x(6)
      dihedral = x(7)
      u1       = x(8)
      u2       = x(9)
      r3(3)    = zpp*(dr4-2.d0*d0)
      angle1   = dacos(u1)
      angle2   = dacos(u2)
      call pes3cvpd(dr1,dr2,dr4,r3,angle1,angle2,dihedral,v)
      return
   end function pes3c

end module testfunc
