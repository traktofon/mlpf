! vim: set ts=3 sw=3 :
module genpot

   use base
   use dof
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine buildpot(fn,dofs,v)
   !--------------------------------------------------------------------
      implicit none
      interface
         function fn(x)
            use base
            real(dbl),intent(in) :: x(:)
            real(dbl)            :: fn
         end function fn
      end interface
      type(dof_tp),intent(in) :: dofs(:)        
      real(dbl),intent(out)   :: v(:)           
      integer                 :: j,f,ndofs      
      integer                 :: idx(size(dofs))
      real(dbl)               :: x(size(dofs))  
      real(dbl)               :: v2sum          

      v2sum  = 0.d0
      ndofs  = size(dofs)
      j      = 1
      idx(:) = 1
      do f=1,ndofs
         x(f) = dofs(f)%p%x(idx(f))
      enddo
  10  v(j) = fn(x)
      v2sum = v2sum + v(j)**2
      j    = j+1
      do f=1,ndofs
         if (idx(f) == dofs(f)%p%gdim) then
            idx(f) = 1
            x(f)   = dofs(f)%p%x(idx(f))
         else
            idx(f) = idx(f)+1
            x(f)   = dofs(f)%p%x(idx(f))
            goto 10
         endif
      enddo
      print *, sqrt(v2sum)
   end subroutine buildpot

end module genpot
