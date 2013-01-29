! vim: set ts=3 sw=3 :
module genpot

   use base
   use logging
   use dof
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine buildpot(fn,dofs,v,vnorm,vmax,vmin)
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
      real(dbl),intent(out)   :: vnorm,vmax,vmin
      integer,save            :: logid=0
      character(len=160)      :: msg
      integer                 :: j,f,ndofs
      integer                 :: idx(size(dofs))
      real(dbl)               :: x(size(dofs))
      real(dbl)               :: v2sum

      call get_logger(logid,"potential")
      write (msg,'(a,i0)') 'Full potential size = ',size(v)
      call write_log(logid, LOGLEVEL_INFO, msg)

      v2sum  = 0.d0
      vmax   = -1.d99
      vmin   =  1.d99
      ndofs  = size(dofs)
      j      = 1
      idx(:) = 1
      do f=1,ndofs
         x(f) = dofs(f)%p%x(idx(f))
      enddo
  10  v(j) = fn(x)
      vmax = max(v(j),vmax)
      vmin = min(v(j),vmin)
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
      vnorm = sqrt(v2sum)

      write (msg,'(a,g22.15)') '||v|| = ',vnorm
      call write_log(logid, LOGLEVEL_INFO, msg)
      write (msg,'(a,g22.15)') 'v_max = ',vmax
      call write_log(logid, LOGLEVEL_INFO, msg)
      write (msg,'(a,g22.15)') 'v_min = ',vmin
      call write_log(logid, LOGLEVEL_INFO, msg)

   end subroutine buildpot


end module genpot
