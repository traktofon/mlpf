module vtree_dbg_m

   use vtree_m
   use logging_m
   implicit none
   private

   public :: examine_vtree

   contains

   !-------------------------------------------------------------------
   subroutine examine_vtree(t)
   ! Print various information about the tree structure.
   ! Mainly for debugging.
   !--------------------------------------------------------------------
      type(vtree_t),intent(in) :: t
      integer                  :: m
      integer,save             :: logid=0
      character(len=400)       :: msg
      character,parameter      :: NL = ACHAR(10)

      call get_logger(logid,"tree")

      write (msg,'(a,4(a,i4,a))') &
         'Tree has:', &
         NL, t%numnodes, ' nodes', &
         NL, t%numdofs, ' dofs', &
         NL, t%numleaves, ' leaves', &
         NL, t%numlayers, ' layers'
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(1x,i0))') &
         'Levels are:', NL, &
         (t%preorder(m)%p%layer, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(1x,l1))') &
         'Leaves?', NL, &
         (t%preorder(m)%p%isleaf, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(1x,i0))') &
         'Leaves are:', NL, &
         (t%leaves(m)%p%num, m=1,t%numleaves)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(1x,i0))') &
         'Pre-order is:', NL, &
         (t%preorder(m)%p%num, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

      write (msg,'(2a,3x,99(1x,i0))') &
         'Post-order is:', NL, &
         (t%postorder(m)%p%num, m=1,t%numnodes)
      call write_log(logid, LOGLEVEL_DEBUG, msg)

   end subroutine examine_vtree

end module vtree_dbg_m
