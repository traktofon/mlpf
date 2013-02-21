module parsetree_m

   use tokenize_m
   use linkedlist_m
   implicit none

   contains

   subroutine xxx(lun)
      integer,intent(in)       :: lun
      character(len=maxtoklen) :: token
      logical                  :: lend,lempty
      type(llist_int_t)        :: fstack
      integer :: f

      do
         call next_token(lun,token,lend)
         if (lend) exit
         if (token == "(") then
            call fstack%push(-1)
         elseif (token == ")") then
            continue
         else
            continue
         endif
      enddo

   end subroutine

end module parsetree_m
