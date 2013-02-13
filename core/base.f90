! vim: set ts=3 sw=3 :
module base

   implicit none

   integer,parameter :: dbl = selected_real_kind(15,300) ! IEEE754 double precision

   contains

   subroutine stopnow(msg)
      character(len=*),intent(in) :: msg
      write (*,'(2a)') "FATAL: ", msg
      stop 1
   end subroutine stopnow

end module base
