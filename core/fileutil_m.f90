module fileutil_m

   use base_m
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine mkdir(path)
   !--------------------------------------------------------------------
      character(len=*),intent(in) :: path
      character(len=1000)         :: cmd
      integer                     :: ierr

      write(cmd,'(3a)') &
         "mkdir -p '", trim(path), "'"
      call EXECUTE_COMMAND_LINE(cmd, WAIT=.true., EXITSTAT=ierr)
      if (ierr /= 0) then
         call stopnow("cannot create directory: "//trim(path))
      endif

   end subroutine mkdir

end module fileutil_m
