module cmdline_m

   use base_m
   implicit none
   private

   public :: cmdline_init, cmdline_next_arg

   integer,save :: argp, argc

   contains

   subroutine cmdline_init
      argc = command_argument_count()
      argp = 1
   end subroutine cmdline_init


   function cmdline_next_arg() result(arg)
      character(len=c5) :: arg, errmsg
      integer           :: ierr
      if (argp > argc) then
         arg = ""
         return
      else
         call get_command_argument(argp, arg, status=ierr)
         if (ierr < 0) then
            write (errmsg, '(a,i0,a)') &
               'command line argument ',argp,' too long'
            call stopnow(errmsg)
         endif
         argp = argp + 1
      endif
   end function cmdline_next_arg

end module cmdline_m
