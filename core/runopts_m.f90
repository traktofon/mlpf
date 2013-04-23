module runopts_m

   use base_m
   implicit none

   character(len=*),parameter :: NOFILE = ""

   type :: runopts_t
      logical           :: lgendvr = .true.
      logical           :: lgenpot = .true.
      logical           :: lgendot = .false.
      character(len=c5) :: namedir = NOFILE
      character(len=c5) :: dvrfile = NOFILE
      character(len=c5) :: potfile = NOFILE
      character(len=c5) :: dotfile = NOFILE
      integer           :: vpotfmt = 1
      real(dbl)         :: rmse    = 0.d0
   end type runopts_t

end module runopts_m
