module runopts_m

   use base_m
   implicit none

   type :: runopts_t
      logical           :: lgendvr
      logical           :: lgenpot
      logical           :: lgendot
      character(len=c5) :: namedir
      character(len=c5) :: dvrfile
      character(len=c5) :: potfile
      character(len=c5) :: dotfile
      integer           :: vpotfmt
      real(dbl)         :: rmse
   end type runopts_t

   character(len=*),parameter :: NOFILE = "(NONE)"

end module runopts_m
