module runopts_m

   use base_m
   implicit none

   type :: runopts_t
      logical           :: lgendvr
      logical           :: lgenpot
      character(len=c5) :: dvrfile
      character(len=c5) :: potfile
      integer           :: vpotfmt
   end type runopts_t

   character(len=*),parameter :: NOFILE = "(NONE)"

end module runopts_m
