module version_m

   use base_m
   implicit none

   integer,parameter :: v_major = 0
   integer,parameter :: v_minor = 9
   integer,parameter :: v_patch = 0

   real(dbl),parameter :: mctdh_compat_versnum = 8.5001d0

   contains

   pure function versnum()
      real(dbl) :: versnum
      versnum = v_major*1.d0 + v_minor*1.d-2 + v_patch*1.d-5
   end function versnum


   pure function verstring()
      character(len=80) :: verstring
      write(verstring,'(i0,a,i0,a,i0,a)') &
         v_major, '.', v_minor, '.', v_patch, ' (hg id = #HGID#)'
   end function verstring
   

   pure function compdate()
      character(len=80) :: compdate
      write(compdate,'(a)') 'compiled on: #COMPDATE#'
   end function compdate


end module version_m
! vim: syntax=fortran :
