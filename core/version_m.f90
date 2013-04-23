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
      character(len=16) :: verstring
      write(verstring,'(i0,a,i0,a,i0)') &
         v_major, '.', v_minor, '.', v_patch
   end function verstring
   
end module version_m
