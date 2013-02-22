module strutil_m

   implicit none
   private

   public :: lcase,ucase,strcmpci

   integer,parameter :: ilca = ichar('a')
   integer,parameter :: ilcz = ichar('z')
   integer,parameter :: iuca = ichar('A')
   integer,parameter :: iucz = ichar('Z')

   contains

   pure function lc(c)
      character,intent(in) :: c
      character            :: lc
      integer              :: i
      i = ichar(c)
      if (i >= iuca .and. i <= iucz) &
         i = i - iuca + ilca
      lc = char(i)
   end function lc


   pure function uc(c)
      character,intent(in) :: c
      character            :: uc
      integer              :: i
      i = ichar(c)
      if (i >= ilca .and. i <= ilcz) &
         i = i - ilca + iuca
      uc = char(i)
   end function uc


   subroutine lcase(str)
      character(len=*),intent(inout) :: str
      integer                        :: i
      do i=1,len(str)
         str(i:i) = lc(str(i:i))
      enddo
   end subroutine lcase


   subroutine ucase(str)
      character(len=*),intent(inout) :: str
      integer                        :: i
      do i=1,len(str)
         str(i:i) = uc(str(i:i))
      enddo
   end subroutine ucase


   pure function strcmpci(a,b) result(cmp)
      character(len=*),intent(in) :: a,b
      integer                     :: cmp
      character                   :: ca,cb
      integer                     :: la,lb,i
      la = len_trim(a)
      lb = len_trim(b)
      cmp = 0
      do i=1,min(la,lb)
         ca = lc(a(i:i))
         cb = lc(b(i:i))
         if (ca<cb) then
            cmp = -1
            return
         elseif (ca>cb) then
            cmp = 1
            return
         endif
      enddo
      if (la<lb) then
         cmp = -1
      elseif (la>lb) then
         cmp = 1
      endif
   end function strcmpci

end module strutil_m
