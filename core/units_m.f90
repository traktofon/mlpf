module units_m

   use base_m
   use tokenize_m
   implicit none

   contains


   function parse_length(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,dflt)
      else
         val = parse_with_unit(t)
      endif
   end function parse_length


   function parse_energy(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,dflt)
      else
         val = parse_with_unit(t)
      endif
   end function parse_energy


   function parse_mass(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,dflt)
      else
         val = parse_with_unit(t)
      endif
   end function parse_mass


   function parse_with_unit(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      character(len=maxtoklen) :: token
      integer :: ierr
      token = t%get()
      read(token,*,iostat=ierr) val
      if (ierr/=0) then
         if (present(dflt)) then
            val = dflt
            return
         else
            call t%error("expected real number")
         endif
      endif
      call t%gofwd
      token = t%get()
      if (token /= ",") return
      call t%gofwd
      token = t%get()
      call stopnow("units_m::parse_with_unit: units not implemented yet")
   end function parse_with_unit

end module units_m
