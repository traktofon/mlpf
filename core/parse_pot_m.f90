module parse_pot_m

   use tokenize_m
   use strutil_m
   implicit none

   contains

   function parse_pot(tkner) result (flag)
      type(tokenizer_t),intent(inout) :: tkner
      logical                         :: flag
      character(len=maxtoklen)        :: token

      token = tkner%get()
      call ucase(token)
      if (token /= "POTENTIAL-SECTION" .and. &
          token /= "POT-SECTION") then
         flag = .false.
         return
      endif
      call tkner%clear_stop
      call tkner%add_stop("END-POTENTIAL-SECTION")
      call tkner%add_stop("END-POT-SECTION")
      call tkner%gofwd()
      flag = .true.

      do
         token = tkner%get()
         if (token == "(END)") exit
         call lcase(token)
         ! TODO
         print *, "POT:", trim(token)
         call tkner%gofwd()
      enddo

      if (tkner%stopreason() /= STOPREASON_STOPWORD) &
         call tkner%error("expected end of section")
      call tkner%clear_stop
      call tkner%gofwd()

   end function parse_pot

end module parse_pot_m
