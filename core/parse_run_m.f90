module parse_run_m

   use tokenize_m
   use strutil_m
   implicit none

   contains

   function parse_run(tkner) result (flag)
      type(tokenizer_t),intent(inout) :: tkner
      logical                         :: flag
      character(len=maxtoklen)        :: token

      token = tkner%get()
      if (strcmpci(token,"RUN-SECTION") /= 0) then
         flag = .false.
         return
      endif
      call tkner%clear_stop
      call tkner%add_stop("END-RUN-SECTION")
      call tkner%gofwd()
      flag = .true.

      write (*,'(a)') '(RUN-SECTION currently ignored)'
      do
         token = tkner%get()
         if (token == "(END)") exit
         call lcase(token)
         ! TODO
         call tkner%gofwd()
      enddo

      if (tkner%stopreason() /= STOPREASON_STOPWORD) &
         call tkner%error("expected end of section")
      call tkner%clear_stop
      call tkner%gofwd()

   end function parse_run

end module parse_run_m
