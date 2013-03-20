module parse_run_m

   use tokenize_m
   use strutil_m
   use runopts_m
   implicit none

   contains

   function parse_run(tkner,runopts) result (flag)
      type(tokenizer_t),intent(inout) :: tkner
      type(runopts_t),intent(out)     :: runopts
      logical                         :: flag
      character(len=maxtoklen)        :: token

      token = tkner%get()
      if (strcmpci(token,"RUN-SECTION") /= 0) then
         flag = .false.
         return
      endif
      call tkner%clear_stop
      call tkner%add_stop("END-RUN-SECTION")
      call tkner%gofwd
      flag = .true.

      ! Set defaults.
      runopts%lgendvr = .true.
      runopts%lgenpot = .true.
      runopts%dvrfile = NOFILE
      runopts%potfile = NOFILE
      runopts%vpotfmt = 1

      do
         token = tkner%get()
         if (token == "(END)") exit
         call lcase(token)
      
         if (token == "gendvr") then
            runopts%lgendvr = .true.
            call tkner%gofwd

         elseif (token == "readdvr") then
            runopts%lgendvr = .false.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               runopts%dvrfile = trim(token)
               call tkner%gofwd
            endif

         elseif (token == "genpot") then
            runopts%lgenpot = .true.
            call tkner%gofwd

         elseif (token == "readpot") then
            runopts%lgenpot = .false.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               runopts%potfile = trim(token)
               call tkner%gofwd
            endif
          
         elseif (token == "vpot-format") then
            call tkner%gofwd
            if (have_option1(tkner)) then
               runopts%vpotfmt = parse_int(tkner)
            else
               call tkner%error("keyword needs an option")
            endif

         elseif (token == "=" .or. token == ",") then
            call tkner%error("expected keyword instead of option")
         else
            call tkner%error("unknown keyword")
         endif
      enddo

      if (tkner%stopreason() /= STOPREASON_STOPWORD) &
         call tkner%error("expected end of section")
      call tkner%clear_stop
      call tkner%gofwd

   end function parse_run

end module parse_run_m
