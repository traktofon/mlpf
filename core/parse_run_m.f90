module parse_run_m

   use tokenize_m
   use units_m
   use strutil_m
   use runopts_m
   implicit none

   contains

   function parse_run(tkner,opts) result (flag)
      type(tokenizer_t),intent(inout) :: tkner
      type(runopts_t),intent(out)     :: opts
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
      opts%lgendvr = .true.
      opts%lgenpot = .true.
      opts%lgendot = .false.
      opts%dvrfile = NOFILE
      opts%potfile = NOFILE
      opts%dotfile = NOFILE
      opts%vpotfmt = 1
      opts%rmse = 0.d0

      do
         token = tkner%get()
         if (token == "(END)") exit
         call lcase(token)
      
         if (token == "gendvr") then
            opts%lgendvr = .true.
            call tkner%gofwd

         elseif (token == "readdvr") then
            opts%lgendvr = .false.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%dvrfile = trim(token)
               call tkner%gofwd
            endif

         elseif (token == "genpot") then
            opts%lgenpot = .true.
            call tkner%gofwd

         elseif (token == "readpot") then
            opts%lgenpot = .false.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%potfile = trim(token)
               call tkner%gofwd
            endif
          
         elseif (token == "vpot-format") then
            call tkner%gofwd
            if (have_option1(tkner)) then
               opts%vpotfmt = parse_int(tkner)
            else
               call tkner%error("keyword needs an option: "//trim(token))
            endif

         elseif (token == "graphviz") then
            opts%lgendot = .true.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%dotfile = trim(token)
               call tkner%gofwd
            endif

         elseif (token == "rmse") then
            call tkner%gofwd
            if (have_option1(tkner)) then
               opts%rmse = parse_energy(tkner)
            else
               call tkner%error("keyword needs an option: "//trim(token))
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
