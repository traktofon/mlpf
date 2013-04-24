module parse_run_m

   use tokenize_m
   use units_m
   use strutil_m
   use runopts_m
   implicit none

   contains

   !--------------------------------------------------------------------
   function parse_run(tkner,opts) result (flag)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: tkner
      type(runopts_t),intent(inout)   :: opts
      logical                         :: flag
      character(len=maxtoklen)        :: token

      ! Check if the RUN-SECTION is starting.
      token = tkner%get()
      if (strcmpci(token,"RUN-SECTION") /= 0) then
         flag = .false.
         return
      endif

      ! Watch out for the end of this section.
      call tkner%clear_stop
      call tkner%add_stop("END-RUN-SECTION")
      call tkner%gofwd
      flag = .true.

      ! Defaults for options are already set.

      ! Loop until the end of the section.
      do
         token = tkner%get()
         if (token == "(END)") exit
         call lcase(token)
      
         if (token == "name") then
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%namedir = trim(token)
               call tkner%gofwd
            else
               call tkner%error("keyword needs an option: "//trim(token))
            endif

         elseif (token == "gendvr") then
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
            opts%lgenpf  = .true.
            call tkner%gofwd

         elseif (token == "readvpot") then
            opts%lgenpot = .false.
            opts%lgenpf  = .true.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%vpotfile = trim(token)
               call tkner%gofwd
            endif

         elseif (token == "readpf") then
            opts%lgenpot = .false.
            opts%lgenpf  = .false.
            call tkner%gofwd
            if (have_option1(tkner)) then
               token = tkner%get()
               opts%pfdir = trim(token)
               call tkner%gofwd
            else
               call tkner%error("keyword needs an option: "//trim(token))
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

      ! Check if the section ended or if there was sth unexpected.
      if (tkner%stopreason() /= STOPREASON_STOPWORD) &
         call tkner%error("expected end of section")
      call tkner%clear_stop
      call tkner%gofwd

   end function parse_run

end module parse_run_m
