module tokenize_m

   use map_str2int_m
   use strutil_m
   use base_m
   implicit none

   integer,parameter :: maxlinlen = 240
   integer,parameter :: maxtoklen = 200
   integer,parameter :: maxtok = 48

   integer,parameter :: STOPREASON_NONE     = 0
   integer,parameter :: STOPREASON_STOPWORD = 1
   integer,parameter :: STOPREASON_EOF      = 2

   type :: tokenizer_t
      character(len=maxtoklen) :: tokbuf(maxtok)
      integer                  :: rpos
      integer                  :: wpos
      type(map_str2int_t)      :: stops
      type(map_str2int_t)      :: ignos
      character(len=c5)        :: filename
      integer                  :: lun
      integer                  :: linenumber
      character(len=maxlinlen) :: currentline
      logical                  :: stopped
      contains
      procedure :: init
      procedure :: add_stop
      procedure :: add_igno
      procedure :: del_stop
      procedure :: del_igno
      procedure :: clear_stop
      procedure :: clear_igno
      procedure :: is_stopped
      procedure :: stopreason
      procedure :: get
      procedure :: put
      procedure :: gofwd
     !procedure :: goback
      procedure :: produce
      procedure :: error
   end type tokenizer_t
 

   contains


   subroutine init(t,filename)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: filename
      integer :: ierr
      t%filename = trim(filename)
      t%linenumber = 0
      t%rpos = 1
      t%wpos = 1
      t%stopped = .false.
      call t%clear_stop
      call t%clear_igno
      open(newunit=t%lun, file=trim(filename), form="formatted", status="old", iostat=ierr)
      if (ierr/=0) call stopnow('cannot open file "'//trim(filename)//'"')
      call t%produce
   end subroutine init


   subroutine clear_stop(t)
      class(tokenizer_t),intent(inout) :: t
      call map_destroy(t%stops)
      t%stopped = .false.
   end subroutine clear_stop


   subroutine add_stop(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      call map_put(t%stops,word,1)
   end subroutine add_stop


   subroutine del_stop(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      call map_del_safe(t%stops,word)
   end subroutine del_stop


   subroutine clear_igno(t)
      class(tokenizer_t),intent(inout) :: t
      call map_destroy(t%ignos)
   end subroutine clear_igno


   subroutine add_igno(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      call map_put(t%ignos,word,1)
   end subroutine add_igno


   subroutine del_igno(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      call map_del_safe(t%ignos,word)
   end subroutine del_igno


   function is_stopped(t)
      class(tokenizer_t),intent(inout) :: t
      logical                          :: is_stopped
      integer                          :: idum
      if (.not. t%stopped) then
         ! check again
         t%stopped = t%rpos==t%wpos .or. &
                     map_has(t%stops, t%tokbuf(t%rpos))
      endif
      is_stopped = t%stopped
   end function is_stopped


   function stopreason(t)
      class(tokenizer_t),intent(in) :: t
      integer                       :: stopreason
      integer                       :: idum
      if (t%rpos == t%wpos) then
         stopreason = STOPREASON_EOF
      elseif (map_has(t%stops, t%tokbuf(t%rpos))) then
         stopreason = STOPREASON_STOPWORD
      else
         stopreason = STOPREASON_NONE
      endif
   end function stopreason


   function get(t) result (token)
      class(tokenizer_t),intent(inout) :: t
      character(len=maxtoklen)         :: token
      do
         if (t%is_stopped()) then
            token = "(END)"
            exit
         else
            token = t%tokbuf(t%rpos)
            if (.not.map_has(t%ignos,token)) exit
            call t%gofwd
         endif
      enddo
   end function get


   subroutine put(t,token)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: token
      t%tokbuf(t%wpos) = token
      t%wpos = mod(t%wpos,maxtok)+1
      if (t%wpos==t%rpos) call t%error("too many tokens")
      t%stopped = .false.
   end subroutine put


   subroutine gofwd(t)
      class(tokenizer_t),intent(inout) :: t
      if (.not.t%is_stopped()) then
         t%rpos = mod(t%rpos,maxtok)+1
         if (t%rpos==t%wpos) then
            call t%produce
         endif
      endif
   end subroutine gofwd


   subroutine produce(t)
      class(tokenizer_t),intent(inout) :: t
      character(len=maxlinlen)         :: line
      character(len=maxtoklen)         :: tbuf(maxtok)
      integer                          :: i0,i1,ntok,k
      character(len=*),parameter       :: ws = " "//ACHAR(10)//ACHAR(13)
      character(len=*),parameter       :: sc = "=,#{}[]()<>"

      do while (t%rpos==t%wpos)
         ! nothing in the buffer, so read next line
         read(t%lun,'(a)',end=500) line
         ! store info (for parse_error)
         t%linenumber = t%linenumber + 1
         t%currentline = line
         ! process the line
         i0 = 1
         lineloop: do while (i0 <= maxlinlen)
            ! find next non-whitespace
            i1 = antiscan(line(i0:maxlinlen),ws)
            if (i1==0) then
               ! no non-whitespace found
               i0 = maxlinlen+1
               cycle
            endif
            i0 = i1+i0-1
            ! find next whitespace
            i1 = scan(line(i0:maxlinlen),ws)
            if (i1==0) then
               i1=maxlinlen
            else
               i1=i1+i0-2
            endif
            ! something is at line(i0:i1)
            call split_word(line(i0:i1), sc, tbuf, ntok)
            ! insert tokens into ringbuffer
            do k=1,ntok
               if (is_comment(tbuf(k))) exit lineloop ! skip rest of line
               call t%put(tbuf(k))
            enddo
            ! continue processing the line
            i0 = i1+1
         enddo lineloop
         ! insert EOL marker
         call t%put("(EOL)")
      enddo
 500  return
   end subroutine produce


   subroutine error(t,msg)
      class(tokenizer_t),intent(in) :: t
      character(len=*),intent(in)   :: msg
      character(len=maxtoklen)      :: token
      if (t%stopreason() == STOPREASON_EOF) then
         token = "(EOF)"
      else
         token = t%tokbuf(t%rpos)
      endif
      write (*,'(4a,/,3x,a,i0,3a,/,3x,3a)') &
         'Parse Error at "',trim(token),'": ',trim(msg), &
         'encountered at line ',t%linenumber,&
         ' of file "',trim(t%filename),'", reading:', &
         '"',trim(t%currentline),'"'
      call stopnow("Input Error")
   end subroutine error


   subroutine split_word(word,splitat,toks,ntok)
      character(len=*),intent(in)          :: word
      character(len=*),intent(in)          :: splitat
      character(len=maxtoklen),intent(out) :: toks(:)
      integer,intent(out)                  :: ntok
      integer                              :: wlen,tlen,i0,i1
      wlen = len(word)
      tlen = len(toks)
      ntok = 0
      i0 = 1
      do while (i0 <= wlen)
         i1 = scan(word(i0:wlen),splitat)
         if (i1==0) then
            i1 = wlen
         elseif (i1==1) then
            i1 = i0
         else
            i1 = i0+i1-2
         endif
         ntok = ntok+1
         if (ntok>tlen) call stopnow("split_word: toks overflow")
         toks(ntok) = word(i0:i1)
         i0 = i1+1
      enddo
      return
   end subroutine split_word


   pure function is_comment(token)
      logical                     :: is_comment
      character(len=*),intent(in) :: token
      is_comment = .false.
      if (token(1:1) == "#") is_comment = .true.
   end function is_comment


   function parse_int(t) result(val)
      type(tokenizer_t),intent(inout) :: t
      integer                         :: val
      character(len=maxtoklen)        :: token
      integer                         :: ierr
      token = t%get()
      read(token,*,iostat=ierr) val
      if (ierr/=0) call t%error("expected integer number")
      call t%gofwd
      return
   end function parse_int


   function parse_real(t) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl)                       :: val
      character(len=maxtoklen)        :: token
      integer                         :: ierr
      token = t%get()
      read(token,*,iostat=ierr) val
      if (ierr/=0) call t%error("expected real number")
      call t%gofwd
      return
   end function parse_real


   function parse_angle(t,flag) result(val)
      type(tokenizer_t),intent(inout) :: t
      logical,intent(out),optional    :: flag
      real(dbl)                       :: val
      character(len=maxtoklen)        :: token
      integer                         :: pord
      token = t%get()
      if (strcmpci(token(1:3),"2pi")==0) then
         if (token(4:4) == " ") then
            pord=1
         elseif (token(4:4) == "/") then
            read (token(5:),*,err=500,end=500) pord
         else
            goto 500
         endif
         if (present(flag)) flag = .true.
         val = 2.d0*pi/pord
         call t%gofwd
      else
         if (.not.present(flag)) goto 500
         flag = .false.
         val  = 0.d0
      endif
      return
 500  call t%error('expected 2pi or 2pi/m')
   end function parse_angle

end module tokenize_m
