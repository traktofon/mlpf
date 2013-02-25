module tokenize_m

   use strutil_m
   implicit none

   integer,parameter :: maxlinlen = 240
   integer,parameter :: maxtoklen = 200
   integer,parameter :: maxtok = 48

   type :: tokenizer_t
      character(len=maxtoklen) :: tokbuf(maxtok)
      integer                  :: rpos
      integer                  :: wpos
      character(len=maxtoklen) :: stoptoks(maxtok)
      integer                  :: nstop
      character(len=maxtoklen) :: ignotoks(maxtok)
      integer                  :: nigno
      integer                  :: lun
      integer                  :: linenumber
      character(len=maxlinlen) :: currentline
      logical                  :: stopped
      contains
      procedure :: init
      procedure :: add_stop
      procedure :: add_igno
      procedure :: clear_stop
      procedure :: clear_igno
      procedure :: is_stopped
      procedure :: get
      procedure :: put
      procedure :: gofwd
     !procedure :: goback
      procedure :: produce
      procedure :: error
   end type tokenizer_t
 

   contains


   subroutine init(t,lun)
      class(tokenizer_t),intent(inout) :: t
      integer,intent(in)               :: lun
      t%lun = lun
      t%rpos = 1
      t%wpos = 1
      t%linenumber = 0
      t%stopped = .false.
      call t%clear_stop
      call t%clear_igno
      call t%produce
   end subroutine init


   subroutine clear_stop(t)
      class(tokenizer_t),intent(inout) :: t
      t%nstop = 0
   end subroutine clear_stop


   subroutine add_stop(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      if (t%nstop == maxtok) &
         call t%error("too many stop-tokens")
      t%nstop = t%nstop + 1
      t%stoptoks(t%nstop) = word
   end subroutine add_stop


   subroutine clear_igno(t)
      class(tokenizer_t),intent(inout) :: t
      t%nigno = 0
   end subroutine clear_igno


   subroutine add_igno(t,word)
      class(tokenizer_t),intent(inout) :: t
      character(len=*),intent(in)      :: word
      if (t%nigno == maxtok) &
         call t%error("too many ignore-tokens")
      t%nigno = t%nigno + 1
      t%ignotoks(t%nigno) = word
   end subroutine add_igno


   function is_stopped(t)
      class(tokenizer_t),intent(inout) :: t
      logical                          :: is_stopped
      if (.not. t%stopped) then
         ! check again
         t%stopped = t%rpos==t%wpos .or. &
                     is_element_of(t%tokbuf(t%rpos), t%stoptoks(1:t%nstop))
      endif
      is_stopped = t%stopped
   end function is_stopped


   function get(t) result (token)
      class(tokenizer_t),intent(inout) :: t
      character(len=maxtoklen)         :: token
      do
         if (t%is_stopped()) then
            token = "(END)"
            exit
         else
            token = t%tokbuf(t%rpos)
            if (.not. is_element_of(token, t%ignotoks(1:t%nigno))) exit
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
      character(len=3),parameter       :: ws = " "//ACHAR(10)//ACHAR(13)

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
            call split_word(line(i0:i1), tbuf, ntok)
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
      token = t%tokbuf(t%rpos)
      write (*,'(4a,/,3x,a,i0,a,/,3x,3a)') &
         'Parse Error at "',trim(token),'": ',trim(msg), &
         'encountered at line ',t%linenumber,', which reads:', &
         '"',trim(t%currentline),'"'
      stop 1
   end subroutine error


   subroutine split_word(word,toks,ntok)
      character(len=*),intent(in)          :: word
      character(len=maxtoklen),intent(out) :: toks(:)
      integer,intent(out)                  :: ntok
      integer                              :: wlen,tlen,i0,i1
      character(len=11),parameter          :: sc = "()=<>,{}[]#"

      wlen = len(word)
      tlen = len(toks)
      ntok = 0
      i0 = 1
      do while (i0 <= wlen)
         i1 = scan(word(i0:wlen),sc)
         if (i1==0) then
            i1 = wlen
         elseif (i1==1) then
            i1 = i0
         else
            i1 = i0+i1-2
         endif
         ntok = ntok+1
         if (ntok>tlen) goto 600
         toks(ntok) = word(i0:i1)
         i0 = i1+1
      enddo
      return

 600  print *, "tbuf overflow"
      stop 1
   end subroutine split_word


   pure function is_comment(token)
      logical                     :: is_comment
      character(len=*),intent(in) :: token
      is_comment = .false.
      if (token(1:1) == "#") is_comment = .true.
   end function is_comment


end module tokenize_m
