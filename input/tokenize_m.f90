module tokenize_m

   use strutil_m
   implicit none

   integer,parameter :: maxlinlen = 240
   integer,parameter :: maxtoklen = 200
   integer,parameter :: maxtok = 48

   character(len=maxtoklen),private :: tokbuf(maxtok)
   integer,private                  :: rpos
   integer,private                  :: wpos
   character(len=maxtoklen),private :: stoptoks(maxtok)
   integer,private                  :: nstop
   character(len=maxtoklen),private :: ignotoks(maxtok)
   integer,private                  :: nigno
   integer,private                  :: lun
   integer,private                  :: currentlinenumber
   character(len=maxlinlen),private :: currentline
 

   contains

   pure function is_element_of(str,strlist)
      logical                     :: is_element_of
      character(len=*),intent(in) :: str
      character(len=*),intent(in) :: strlist(:)
      integer                     :: w
      is_element_of = .false.
      do w=1,size(strlist)
         if (strcmpci(str,strlist(w))==0) then
            is_element_of = .true.
            exit
         endif
      enddo
   end function is_element_of


   subroutine init_tokenize(lun1,lend)
      integer,intent(in)  :: lun1
      logical,intent(out) :: lend
      lun = lun1
      rpos = 1
      wpos = 1
      nstop = 0
      nigno = 0
      currentlinenumber = 0
      call produce_tokens(lend)
   end subroutine init_tokenize


   subroutine set_stoptoks(words)
      character(len=*),intent(in) :: words(:)
      nstop = size(words)
      if (nstop>0) stoptoks(1:nstop) = words(:)
   end subroutine set_stoptoks


   subroutine set_ignotoks(words)
      character(len=*),intent(in) :: words(:)
      nigno = size(words)
      if (nigno>0) ignotoks(1:nigno) = words(:)
   end subroutine set_ignotoks


   function is_stopped()
      logical :: is_stopped
      is_stopped = (rpos==wpos) .or. is_element_of(tokbuf(rpos), stoptoks(1:nstop))
   end function is_stopped


   subroutine get_token(token)
      character(len=maxtoklen),intent(out) :: token
      logical                              :: lend
      do
         if (rpos==wpos) then
            token = "(EOF)"
            exit
         else
            token = tokbuf(rpos)
            if (is_stopped()) exit
            if (.not. is_element_of(token, ignotoks(1:nigno))) exit
            call next_token(lend)
         endif
      enddo
   end subroutine get_token


   subroutine next_token(lend)
      logical,intent(out) :: lend
      if (is_stopped()) then
         lend = .true.
      else
         lend = .false.
         rpos = mod(rpos,maxtok)+1
         if (rpos==wpos) then
            call produce_tokens(lend)
         endif
      endif
   end subroutine next_token


   subroutine produce_tokens(lend)
      logical,intent(out)        :: lend
      character(len=maxlinlen)   :: line
      character(len=maxtoklen)   :: tbuf(maxtok)
      integer                    :: i0,i1,ntok,k
      character(len=3),parameter :: ws = " "//ACHAR(10)//ACHAR(13)

      lend = .false.
      do while (rpos==wpos)
         ! nothing in the buffer, so read next line
         read(lun,'(a)',end=500) line
         ! store info (for parse_error)
         currentlinenumber = currentlinenumber + 1
         currentline = line
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
               tokbuf(wpos) = tbuf(k)
               wpos = mod(wpos,maxtok)+1
               if (wpos==rpos) goto 600
               enddo
            ! continue processing the line
            i0 = i1+1
         enddo lineloop
         ! insert EOL marker
         tokbuf(wpos) = "(EOL)"
         wpos = mod(wpos,maxtok)+1
         if (wpos==rpos) goto 600
      enddo
      return

 500  lend = .true.
      return
 
 600  call parse_error("too many tokens")

   end subroutine produce_tokens



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


   function is_comment(token)
      logical          :: is_comment
      character(len=*) :: token
      is_comment = .false.
      if (token(1:1) == "#") is_comment = .true.
   end function is_comment


   subroutine parse_error(msg)
      character(len=*),intent(in) :: msg
      character(len=maxtoklen)    :: token
      call get_token(token)
      write (*,'(4a,/,3x,a,i0,a,/,3x,3a)') &
         'Parse Error at "',trim(token),'": ',trim(msg), &
         'encountered at line ',currentlinenumber,', which reads:', &
         '"',trim(currentline),'"'
      stop 1
   end subroutine parse_error


end module tokenize_m
