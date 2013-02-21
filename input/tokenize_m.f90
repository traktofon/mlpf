module tokenize_m

   implicit none

   integer,parameter :: maxlinlen = 240
   integer,parameter :: maxtoklen = 200
   integer,parameter :: maxtok = 48

   character(len=maxtoklen),private :: tokbuf(maxtok)
   integer,private                  :: rpos = 1
   integer,private                  :: wpos = 1
   integer,private                  :: lun = 0
 

   contains


   function antiscan(string,set) result (pos)
   ! return position of first character in string that does *not* belong to set
   ! or zero if there is no such character
      character(len=*),intent(in) :: string,set
      integer                     :: pos
      integer                     :: i,k
      pos = 0
      strloop: do i=1,len(string)
         do k=1,len(set)
            if (string(i:i) == set(k:k)) cycle strloop
         enddo
         pos = i
         exit
      enddo strloop
   end function antiscan


   subroutine init_tokenize(lun1,lend)
      integer,intent(in)  :: lun1
      logical,intent(out) :: lend
      lun = lun1
      rpos = 1
      wpos = 1
      call produce_tokens(lend)
   end subroutine init_tokenize


   subroutine get_token(token)
      character(len=maxtoklen),intent(out) :: token
      if (rpos==wpos) then
         token = "(EOF)"
      else
         token = tokbuf(rpos)
      endif
   end subroutine get_token


   subroutine next_token(lend)
      logical,intent(out) :: lend
      rpos = mod(rpos,maxtok)+1
      if (rpos==wpos) then
         call produce_tokens(lend)
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
            do k=1,ntok
               if (is_comment(tbuf(k))) exit lineloop ! skip rest of line
               tokbuf(wpos) = tbuf(k)
               wpos = mod(wpos,maxtok)+1
               if (wpos==rpos) goto 600
            enddo
            ! continue processing the line
            i0 = i1+1
         enddo lineloop
      enddo
      return

 500  lend = .true.
      return

 600  print *, "tokbuf overflow"
      stop 1
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


end module tokenize_m
