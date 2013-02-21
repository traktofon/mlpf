module tokenize_m

   implicit none

   integer,parameter :: maxlinlen = 240
   integer,parameter :: maxtoklen = 200
   integer,parameter :: maxtok = 48


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


   subroutine next_token(lun,token,lend)
      integer,intent(in)                   :: lun
      character(len=maxtoklen),intent(out) :: token
      logical,intent(out)                  :: lend
      character(len=maxlinlen)             :: line
      character(len=maxtoklen)             :: tbuf(maxtok)
      integer                              :: i0,i1,ntok,k
      character(len=maxtoklen),save        :: tokbuf(maxtok)
      integer,save                         :: rpos = 1
      integer,save                         :: wpos = 1
      character(len=3),parameter           :: ws = " "//ACHAR(10)//ACHAR(13)

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

      ! get token from buffer
      lend = .false.
      token = tokbuf(rpos)
      rpos  = mod(rpos,maxtok)+1
      return

 500  lend = .true.
      token = " "
      return

 600  print *, "tokbuf overflow"
      stop 1
   end subroutine next_token



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
