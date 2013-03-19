! vim: set ts=3 sw=3 :
module genpot_m

   use base_m
   use logging_m
   use dof_m
   use dof_io_m
   use mmap_m
   implicit none

   contains

   !--------------------------------------------------------------------
   subroutine buildpot(fn,dofs,v,vnorm,vmax,vmin)
   !--------------------------------------------------------------------
      interface
         function fn(x)
            use base_m
            real(dbl),intent(in) :: x(:)
            real(dbl)            :: fn
         end function fn
      end interface
      type(dof_tp),intent(in) :: dofs(:)
      real(dbl),intent(out)   :: v(:)
      real(dbl),intent(out)   :: vnorm,vmax,vmin
      integer,save            :: logid=0
      character(len=160)      :: msg
      integer                 :: j,f,ndofs
      integer                 :: idx(size(dofs))
      real(dbl)               :: x(size(dofs))
      real(dbl)               :: v2sum

      call get_logger(logid,"data")
      write (msg,'(a,i0)') 'Full potential size = ',size(v)
      call write_log(logid, LOGLEVEL_INFO, msg)

      v2sum  = 0.d0
      vmax   = -1.d99
      vmin   =  1.d99
      ndofs  = size(dofs)
      j      = 1
      idx(:) = 1
      do f=1,ndofs
         x(f) = dofs(f)%p%x(idx(f))
      enddo
  10  v(j) = fn(x)
      vmax = max(v(j),vmax)
      vmin = min(v(j),vmin)
      v2sum = v2sum + v(j)**2
      j    = j+1
      do f=1,ndofs
         if (idx(f) == dofs(f)%p%gdim) then
            idx(f) = 1
            x(f)   = dofs(f)%p%x(idx(f))
         else
            idx(f) = idx(f)+1
            x(f)   = dofs(f)%p%x(idx(f))
            goto 10
         endif
      enddo
      vnorm = sqrt(v2sum)

      write (msg,'(a,g22.15)') '||v|| = ',vnorm
      call write_log(logid, LOGLEVEL_INFO, msg)
      write (msg,'(a,g22.15)') 'v_max = ',vmax
      call write_log(logid, LOGLEVEL_INFO, msg)
      write (msg,'(a,g22.15)') 'v_min = ',vmin
      call write_log(logid, LOGLEVEL_INFO, msg)

   end subroutine buildpot



   !--------------------------------------------------------------------
   subroutine loadpot(fname,vfmt,dofs,v)
   !--------------------------------------------------------------------
      character(len=*),intent(in) :: fname
      integer,intent(in)          :: vfmt
      type(dof_tp),pointer        :: dofs(:)
      real(dbl),pointer           :: v(:)
      integer                     :: lun,ierr,f,j
      real(dbl)                   :: fver
      integer*8                   :: vlen
      integer                     :: nitems,niobuf,vpos,nleft,nread

      ! Open the vpot file.
      open(newunit=lun, file=trim(fname), status="old", form="unformatted", iostat=ierr)
      if (ierr /= 0) &
         call stopnow("cannot open file: "//trim(fname))

      ! Read Headers...
      ! File version
      read(unit=lun,iostat=ierr) fver
      if (ierr /=0 ) &
         call stopnow("cannot read file version: "//trim(fname))
      ! DVR information
      dofs => rddvrdef(lun,fver)
      ! vpot parameters -- TODO: which of those information do we need?
      call rdvpotpars(lun)
      ! ...done.

      ! How much data is there?
      vlen = 1
      do f=1,size(dofs)
         vlen = vlen * dofs(f)%p%gdim
      enddo

      ! Get the data.
      if (vfmt==1) then
         ! Traditional vpot file. No mapping.
         allocate(v(vlen),stat=ierr)
         if (ierr /=0 ) &
            call stopnow("loadpot: cannot allocate memory")
         read(unit=lun,err=500) nitems,niobuf
         if (nitems /= vlen) &
            call stopnow("loadpot: inconsistent data sizes")
         vpos = 0
         nleft = nitems
         do while (nleft>0)
            nread = min(nleft,niobuf)
            read(unit=lun,err=500) (v(vpos+j), j=1,nread)
            nleft = nleft - nread
            vpos = vpos + nread
         enddo
         close(lun)

      else
         ! Separate file with vpot data (no record markers).
         close(lun)
         v => mmap_dbl(trim(fname)//".raw", vlen, 0)
      endif
      return
 
 500  call stopnow("error reading file: "//trim(fname))

   end subroutine loadpot



   !--------------------------------------------------------------------
   subroutine rdvpotpars(lun)
   !--------------------------------------------------------------------
      integer,intent(in) :: lun
      character(len=80)  :: datetime
      character(len=c2)  :: peslabel
      character(len=c5)  :: pesopts
      integer            :: nlines,i,lsp
      read(unit=lun,err=500) datetime
      ! skip blabla
      read(unit=lun,err=500) nlines
      do i=1,nlines
         read(unit=lun,err=500)
      enddo
      ! PES parameters
      read(unit=lun,err=500) peslabel
      read(unit=lun,err=500) pesopts
      ! Data type
      read(unit=lun,err=500) lsp
      ! Other data, ignore.
      read(unit=lun,err=500) ! lftp,pdofft
      read(unit=lun,err=500) ! pesmin,pesmax
      read(unit=lun,err=500) ! lcspot,jtot,jbf,csmass
      return
 500  call stopnow("error reading vpot parameters")
   end subroutine rdvpotpars

end module genpot_m
