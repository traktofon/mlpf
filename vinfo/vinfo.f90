program vinfo

   use base_m
   use cmdline_m
   use strutil_m
   use dof_m
   use dof_io_m
   use meta_dof_m
   use genpot_m
   implicit none

   character(len=240) :: fname
   character(len=c5) :: arg
   integer :: lun,ierr,vfmt
   integer*4 :: nitems1,niobuf1
   integer :: nitems,niobuf,nleft,nread
   real(dbl) :: fver
   type(dof_tp),pointer :: dofs(:)
   integer :: ndof,f,j
   integer*8 :: vlen
   real(dbl),allocatable :: buf(:)
   real(dbl) :: v,vmin,vmax,vsum,v2sum
   

   call cmdline_init
   arg = cmdline_next_arg()
   if (arg == "") &
      call stopnow("missing filename")
   fname = trim(arg)

   if (endswith(fname,"vpot2")) then
      vfmt = 2
   else
      vfmt = 1
   endif

   call init_doftyps

   open(newunit=lun, file=trim(fname), status="old", form="unformatted", err=510)
   read(unit=lun,iostat=ierr) fver
   if (ierr /= 0) &
      call stopnow("cannot read file version: "//trim(fname))
   dofs => rddvrdef(lun,fver)
   call rdvpotpars(lun)

   vlen = 1
   ndof = size(dofs)
   do f=1,ndof
      vlen = vlen * dofs(f)%p%gdim
   enddo

   if (vfmt==1) then
      read(unit=lun,err=500) nitems1,niobuf1
      nitems = nitems1
      niobuf = niobuf1
      if (nitems /= vlen) &
         call stopnow("inconsistent data sizes: nitems /= vlen")
   else
      close(lun)
      fname = trim(fname)//".raw"
      open(newunit=lun, file=trim(fname), status="old", form="unformatted", access="stream", err=510)
      nitems = vlen
      niobuf = min(128*1024,nitems)
   endif

   allocate(buf(niobuf))
   vmin = 1.d99
   vmax = -1.d99
   vsum = 0.d0
   v2sum = 0.d0

   nleft = nitems
   do while (nleft>0)
      nread = min(nleft,niobuf)
      read(unit=lun,err=500) (buf(j), j=1,nread)
      do j=1,nread
         v = buf(j)
         vmax = max(vmax,v)
         vmin = min(vmin,v)
         vsum = vsum + v
         v2sum = v2sum + v**2
      enddo
      nleft = nleft - nread
   enddo

   write (*,'(a,es15.6)') "Vmax  = ",vmax
   write (*,'(a,es15.6)') "Vmin  = ",vmin
   write (*,'(a,es15.6)') "Vavg  = ",vsum/vlen
   write (*,'(a,es15.6)') "Vstd  = ",sqrt(v2sum/vlen - (vsum/vlen)**2)
   write (*,'(a,es15.6)') "Vnorm = ",sqrt(v2sum/vlen)

   close(lun)
   deallocate(buf)
   return

 500  call stopnow("error reading data from file: "//trim(fname))
 510  call stopnow("cannot open file: "//trim(fname))

end program vinfo
