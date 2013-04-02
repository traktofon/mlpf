! vim: set ts=3 sw=3 :
module genpot_m

   use base_m
   use logging_m
   use dof_m
   use dof_io_m
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
      real(dbl),pointer       :: v(:)
      real(dbl),intent(out)   :: vnorm,vmax,vmin
      integer,save            :: logid=0
      character(len=160)      :: msg
      integer                 :: vlen,ierr,j,f,ndofs
      integer                 :: idx(size(dofs))
      real(dbl)               :: x(size(dofs))
      real(dbl)               :: v2sum

      call get_logger(logid,"data")
      write (msg,'(a,i0)') 'Full potential size = ',size(v)
      call write_log(logid, LOGLEVEL_INFO, msg)

      vlen = gridsize(dofs)
      allocate(v(vlen),stat=ierr)
      if (ierr /= 0) &
         call stopnow("buildpot: cannot allocate memory")

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
      integer                     :: lun,ierr,ndof,f,j,k
      type(dof_tp),pointer        :: vdofs(:)
      real(dbl),allocatable       :: buf(:)
      real(dbl)                   :: fver
      integer*8                   :: vlen
      integer                     :: nitems,niobuf,nleft,nread
      character(len=c5)           :: fname1
      integer                     :: vdim(size(dofs))
      integer                     :: vidx(size(dofs))
      integer                     :: fmap(size(dofs))

      ! Open the vpot file.
      fname1 = fname
      open(newunit=lun, file=trim(fname1), status="old", form="unformatted", err=510)

      ! Read Headers...
      ! File version
      read(unit=lun,iostat=ierr) fver
      if (ierr /=0 ) &
         call stopnow("cannot read file version: "//trim(fname1))
      ! DVR information
      vdofs => rddvrdef(lun,fver)
      ! vpot parameters -- TODO: which of those information do we need?
      call rdvpotpars(lun)
      ! ...done.

      ! map vpot-DOFs to system-DOFs
      ndof = size(dofs)
      do f=1,ndof
         fmap(f) = find_dofnum_by_label(dofs(f)%p%label, vdofs)
         if (fmap(f) == 0) &
            call stopnow("DOF not present in vpot: "//trim(dofs(f)%p%label))
      enddo

      ! check vpot DOFs and system DOFs for consistency
      do f=1,ndof
         if (.not. are_dofs_compatible(dofs(f)%p, vdofs(fmap(f))%p)) then
            call stopnow('Definition of DOF "' // trim(dofs(f)%p%label) // &
               '" from vpot-file differs from system definition.')
         endif
      enddo

      ! How much data is there?
      vlen = 1
      do f=1,size(vdofs)
         vdim(f) = vdofs(f)%p%gdim
         vlen = vlen * vdim(f)
      enddo
      allocate(v(vlen),stat=ierr)
      if (ierr /=0 ) &
         call stopnow("loadpot: cannot allocate memory")

      ! Prepare for reading the data.
      if (vfmt==1) then
         ! Traditional vpot file.
         read(unit=lun,err=500) nitems,niobuf
         if (nitems /= vlen) &
            call stopnow("loadpot: inconsistent data sizes")
      else
         ! Separate file with vpot data (no record markers).
         close(lun)
         fname1 = trim(fname)//".raw"
         open(newunit=lun, file=trim(fname1), status="old", &
              form="unformatted", access="stream", iostat=ierr)
         nitems = vlen
         niobuf = min(128*1024, nitems)
      endif

      ! Read all data.
      allocate(buf(niobuf))
      vidx = 0
      nleft = nitems
      do while (nleft>0)
         nread = min(nleft,niobuf)
         read(unit=lun,err=500) (buf(j), j=1,nread)
         do j=1,nread
            ! remap the grid point from vpot-order to system-order
            k = 0
            do f=ndof,1,-1
               k = k*vdim(fmap(f)) + vidx(fmap(f))
            enddo
            ! store data of this grid point
            v(k+1) = buf(j)
            ! advance to next grid point
            do f=1,ndof
               vidx(f) = vidx(f)+1
               if (vidx(f) < vdim(f)) exit
               vidx(f) = 0
            enddo
         enddo
         nleft = nleft - nread
      enddo

      close(lun)
      deallocate(buf)
      return
 
 500  call stopnow("error reading file: "//trim(fname1))
 510  call stopnow("cannot open file: "//trim(fname1))

      contains

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
