! vim: set ts=3 sw=3 :
module logging_m

   implicit none
   private

   public :: open_logfile, close_logfile, &
             get_logger, set_logger, write_log

   integer,parameter,public :: &
      LOGLEVEL_DEBUG = 10, &
      LOGLEVEL_INFO  = 20, &
      LOGLEVEL_WARN  = 30, &
      LOGLEVEL_ERROR = 40


   integer,parameter :: max_loggers=50
   integer,parameter :: logger_facility_len=40
   integer,parameter :: maxunit=999


   type :: logger_t
      character(len=logger_facility_len) :: facility
      integer                            :: level
      integer                            :: lun
   end type logger_t


   integer,save        :: current_logger_id=0
   type(logger_t),save :: loggers(max_loggers)


   contains


   subroutine get_logger(id,facility)
      implicit none
      integer,intent(inout)       :: id
      character(len=*),intent(in) :: facility
      integer                     :: i
      ! if id is already set, do nothing
      if (id /= 0) return
      ! search existing loggers for the facility
      do i=1,current_logger_id
         if (loggers(i)%facility == facility) then
            id = i
            return
         endif
      enddo
      ! not found, so create new logger
      if (current_logger_id == max_loggers) then
         write (*,'(a)') 'get_logger: too many loggers'
         stop 1
      endif
      current_logger_id = current_logger_id + 1
      id = current_logger_id
      loggers(id)%facility = facility
      ! by default, log everything to stdout
      loggers(id)%level = 0 
      loggers(id)%lun   = 6
   end subroutine get_logger


   subroutine set_logger(facility,level,lun)
      implicit none
      character(len=*),intent(in) :: facility
      integer,intent(in)          :: level
      integer,intent(in),optional :: lun
      integer                     :: id
      id = 0
      call get_logger(id,facility)
      loggers(id)%level = level
      if (present(lun))  loggers(id)%lun = lun
   end subroutine set_logger


   subroutine write_log(id,level,msg)
   ! log a message of the given level to the logger with the given id
      implicit none
      integer,intent(in)          :: id,level
      character(len=*),intent(in) :: msg
      integer                     :: lun
      if (level >= loggers(id)%level) then
         lun = loggers(id)%lun
         write (lun,'(a)') trim(msg)
         call flush(lun)
      endif
   end subroutine write_log


   subroutine open_logfile(lun,filename)
      implicit none
      character(len=*),intent(in) :: filename
      integer,intent(out)         :: lun
      integer                     :: ierr
      lun = getlun()
      if (lun < 0) then
         write (*,'(a)') 'open_logfile: no I/O unit available'
      endif
      open (unit=lun,iostat=ierr,file=filename,form='formatted',status='unknown')
      if (ierr /= 0) then
         write (*,'(a,i0,2a)') &
            'open_logfile: cannot open, status=', ierr, ', file=', trim(filename)
         stop 1
      endif
   end subroutine open_logfile


   subroutine close_logfile(lun)
      implicit none
      integer,intent(in) :: lun
      close(lun)
   end subroutine close_logfile


   function getlun() result (lun)
   ! finds an unused unit number
      implicit none      
      integer :: lun
      logical :: isopen,exists
      do lun=maxunit,1,-1
         inquire(unit=lun,exist=exists,opened=isopen)
         if (exists.and..not.isopen) return
      enddo
      lun = -1 ! signal error
   end function getlun


end module logging_m
