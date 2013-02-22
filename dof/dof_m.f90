module dof_m

   use base_m
   implicit none

   type,abstract :: dof_t
      character(len=c1)     :: label
      integer               :: gdim
      real(dbl),allocatable :: x(:)
      real(dbl),allocatable :: w(:)
      contains
      procedure(init_dof),deferred :: init
      procedure                    :: set_label
   end type dof_t

   type :: dof_tp
      class(dof_t),pointer :: p => null()
   end type dof_tp


   abstract interface
      subroutine init_dof(dof)
         import :: dof_t
         class(dof_t),intent(inout) :: dof
      end subroutine init_dof

      subroutine parse_dof(dof,tokens)
         import :: dof_t,c1
         class(dof_t),pointer         :: dof
         character(len=c1),intent(in) :: tokens(:)
      end subroutine parse_dof

      subroutine unpickle_dof(dof,gdim,ipar,rpar)
         import :: dof_t,dbl
         class(dof_t),pointer        :: dof
         integer,intent(in)          :: gdim
         integer,intent(in)          :: ipar(:)
         real(dbl),intent(in)        :: rpar(:)
      end subroutine unpickle_dof
   end interface


   type :: doftyp_t
      character(len=c1)                      :: sname
      integer                                :: id
      procedure(parse_dof),nopass,pointer    :: parse
      procedure(unpickle_dof),nopass,pointer :: unpickle
   end type doftyp_t

   integer               :: num_doftyps = 0
   integer,parameter     :: max_doftyps = 32
   type(doftyp_t),target :: doftyps(max_doftyps)


   contains


   subroutine set_label(dof,label)
      class(dof_t),intent(inout)  :: dof
      character(len=*),intent(in) :: label
      dof%label = label
   end subroutine set_label


   subroutine register_doftyp(sname,id,parse,unpickle)
      character(len=*),intent(in)     :: sname
      integer                         :: id
      procedure(parse_dof),pointer    :: parse
      procedure(unpickle_dof),pointer :: unpickle
      integer                         :: i
      if (num_doftyps == max_doftyps) &
         call stopnow("dof_m::register_doftyp : too many doftyps")
      i = num_doftyps + 1
      doftyps(i)%sname = sname
      doftyps(i)%id = id
      doftyps(i)%parse => parse
      doftyps(i)%unpickle => unpickle
      num_doftyps = i
   end subroutine register_doftyp


   function find_doftyp_by_sname(sname) result(doftyp)
      type(doftyp_t),pointer      :: doftyp
      character(len=*),intent(in) :: sname
      integer                     :: i
      nullify(doftyp)
      do i=1,num_doftyps
         if (doftyps(i)%sname == sname) then
            doftyp => doftyps(i)
            return
         endif
      enddo
   end function find_doftyp_by_sname


   function find_doftyp_by_id(id) result(doftyp)
      type(doftyp_t),pointer :: doftyp
      integer,intent(in)     :: id
      integer                :: i
      nullify(doftyp)
      do i=1,num_doftyps
         if (doftyps(i)%id == id) then
            doftyp => doftyps(i)
            return
         endif
      enddo
   end function find_doftyp_by_id

end module dof_m
