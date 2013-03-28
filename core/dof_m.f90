module dof_m

   use base_m
   use tokenize_m
   use strutil_m
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

      subroutine parse_dof(dof,tkner)
         import :: dof_t,tokenizer_t
         class(dof_t),pointer            :: dof
         type(tokenizer_t),intent(inout) :: tkner
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

   integer,save               :: num_doftyps = 0
   integer,parameter          :: max_doftyps = 32
   type(doftyp_t),target,save :: doftyps(max_doftyps)


   contains


   subroutine set_label(dof,label)
      class(dof_t),intent(inout)  :: dof
      character(len=*),intent(in) :: label
      dof%label = label
   end subroutine set_label


   function find_dofnum_by_label(label,dofs) result(dofnum)
      character(len=*),intent(in) :: label
      type(dof_tp),intent(in)     :: dofs(:)
      integer                     :: dofnum
      integer                     :: f
      dofnum = 0
      do f=1,size(dofs)
         if (label == dofs(f)%p%label) then
            dofnum = f
            return
         endif
      enddo
   end function find_dofnum_by_label


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


   function find_doftyp_by_name(sname) result(doftyp)
      type(doftyp_t),pointer      :: doftyp
      character(len=*),intent(in) :: sname
      integer                     :: i
      nullify(doftyp)
      do i=1,num_doftyps
         if (strcmpci(doftyps(i)%sname,sname) == 0) then
            doftyp => doftyps(i)
            return
         endif
      enddo
   end function find_doftyp_by_name


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


   function gridsize(dofs)
      integer*8               :: gridsize
      type(dof_tp),intent(in) :: dofs(:)
      integer                 :: f
      gridsize = 1
      do f=1,size(dofs)
         gridsize = gridsize * dofs(f)%p%gdim
      enddo
   end function gridsize


   subroutine dofs_shape(dofs,vdim)
      type(dof_tp),intent(in) :: dofs(:)
      integer,intent(out)     :: vdim(:)
      integer                 :: f
      do f=1,size(dofs)
         vdim(f) = dofs(f)%p%gdim
      enddo
   end subroutine dofs_shape

end module dof_m
