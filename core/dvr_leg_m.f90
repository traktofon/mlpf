module dvr_leg_m

   use dvr_m
   use dof_m
   use tokenize_m
   use base_m
   implicit none

   type,extends(dvr_t) :: dvr_leg_t
      integer :: blz ! magnetic quantum number
      integer :: sym ! 0:all  1:odd  2:even
      contains
      procedure :: init   => init_leg
      procedure :: pickle => pickle_leg
   end type dvr_leg_t

   integer,parameter,private :: typid = 2 ! MCTDH basis type

   contains


   !--------------------------------------------------------------------
   subroutine parse_leg(dof,tkner)
   !--------------------------------------------------------------------
   ! legdvr :~ INTEGER INTEGER ( "all" | "odd" | "even" )
   !--------------------------------------------------------------------
      class(dof_t),pointer            :: dof
      type(tokenizer_t),intent(inout) :: tkner
      allocate(dvr_leg_t::dof)
      select type (dof)
      type is (dvr_leg_t)
      dof%gdim = parse_int(tkner)
      dof%blz  = parse_int(tkner)
      dof%sym  = parse_sym(tkner)
      end select
   end subroutine parse_leg


   !--------------------------------------------------------------------
   subroutine unpickle_leg(dof,gdim,ipar,rpar)
   !--------------------------------------------------------------------
      class(dof_t),pointer        :: dof
      integer,intent(in)          :: gdim
      integer,intent(in)          :: ipar(:)
      real(dbl),intent(in)        :: rpar(:)
      allocate(dvr_leg_t::dof)
      dof%gdim  = gdim
      select type(dof)
      type is (dvr_leg_t)
      dof%blz = ipar(1)
      dof%sym = ipar(2)
      end select
   end subroutine unpickle_leg


   !--------------------------------------------------------------------
   subroutine pickle_leg(dof,id,ipar,rpar)
   !--------------------------------------------------------------------
      class(dvr_leg_t),intent(inout) :: dof
      integer,intent(out)            :: id
      integer,intent(out)            :: ipar(:)
      real(dbl),intent(out)          :: rpar(:)
      id = typid
      ipar = 0
      rpar = 0.d0
      ipar(1) = dof%blz
      ipar(2) = dof%sym
   end subroutine pickle_leg


   !--------------------------------------------------------------------
   subroutine init_leg(dof)
   !--------------------------------------------------------------------
      class(dvr_leg_t),intent(inout) :: dof
      ! call general DVR constructor
      call dof%init_dvr
      ! TODO
      call stopnow("TODO init_leg")
   end subroutine init_leg


   !--------------------------------------------------------------------
   subroutine init_doftyp_leg
   !--------------------------------------------------------------------
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      ! using pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse_leg
      u => unpickle_leg
      call register_doftyp("Leg", typid, p, u)
   end subroutine init_doftyp_leg

end module dvr_leg_m
