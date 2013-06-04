module dof_k_m

   use dof_m
   use tokenize_m
   use base_m
   implicit none

   type,extends(dof_t) :: dof_k_t
      integer :: kmin  ! minimum value of magnetic quantum number
      integer :: kmax  ! maximum value of magnetic quantum number
      integer :: dk    ! increment in k
      contains
      procedure :: init   => init_k
      procedure :: pickle => pickle_k
   end type dof_k_t

   integer,parameter,private :: typid = 8 ! MCTDH basis type

   contains


   !--------------------------------------------------------------------
   subroutine parse_k(dof,tkner)
   !--------------------------------------------------------------------
   ! kdvr :~ INTEGER INTEGER INTEGER?
   !--------------------------------------------------------------------
      class(dof_t),pointer            :: dof
      type(tokenizer_t),intent(inout) :: tkner
      character(len=maxtoklen)        :: token
      allocate(dof_k_t::dof)
      select type (dof)
      type is (dof_k_t)
      dof%kmin = parse_int(tkner)
      dof%kmax = parse_int(tkner)
      token = tkner%get()
      if (token == "(EOL)") then
         dof%dk = 1
      else
         dof%dk = parse_int(tkner)
      endif
      dof%gdim = (dof%kmax-dof%kmin)/dof%dk + 1
      end select
   end subroutine parse_k


   !--------------------------------------------------------------------
   subroutine unpickle_k(dof,gdim,ipar,rpar)
   !--------------------------------------------------------------------
      class(dof_t),pointer :: dof
      integer,intent(in)   :: gdim
      integer,intent(in)   :: ipar(:)
      real(dbl),intent(in) :: rpar(:)
      allocate(dof_k_t::dof)
      dof%gdim  = gdim
      select type(dof)
      type is (dof_k_t)
      dof%kmin = ipar(1)
      dof%kmax = ipar(2)
      dof%dk   = ipar(5)
      end select
   end subroutine unpickle_k


   !--------------------------------------------------------------------
   subroutine pickle_k(dof,id,ipar,rpar)
   !--------------------------------------------------------------------
      class(dof_k_t),intent(inout) :: dof
      integer,intent(out)          :: id
      integer,intent(out)          :: ipar(:)
      real(dbl),intent(out)        :: rpar(:)
      id = typid
      ipar = 0
      rpar = 0.d0
      ipar(1) = dof%kmin
      ipar(2) = dof%kmax
      if (dof%kmin * dof%kmax < 0) then
         ipar(3) = 0
      else
         ipar(3) = min(abs(dof%kmin),abs(dof%kmax))
      endif
      ipar(4) = max(abs(dof%kmin),abs(dof%kmax))
      ipar(5) = dof%dk
   end subroutine pickle_k


   !--------------------------------------------------------------------
   subroutine init_k(dof)
   !--------------------------------------------------------------------
      class(dof_k_t),intent(inout) :: dof
      real(dbl)                    :: dx,w
      integer                      :: g
      if (dof%initialized) return
      dof%initialized = .true.
      ! allocate space for gridpoints and weights
      allocate(dof%x(dof%gdim))
      allocate(dof%w(dof%gdim))
      ! gridpoints are simply equidistant
      ! weights are constant
      dx = dof%dk
      w  = sqrt(dx) ! TODO: check
      do g = 1, dof%gdim
         dof%x(g) = dof%kmin + (g-1)*dof%dk
         dof%w(g) = w
      enddo
   end subroutine init_k


   !--------------------------------------------------------------------
   subroutine init_doftyp_k
   !--------------------------------------------------------------------
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      ! using pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse_k
      u => unpickle_k
      call register_doftyp("k", typid, p, u)
   end subroutine init_doftyp_k

end module dof_k_m
