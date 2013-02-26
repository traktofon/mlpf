module dvr_sin_m

   use dvr_m
   use dof_m
   use tokenize_m
   use units_m
   use base_m
   implicit none

   type,extends(dvr_t) :: dvr_sin_t
      real(dbl) :: xi  ! first grid point
      real(dbl) :: xf  ! last grid point
      integer   :: typ ! 0=normal 1=sdq 2=spin
      contains
      procedure :: init => init_sin
   end type dvr_sin_t

   contains


   !--------------------------------------------------------------------
   subroutine parse_sin(dof,tkner)
   !--------------------------------------------------------------------
   ! sindvr :~ INTEGER ( length length ( "short" | "long )? ) | angle ) "sdq"?
   !--------------------------------------------------------------------
      class(dof_t),pointer            :: dof
      type(tokenizer_t),intent(inout) :: tkner
      character(len=maxtoklen)        :: token
      real(dbl) :: r1,dx
      logical :: have2pi
      allocate(dvr_sin_t::dof)
      select type (dof)
      type is (dvr_sin_t)

      dof%gdim = parse_int(tkner)
      dof%typ = 0
      r1 = parse_angle(tkner,have2pi)
      if (have2pi) then
         dx = r1/(dof%gdim + 1)
         dof%xi = dx
         dof%xf = r1-dx
      else
         dof%xi = parse_length(tkner)
         dof%xf = parse_length(tkner)
         token = tkner%get()
         if (token=="short") then
            call tkner%gofwd
         elseif (token=="long") then
            dx = (dof%xf - dof%xi)/(dof%gdim + 1)
            dof%xi = dof%xi + dx
            dof%xf = dof%xf - dx
            call tkner%gofwd
         endif
      endif
      token = tkner%get()
      if (token=="sdq") then
         dof%typ = 1
         call tkner%gofwd
      endif

      end select
   end subroutine parse_sin


   subroutine unpickle_sin(dof,gdim,ipar,rpar)
      class(dof_t),pointer        :: dof
      integer,intent(in)          :: gdim
      integer,intent(in)          :: ipar(:)
      real(dbl),intent(in)        :: rpar(:)
      allocate(dvr_sin_t::dof)
      dof%gdim  = gdim
      select type(dof)
      type is (dvr_sin_t)
      dof%xi   = rpar(1)
      dof%xf   = rpar(2)
      dof%typ  = ipar(1)
      end select
   end subroutine unpickle_sin


   subroutine init_sin(dof)
      class(dvr_sin_t),intent(inout) :: dof
      real(dbl)                      :: dx,w
      integer                        :: g
      ! call general DVR constructor
      call dof%init_dvr
      ! gridpoints are simply equidistant
      ! weights are constant
      dx = (dof%xf - dof%xi) / (dof%gdim - 1)
      w  = sqrt(dx)
      do g = 1, dof%gdim
         dof%x(g) = dof%xi + (g-1)*dx
         dof%w(g) = w 
      enddo
      ! TODO: set trafo,d1mat,d2mat
      dof%trafo = 0.d0
      dof%d1mat = 0.d0
      dof%d2mat = 0.d0
   end subroutine init_sin


   subroutine init_doftyp_sin
      integer                         :: id
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      id = 3 ! MCTDH basis type
      ! using pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse_sin
      u => unpickle_sin
      call register_doftyp("sin", id, p, u)
   end subroutine init_doftyp_sin

end module dvr_sin_m
