module dvr_exp_m

   use dvr_m
   use dof_m
   use tokenize_m
   use units_m
   use base_m
   implicit none

   type,extends(dvr_t) :: dvr_exp_t
      real(dbl) :: xi  ! first grid point
      real(dbl) :: xf  ! last grid point
      contains
      procedure :: init   => init_exp
      procedure :: pickle => pickle_exp
   end type dvr_exp_t

   integer,parameter,private :: typid = 5 ! MCTDH basis type

   contains


   !--------------------------------------------------------------------
   subroutine parse_fft_bounds(tkner,gdim,xi,xf)
   !--------------------------------------------------------------------
   ! fftdvr :~ INTEGER ( length length | angle ) ( "linear" | "periodic" | "s-periodic" )?
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: tkner
      integer,intent(out)             :: gdim
      real(dbl),intent(out)           :: xi,xf
      character(len=maxtoklen)        :: token
      logical                         :: have2pi
      real(dbl)                       :: r1,dx
      integer                         :: ptyp
      gdim = parse_int(tkner)
      r1 = parse_angle(tkner,have2pi)
      if (have2pi) then
         xi = 0.d0
         xf = r1
         ptyp = 1 ! default periodic
      else
         xi = parse_length(tkner)
         xf = parse_length(tkner)
         ptyp = 0 ! default linear
      endif
      token = tkner%get()
      if (token=="linear") then
         ptyp = 0
         call tkner%gofwd
      elseif (token=="periodic") then
         ptyp = 1
         call tkner%gofwd
      elseif (token=="s-periodic") then
         ptyp = 2
         call tkner%gofwd
      endif
      if (ptyp==1) then
         dx = (xf-xi)/gdim
         xf = xf-dx ! shift last gridpoint
      elseif (ptyp==2) then
         dx = (xf-xi)/(2*gdim)
         xi = xi+dx ! shift first and
         xf = xf-dx ! last gridpoints
      endif
   end subroutine parse_fft_bounds


   !--------------------------------------------------------------------
   subroutine parse_exp(dof,tkner)
   !--------------------------------------------------------------------
      class(dof_t),pointer            :: dof
      type(tokenizer_t),intent(inout) :: tkner
      allocate(dvr_exp_t::dof)
      select type (dof)
      type is (dvr_exp_t)
      call parse_fft_bounds(tkner, dof%gdim, dof%xi, dof%xf)
      end select
   end subroutine parse_exp


   !--------------------------------------------------------------------
   subroutine unpickle_exp(dof,gdim,ipar,rpar)
   !--------------------------------------------------------------------
      class(dof_t),pointer        :: dof
      integer,intent(in)          :: gdim
      integer,intent(in)          :: ipar(:)
      real(dbl),intent(in)        :: rpar(:)
      allocate(dvr_exp_t::dof)
      dof%gdim  = gdim
      select type(dof)
      type is (dvr_exp_t)
      dof%xi   = rpar(1)
      dof%xf   = rpar(2)
      end select
   end subroutine unpickle_exp


   !--------------------------------------------------------------------
   subroutine pickle_exp(dof,id,ipar,rpar)
   !--------------------------------------------------------------------
      class(dvr_exp_t),intent(inout) :: dof
      integer,intent(out)            :: id
      integer,intent(out)            :: ipar(:)
      real(dbl),intent(out)          :: rpar(:)
      id = typid
      ipar = 0
      rpar = 0.d0
      rpar(1) = dof%xi
      rpar(2) = dof%xf
   end subroutine pickle_exp


   !--------------------------------------------------------------------
   subroutine init_exp(dof)
   !--------------------------------------------------------------------
      class(dvr_exp_t),intent(inout) :: dof
      real(dbl)                      :: dx,w
      integer                        :: g
      ! check
      if (mod(dof%gdim,2) == 0) &
         call stopnow("dvr_exp_m::init_exp : number of grid points must be odd for exp-DVR")
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
   end subroutine init_exp


   !--------------------------------------------------------------------
   subroutine init_doftyp_exp
   !--------------------------------------------------------------------
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      ! using pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse_exp
      u => unpickle_exp
      call register_doftyp("exp", typid, p, u)
   end subroutine init_doftyp_exp

end module dvr_exp_m
