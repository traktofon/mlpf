module dvr_exp_m

   use dvr_m
   use dof_m
   use base_m
   implicit none

   type,extends(dvr_t) :: dvr_exp_t
      real(dbl) :: xi  ! first grid point
      real(dbl) :: xf  ! last grid point
      contains
      procedure :: init => init_exp
   end type dvr_exp_t

   contains


   subroutine parse_exp(dof,tokens)
      class(dof_t),pointer        :: dof
      character(len=*),intent(in) :: tokens(:)
      allocate(dvr_exp_t::dof)
      call stopnow("dvr_exp_m::parse not implemented")
   end subroutine parse_exp


   subroutine unpickle_exp(dof,gdim,label,ipar,rpar)
      class(dof_t),pointer        :: dof
      integer,intent(in)          :: gdim
      character(len=*),intent(in) :: label
      integer,intent(in)          :: ipar(:)
      real(dbl),intent(in)        :: rpar(:)
      allocate(dvr_exp_t::dof)
      dof%gdim  = gdim
      dof%label = label
      select type(dof)
      type is (dvr_exp_t)
      dof%xi   = rpar(1)
      dof%xf   = rpar(2)
      end select
   end subroutine unpickle_exp


   subroutine init_exp(dof)
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


   subroutine init_doftyp_exp
      integer                         :: id
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      id = 5 ! MCTDH basis type
      ! uexpg pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse_exp
      u => unpickle_exp
      call register_doftyp("exp", id, p, u)
   end subroutine init_doftyp_exp

end module dvr_exp_m
