module dvr_ho_m

   use dvr_m
   use dof_m
   use base
   implicit none

   type,extends(dvr_t) :: dvr_ho_t
      real(dbl) :: xeq ! equilibrium position
      real(dbl) :: fm  ! frequency*mass | xf-xi
      integer   :: typ ! 0=fm 1=xi-xf
      contains
      procedure :: init => init_ho
   end type dvr_ho_t

   contains


   subroutine parse(dof,tokens)
      class(dof_t),pointer        :: dof
      character(len=*),intent(in) :: tokens(:)
      allocate(dvr_ho_t::dof)
      call stopnow("dvr_ho_m::parse not implemented")
   end subroutine parse


   subroutine unpickle(dof,gdim,label,ipar,rpar)
      class(dof_t),pointer        :: dof
      integer,intent(in)          :: gdim
      character(len=*),intent(in) :: label
      integer,intent(in)          :: ipar(:)
      real(dbl),intent(in)        :: rpar(:)
      real(dbl)                   :: xi,xf
      allocate(dvr_ho_t::dof)
      dof%gdim  = gdim
      dof%label = label
      select type(dof)
      type is (dvr_ho_t)
      dof%typ  = ipar(1)
      if (dof%typ == 0) then
         dof%xeq = rpar(1)
         dof%fm  = rpar(2)*rpar(3)
      else
         xi = rpar(1)
         xf = rpar(2)
         dof%xeq = 0.5d0*(xi+xf)
         dof%fm  = xf-xi
      endif
      end select
   end subroutine unpickle


   subroutine init_ho(dof)
      class(dvr_ho_t),intent(inout) :: dof
      integer                       :: g,ierr
      real(dbl)                     :: dx,fac,pi4,ep,w
      ! Call general DVR constructor.
      call dof%init_dvr
      ! Preparations for xi-xf case.
      if (dof%typ == 1) then
         dx = dof%fm
         dof%fm = 1.d0
      endif
      ! Coordinate matrix in FBR is tridiagonal.
      ! Store diagonal in x, subdiagonal in d1mat.
      dof%x = 0.d0
      do g = 1, dof%gdim-1
         dof%d1mat(g,1) = sqrt(g/(2.d0 * dof%fm))
      enddo
      ! Diagonalize it. Use d2mat as workspace.
      call DSTEQR('I', dof%gdim, dof%x, dof%d1mat, dof%trafo, dof%gdim, dof%d2mat, ierr)
      if (ierr /= 0) &
         call stopnow("dvr_ho_m::init_ho : error in diagonalising HO-coordinate-matrix")
      ! Rescale coordinates for xi-xf case.
      if (dof%typ == 1) then
         fac = dx/(dof%x(dof%gdim) - dof%x(1))
         dof%x = fac * dof%x
         dof%fm = fac**(-2)
         dof%typ = 0 ! unset xi-xf flag
      endif
      ! Adjust sign of eigenvectors to ensure positive weights.
      do g = 1, dof%gdim
         if (dof%trafo(1,g) < 0.d0) &
            dof%trafo(:,g) = -dof%trafo(:,g)
      enddo
      ! TODO: set d1mat,d2mat
      dof%d1mat = 0.d0
      dof%d2mat = 0.d0
      ! Compute weights.
      pi4 = PI**0.25d0
      do g = 1, dof%gdim
         ep = 0.5d0 * dof%fm * dof%x(g)**2
         if (ep < 600.d0) then
            w = dof%trafo(1,g) * pi4 * dof%fm**(-0.25d0) * exp(ep)
         else ! fix for large grids
            if (g==1) then
               w = sqrt(dof%x(g+1) - dof%x(g))
            elseif (g==dof%gdim) then
               w = sqrt(dof%x(g) - dof%x(g-1))
            else
               w = sqrt(0.5d0*(dof%x(g+1) - dof%x(g-1)))
            endif
         endif
         dof%w(g) = w
      enddo
      ! Shift coordinates.
      dof%x = dof%x + dof%xeq
   end subroutine init_ho


   subroutine init_doftyp_ho
      integer                         :: id
      procedure(parse_dof),pointer    :: p
      procedure(unpickle_dof),pointer :: u
      id = 1 ! MCTDH basis type
      ! using pointers is not strictly necessary, but this makes it
      ! more likely that the compiler checks for interface mismatch
      p => parse
      u => unpickle
      call register_dof("HO", id, p, u)
   end subroutine init_doftyp_ho

end module dvr_ho_m
