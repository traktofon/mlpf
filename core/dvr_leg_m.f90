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
      integer                        :: m,g,g1,g2,j,j1,jmin,ierr
      real(dbl)                      :: k,fac

      ! call general DVR constructor
      call dof%init_dvr

      !-----------------------------------------------------------------
      ! Compute grid points, DVR matrices.
      !-----------------------------------------------------------------
      m = dof%blz
      if (m<0) &
         call stopnow("init_leg: magnetic quantum number must not be negative")

      if (dof%sym == 0) then ! no symmetry
         !--------------------------------------------------------------
         ! Calculate symmetric tridiagonal matrix representation X_mn of
         ! the cos(theta) operator (note that the diagonal elements are
         ! always zero).
         !--------------------------------------------------------------

         ! set diagonal elements of X_mn to 0
         dof%x = 0.d0
         ! set subdiagonal elements of X_mn
         do g=2,dof%gdim
            j = g-1 + m
            dof%d2mat(g,1) = sqrt((j**2 - m**2)/(4.d0*j**2 - 1.d0))
         enddo
         ! diagonalize cos(theta)-matrix to obtain "ort" and "trafo"
         ! d2mat contains the sub-diagonal and d1mat serves as a work array
         call DSTEQR('I', dof%gdim, dof%x, dof%d2mat, dof%trafo, dof%gdim, dof%d1mat, ierr)
         if (ierr /= 0) &
            call stopnow("error in diagonalising cos(theta)-matrix")
         ! change the sign of the eigenvector if necessary, and thus
         ! ensure that the DVR-weights are all positive
         do g=1,dof%gdim
            if (dof%trafo(1,g) < 0.d0) &
               dof%trafo(:,g) = -dof%trafo(:,g)
         enddo
         ! calculate second derivative matrix in corresponding DVR basis
         dof%d2mat = 0.d0
         do g=1,dof%gdim
            ! -(I/2)*Kinetic-energy in Leg-FBR basis (diagonal)
            k = (g+m)*(g+m - 1)
            ! transform to DVR basis
            do g1=1,dof%gdim
               do g2=1,dof%gdim
                  dof%d2mat(g2,g1) = dof%d2mat(g2,g1) - dof%trafo(g,g2)*k*dof%trafo(g,g1)
               enddo
            enddo
         enddo
         ! calculate DVR matrix representation of of d/dtheta * sin(theta)
         dof%d1mat = 0.d0
         do g=1,dof%gdim-1
            ! super-diagonal matrix element in FBR basis (j is really j+1)
            j = g + m
            k = j*sqrt((j**2 - m**2)/(4.d0*j**2 - 1.d0))
            ! transform from FBR to DVR basis
            do g1=1,dof%gdim
               do g2=1,dof%gdim
                  dof%d1mat(g2,g1) = dof%d1mat(g2,g1) &
                     + dof%trafo(g+1,g2) * dof%trafo(g,g1) * k &
                     - dof%trafo(g,g2) * dof%trafo(g+1,g1) * k
               enddo
            enddo
         enddo
         ! rescale grid points cos(theta) --> theta
         do g=1,dof%gdim
            dof%x(g) = acos(dof%x(g))
         enddo

      else ! symmetry
         !--------------------------------------------------------------
         ! Calculate symmetric tridiagonal matrix representation X_mn of
         ! the cos^2(theta) operator (note that the diagonal elements
         ! are not zero).
         !--------------------------------------------------------------

         ! set minimum value of j
         jmin = m + mod(m + dof%sym, 2)
         ! set diagonal elements for X_mn
         do g=1,dof%gdim
            j = 2*(g-1) + jmin
            j1 = j + 1
            dof%x(g) = (j**2 - m**2)/(4.d0*j**2 - 1.d0) &
                     + (j1**2 - m**2)/(4.d0*j1**2 - 1.d0)
         enddo
         ! set super-diagonal elements for X_mn (j is really j+2)
         do g=2,dof%gdim
            j = 2*(g-1) + jmin
            j1 = j - 1
            dof%d2mat(g-1,1) = sqrt((j**2 - m**2)/(4.d0*j**2 - 1.d0) &
                                  * (j1**2 - m**2)/(4.d0*j1**2 - 1.d0))
         enddo
         ! diagonalize cos^2(theta) matrix to obtain "ort" and "trafo"
         call DSTEQR('I', dof%gdim, dof%x, dof%d2mat, dof%trafo, dof%gdim, dof%d1mat, ierr)
         if (ierr /= 0) &
            call stopnow("error in diagonalising cos^2(theta)-matrix")
         ! change the sign of the eigenvector if necessary, and thus
         ! ensure that the DVR-weights are all positive
         do g=1,dof%gdim
            if (dof%trafo(1,g) < 0.d0) &
               dof%trafo(:,g) = -dof%trafo(:,g)
         enddo
         ! calculate second derivative matrix in corresponding DVR basis
         dof%d2mat = 0.d0
         do g=1,dof%gdim
            ! -(I/2)*Kinetic-energy in Leg-FBR basis (diagonal)
            j = 2*(g-1) + jmin
            k = j*(j+1)
            ! transform to DVR basis
            do g1=1,dof%gdim
               do g2=1,dof%gdim
                  dof%d2mat(g2,g1) = dof%d2mat(g2,g1) - dof%trafo(g,g2)*k*dof%trafo(g,g1)
               enddo
            enddo
         enddo
         ! calculate DVR matrix representation of:
         ! 0.5*( cos(theta)*d/dtheta*sin(theta)  + d/dtheta*sin(theta)*cos(theta) )
         dof%d1mat = 0.d0
         do g=1,dof%gdim-1
            ! super-diagonal matrix element in FBR basis
            j = 2*(g-1) + jmin
            k = (j+1.5d0) * sqrt(((j+1)**2 - m**2)/(4.d0*(j+1)**2 - 1.d0)) &
                          * sqrt(((j+2)**2 - m**2)/(4.d0*(j+2)**2 - 1.d0))
            ! transform from FBR to DVR basis
            do g1=1,dof%gdim
               do g2=1,dof%gdim
                  dof%d1mat(g2,g1) = dof%d1mat(g2,g1) &
                     + dof%trafo(g+1,g2) * dof%trafo(g,g1) * k &
                     - dof%trafo(g,g2) * dof%trafo(g+1,g1) * k
               enddo
            enddo
         enddo
         ! rescale grid points cos^2(theta) --> theta
         do g=1,dof%gdim
            dof%x(g)=acos(sqrt(dof%x(g)))
         enddo

      endif

      !-----------------------------------------------------------------
      ! Compute the weights.
      !-----------------------------------------------------------------
      fac = 1.d0/dble(m + 1)
      do j=1,m+1
         fac = fac*2*j/dble(2*j-1)
      enddo
      fac = sqrt(fac)
      do g=1,dof%gdim
         dof%w(g) = dof%trafo(1,g)*fac
         if (m>0) then
            dof%w(g) = dof%w(g)/sin(dof%x(g))**m
         endif
         if ((dof%sym > 0) .and. (mod(dof%sym + m, 2) == 1)) &
            dof%w(g) = dof%w(g)/(sqrt(dble(2*m+3))*cos(dof%x(g)))
      enddo

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
