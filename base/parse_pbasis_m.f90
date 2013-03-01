module parse_pbasis_m

   use tokenize_m
   use units_m
   use dof_m
   use strutil_m
   use base_m
   implicit none

   integer,parameter,private :: maxdofs   = 100

   contains


   !--------------------------------------------------------------------
   subroutine parse_pbasis_line(tkner,dof)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: tkner
      class(dof_t),pointer            :: dof
      character(len=maxtoklen)        :: token
      character(len=c1)               :: mlabel,blabel
      type(doftyp_t),pointer          :: btyp
      procedure(parse_dof),pointer    :: bparse

      nullify(dof)
      call tkner%add_igno("(EOL)")

      token = tkner%get()
      if (token == "(END)") return
      if (len_trim(token) > c1) &
         call tkner%error("modelabel too long")
      mlabel = trim(token)
      call tkner%gofwd
      call tkner%clear_igno

      token = tkner%get()
      blabel = trim(token)
      btyp => find_doftyp_by_name(blabel)
      if (.not.associated(btyp)) &
         call tkner%error("unknown basis type")
      call tkner%gofwd

      call tkner%clear_igno
      bparse => btyp%parse
      call bparse(dof,tkner)
      call dof%set_label(mlabel)

      token = tkner%get()
      if (token /= "(EOL)") &
         call tkner%error("unexpected extra arguments")

   end subroutine parse_pbasis_line



   !--------------------------------------------------------------------
   function parse_pbasis(tkner,dofs) result(flag)
   !--------------------------------------------------------------------
      type(tokenizer_t),intent(inout) :: tkner
      type(dof_tp),allocatable        :: dofs(:)
      logical                         :: flag
      character(len=maxtoklen)        :: token
      class(dof_t),pointer            :: dof
      integer                         :: ndof
      type(dof_tp)                    :: dofs1(maxdofs)

      token = tkner%get()
      call ucase(token)
      if (token /= "PBASIS-SECTION" .and. &
          token /= "PRIMITIVE-BASIS-SECTION") then
         flag = .false.
         return
      endif
      call tkner%clear_stop
      call tkner%add_stop("END-PBASIS-SECTION")
      call tkner%add_stop("END-PRIMITIVE-BASIS-SECTION")
      call tkner%gofwd

      ndof = 0
      do
         call parse_pbasis_line(tkner,dof)
         if (.not.associated(dof)) exit
         if (ndof==maxdofs) &
            call stopnow("parse_pbasis_m::maxdofs too small")
         ndof = ndof+1
         dofs1(ndof)%p => dof
      enddo
      if (allocated(dofs)) deallocate(dofs)
      allocate(dofs(ndof))
      dofs(:) = dofs1(1:ndof)
      flag = .true.

      if (tkner%stopreason() /= STOPREASON_STOPWORD) &
         call tkner%error("expected end of section")
      call tkner%clear_stop
      call tkner%gofwd

   end function parse_pbasis


end module parse_pbasis_m
