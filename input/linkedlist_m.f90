module linkedlist_m

   implicit none

   type :: llitem_t
      class(llitem_t),pointer :: prev => null()
      class(llitem_t),pointer :: next => null()
   end type llitem_t

   type :: llist_t
      class(llitem_t),pointer :: first => null()
      class(llitem_t),pointer :: last => null()
      contains
      procedure :: llpush
      procedure :: llpop
   end type llist_t


   type,extends(llitem_t) :: llitem_int_t
      integer :: val
   end type llitem_int_t

   type,extends(llist_t) :: llist_int_t
      contains
      procedure :: push => push_int
      procedure :: pop  => pop_int
   end type llist_int_t

   contains


   subroutine llpush(llist,item)
      class(llist_t),intent(inout) :: llist
      class(llitem_t),pointer      :: item
      if (associated(llist%last)) then
         llist%last%next => item
         item%prev => llist%last
      endif
      llist%last => item
      if (.not.associated(llist%first)) then
         llist%first => item
      endif
   end subroutine llpush


   subroutine llpop(llist,item,lempty)
      class(llist_t),intent(inout) :: llist
      class(llitem_t),pointer      :: item
      logical,intent(out)          :: lempty
      if (associated(llist%last)) then
         item => llist%last
         if (associated(item%prev)) then
            llist%last => item%prev
            nullify(llist%last%next)
         else
            nullify(llist%last)
            nullify(llist%first)
         endif
         nullify(item%prev)
         nullify(item%next)
         lempty = .false.
      else
         nullify(item)
         lempty = .true.
      endif
   end subroutine llpop


   subroutine push_int(llist,val)
      class(llist_int_t),intent(inout) :: llist
      integer,intent(in)               :: val
      class(llitem_t),pointer          :: item
      allocate(llitem_int_t::item)
      select type (item)
      class is (llitem_int_t)
      item%val = val
      end select
      call llist%llpush(item)
   end subroutine push_int


   subroutine pop_int(llist,val,lempty)
      class(llist_int_t),intent(inout) :: llist
      integer,intent(out)              :: val
      logical,intent(out)              :: lempty
      class(llitem_t),pointer          :: item
      call llist%llpop(item,lempty)
      if (.not.lempty) then
         select type (item)
         class is (llitem_int_t)
         val = item%val
         end select
         deallocate(item)
      endif
   end subroutine pop_int

end module linkedlist_m
