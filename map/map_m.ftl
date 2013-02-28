!=======================================================================
module #mod#_m
!=======================================================================
!
! This module implements a map (aka dictionary, aka associative list).
! 
! While hash tables are usually considered to be the asymptotically best
! implementation for maps, here we use a variant of self-balancing
! binary trees, namely a left-leaning red-black tree (LLRBT), described
! by Robert Sedgewick in 2008:
!   http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf
!
! Supported operations are:
! * map_get: return the value associated with a key
! * map_put: store a key and its associated value
! * map_del: remove a key/value pair
! * map_destroy: delete all items from the map
! * map_iter_start: initialize iteration over all items
! * map_iter_next: get next key/value pair
!
! Some advantages of LLRBTs are:
! - no memory overhead for maps with few items
! - fairly easy to implement (due to existing example code)
! - no tuning required (nor possible)
!
! The Fortran code here is a direct translation of the Java example code
! from Sedgewick's paper. (The 2-3 tree variant is implemented here.)
! Some small helper routines are not given in the paper, and their
! implementation has been gleaned from a C implementation by Lee Stanza:
!   http://www.teachsolaisgames.com/articles/balanced_left_leaning.html
!
! FO 02/2013
!
!=======================================================================

   #use#
   implicit none
   private
   
   public :: map_get, map_put, map_del, map_destroy, map_iter_start, map_iter_next

   type,public :: #mod#_t
      type(mapnode_t),private,pointer   :: root => null()
      type(iterstack_t),private,pointer :: iter => null()
   end type #mod#_t


   interface map_get
      module procedure map1_get
   end interface map_get
   
   interface map_put
      module procedure map1_put
   end interface map_put
   
   interface map_del
      module procedure map1_del
   end interface map_del
   
   interface map_destroy
      module procedure map1_destroy
   end interface map_destroy
   
   interface map_iter_start
      module procedure map1_iter_start
   end interface map_iter_start
   
   interface map_iter_next
      module procedure map1_iter_next
   end interface map_iter_next


   logical,parameter :: BLACK  = .false.
   logical,parameter :: RED    = .true.

   type :: mapnode_t
      #key1# :: key
      #val# :: val
      type(mapnode_t),pointer :: left => null()
      type(mapnode_t),pointer :: right => null()
      logical :: color
   end type mapnode_t

   type :: iterstack_t
      type(mapnode_t),pointer   :: node => null()
      type(iterstack_t),pointer :: next => null()
      type(iterstack_t),pointer :: prev => null()
   end type


   contains


   !--------------------------------------------------------------------
   function map1_get(map,key,val) result(found)
   !--------------------------------------------------------------------
   ! If the map has an item with the given key, return its value in val,
   ! and set found to .true., otherwise found is .false. and val is
   ! undefined.
   ! * runtime ~ O(log n)
   !--------------------------------------------------------------------
      type(#mod#_t),intent(in) :: map
      #key#,intent(in) :: key
      #val#,intent(out) :: val
      logical :: found
      type(mapnode_t),pointer :: node
      integer :: cmp
      node => map%root
      do while (associated(node))
         cmp = #cmp#(key,node%key)
         if (cmp<0) then
            node => node%left
         elseif (cmp>0) then
            node => node%right
         else
            val = node%val
            found = .true.
            return
         endif
      enddo
      found = .false.
   end function map1_get


   !--------------------------------------------------------------------
   subroutine map1_put(map,key,val)
   !--------------------------------------------------------------------
   ! Store the given key/value-pair in the map. If there already is an
   ! item with the given key present, it will be overwritten.
   ! * runtime ~ O(log n)
   !--------------------------------------------------------------------
      type(#mod#_t),intent(inout) :: map
      #key#,intent(in) :: key
      #val#,intent(in) :: val
      call insert(map%root,key,val)
      map%root%color = BLACK
   end subroutine map1_put


   !--------------------------------------------------------------------
   subroutine map1_del(map,key)
   !--------------------------------------------------------------------
   ! Deletes the item with given key from the map.
   ! * runtime ~ O(log n)
   !--------------------------------------------------------------------
      type(#mod#_t),intent(inout) :: map
      #key#,intent(in) :: key
      if (associated(map%root)) then
         call delete(map%root,key)
         if (associated(map%root)) &
            map%root%color = BLACK
      endif
   end subroutine map1_del


   !--------------------------------------------------------------------
   subroutine map1_destroy(map)
   !--------------------------------------------------------------------
   ! Deletes all items from the map.
   ! * runtime ~ Theta(n)
   !--------------------------------------------------------------------
      type(#mod#_t),intent(inout) :: map
      call iter_empty(map%iter)
      call dispose(map%root)
   end subroutine map1_destroy


   !--------------------------------------------------------------------
   subroutine map1_iter_start(map)
   !--------------------------------------------------------------------
   ! Initializes iteration over all items in the map. A following call
   ! to map_iter_next will return the item with the smallest key.
   !--------------------------------------------------------------------
      type(#mod#_t),intent(inout) :: map
      call iter_empty(map%iter)
      call iter_fill(map%iter,map%root)
   end subroutine map1_iter_start


   !--------------------------------------------------------------------
   function map1_iter_next(map,key,val) result(found)
   !--------------------------------------------------------------------
   ! In an iteration over all items (started with map_iter_start),
   ! returns key/value of the next item.  If the iteration is finished,
   ! found will be .false., otherwise .true.
   !--------------------------------------------------------------------
      type(#mod#_t),intent(inout) :: map
      #key#,intent(out) :: key
      #val#,intent(out) :: val
      logical :: found
      type(mapnode_t),pointer :: node
      if (.not.associated(map%iter)) then
         found = .false.
         return
      endif
      node => map%iter%node
      call iter_pop(map%iter)
      call iter_fill(map%iter,node%right)
      key = node%key
      val = node%val
      found = .true.
   end function map1_iter_next


   !====================================================================
   ! End of public methods. Now follows the actual implementation of
   ! the LLRBT. If you want to understand this, see the links cited at
   ! the top.
   !====================================================================


   !--------------------------------------------------------------------
   recursive subroutine inorder(node,func)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      interface
         subroutine func(key,val)
            #use#
            #key#,intent(in) :: key
            #val#,intent(inout) :: val
         end subroutine func
      end interface
      if (associated(node)) then
         call inorder(node%left,func)
         call func(node%key,node%val)
         call inorder(node%right,func)
      endif
   end subroutine inorder


   !--------------------------------------------------------------------
   function new_node(key,val) result(node)
   !--------------------------------------------------------------------
      #key#,intent(in) :: key
      #val#,intent(in) :: val
      type(mapnode_t),pointer :: node
      allocate(node)
      node%key = key
      node%val = val
      node%color = RED
   end function new_node


   !--------------------------------------------------------------------
   subroutine rot_left(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      type(mapnode_t),pointer :: tmp
      tmp => node%right
      node%right => tmp%left
      tmp%left => node
      tmp%color = node%color
      node%color = RED
      node => tmp
   end subroutine rot_left


   !--------------------------------------------------------------------
   subroutine rot_right(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      type(mapnode_t),pointer :: tmp
      tmp => node%left
      node%left => tmp%right
      tmp%right => node
      tmp%color = node%color
      node%color = RED
      node => tmp
   end subroutine rot_right


   !--------------------------------------------------------------------
   subroutine flipcolor(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      node%color = .not. node%color
      node%left%color = .not. node%left%color
      node%right%color = .not. node%right%color
   end subroutine flipcolor


   !--------------------------------------------------------------------
   function is_red(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      logical                 :: is_red
      if (associated(node)) then
         is_red = (node%color .eqv. RED)
      else
         is_red = .false.
      endif
   end function is_red


   !--------------------------------------------------------------------
   recursive subroutine insert(node,key,val)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      #key#,intent(in) :: key
      #val#,intent(in) :: val
      integer :: cmp
      if (.not.associated(node)) then
         node => new_node(key,val)
         return
      endif
      cmp = #cmp#(key,node%key)
      if (cmp<0) then
         call insert(node%left,key,val)
      elseif (cmp>0) then
         call insert(node%right,key,val)
      else
         node%val = val
      endif
      if (is_red(node%right).and..not.is_red(node%left)) then
         call rot_left(node)
      endif
      if (is_red(node%left)) then
         if (is_red(node%left%left)) call rot_right(node)
      endif
      if (is_red(node%left).and.is_red(node%right)) then
         call flipcolor(node)
      endif
   end subroutine insert


   !--------------------------------------------------------------------
   subroutine mv_red_left(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      call flipcolor(node)
      if (is_red(node%right%left)) then
         call rot_right(node%right)
         call rot_left(node)
         call flipcolor(node)
      endif
   end subroutine mv_red_left
   

   !--------------------------------------------------------------------
   subroutine mv_red_right(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      call flipcolor(node)
      if (is_red(node%left%left)) then
         call rot_right(node)
         call flipcolor(node)
      endif
   end subroutine mv_red_right


   !--------------------------------------------------------------------
   subroutine fixup(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      if (is_red(node%right)) call rot_left(node)
      if (is_red(node%left)) then
         if (is_red(node%left%left)) call rot_right(node)
      endif
      if (is_red(node%left).and.is_red(node%right)) then
         call flipcolor(node)
      endif
   end subroutine fixup


   !--------------------------------------------------------------------
   function find_min(node) result(node1)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node,node1
      node1 => node
      do while (associated(node1%left))
         node1 => node1%left
      enddo
   end function find_min


   !--------------------------------------------------------------------
   recursive subroutine del_min(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      if (.not.associated(node%left)) then
         deallocate(node)
         return
      endif
      if (.not.is_red(node%left) .and. .not.is_red(node%left%left)) then
         call mv_red_left(node)
      endif
      call del_min(node%left)
      call fixup(node)
   end subroutine del_min


   !--------------------------------------------------------------------
   recursive subroutine delete(n,key)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: n
      #key#,intent(in) :: key
      type(mapnode_t),pointer :: m
      if (#cmp#(key,n%key)<0) then
         if (.not.is_red(n%left) .and. .not.is_red(n%left%left)) &
            call mv_red_left(n)
         call delete(n%left,key)
      else
         if (is_red(n%left)) call rot_right(n)
         if (#cmp#(key,n%key)==0 .and. .not.associated(n%right)) then
            deallocate(n)
            return
         endif
         if (.not.is_red(n%right) .and. .not.is_red(n%right%left)) &
            call mv_red_right(n)
         if (#cmp#(key,n%key)==0) then
            m => find_min(n%right)
            n%key = m%key
            n%val = m%val
            call del_min(n%right)
         else
            call delete(n%right,key)
         endif
      endif
      call fixup(n)
   end subroutine delete


   !--------------------------------------------------------------------
   recursive subroutine dispose(node)
   !--------------------------------------------------------------------
      type(mapnode_t),pointer :: node
      if (associated(node)) then
         call dispose(node%left)
         call dispose(node%right)
         deallocate(node)
      endif
   end subroutine dispose


   !====================================================================
   ! Now follows the code for iterating over all items in the map.
   ! This could be done much more naturally with callback functions,
   ! but Fortran's lack of closures makes that approach cumbersome to
   ! use. Also, our nodes don't have parent pointers. So we use a simple
   ! stack to keep track of the nodes encountered on the way down.
   !====================================================================


   !--------------------------------------------------------------------
   subroutine iter_pop(iter)
   !--------------------------------------------------------------------
      type(iterstack_t),pointer :: iter
      type(iterstack_t),pointer :: prev
      if (.not.associated(iter)) return
      if (associated(iter%prev)) then
         prev => iter%prev
         nullify(prev%next)
         deallocate(iter)
         iter => prev
      else
         deallocate(iter)
      endif
   end subroutine iter_pop


   !--------------------------------------------------------------------
   subroutine iter_push(iter,node)
   !--------------------------------------------------------------------
      type(iterstack_t),pointer :: iter
      type(mapnode_t),pointer   :: node
      type(iterstack_t),pointer :: new
      allocate(new)
      new%node => node
      if (associated(iter)) then
         new%prev => iter
         iter%next => new
      endif
      iter => new
   end subroutine iter_push


   !--------------------------------------------------------------------
   subroutine iter_fill(iter,node)
   !--------------------------------------------------------------------
      type(iterstack_t),pointer :: iter
      type(mapnode_t),pointer   :: node,node1
      node1 => node
      do while (associated(node1))
         call iter_push(iter,node1)
         node1 => node1%left
      enddo
   end subroutine iter_fill


   !--------------------------------------------------------------------
   subroutine iter_empty(iter)
   !--------------------------------------------------------------------
      type(iterstack_t),pointer :: iter
      do while (associated(iter))
         call iter_pop(iter)
      enddo
   end subroutine iter_empty


end module #mod#_m
! vim: set syntax=fortran ts=3 sw=3 expandtab :
