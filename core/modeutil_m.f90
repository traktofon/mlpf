! vim: set ts=3 sw=3 :
module modeutil_m

   use dof_m
   use vtree_m
   implicit none

   contains

   !--------------------------------------------------------------------
   function map_dofmod(tree) result(map)
   ! Builds a mapping from DOFs to leaves, such that
   ! map(i) = m  where  i \in tree.leaves[m].dofnums
   !--------------------------------------------------------------------
      type(vtree_t),intent(in) :: tree
      integer,pointer          :: map(:)
      integer                  :: ndof,m,f,i
      type(vnode_t),pointer    :: no
      character(len=80)        :: errmsg
      ndof = tree%numdofs
      allocate(map(ndof))
      map = 0
      do m=1,tree%numleaves
         no => tree%leaves(m)%p
         do i=1,no%nmodes
            f = no%dofnums(i)
            map(f) = m
         enddo
      enddo
      ! Check that all DOFs have been mapped.
      do f=1,ndof
         if (map(f)==0) then
            write(errmsg,'(a,i0)') &
               "map_dofmod: failed to find leaf-node for DOF #",f
            call stopnow(errmsg)
         endif
      enddo
   end function map_dofmod


   !--------------------------------------------------------------------
   function map_modes(dofs1,tree1,dofs2,tree2,ierr) result(modmap)
   ! Tries to match the primitive modes (stored in the leaves) of tree1
   ! with those of tree2, such that
   !   modmap(m1) = m2
   ! where the m2-th tree2-leaf corresponds to the m1-th tree1-leaf.
   !--------------------------------------------------------------------
   ! ierr =
   !   0      : everything OK
   !   1      : number of dofs is different
   !   2      : number of modes is different
   !   1000+f : dof f from dofs1 has no corresponding dof in dofs2
   !   2000+f : dof f from dofs1 has a corresponding dof, but the grids don't match
   !   3000+m : mode m from tree1 has no corresponding mode in tree2
   !--------------------------------------------------------------------
      type(dof_tp),intent(inout) :: dofs1(:),dofs2(:)
      type(vtree_t),intent(in)   :: tree1,tree2
      integer,intent(out)        :: ierr
      integer,pointer            :: modmap(:)
      integer                    :: ndof1,ndof2,nmode1,nmode2
      integer                    :: m1,m2,f1,f2,i
      integer                    :: fmap(size(dofs1))
      integer,pointer            :: dofmod2(:)
      type(vnode_t),pointer      :: no1,no2

      ierr = 0
      nullify(modmap)

      ! Check number of DOFs.
      ndof1 = tree1%numdofs
      ndof2 = tree2%numdofs
      if (ndof1 /= ndof2) then
         ierr = 1
         return
      endif

      ! Check number of modes.
      nmode1 = tree1%numleaves
      nmode2 = tree2%numleaves
      if (nmode1 /= nmode2) then
         ierr = 2
         return
      endif

      ! Build the map from dofs1 to dofs2.
      call map_dofs(dofs1,dofs2,fmap,ierr)
      if (ierr /= 0) then
         ierr = 1000 + ierr
         return
      endif

      ! Check the DOFs for compatibility
      do f1=1,ndof1
         f2 = fmap(f1)
         if (.not.are_dofs_compatible(dofs1(f1)%p, dofs2(f2)%p)) then
            ierr = 2000 + f1
            return
         endif
      enddo

      ! Make it easy to find the tree2-leaf that a certain
      ! DOF belongs too.
      dofmod2 => map_dofmod(tree2)

      ! Make space for the result.
      allocate(modmap(nmode1))
      modmap = 0

      ! Go through all leaves in tree1.
      do m1=1,nmode1
         no1 => tree1%leaves(m1)%p
         ! Find the corresponding leaf in tree2.
         ! We do this by looking at the first DOF.
         f1 = no1%dofnums(1)
         f2 = fmap(f1)
         m2 = dofmod2(f2)
         no2 => tree2%leaves(m2)%p
         ! Check that the leaves contain the same DOFs, in the same order.
         if (no1%nmodes /= no2%nmodes) then
            ierr = 3000 + m1
            goto 500
         endif
         do i=1,no1%nmodes
            f1 = no1%dofnums(i)
            f2 = no2%dofnums(i)
            if (f2 /= fmap(f1)) then
               ierr = 3000 + m1
               goto 500
            endif
         enddo
         ! Modes match, store the result.
         modmap(m1) = m2
         ! DEBUG
         write(*,'(a,i0,a,i0)') &
            "node match: tree1@",no1%num," ~= tree2@",no2%num
      enddo

      ! Clean up.
 500  deallocate(dofmod2)
      return

   end function map_modes

end module modeutil_m
