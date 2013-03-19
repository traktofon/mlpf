program mlpf

   use tokenize_m
   use units_m
   use parse_run_m
   use parse_pot_m
   use parse_pbasis_m
   use parse_tree_m
   use tree_m
   use hiertuck_m
   use graphviz_m
   implicit none

   type(node_t),pointer     :: topnode
   type(tree_t),pointer     :: tree
   type(tokenizer_t)        :: tkner
   character(len=maxtoklen) :: token
   integer                  :: m,idot
   class(dof_t),pointer     :: dof
   type(dof_tp),allocatable :: dofs(:)
   logical                  :: have_run
   logical                  :: have_pot
   logical                  :: have_pbasis
   logical                  :: have_tree
   integer                  :: f,i


   call init_doftyps
   call init_units
   call tkner%init("test.inp")

   allocate(dofs(0))
   have_run    = .false.
   have_pot    = .false.
   have_pbasis = .false.
   have_tree   = .false.

   do
      call tkner%add_stop("END-INPUT")
      call tkner%add_igno("(EOL)")

      token = tkner%get()
      if (tkner%stopreason() == STOPREASON_STOPWORD) exit

      if (parse_run(tkner)) then
         if (have_run) call tkner%error("duplicate RUN-SECTION")
         have_run = .true.
   
      elseif (parse_pot(tkner)) then
         if (have_pot) call tkner%error("duplicate POT-SECTION")
         have_pot = .true.
      
      elseif (parse_pbasis(tkner,dofs)) then
         if (have_pbasis) call tkner%error("duplicate PBASIS-SECTION")
         have_pbasis = .true.

      elseif (parse_tree(tkner,dofs,topnode)) then
         if (have_tree) call tkner%error("duplicate TREE-SECTION")
         have_tree = .true.

      else
         call tkner%error("expected END-INPUT or start of section")
      endif

      call tkner%clear_igno
      call tkner%clear_stop
   enddo

   if (.not.have_run) &
      call stopnow("missing RUN-SECTION")

   if (.not.have_pot) &
      call stopnow("missing POTENTIAL-SECTION")
   ! several possibilities:
   ! - have to build potential
   ! - have to load potential (vpot)
   ! - have to load 1st stage potfit (natpot)
   ! determined by input file!

   if (.not.have_pbasis) &
      call stopnow("missing PBASIS-SECTION")
   do f=1,size(dofs)
      dof => dofs(f)%p
      call dof%init
      write(*,'(i2,99(2x,f5.2))') f, (dof%x(i), i=1,size(dof%x))
   enddo

   if (.not.have_tree) &
      call stopnow("missing TREE-SECTION")
   tree => make_tree(topnode)
   do m=1,tree%numleaves
      call init_leaf(tree%leaves(m)%p, dofs)
   enddo
   call open_logfile(idot,"tree.dot")
   call mkdot(idot,tree,dofs,2)
   call flush(idot)
   call close_logfile(idot)

end program mlpf
