program mlpf

   use tokenize_m
   use parse_run_m
   use parse_pot_m
   use parse_pbasis_m
   use parse_tree_m
   use inp_tree_m
   use tree_m
   use runopts_m
   use meta_dof_m
   use hiertuck_m
   use graphviz_m
   use units_m
   use base_m
   implicit none

   character(len=c5)        :: inpfile
   type(inp_node_t),pointer :: inptree
   type(node_t),pointer     :: topnode
   type(tree_t),pointer     :: tree
   type(tokenizer_t)        :: tkner
   character(len=maxtoklen) :: token
   type(runopts_t)          :: runopts
   integer                  :: m,idot
   class(dof_t),pointer     :: dof
   type(dof_tp),allocatable :: dofs(:)
   logical                  :: have_run
   logical                  :: have_pot
   logical                  :: have_pbasis
   logical                  :: have_tree
   integer                  :: f,i


   call get_command_argument(1,inpfile)
   if (inpfile == "") then
      call usage
      stop 1
   endif

   call init_doftyps
   call init_units
   call tkner%init(trim(inpfile))

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

      if (parse_run(tkner,runopts)) then
         if (have_run) call tkner%error("duplicate RUN-SECTION")
         have_run = .true.
   
      elseif (parse_pot(tkner)) then
         if (have_pot) call tkner%error("duplicate POT-SECTION")
         have_pot = .true.
      
      elseif (parse_pbasis(tkner,dofs)) then
         if (have_pbasis) call tkner%error("duplicate PBASIS-SECTION")
         have_pbasis = .true.

      elseif (parse_tree(tkner,inptree)) then
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

   topnode => inp2node(inptree,dofs)
   tree => make_tree(topnode)
   do m=1,tree%numleaves
      call init_leaf(tree%leaves(m)%p, dofs)
   enddo
   call open_logfile(idot,"tree.dot")
   call mkdot(idot,tree,dofs,2)
   call flush(idot)
   call close_logfile(idot)

   call dispose_inp_node(inptree)

   contains

   subroutine usage
      write (*,'(a,4(/a))') &
         "------------------------------------------------------------",&
         "Usage: mlpf [ options ] inputfile",&
         "Options:",&
         "  (none)",&
         "------------------------------------------------------------"
   end subroutine usage

end program mlpf

