!=======================================================================
program mlpf
!=======================================================================

   use tokenize_m
   use parse_run_m
   use parse_pot_m
   use parse_pbasis_m
   use parse_tree_m
   use inp_tree_m
   use tree_m
   use runopts_m
   use meta_dof_m
   use dof_io_m
   use genpot_m
   use hiertuck_m
   use graphviz_m
   use units_m
   use strutil_m
   use fileutil_m
   use base_m
   implicit none

   character(len=c5)        :: inpfile
   type(inp_node_t),pointer :: inptree => null()
   type(tree_t),pointer     :: tree
   type(runopts_t)          :: runopts
   type(dof_tp),pointer     :: dofs(:)  => null()
   type(dof_tp),pointer     :: pdofs(:) => null()
   logical                  :: have_run
   logical                  :: have_pot
   logical                  :: have_pbasis
   logical                  :: have_tree
   real(dbl),pointer        :: v(:)
   integer,pointer          :: vdim(:)
   real(dbl)                :: accerr2,err2limit

   call get_command_argument(1,inpfile)
   if (inpfile == "") then
      call usage
      stop 1
   endif

   ! initialize global data structures
   call init_doftyps
   call init_units

   ! parse the input file
   call runinp

   ! do some elementary checks
   if (.not.have_run) &
      call stopnow("missing RUN-SECTION")
   if (runopts%namedir == NOFILE) &
      call stopnow("no name-directory specified")
   call mkdir(runopts%namedir)

   ! create the DVR definition -> dofs
   call rundvr

   ! create the potential data -> v
   call runpot

   ! create the MLPF tree -> tree
   call runtree
   call dispose_inp_node(inptree)

   ! create an initial potfit
   call runpf

   ! create the multi-layer potfit
   call runmlpf

   ! create the graphviz input file
   call rundot

!=======================================================================
   contains
!=======================================================================

   !--------------------------------------------------------------------
   subroutine usage
   !--------------------------------------------------------------------
      write (*,'(a,2(/a))') &
         "Usage: mlpf [ options ] inputfile",&
         "Options:",&
         "  (currently none)"
   end subroutine usage


   !--------------------------------------------------------------------
   subroutine runinp
   !--------------------------------------------------------------------
      type(tokenizer_t)        :: tkner
      character(len=maxtoklen) :: token

      call tkner%init(trim(inpfile))

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
         
         elseif (parse_pbasis(tkner,pdofs)) then
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

      call tkner%dispose
      
   end subroutine runinp


   !--------------------------------------------------------------------
   subroutine rundvr
   !--------------------------------------------------------------------
      character(len=c5) :: fname
      integer           :: lun,ierr
      real(dbl)         :: fver

      if (runopts%lgendvr) then
         ! build the DVR
         if (.not.have_pbasis) &
            call stopnow("missing PBASIS-SECTION")
         dofs => pdofs

      else
         ! read the DVR
         fname = runopts%dvrfile
         if (fname == NOFILE) &
            return ! DVR info from vpot file will be used
         ! open the DVR file
         open(newunit=lun, file=trim(fname), status="old", form="unformatted", iostat=ierr)
         if (ierr /= 0) &
            call stopnow("cannot open file: "//trim(fname))
         ! read headers
         read(unit=lun,iostat=ierr) fver
         if (ierr /=0 ) &
            call stopnow("cannot read file version: "//trim(fname))
         dofs => rddvrdef(lun,fver)
         close(lun)

      endif
   end subroutine rundvr


   !--------------------------------------------------------------------
   subroutine initdvr(dofs)
   !--------------------------------------------------------------------
      type(dof_tp),pointer :: dofs(:)
      class(dof_t),pointer :: dof
      integer              :: f
      do f=1,size(dofs)
         dof => dofs(f)%p
         call dof%init
      enddo
   end subroutine initdvr


   !--------------------------------------------------------------------
   subroutine runpot
   !--------------------------------------------------------------------
      character(len=c5)    :: fname
      type(dof_tp),pointer :: vdofs(:)

      ! several possibilities:
      ! - have to build potential
      ! - have to load potential (vpot)
      ! - have to load 1st stage potfit (natpot)
      ! determined by input file!

      if (runopts%lgenpot) then
         ! build potential
         if (.not.have_pot) &
            call stopnow("missing POTENTIAL-SECTION")
         if (.not.associated(dofs)) &
            call stopnow("cannot build potential without DVR definition")
         call initdvr(dofs)
         ! TODO:
         ! implement parsing of POT-SECTION
         ! map defined DOFs to the DOFs required by the pot.func.
         ! call buildpot
         call stopnow("TODO buildpot")

      else
         ! read potential
         fname = runopts%potfile
         if (fname == NOFILE) fname = ""
         if (fname == "" .or. endswith(fname,"/")) then
            ! if no filename was specified, set default
            ! if pathname was specified, add default vpot filename
            if (runopts%vpotfmt == 1) then
               fname = trim(fname)//"vpot"
            else
               fname = trim(fname)//"vpot2"
            endif
         endif
         call loadpot(fname, runopts%vpotfmt, vdofs, v)
         ! TODO: check DVR consistency
         if (have_pbasis) &
            write(*,*) 'WARNING: using DVR from vpot file, ignoring PBASIS-SECTION'
         dofs => vdofs

      endif
   end subroutine runpot


   !--------------------------------------------------------------------
   subroutine runtree
   !--------------------------------------------------------------------
      type(node_t),pointer :: topnode
      integer              :: m

      if (.not.have_tree) &
         call stopnow("missing TREE-SECTION")

      topnode => inp2node(inptree,dofs)
      tree => make_tree(topnode)
      do m=1,tree%numleaves
         call init_leaf(tree%leaves(m)%p, dofs)
      enddo

   end subroutine runtree


   !--------------------------------------------------------------------
   subroutine rundot
   !--------------------------------------------------------------------
      integer           :: lun,ierr
      character(len=c5) :: dotfile

      if (runopts%lgendot) then
         dotfile = runopts%dotfile
         if (dotfile == NOFILE)  dotfile="tree.dot"
         dotfile = trim(runopts%namedir) // "/" // trim(dotfile)
         open(newunit=lun, file=trim(dotfile), status="unknown", form="formatted", iostat=ierr)
         if (ierr /= 0) &
            call stopnow("cannot create file: "//trim(dotfile))
         call mkdot(lun, tree, dofs)
         call flush(lun)
         close(lun)
      endif

   end subroutine rundot


   !--------------------------------------------------------------------
   subroutine runpf
   !--------------------------------------------------------------------
      integer   :: vlen
      real(dbl) :: err2

      allocate(vdim(tree%numleaves))
      call leaf_shape(tree,vdim,vlen)
      err2limit = vlen * (runopts%rmse)**2
      call potfit_from_v(tree, v, vdim, err2limit, err2)
      accerr2 = accerr2 + err2

   end subroutine runpf
   

   !--------------------------------------------------------------------
   subroutine runmlpf
   !--------------------------------------------------------------------
      integer   :: vlen
      real(dbl) :: err2

      vlen = product(vdim)
      call compute_ht(tree, v(1:vlen), vdim, err2limit-accerr2, err2)
      accerr2 = accerr2 + err2

   end subroutine runmlpf


end program mlpf

