obj/base.o obj/base.mod: core/base.f90
obj/logging.o obj/logging.mod: core/logging.f90
obj/dof.o obj/dof.mod: core/dof.f90 obj/base.mod
obj/genpot.o obj/genpot.mod: core/genpot.f90 obj/base.mod obj/logging.mod \
 obj/dof.mod obj/base.mod
obj/linear.o obj/linear.mod: core/linear.f90 obj/base.mod
obj/modeutil.o obj/modeutil.mod: core/modeutil.f90
obj/tree.o obj/tree.mod: core/tree.f90 obj/base.mod obj/logging.mod
obj/graphviz.o obj/graphviz.mod: core/graphviz.f90 obj/dof.mod obj/tree.mod
obj/tuckerdecomp.o obj/tuckerdecomp.mod: core/tuckerdecomp.f90 obj/base.mod \
 obj/modeutil.mod obj/linear.mod
obj/hiertuck.o obj/hiertuck.mod: core/hiertuck.f90 obj/logging.mod \
 obj/dof.mod obj/tree.mod obj/tuckerdecomp.mod
