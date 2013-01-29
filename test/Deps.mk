obj/testfunc.o obj/testfunc.mod: test/testfunc.f90 obj/base.mod
obj/test.o: test/test.f90 obj/dof.mod obj/tree.mod obj/graphviz.mod \
 obj/genpot.mod obj/tuckerdecomp.mod obj/hiertuck.mod obj/linear.mod \
 obj/testfunc.mod
