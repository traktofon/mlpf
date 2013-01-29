obj/zundelmod.o obj/zundelmod.mod: zundel/zundelmod.f90 obj/base.mod
obj/zundeltest.o: zundel/zundeltest.f90 obj/dof.mod obj/tree.mod \
 obj/graphviz.mod obj/genpot.mod obj/zundelmod.mod obj/tuckerdecomp.mod \
 obj/hiertuck.mod obj/linear.mod
