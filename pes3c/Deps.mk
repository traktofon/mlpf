obj/pes3cvpd.o: pes3c/pes3cvpd.f
obj/pes3cmod.o obj/pes3cmod.mod: pes3c/pes3cmod.f90 obj/base.mod
obj/pes3ctest.o: pes3c/pes3ctest.f90 obj/dof.mod obj/tree.mod \
 obj/graphviz.mod obj/genpot.mod obj/pes3cmod.mod obj/tuckerdecomp.mod \
 obj/hiertuck.mod obj/linear.mod
