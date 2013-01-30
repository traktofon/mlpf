obj/hfco.o: hfco/hfco.f
obj/hfcomod.o obj/hfcomod.mod: hfco/hfcomod.f90 obj/base.mod
obj/hfcotest.o: hfco/hfcotest.f90 obj/dof.mod obj/tree.mod obj/graphviz.mod \
 obj/genpot.mod obj/hfcomod.mod obj/tuckerdecomp.mod obj/hiertuck.mod \
 obj/linear.mod
