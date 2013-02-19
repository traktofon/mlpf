obj/testfunc_m.o obj/testfunc_m.mod: test/testfunc_m.f90 obj/base_m.mod
obj/test.o: test/test.f90 obj/logging_m.mod obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/tuckerdecomp_m.mod \
 obj/hiertuck_m.mod obj/linear_m.mod obj/testfunc_m.mod
