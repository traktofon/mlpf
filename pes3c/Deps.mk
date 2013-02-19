obj/pes3cvpd.o: pes3c/pes3cvpd.f
obj/pes3c_m.o obj/pes3c_m.mod: pes3c/pes3c_m.f90 obj/base_m.mod
obj/pes3ctest.o: pes3c/pes3ctest.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/pes3c_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
