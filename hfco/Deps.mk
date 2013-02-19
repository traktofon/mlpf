obj/hfco.o: hfco/hfco.f
obj/hfco_m.o obj/hfco_m.mod: hfco/hfco_m.f90 obj/base_m.mod
obj/hfcotest.o: hfco/hfcotest.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/hfco_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
