obj/genpot_m.o obj/genpot_m.mod: mlpf/genpot_m.f90 obj/base_m.mod \
 obj/logging_m.mod obj/dof_m.mod obj/base_m.mod
obj/graphviz_m.o obj/graphviz_m.mod: mlpf/graphviz_m.f90 obj/dof_m.mod \
 obj/tree_m.mod
obj/hiertuck_m.o obj/hiertuck_m.mod: mlpf/hiertuck_m.f90 obj/logging_m.mod \
 obj/dof_m.mod obj/tree_m.mod obj/tuckerdecomp_m.mod
