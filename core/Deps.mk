obj/base_m.o obj/base_m.mod: core/base_m.f90
obj/logging_m.o obj/logging_m.mod: core/logging_m.f90
obj/linear_m.o obj/linear_m.mod: core/linear_m.f90 obj/base_m.mod
obj/modeutil_m.o obj/modeutil_m.mod: core/modeutil_m.f90
obj/tree_m.o obj/tree_m.mod: core/tree_m.f90 obj/base_m.mod obj/logging_m.mod
obj/tuckerdecomp_m.o obj/tuckerdecomp_m.mod: core/tuckerdecomp_m.f90 \
 obj/base_m.mod obj/modeutil_m.mod obj/linear_m.mod
