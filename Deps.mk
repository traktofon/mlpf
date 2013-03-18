obj/base_m.o obj/base_m.mod: core/base_m.f90
obj/dof_io_m.o obj/dof_io_m.mod: core/dof_io_m.f90 obj/dof_m.mod \
 obj/base_m.mod
obj/dof_m.o obj/dof_m.mod: core/dof_m.f90 obj/base_m.mod obj/tokenize_m.mod
obj/dvr_exp_m.o obj/dvr_exp_m.mod: core/dvr_exp_m.f90 obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/dvr_ho_m.o obj/dvr_ho_m.mod: core/dvr_ho_m.f90 obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/dvr_m.o obj/dvr_m.mod: core/dvr_m.f90 obj/dof_m.mod obj/base_m.mod
obj/dvr_sin_m.o obj/dvr_sin_m.mod: core/dvr_sin_m.f90 obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/genpot_m.o obj/genpot_m.mod: core/genpot_m.f90 obj/base_m.mod \
 obj/logging_m.mod obj/dof_m.mod obj/dof_io_m.mod obj/mmap_m.mod \
 obj/base_m.mod
obj/graphviz_m.o obj/graphviz_m.mod: core/graphviz_m.f90 obj/dof_m.mod \
 obj/tree_m.mod
obj/hiertuck_m.o obj/hiertuck_m.mod: core/hiertuck_m.f90 obj/logging_m.mod \
 obj/dof_m.mod obj/tree_m.mod obj/tuckerdecomp_m.mod
obj/linear_m.o obj/linear_m.mod: core/linear_m.f90 obj/base_m.mod
obj/logging_m.o obj/logging_m.mod: core/logging_m.f90
obj/map_str2dbl_m.o obj/map_str2dbl_m.mod: core/map_str2dbl_m.f90 \
 obj/base_m.mod obj/strutil_m.mod
obj/map_str2int_m.o obj/map_str2int_m.mod: core/map_str2int_m.f90 \
 obj/strutil_m.mod
obj/meta_dof_m.o obj/meta_dof_m.mod: core/meta_dof_m.f90 obj/dvr_ho_m.mod \
 obj/dvr_sin_m.mod obj/dvr_exp_m.mod obj/dof_m.mod obj/base_m.mod
obj/mmap_m.o obj/mmap_m.mod: core/mmap_m.f90 obj/base_m.mod
obj/modeutil_m.o obj/modeutil_m.mod: core/modeutil_m.f90
obj/parse_pbasis_m.o obj/parse_pbasis_m.mod: core/parse_pbasis_m.f90 \
 obj/tokenize_m.mod obj/units_m.mod obj/dof_m.mod obj/strutil_m.mod \
 obj/base_m.mod
obj/parse_pot_m.o obj/parse_pot_m.mod: core/parse_pot_m.f90 \
 obj/tokenize_m.mod obj/strutil_m.mod
obj/parse_run_m.o obj/parse_run_m.mod: core/parse_run_m.f90 \
 obj/tokenize_m.mod obj/strutil_m.mod
obj/parse_tree_m.o obj/parse_tree_m.mod: core/parse_tree_m.f90 \
 obj/tokenize_m.mod obj/tree_m.mod obj/meta_dof_m.mod obj/dof_m.mod \
 obj/strutil_m.mod obj/base_m.mod
obj/strutil_m.o obj/strutil_m.mod: core/strutil_m.f90
obj/tokenize_m.o obj/tokenize_m.mod: core/tokenize_m.f90 \
 obj/map_str2int_m.mod obj/strutil_m.mod obj/base_m.mod
obj/tree_m.o obj/tree_m.mod: core/tree_m.f90 obj/base_m.mod obj/logging_m.mod
obj/tuckerdecomp_m.o obj/tuckerdecomp_m.mod: core/tuckerdecomp_m.f90 \
 obj/base_m.mod obj/modeutil_m.mod obj/linear_m.mod
obj/units_m.o obj/units_m.mod: core/units_m.f90 obj/base_m.mod \
 obj/tokenize_m.mod obj/map_str2dbl_m.mod
obj/c_io.o: core/c_io.c
obj/testfunc_m.o obj/testfunc_m.mod: test/testfunc_m.f90 obj/base_m.mod
obj/test.o: test/test.f90 obj/logging_m.mod obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/tuckerdecomp_m.mod \
 obj/hiertuck_m.mod obj/linear_m.mod obj/testfunc_m.mod
obj/pes3cvpd.o: pes3c/pes3cvpd.f
obj/pes3c_m.o obj/pes3c_m.mod: pes3c/pes3c_m.f90 obj/base_m.mod
obj/pes3ctest.o: pes3c/pes3ctest.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/pes3c_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
obj/zundel_m.o obj/zundel_m.mod: zundel/zundel_m.f90 obj/base_m.mod
obj/zundeltest.o: zundel/zundeltest.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/zundel_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
obj/hfco.o: hfco/hfco.f
obj/hfco_m.o obj/hfco_m.mod: hfco/hfco_m.f90 obj/base_m.mod
obj/hfcotest.o: hfco/hfcotest.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/hfco_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
obj/coul4_m.o obj/coul4_m.mod: coul4/coul4_m.f90 obj/base_m.mod
obj/coul4test.o: coul4/coul4test.f90 obj/meta_dof_m.mod obj/tree_m.mod \
 obj/graphviz_m.mod obj/genpot_m.mod obj/coul4_m.mod \
 obj/tuckerdecomp_m.mod obj/hiertuck_m.mod obj/linear_m.mod
