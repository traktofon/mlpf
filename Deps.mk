obj/base_m.o obj/base_m.mod: 
obj/dof_io_m.o obj/dof_io_m.mod:  obj/dof_m.mod \
 obj/base_m.mod
obj/dof_m.o obj/dof_m.mod:  obj/base_m.mod obj/tokenize_m.mod \
 obj/strutil_m.mod
obj/dvr_exp_m.o obj/dvr_exp_m.mod:  obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/dvr_fft_m.o obj/dvr_fft_m.mod:  obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/dvr_exp_m.mod obj/base_m.mod
obj/dvr_ho_m.o obj/dvr_ho_m.mod:  obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/dvr_leg_m.o obj/dvr_leg_m.mod:  obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/base_m.mod
obj/dvr_m.o obj/dvr_m.mod:  obj/dof_m.mod obj/base_m.mod
obj/dvr_sin_m.o obj/dvr_sin_m.mod:  obj/dvr_m.mod \
 obj/dof_m.mod obj/tokenize_m.mod obj/units_m.mod obj/base_m.mod
obj/fileutil_m.o obj/fileutil_m.mod:  obj/base_m.mod
obj/genpot_m.o obj/genpot_m.mod:  obj/base_m.mod \
 obj/logging_m.mod obj/dof_m.mod obj/dof_io_m.mod obj/base_m.mod
obj/graphviz_m.o obj/graphviz_m.mod:  obj/dof_m.mod \
 obj/vtree_m.mod
obj/hiertuck_m.o obj/hiertuck_m.mod:  obj/logging_m.mod \
 obj/dof_m.mod obj/vtree_m.mod obj/tuckerdecomp_m.mod
obj/itree_m.o obj/itree_m.mod:  obj/vtree_m.mod obj/dof_m.mod \
 obj/base_m.mod
obj/linear_m.o obj/linear_m.mod:  obj/base_m.mod
obj/logging_m.o obj/logging_m.mod: 
obj/map_str2dbl_m.o obj/map_str2dbl_m.mod:  \
 obj/base_m.mod obj/strutil_m.mod
obj/map_str2int_m.o obj/map_str2int_m.mod:  \
 obj/strutil_m.mod
obj/meta_dof_m.o obj/meta_dof_m.mod:  obj/dvr_ho_m.mod \
 obj/dvr_leg_m.mod obj/dvr_sin_m.mod obj/dvr_fft_m.mod obj/dvr_exp_m.mod \
 obj/dof_m.mod obj/base_m.mod
obj/mmap_m.o obj/mmap_m.mod:  obj/base_m.mod
obj/modeutil_m.o obj/modeutil_m.mod: 
obj/parse_pbasis_m.o obj/parse_pbasis_m.mod:  \
 obj/tokenize_m.mod obj/units_m.mod obj/dof_m.mod obj/strutil_m.mod \
 obj/base_m.mod
obj/parse_pot_m.o obj/parse_pot_m.mod:  \
 obj/tokenize_m.mod obj/strutil_m.mod
obj/parse_run_m.o obj/parse_run_m.mod:  \
 obj/tokenize_m.mod obj/units_m.mod obj/strutil_m.mod obj/runopts_m.mod
obj/parse_tree_m.o obj/parse_tree_m.mod:  \
 obj/tokenize_m.mod obj/itree_m.mod obj/base_m.mod
obj/runopts_m.o obj/runopts_m.mod:  obj/base_m.mod
obj/strutil_m.o obj/strutil_m.mod: 
obj/tokenize_m.o obj/tokenize_m.mod:  \
 obj/map_str2int_m.mod obj/strutil_m.mod obj/base_m.mod
obj/tuckerdecomp_m.o obj/tuckerdecomp_m.mod:  \
 obj/base_m.mod obj/modeutil_m.mod obj/linear_m.mod
obj/units_m.o obj/units_m.mod:  obj/base_m.mod \
 obj/tokenize_m.mod obj/map_str2dbl_m.mod
obj/vtree_m.o obj/vtree_m.mod:  obj/base_m.mod \
 obj/logging_m.mod
obj/c_io.o: core/c_io.c
obj/mlpf.o:  obj/tokenize_m.mod obj/parse_run_m.mod \
 obj/parse_pot_m.mod obj/parse_pbasis_m.mod obj/parse_tree_m.mod \
 obj/itree_m.mod obj/vtree_m.mod obj/runopts_m.mod obj/meta_dof_m.mod \
 obj/dof_io_m.mod obj/genpot_m.mod obj/hiertuck_m.mod obj/graphviz_m.mod \
 obj/units_m.mod obj/strutil_m.mod obj/fileutil_m.mod obj/base_m.mod
