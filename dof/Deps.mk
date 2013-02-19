obj/dof_m.o obj/dof_m.mod: dof/dof_m.f90 obj/base.mod
obj/dof_io_m.o obj/dof_io_m.mod: dof/dof_io_m.f90 obj/dof_m.mod obj/base.mod
obj/dvr_m.o obj/dvr_m.mod: dof/dvr_m.f90 obj/dof_m.mod obj/base.mod
obj/dvr_ho_m.o obj/dvr_ho_m.mod: dof/dvr_ho_m.f90 obj/dvr_m.mod obj/dof_m.mod \
 obj/base.mod
obj/dvr_sin_m.o obj/dvr_sin_m.mod: dof/dvr_sin_m.f90 obj/dvr_m.mod \
 obj/dof_m.mod obj/base.mod
obj/dvr_exp_m.o obj/dvr_exp_m.mod: dof/dvr_exp_m.f90 obj/dvr_m.mod \
 obj/dof_m.mod obj/base.mod
obj/meta_dof_m.o obj/meta_dof_m.mod: dof/meta_dof_m.f90 obj/dvr_ho_m.mod \
 obj/dvr_sin_m.mod obj/dvr_exp_m.mod
