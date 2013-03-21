# vim: set ts=8 noexpandtab :

FSOURCES := \
   base_m.f90         \
   dof_io_m.f90       \
   dof_m.f90          \
   dvr_exp_m.f90      \
   dvr_ho_m.f90       \
   dvr_leg_m.f90      \
   dvr_m.f90          \
   dvr_sin_m.f90      \
   genpot_m.f90       \
   graphviz_m.f90     \
   hiertuck_m.f90     \
   inp_tree_m.f90     \
   linear_m.f90       \
   logging_m.f90      \
   map_str2dbl_m.f90  \
   map_str2int_m.f90  \
   meta_dof_m.f90     \
   mmap_m.f90         \
   modeutil_m.f90     \
   parse_pbasis_m.f90 \
   parse_pot_m.f90    \
   parse_run_m.f90    \
   parse_tree_m.f90   \
   runopts_m.f90      \
   strutil_m.f90      \
   tokenize_m.f90     \
   tree_m.f90         \
   tuckerdecomp_m.f90 \
   units_m.f90

CSOURCES := \
   c_io.c

SRC_$(dd) := \
   $(addprefix $(dd)/,$(FSOURCES)) \
   $(addprefix $(dd)/,$(CSOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(FSOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(FSOURCES:.f90=.o)) \
   $(addprefix $(OBJDIR)/,$(CSOURCES:.c=.o))

ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))
