# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   dof_m.f90 \
   dof_io_m.f90 \
   dvr_m.f90 \
   dvr_ho_m.f90 \
   dvr_sin_m.f90 \
   dvr_exp_m.f90 \
   meta_dof_m.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(MODS_parse) $(MODS_dof) $(SRC_dof)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -Jtmp $(SRC_dof) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" -e "s@tmp/@obj/@" > $@

-include $(DEP)

