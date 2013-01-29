# vim: set ts=8 noexpandtab :

DEP := $(dd)/Deps.mk

SOURCES := \
   base.f90 \
   logging.f90 \
   dof.f90 \
   genpot.f90 \
   linear.f90 \
   modeutil.f90 \
   tree.f90 \
   graphviz.f90 \
   tuckerdecomp.f90 \
   hiertuck.f90 \

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.mod))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))


# module dependencies

$(DEP): $(MODS_core) $(SRC_core)
	$(FC) $(DEPFLAGS) -I$(OBJDIR) -J$(OBJDIR) $(SRC_core) | sed -e "s@^\(\S\)@$(OBJDIR)/\1@" > $@

-include $(DEP)

