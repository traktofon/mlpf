# vim: set ts=8 noexpandtab :

TGT := $(dd)/coul4test

SOURCES := \
   coul4_m.f90 \
   coul4test.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/coul4_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_coul4) $(OBJS_core)
	$(LINK) -o $@ $+ $(LIBS)

