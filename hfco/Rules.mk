# vim: set ts=8 noexpandtab :

TGT := $(dd)/hfcotest

SOURCES90 := \
   hfco_m.f90 \
   hfcotest.f90
SOURCES77 := \
   hfco.f

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES77)) \
   $(addprefix $(dd)/,$(SOURCES90))

MODS_$(dd) := \
   $(OBJDIR)/hfco_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES77:.f=.o)) \
   $(addprefix $(OBJDIR)/,$(SOURCES90:.f90=.o))

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_hfco) $(OBJS_core)
	$(LINK) -o $@ $+ $(LIBS)

