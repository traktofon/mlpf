# vim: set ts=8 noexpandtab :

TGT := $(dd)/mlpf

SOURCES := \
   mlpf.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_main) $(OBJS_core)
	$(LINK) -o $@ $+ $(LIBS)

