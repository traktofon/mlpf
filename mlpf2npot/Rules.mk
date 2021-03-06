# vim: set ts=8 noexpandtab :

TGT := $(BINDIR)/mlpf2npot

vpath %.f90 $(dd)

SOURCES := \
   mlpf2npot.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o))

PROGS += $(TGT)
ALLSOURCES += $(SRC_$(dd))

# build rule

$(TGT): $(OBJS_mlpf2npot) $(OBJS_core)
	$(dir_guard)
	$(LINK) -o $@ $+ $(LIBS)

