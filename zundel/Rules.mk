# vim: set ts=8 noexpandtab :

TGT := $(dd)/zundeltest

SOURCES := \
   zundel_m.f90 \
   zundeltest.f90

SRC_$(dd) := \
   $(addprefix $(dd)/,$(SOURCES))

MODS_$(dd) := \
   $(OBJDIR)/zundel_m.mod

OBJS_$(dd) := \
   $(addprefix $(OBJDIR)/,$(SOURCES:.f90=.o)) \
   $(OBJDIR)/h5o2.pes.o

PROGS := $(PROGS) $(TGT)
ALLSOURCES := $(ALLSOURCES) $(SRC_$(dd))

# build rules

$(TGT): $(OBJS_zundel) $(OBJS_COMMON)
	$(FC) -o $@ $+ $(LIBS)

$(OBJDIR)/h5o2.pes.o: h5o2.pes.f90
	$(FC) $(FFLAGS) -ffixed-form -c -o $@ $<

