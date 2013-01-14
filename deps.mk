base.o base.mod: base.f90
dof.o dof.mod: dof.f90 base.mod
genpot.o genpot.mod: genpot.f90 base.mod dof.mod base.mod
graphviz.o graphviz.mod: graphviz.f90 dof.mod tree.mod
linear.o linear.mod: linear.f90 base.mod
modeutil.o modeutil.mod: modeutil.f90
tree.o tree.mod: tree.f90 base.mod
tuckerdecomp.o tuckerdecomp.mod: tuckerdecomp.f90 base.mod modeutil.mod \
 linear.mod
hiertuck.o hiertuck.mod: hiertuck.f90 dof.mod tree.mod tuckerdecomp.mod
testfunc.o testfunc.mod: testfunc.f90 base.mod
