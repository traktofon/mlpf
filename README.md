MLPF: Multi-Layer Potfit
========================

This is an auxiliary program for the [Heidelberg MCTDH Package](http://mctdh.uni-hd.de/).
It converts a potential energy surface into MLPF format, which can be used by MCTDH as
a multi-layer operator in ML-MCTDH runs.  Compared to potentials fitted by regular Potfit,
MLPF potentials are much more compact, so that ML-MCTDH will run much more efficiently,
especially if high-accuracy fits are used.

Support for MLPF in MCTDH exists currently in a development branch, and is planned to be
released with version 8.6.


Prerequisites
-------------

If you are using MCTDH, chances are very high that all of the prerequisites are already installed.
For reference, this is the software required to build MLPF:

* Compilers for C and Modern Fortran. Versions of [GCC](http://gcc.gnu.org/) since 4.8 should work fine.
* An installation of BLAS (e.g. [OpenBLAS](http://www.openblas.net/)) and LAPACK.
* GNU Make.
* Perl and Bash.


Compilation
-----------

Before first compiling, copy the file `local.mk.sample` into `local.mk` and edit it, following
the comments inside. This is for specifying the compiler and its options, and how to link with
BLAS/LAPACK on your system.

To compile, simply run `make`. To compile faster, use `make -j N`, where N is the number of
CPU cores you want to use.  This will create a number of executables in the subdirectory
`bin/`, namely:

* **mlpf** -- the MLPF fitting program
* **mlpf2npot** -- a program to convert an MLPF potential into a Potfit potential
* **vinfo** -- a simple test/example program which prints out information about _vpot_ or _vpot2_ files, which are used to store full-grid potentials.


Usage
-----



