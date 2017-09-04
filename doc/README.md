MLPF: Multi-Layer Potfit
========================

This is an auxiliary program for the [Heidelberg MCTDH Package](http://mctdh.uni-hd.de/).
It converts a potential energy surface into MLPF format, which can be used by MCTDH as
a multi-layer operator in ML-MCTDH runs.  Compared to potentials fitted by regular Potfit,
MLPF potentials are much more compact, so that ML-MCTDH will run much more efficiently,
especially if high-accuracy fits are used.  For more information and details, see
[my paper](http://dx.doi.org/10.1063/1.4856135), also available on [arXiv](http://arxiv.org/abs/1309.5060).

MCTDH version 8.5.7 or higher is required to make use of potentials in MLPF format.

In case of problems with this software, please check the [Bitbucket repository](https://bitbucket.org/frankotto/mlpf)
for a potentially updated version, or to [report issues](https://bitbucket.org/frankotto/mlpf/issues).

If you use MLPF for your research, you are requested to cite the following paper:

> Frank Otto,
> "Multi-Layer Potfit: An Accurate Potential Representation for Efficient High-Dimensional Quantum Dynamics."
> J. Chem. Phys. 140, 014106, (2014)



Documentation
-------------

This is the README for MLPF.  It and the rest of the documentation are stored in the `doc/`
subdirectory in Markdown format.  For convenience, the documentation is also available in
HTML format.  Please point your web browser to the `README.html` inside `doc/`.


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

### mlpf ###

To run the MLPF fitting program, use

    mlpf inputfile

where the syntax of the _inputfile_ is described in the [input file documentation](input.md).
The generated output files are:

* `mlpf.dat` -- a binary file containing the MLPF-format potential, which can
  be read by MCTDH using the `mlpf`-statement in the operator file.
* `log` -- a log file, containing information about estimated errors and
  the number of single-particle potentials for each mode.
* `tree.dot` -- (optional) an input file for Graphviz's dot utility, to
  display the mode hierarchy.


### mlpf2npot ###

To convert a potential in MLPF format into Natpot format, use

    mlpf2npot mlpf.dat modc

where _mlpf.dat_ is the output file of `mlpf`, and _modc_ is an integer which specifies
the number of the contracted mode.  It generates a file called _natpot_ which can be
used by several utilities in the MCTDH suite.


### vinfo ###

To get some information about a full-grid potential, stored in a `vpot` or `vpot2` file, use

    vinfo vpotfile

where _vpotfile_ is the filename of the full-grid potential. It prints out some statistical
information about the potential: maximum, minimum, mean, standard deviation, and the L2 norm.

