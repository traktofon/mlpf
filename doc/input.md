MLPF Input File Syntax
======================

Input files for MLPF use a syntax similar to MCTDH.
Parameters are organized into sections, which are enclosed in directives like
```
XXX-SECTION
END-XXX-SECTION
```
The file must end with the directive `END-INPUT`.

The following describes the individual sections.


RUN-SECTION
-----------

This section controls parameters for input, output, and what exactly to compute.
Directives are simply keywords, or of the form `keyword = parameter(s)`.
The following table lists the available directives.
_S_ specifies a string, _I_ an integer, and _R_ a real-valued parameter, which might
include a unit (in the same manner as MCTDH, e.g. "1.0,eV"). Square brackets indicate
optional parameters.

| Directive         | Description                                                                                                                 |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------- |
| name = _S_        | Output will be put into the directory _S_.                                                                                  |
| gendvr            | DVR information will be generated and stored it in the file `dvr`.  This requires a PRIMITIVE-BASIS-SECTION to be present   |
| readdvr [= _S_]   | DVR information will be read from the file _S_. If unspecified, _S_ is `dvr`.                                               |
| genpot            | A full-grid potential will be computed and stored, in a format depending on _vpot-format_.                                  |
| readvpot [= _S_]  | A full-grid potential will be read from the file _S_. If unspecified, _S_ is `vpot` or `vpot2`, depending on _vpot-format_. |
| readnpot [= _S_]  | A Potfit potential will be read from the file _S_. If unspecified, _S_ is `natpot`.                                         |
| vpot-format = _I_ | Format for full-grid potential. _I_ must be 1 or 2.                                                                         |
| rmse = _R_        | Set the maximum allowed root-mean-square error of the fitted potential to _R_.                                              |
| graphviz          | Produce a file `tree.dot` which can be used to visualize the multi-layer tree using [Graphviz](http://www.graphviz.org/).   |


PRIMITIVE-BASIS-SECTION
-----------------------
(may be abbreviated to PBASIS-SECTION)

This section defines the coordinate system, and which grid points to use for each coordinate.
Each line describes one primitive coordinate, using the syntax

    Modelabel   BasisType   Parameters...

This is the same format as the corresponding section in MCTDH.  Therefore, please see the
[corresponding MCTDH documentation](http://www.pci.uni-heidelberg.de/tc/usr/mctdh/doc/mctdh/input_docu.html#pbasiskey)
for details of specifying the basis and parameters.

The idea is that you can copy & paste your PBASIS-SECTION from your MCTDH input file.
However, there are currently two caveats:

* MLPF was written from scratch, so the input file parser may not reproduce all MCTDH quirks and bugs accurately.
* At the moment, only a subset of primitive basis types are implemented by MLPF, namely:
    * **HO** -- Harmonic oscillator (Hermite) DVR
    * **Leg** -- Legendre DVR
    * **sin** -- Sine (Chebychev) DVR
    * **FFT** -- Fast Fourier transform collocation
    * **exp** -- Exponential DVR (periodic boundary conditions)
    * **K** -- magnetic quantum number basis


TREE-SECTION
------------

This section defines how the primitive modes are hierarchically combined
into larger and larger modes, and (optionally) specifies how many "single-particle potentials"
should be used for each mode.  So this section fulfills a similar role to the
`ML-BASIS-SECTION` used in MCTDH for ML-MCTDH runs, but the syntax is rather different.

