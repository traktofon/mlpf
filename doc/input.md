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

TREE-SECTION
------------
