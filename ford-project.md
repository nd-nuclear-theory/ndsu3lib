---
src_dir: ./
output_dir: ./doc/ford
project: ndsu3lib
summary: Notre Dame SU(3) Coupling and Recoupling Coefficient Library
author: Jakub Herko
email: jherko@triumf.ca
exclude_dir: doc
exclude_dir: build
graph: true
display: public
         private
         protected
sort: permission-alpha
dbg: true
coloured_edges: true
extra_filetypes: h   //
                 cpp //
---

ndsu3lib is a Fortran library for calculation of SU(3)-SU(2)xU(1) reduced coupling coefficients, SU(3)-SO(3) reduced coupling coefficients, and SU(3) recoupling coefficients, namely U, Z, and 9-(lambda,mu) coefficients.

A simple program ndsu3lib_example demonstrating usage of the library is provided.

C++ wrappers of Fortran subroutines are provided. The C++ functions are in the file ndsu3lib.h along with interface documentation.

A simple program ndsu3lib_example_cpp demonstrating usage of the C++wrappers is provided.

For installation instructions see INSTALL.md.
