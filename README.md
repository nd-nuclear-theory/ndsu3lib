ndsu3lib is a Fortran library for calculation of SU(3)-SU(2)xU(1) reduced coupling coefficients, SU(3)-SO(3) reduced coupling coefficients, and SU(3) recoupling coefficients, namely U, Z, and 9-(lambda,mu) coefficients.

A simple program ndsu3lib_example demonstrating usage of the library is provided.

C++ wrappers of Fortran subroutines are provided. The C++ functions are in the file ndsu3lib.h along with interface documentation.

A simple program ndsu3lib_example_cpp demonstrating usage of the C++wrappers is provided.

The following libraries are required:
 - Lapack (unless recoupling coefficients are not going to be calculated)
 - WIGXJPF or GNU Scientific Library (unless only SU(3)-SU(2)xU(1) reduced Wigner coefficients are going to be calculated)
 - optionally MPFUN2020-Fort (recommended for SU(3)-SO(3) reduced coupling coefficients)
 
Download links and documentation for WIGXJPF and MPFUN2020-Fort can be found here: 

http://fy.chalmers.se/subatom/wigxjpf/

https://www.davidhbailey.com/dhbsoftware/

For installation instructions see cmake_install_notes.txt.

Contact:

Jakub Herko
TRIUMF
4004 Wesbrook Mall
Vancouver, British Columbia V6T 2A3
Canada
e-mail: jherko@triumf.ca
