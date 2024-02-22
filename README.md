ndsu3lib is a Fortran library for evaluation of SU(3)-SU(2)xU(1) reduced Wigner coefficients, SU(3)-SO(3) reduced Wigner coefficients and SU(3) recoupling coefficients, namely U, Z and 9-(lambda,mu) coefficients.

A simple program ndsu3lib_example demonstrating usage of the library is provided.

C++ wrappers of Fortran subroutines are provided. The C++ functions are in the file ndsu3lib.h along with interface documentation.

The following libraries must be installed:
 - Lapack (unless recoupling coefficients are not going to be calculated)
 - WIGXJPF or GNU Scientific Library (unless only SU(3)-SU(2)xU(1) reduced Wigner coefficients are going to be calculated)
 - optionally MPFUN2020-Fort (recommended for SU(3)-SO(3) reduced Wigner coefficients)
 
Download links and documentation for WIGXJPF and MPFUN2020-Fort can be found here: 

http://fy.chalmers.se/subatom/wigxjpf/

https://www.davidhbailey.com/dhbsoftware/

The Variant 2 of MPFUN2020-Fort is used, therefore use the script gnu-complib2.scr (if using gfortan) or intel-complib2.scr (if using ifort). 

Contact:

Jakub Herko
Department of Physics
University of Notre Dame
Notre Dame, Indiana 46556-5670
USA
e-mail: jherko@nd.edu
