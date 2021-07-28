Contents
========

1. Purpose
2. License
3. Usage
4. Compilation
5. Contact


1. Purpose
==========

ndsu3lib is a Fortran library for evaluation of SU(3)-SU(2)xU(1) reduced Wigner coefficients, SU(3)-SO(3) reduced Wigner coefficients and SU(3) recoupling coefficients, namely U, Z and 9-(lambda,mu) coefficients.

2. License
==========

See LICENSE.md.

3. Usage
========

The SU(3)-SU(2)xU(1) reduced Wigner coefficients are obtained by calling the subroutine calculate_wigner_canonical.
The SU(3)-SO(3) reduced Wigner coefficients are obtained by calling the subroutine calculate_wigner_su3so3.
The U recoupling coefficients are obtained by calling the subroutine calculate_u_coeff.
The Z recoupling coefficients are obtained by calling the subroutine calculate_z_coeff.
The 9-(lambda,mu) coefficients are obtained by calling the subroutine calculate_9_lambda_mu.

Code of each subroutine is in a separate text file with the same name. At the beginning of each text file there are comments documenting the interface of the subroutine. All the subroutines except the first one require explicit interface. A simple program ndsu3lib_example demonstrating usage of the subroutines is provided.

C++ wrappers of the Fortran subroutines are provided. Prototypes of the C++ functions are in the file c++wrappers.h along with interface documentation.

4. Compilation
==============

First, the following libraries must be installed:
 - Intel Math Kernel Library (unless recoupling coefficients are not going to be calculated)
 - WIGXJPF or GNU Scientific Library (unless only SU(3)-SU(2)xU(1) reduced Wigner coefficients are going to be calculated)
 - optionally MPFUN2020-Fort (recommended for SU(3)-SO(3) reduced Wigner coefficients)
 
Download links and documentation for WIGXJPF and MPFUN2020-Fort can be found here: 

http://fy.chalmers.se/subatom/wigxjpf/

https://www.davidhbailey.com/dhbsoftware/

The Variant 2 of MPFUN2020-Fort is used, therefore use the script ./gnu-complib2.scr (if using gfortan) or ./intel-complib2.scr (if using ifort). This generates the following object files:

mpmodule.o
mpfuna.o
mpfunbq.o
mpfunc.o
mpfund.o
mpfune.o
mpfunf.o
mpfungq2.o
mpfunhq2.o
mpmask13.o
second.o

For SU(3)-SU(2)xU(1) reduced Wigner coefficients calculated by the subroutine calculate_wigner_canonical the following files must be compiled:

binomial_coeff.F90
calculate_wigner_canonical.f90
wigner_canonical_extremal.f90
wigner_canonical.f90
outer_multiplicity.f90

For SU(3)-SU(2)xU(1) reduced Wigner coefficients calculated by the subroutine calculate_wigner_su3so3 the following files must be compiled:

binomial_coeff.F90
I_S_module.F90
calculate_wigner_su3so3.F90
wigner_canonical_extremal.f90
wigner_canonical.f90
outer_multiplicity.f90
wigner_su3so3.F90

5. Contact
==========

Jakub Herko
Department of Physics
University of Notre Dame
Notre Dame, Indiana 46556-5670
USA
e-mail: jherko@nd.edu
