Instructions for installing ndsu3lib with cmake. 

Requires cmake version 3.15 or higher 

The following libraries are required:
 - Lapack (unless recoupling coefficients are not going to be calculated)
 - WIGXJPF or GNU Scientific Library (unless only SU(3)-U(1)xSU(2) reduced Wigner coefficients are going to be calculated)
 - optionally MPFUN2020-Fort (recommended for SU(3)-SO(3) reduced coupling coefficients)

Download links and documentation for WIGXJPF and MPFUN2020-Fort can be found here:

http://fy.chalmers.se/subatom/wigxjpf/

https://www.davidhbailey.com/dhbsoftware/

To do a basic installation (with default values for compile flags) run  
	cmake -B <build-dir> 

Then to compile 
	cmake --build <build-dir>

or, if you get impatient
	cmake --build <build-dir> -j<N>

To turn on debug mode.  Include -DCMAKE_BUILD_TYPE=Debug when running cmake -B.  I.e., 
	cmake -B <build-dir> -DCMAKE_BUILD_TYPE=Debug

To build with debug off, include -DCMAKE_BUILD_TYPE=Release

If install at NERSC, you should include -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment.

To choose which SO(3) coefficient library, set -DSO3COEF_LIBRARY to either gsl or wigxjpf, 
or run ccmake <build-dir> and select the appropriate option (arrow down to SO3COEFF_LIBRARY and hit enter until it shows desired libary). Default is -DSO3COEF_LIBRARY=gsl

To choose precision, set -DNDSU3LIB_PRECISION to double, quad, multi or multiquad. Or run ccmake and select desired precision.  

Default values are 
	NDSU3LIB_PRECISION=double
	SO3COEF_LIBRARY=gsl	

For example: 

If you're compiling with the gnu compiler and run

	cmake -B build -DSO3COEF_LIBRARY=wigxjpf -DNDSU3LIB_PRECISION=quad
	cmake --build build

Then ndsu3lib will compile with 
	-DNDSU3LIB_QUAD_GNU -DNDSU3LIB_RACAH_WIGXJPF -DNDSU3LIB_WSO3_WIGXJPF


If -DSO3COEF_LIBRARY=wigxjpf or -DNDSU3LIB_PRECISION=multi or -DNDSU3LIB_PRECISION=multiquad then cmake will look for
libraries wigxjpf and mpfun20-fort, respectively.  The library is found, then it is simply included.  If the library is not found, then cmake will make a shallow git clone and copy the necessary sources to _dep dir in <build-dir>. The library is then compiled along with ndsu3lib. 

The git repository that cmake clones from is set in retrieve.cmake.  
