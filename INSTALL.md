# 0. Retrieving and installing source

  Change to the directory where you want the repository to be installed, e.g.,

  ~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~

  Clone the `ndsu3lib` repository.

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/ndsu3lib.git
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Supporting libraries

  First, some supporting libraries must be installed, namely:

  LAPACK

  GSL or Fortran version of WIGXJPF

  MPFUN2020-Fort, version 2 (optional)

  Download links and documentation for these libraries can be found here:

  LAPACK: https://www.netlib.org/lapack/#_lapack_version_3_12_1

  GSL: https://www.gnu.org/software/gsl/

  WIGXJPF: http://fy.chalmers.se/subatom/wigxjpf/

  MPFUN2020-Fort: https://www.davidhbailey.com/dhbsoftware/

# 2. Configuration

  In the `ndsu3lib` code directory, configure the project using CMake. If using GSL and only double precision arithmetic, use

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cmake -B build \
      -DBLAS_LIBRARIES=<directory with BLAS and LAPACK library files>/librefblas.a \
      -DGSL_INCLUDE_DIR=<directory with GSL header files> \
      -DGSL_LIBRARY=<directory with GSL library files>/libgsl.a \
      -DGSL_CBLAS_LIBRARY=<directory with GSL library files>/libgslcblas.a
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  On clusters, LAPACK and GSL may by available by loading appropriate modules. Then, these libraries do not need to be installed, and the variables BLAS_LIBRARIES, GSL_INCLUDE_DIR, GSL_LIBRARY, and GSL_CBLAS_LIBRARY above do not need to be defined.

  To enable quadruple precision arithmetic without multiprecision arithmetic, include

  ~~~~~~~~~~~~~~~~
  -DPRECISION=quad
  ~~~~~~~~~~~~~~~~

  To enable multiprecision arithmetic without quadruple precision arithmetic, include

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DPRECISION=multi \
  -DMPFUN20_DIR=<directory containing MPFUN20 *.o and *.mod files>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  To enable both quadruple precision and multiprecision arithmetic (recommended), include
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DPRECISION=multiquad \
  -DMPFUN20_DIR=<directory containing MPFUN20 *.o and *.mod files>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If using WIGXJPF (recommended), include

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DSU2COEF_LIBRARY=wigxjpf \
  -DWIGXJPF_INC_DIR=<directory with fwigxjpf.mod file> \
  -DWIGXJPF_LIB_DIR=<directory with libwigxjpf.a library file>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  For OpenMP multithreaded applications, include
  
  ~~~~~~~~~~~~~~~~
  -DOPENMP=ON
  ~~~~~~~~~~~~~~~~

  If CMake does not find the WIGXJPF or MPFUN2020-Fort library, it can make a shallow git clone and copy the necessary sources to _dep directory in the build directory. The library is then compiled along with ndsu3lib. To enable this, include

  ~~~~~~~~~~~~~~~~
  -DFETCH=ON
  ~~~~~~~~~~~~~~~~

  The git repositories that CMake clones from are set in fetch_declarations.cmake.

# 3. Compilation and installation

  To compile the ndsu3lib library itself along with the example programs:

  ~~~~~~~~~~~~~~~~
  % cmake --build build
  ~~~~~~~~~~~~~~~~

  The example programs can then be run:

  ~~~~~~~~~~~~~~~~
  % ./build/ndsu3lib_example
  % ./build/ndsu3lib_example_cpp
  ~~~~~~~~~~~~~~~~

  The output should match the example_output.txt file.

  To install the ndsu3lib library:

  ~~~~~~~~~~~~~~~~
  % cmake --install build --prefix <prefix>
  ~~~~~~~~~~~~~~~~

# 4. Linking

  A Fortran code program.f90 using ndsu3lib can be compiled by (using GCC and GSL)
 
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % gfortran program.f90 -I <prefix>/include/ndsu3lib/mod \
      -L <prefix>/lib -lndsu3lib \
      -L <directory with BLAS and LAPACK library files> -llapack -lrefblas \
      -L <directory with GSL library files> -lgsl -lgslcblas
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If using WIGXJPF, include

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -L <directory with libwigxjpf.a library file> -lwigxjpf
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If using multiprecision arithmetic (supported by the MPFUN2020-Fort library), include

  ~~~~~~~~~~~~~~~~
  -lmpfun20.
  ~~~~~~~~~~~~~~~~

  A C++ code using ndsu3lib must contain

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #include "<ndsu3lib/ndsu3lib.h"
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  and can be compiled just like the Fortran program above, with

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  g++ -lgfortran -I <prefix>/include
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  instead of

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gfortran -I <prefix>/include/ndsu3lib/mod
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  and including

  ~~~~~~~~~~~~~~~~
  -lquadmath
  ~~~~~~~~~~~~~~~~

  to enable quadruple precision arithmetic.
