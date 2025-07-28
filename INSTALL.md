# 0. Retrieving and installing source

  Change to the directory where you want the repository to be installed, e.g.,

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Clone the `ndsu3lib` repository:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/ndsu3lib.git
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Supporting libraries

  First, some supporting libraries must be installed, if they are not already
  available on your system, namely:

  - LAPACK (or a vendor optimized implementation)

  - WIGXJPF (installed for use from Fortran) or GNU Scientific Library (GSL)

  - MPFUN20-Fort (version 2) (optional, used only for SU(3)-SO(3) reduced
    coupling coefficients)

  Download links and documentation for these libraries can be found here:

  - LAPACK: https://www.netlib.org/lapack/#_lapack_version_3_12_1

  - WIGXJPF: http://fy.chalmers.se/subatom/wigxjpf/

  - GSL: https://www.gnu.org/software/gsl/

  - MPFUN20-Fort: https://www.davidhbailey.com/dhbsoftware/

  On clusters, LAPACK and GSL may by available by loading appropriate modules.

# 2. Configuration

  In the `ndsu3lib` code directory, configure the project using CMake.  The
  basic command is

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cmake -B <build_dir>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Here `<build_dir>` is whatever subdirectory name you choose for CMake to use
  as the build directory, e.g., `build`.
  
  This will default to using GSL for angular momentum coefficients, and only
  enabling double precision arithmetic.  Other choices are documented below.
  
  This basic form of the configuration command also assumes that any installed
  libraries are already in the compiler library search path, either because they
  were loaded as a module or since they are installed under one of the standard
  installation prefixes on your system (e.g., in `/usr`).  If you have installed
  these libraries to custom locations, these needs to be specified via
  additional flags to CMake, also as documented below.
  
## 2.1. Configuration mode flags
  
  To use WIGXJPF (recommended), rather than GSL, for angular momentum
  coefficients, specify

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DSU2COEF_LIBRARY=wigxjpf
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  For SU(3)-SO(3) reduced coupling coefficients, floating-point calculations can
  be done in either double or quadruple or arbitrary precision (multiprecision).
  If you do not need to compute these coefficients, you do not need the
  MPFUN20-Fort library and you do not need to enable the quadruple precision or
  multiprecision arithmetic.

  To enable both quadruple precision and multiprecision arithmetic
  (recommended for SU(3)-SO(3) reduced coupling coefficients), specify
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DPRECISION=multiquad
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  However, if you wish to enable just quadruple precision arithmetic (without
  multiprecision arithmetic), specify

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DPRECISION=quad
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Or, to enable just multiprecision arithmetic (without quadruple precision
  arithmetic), include

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DPRECISION=multi
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  For OpenMP multithreaded applications, include
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DOPENMP=ON
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2.2. Library path configuration variables

  If you have installed any of the aforementioned libraries to a custom
  location, you can use configuration variables to inform CMake of this
  location.
  
  CMake has a built-in capability to find various vendor optimized LAPACK
  implemenatations (recommended).  If, however, you have installed the reference
  implementation of LAPACK, a custom path to the installation can be specified:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DBLAS_LIBRARIES=<directory with BLAS/LAPACK library files>/librefblas.a
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  See also: https://cmake.org/cmake/help/latest/module/FindBLAS.html
  
  To specify a custom location for WIGXJPF:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DWIGXJPF_INC_DIR=<directory with fwigxjpf.mod file>
  -DWIGXJPF_LIB_DIR=<directory with libwigxjpf.a library file>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  To specify a custom location for GSL:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DGSL_ROOT_DIR=<installation prefix for GSL>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  See also: https://cmake.org/cmake/help/latest/module/FindGSL.html
  
  To specify a custom location for MPFUN20-Fort:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -DMPFUN20_DIR=<directory containing MPFUN20 *.o and *.mod files>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3. Compilation and installation

  To compile the `ndsu3lib` library itself along with the example programs:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cmake --build <build_dir>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  The example programs can then be run:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % ./build/ndsu3lib_example
  % ./build/ndsu3lib_example_c
  % ./build/ndsu3lib_example_cpp
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  These are the example programs implemented in Fortran, C, and C++, respectively.
  The output of the Fortran code should match the `example_output.txt` file.
  The output of the C/C++ example program should be the same, to within formatting
  differences.

  To install the `ndsu3lib` library:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cmake --install <build_dir> --prefix <prefix>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Here `<prefix>` specifies the installation prefix (i.e., the parent directory
  to the `include` and `lib` subdirectories, in which the Fortran module files,
  C++ header files, library binary files, and CMake configuration files will be
  installed).
  
# 4. Linking to `ndsu3lib`

  We provide here simple examples illustrating compiling and linking Fortran or
  C++ programs which use the `ndsu3lib` library.  For purposes of illustration,
  we assume the GNU Compiler Collection (GCC) compilers are being used.

# 4.1. Linking from a Fortran program

  We start from the basic compilation command for a program `program.f90`:
 
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % gfortran program.f90 -o program -lndsu3lib
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Unless you have installed `ndsu3lib` globally, you will also need to specify
  the location of the module and library files:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -I <prefix>/include/ndsu3lib/mod -L <prefix>/lib
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
  You will need to provide flags for the BLAS and LAPACK libraries.  See the
  documentation for your vendor optimized implementation (recommended).
  However, if you built `ndsu3lib` to use the reference implementation of
  LAPACK:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -llapack -lrefblas
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  If you built `ndsu3lib` to use WIGXJPF for angular momentum coefficients:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -lwigxjpf
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If you built `ndsu3lib` to use GSL for angular momentum coefficients:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -lgsl -lgslcblas
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If you built `ndsu3lib` to use multiprecision arithmetic (supported by the
  MPFUN20-Fort library):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -lmpfun20
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  If you built `ndsu3lib` to use OpenMP, linkage to OpenMP libraries will also
  need to be specified:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -fopenmp
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Note that the `-fopenmp` flag here is specific to GCC.  See the relevant
  documentation for other compilers.

  If you have not installed these libraries globally, or in locations specified
  through the `LIBRARY_PATH` environment variable (e.g., after loading a
  module), you will also need to specify their locations with an `-L` flag.
  However, the MPFUN20 library file is generated by CMake when it builds
  `ndsu3lib`, and is installed alongside the `ndsu3lib` library file, so you do
  not need to independently specify a location for it.

# 4.2. Linking from a C++ program

  A C++ code using `ndsu3lib` must contain

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #include <ndsu3lib/ndsu3lib.h>
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  The basic compilation command for a program `program.cpp` is

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % g++ program.cpp -o program -lndsu3lib -lgfortran
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Note that the `gfortran` library name here is specific to GCC.  See the
  relevant documentation for other compilers.
  
  Unless you have installed `ndsu3lib` globally, you will also need to specify
  the location of the header and library files:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -I <prefix>/include -L <prefix>/lib
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  If you built `ndsu3lib` to use quadruple precision arithmetic:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -lquadmath
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Note that the `quadmath` library name here is specific to GCC.  See the
  relevant documentation for other compilers.

  Then you will need to provide flags for linkage to external libraries as
  described above for linking from a Fortran program.
