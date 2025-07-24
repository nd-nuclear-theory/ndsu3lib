/*****************************************************

 ndsu3lib.h
//! C++ wrappers for Fortran subroutines and functions

 Jakub Herko
 University of Notre Dame

 SPDX-License-Identifier: MIT

******************************************************/
#ifndef NDSU3LIB_H_
#define NDSU3LIB_H_

#include <stdbool.h>

// NDSU3LIB_ALWAYS_INLINE macro for various compilers
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__IBMC__)
  #define NDSU3LIB_ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
  #define NDSU3LIB_ALWAYS_INLINE __forceinline
#elif defined(__NVCC__)
  #define NDSU3LIB_ALWAYS_INLINE __forceinline__ inline
#elif defined(_CRAYC)
  #define NDSU3LIB_ALWAYS_INLINE _Pragma("inline") inline
#else
  #define NDSU3LIB_ALWAYS_INLINE inline
#endif

#ifdef __cplusplus
namespace ndsu3lib
{
   extern "C"
   {
#endif
   typedef struct {
      int lambda, mu;
   } su3irrep;

   extern int outer_multiplicity(su3irrep irrep1, su3irrep irrep2, su3irrep irrep3);
      // Multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
      // Input arguments: irrep1,irrep2,irrep3

   extern int inner_multiplicity(su3irrep irrep, int L);
      // Inner multiplicity of L within SU(3) irrep (lambda,mu)
      // Input arguments: irrep,L

   extern void initialize_ndsu3lib(bool wso3, bool openmp, int lmpmu);
      // ndsu3lib initialization subroutine
      // Must be called by the main program before calling
      // ndsu3lib subroutines for SU(3) coupling or recoupling coefficients.
      // In OpenMP parallelized programs it should be called by each thread.
      //
      // Input arguments: wso3,openmp,lmpmu
      //
      // wso3 must be true if SU(3)-SO(3) coupling coefficients are going to be
      // calculated.
      // openmp must be .TRUE. if OpenMP is used, otherwise it must be .FALSE.
      // lmpmu should be greater than or equal to the maximal expected value of lambda+mu.

   extern void finalize_ndsu3lib(bool wso3);
      // This subroutine should be called by each thread once SU(3) coupling or
      // recoupling coefficients are not going to be calculated anymore to free memory.
      //
      // Input argument: wso3
      //
      // wso3 should be true if initialize_ndsu3lib or initialize_ndsu3lib_thread
      // was called with the first argument being true.

   extern void calculate_coupling_canonical(
         su3irrep irrep1, su3irrep irrep2, su3irrep irrep3,
         int epsilon3, int Lambda32, int dimpq, int dimw, int rhomax,
         int* numb, double wigner[], int p1a[], int p2a[], int q2a[]
      );
      // Calculate SU(3)-SU(2)xU(1) reduced coupling coefficients
      // / (lambda1,mu1)    (lambda2,mu2)   ||  (lambda3,mu3)  \
      // \epsilon1,Lambda1 epsilon2,Lambda2 || epsilon3,Lambda3/rho
      // for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3=Lambda32/2
      //
      // Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,dimpq,dimw,rhomax
      // Output arguments: numb,wigner,p1a,p2a,q2a
      //
      // rmomax is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
      // dimpq is size of arrays p1a,p2a,q2a, which should be at least (max(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
      // dimw is size of array wigner, which should be at least rhomax*dimpq
      //
      // Reduced Wigner coefficient is wigner[ind], where
      //   ind=i+numb*(rho-1)
      //   0<=i<=numb-1
      //   p1=p1a(i)
      //   p2=p2a(i)
      //   q2=q2a(i)
      //   q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
      //   epsilon2=2*lambda2+mu2-3*(p2+q2)
      //   epsilon1=epsilon3-epsilon2
      //   Lambda1=(mu1+p1-q1)/2
      //   Lambda2=(mu2+p2-q2)/2

   extern int calculate_u_coef(
         su3irrep irrep1, su3irrep irrep2, su3irrep irrep,
         su3irrep irrep3, su3irrep irrep12, su3irrep irrep23,
         int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd,
         double rac[]
      );
      // Calculate SU(3) recoupling U coefficients
      // U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
      // for given lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
      // calling Lapack subroutine dgesv to solve system of linear equations
      // Returns is 0 iff Lapack subroutine dgesv ran without errors.
      //
      // Input arguments: irrep1,irrep2,irrep,irrep3,irrep12,irrep23,rhomaxa,rhomaxb,rhomaxc,rhomaxd
      // Output argument: rac
      //
      // rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      // rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      // rhomaxc is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23),
      // rhomaxd is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu),
      // rac(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1)
      //   is U coefficient for given rhoa,rhob,rhoc,rhod,

   extern int calculate_z_coef(
         su3irrep irrep2, su3irrep irrep1, su3irrep irrep,
         su3irrep irrep3, su3irrep irrep12, su3irrep irrep13,
         int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd,
         double Zcoeff[]
      );
      // Calculate SU(3) recoupling Z coefficients
      // Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
      // for given lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
      // calling Lapack subroutine dgesv to solve system of linear equations
      // Returns is 0 iff Lapack subroutine dgesv ran without errors.
      //
      // Input arguments: irrep2,irrep1,irrep,irrep3,irrep12,irrep13,rhomaxa,rhomaxb,rhomaxc,rhomaxd
      // Output argument: Zcoeff
      //
      // rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      // rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      // rhomaxc is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      // rhomaxd is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu),
      // Zcoeff(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1)
      //   is Z coefficient for given rhoa,rhob,rhoc,rhod,

   extern int calculate_9_lambda_mu(
         su3irrep irrep1, su3irrep irrep2, su3irrep irrep12,
         su3irrep irrep3, su3irrep irrep4, su3irrep irrep34,
         su3irrep irrep13, su3irrep irrep24, su3irrep irrep,
         int rhomax12, int rhomax34, int rhomax1234, int rhomax13, int rhomax24, int rhomax1324,
         double ninelm[]
      );
      // Calculate 9-(lambda,mu) coefficients
      //
      // | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
      // | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
      // |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
      // |     rho13           rho24          rho1234            |
      //
      // for given lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,
      // lambda34,mu34,lambda13,mu13,lambda24,mu24,lambda,mu
      // from U and Z coefficients obtained by calling calculate_u_coef and calculate_z_coef
      // Returns is 0 iff Lapack subroutine dgesv ran without errors.
      //
      // Input arguments: irrep1,irrep2,irrep12,irrep3,irrep4,irrep34,irrep13,irrep24,irrep,
      //                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
      // Output argument: ninelm
      //
      // rhomax12 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      // rhomax34 is multiplicity of SU(3) coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34),
      // rhomax1234 is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu),
      // rhomax13 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      // rhomax24 is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24),
      // rhomax1324 is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu),
      // ninelm(ind), where ind=rho12+rhomax12*(rho34-1)+rhomax12*rhomax34*(rho1234-1)+rhomax12*rhomax34*rhomax1234*(rho13-1)
      //   +rhomax12*rhomax34*rhomax1234*rhomax13*(rho24-1)+rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*(rho1324-1)-1,
      //   is 9-(lambda,mu) coefficient for given rho12,rho34,rho1234,rho13,rho24,rho1324,

   extern void calculate_coupling_su3so3(
         su3irrep irrep1, int L1, su3irrep irrep2, int L2, su3irrep irrep3, int L3,
         int kappa1max, int kappa2max, int kappa3max, int rhomax,
         double wigner[]
      );
      // Calculate SU(3)-SO(3) reduced coupling coefficients
      // /(lambda1,mu1) (lambda2,mu2) || (lambda3,mu3)\
      // \  kappa1,L1     kappa2,L2   ||   kappa3,L3  /rho
      // for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3
      //
      // Input arguments: irrep1,L1,irrep2,L2,irrep3,L3,kappa1max,kappa2max,kappa3max,rhomax
      // Output argument: wigner
      //
      // kappa1max is inner multiplicity of L1 within SU(3) irrep (lambda1,mu1),
      // kappa2max is inner multiplicity of L2 within SU(3) irrep (lambda2,mu2),
      // kappa3max is inner multiplicity of L3 within SU(3) irrep (lambda3,mu3),
      // rhomax is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3),
      // wigner(ind), where
      //   ind=kappa1-1+kappa1max*(kappa2-1)+kappa1max*kappa2max*(kappa3-1)+kappa1max*kappa2max*kappa3max*(rho-1),
      //   is reduced coupling coefficient for given kappa1,kappa2,kappa3,rho.

#ifdef __cplusplus
}  // extern "C"

   // C++ pass-by-reference wrappers

   NDSU3LIB_ALWAYS_INLINE
   void calculate_coupling_canonical(
         const su3irrep& irrep1, const su3irrep& irrep2, const su3irrep& irrep3,
         const int& epsilon3, const int& Lambda32, const int& dimpq, const int& dimw, const int& rhomax,
         int& numb, double wigner[], int p1a[], int p2a[], int q2a[]
      )
   {
      calculate_coupling_canonical(
            irrep1, irrep2, irrep3,
            epsilon3, Lambda32, dimpq, dimw, rhomax,
            &numb, wigner, p1a, p2a, q2a
         );
   }

}  // namespace ndsu3lib
#endif  // __cplusplus

#endif  // NDSU3LIB_H_
