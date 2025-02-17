/*****************************************************

 ndsu3lib.h
//! C++ wrappers for Fortran subroutines and functions

 Jakub Herko
 University of Notre Dame

 SPDX-License-Identifier: MIT

******************************************************/
#ifndef NDSU3LIB_H_
#define NDSU3LIB_H_

namespace ndsu3lib
{
   typedef struct {
      int lambda, mu;
   } SU3Irrep;

   extern "C"
   {
      namespace fortran
      {
	 // for internal use, not intended to be called by the user

         extern int outer_multiplicity(SU3Irrep irrep1, SU3Irrep irrep2, SU3Irrep irrep3);
   
         extern int inner_multiplicity(SU3Irrep irrep, int L);

         extern void initialize_ndsu3lib(bool wso3, int lmpmu);

         extern void initialize_ndsu3lib_thread(bool wso3, int lmpmu);

         extern void finalize_ndsu3lib(bool wso3);

         extern void calculate_coupling_canonical(
             const SU3Irrep* irrep1, const SU3Irrep* irrep2, const SU3Irrep* irrep3,
	     const int* epsilon3, const int* Lambda32, const int* dimpq, const int* dimw, const int* rhomax,
  	     int* numb, double wigner[], int p1a[], int p2a[], int q2a[]
           );

         extern void calculate_u_coef(
             SU3Irrep irrep1, SU3Irrep irrep2, SU3Irrep irrep,
             SU3Irrep irrep3, SU3Irrep irrep12, SU3Irrep irrep23,
             int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd,
             double* rac, int* info
           );

         extern void calculate_z_coef(
             SU3Irrep irrep2, SU3Irrep irrep1, SU3Irrep irrep,
             SU3Irrep irrep3, SU3Irrep irrep12, SU3Irrep irrep13,
             int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd,
             double* Zcoeff, int* info
           );

         extern void calculate_9_lambda_mu(
             SU3Irrep irrep1, SU3Irrep irrep2, SU3Irrep irrep12,
             SU3Irrep irrep3, SU3Irrep irrep4, SU3Irrep irrep34,
             SU3Irrep irrep13, SU3Irrep irrep24, SU3Irrep irrep,
             int rhomax12, int rhomax34, int rhomax1234, int rhomax13, int rhomax24, int rhomax1324,
             double* ninelm, int* info
           );

         extern void calculate_coupling_su3so3(
             SU3Irrep irrep1, int L1, SU3Irrep irrep2, int L2, SU3Irrep irrep3, int L3,
             int kappa1max, int kappa2max, int kappa3max, int rhomax,
             double* wigner
           );

      } // namespace fortran
   }

   inline
   int OuterMultiplicity(const SU3Irrep& irrep1, const SU3Irrep& irrep2, const SU3Irrep& irrep3)
   {
      // Multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
      // Input arguments: irrep1,irrep2,irrep3
      return fortran::outer_multiplicity(irrep1, irrep2, irrep3);
   }

   inline
   int InnerMultiplicity(const SU3Irrep& irrep, const int& L)
   {
      // Inner multiplicity of L within SU(3) irrep (lambda,mu)
      // Input arguments: irrep,L
      return fortran::inner_multiplicity(irrep, L);
   }

   inline
   void InitializeNdsu3lib(const bool& wso3, const int& lmpmu)
   {
      // ndsu3lib initialization subroutine
      // This subroutine must be called by the main program before calling ndsu3lib
      // subroutines for SU(3) coupling or recoupling coefficients.
      //
      // Input arguments: wso3,lmpmu
      //
      // wso3 must be true if SU(3)-SO(3) coupling coefficients are going to be
      // calculated.
      // lmpmu should be greater than or equal to the maximal expected value of lambda+mu.
      fortran::initialize_ndsu3lib(wso3, lmpmu);
   }

   inline
   void InitializeNdsu3libThread(const bool& wso3, const int& lmpmu)
   {
      // ndsu3lib initialization subroutine to be called by each thread
      //
      // Input arguments: wso3,lmpmu
      //
      // wso3 must be true if SU(3)-SO(3) coupling coefficients are going to be
      // calculated.
      // lmpmu should be greater than or equal to the maximal expected value of lambda+mu.
      fortran::initialize_ndsu3lib_thread(wso3, lmpmu);
   }

   inline
   void FinalizeNdsu3lib(const bool& wso3)
   {
      // This subroutine should be called by each thread once SU(3) coupling or
      // recoupling coefficients are not going to be calculated anymore to free memory.
      //
      // Input argument: wso3
      // 
      // wso3 should be true if initialize_ndsu3lib or initialize_ndsu3lib_thread
      // was called with the first argument being true.
      fortran::finalize_ndsu3lib(wso3);
   }

   inline
   void CalculateCouplingCanonical(
       const SU3Irrep& irrep1, const SU3Irrep& irrep2, const SU3Irrep& irrep3,
       const int& epsilon3, const int& Lambda32, const int& dimpq, const int& dimw, const int& rhomax,
       int& numb, double wigner[], int p1a[], int p2a[], int q2a[]
     )
   {
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
      fortran::calculate_coupling_canonical(
          &irrep1, &irrep2, &irrep3,
          &epsilon3, &Lambda32, &dimpq, &dimw, &rhomax,
          &numb, wigner, p1a, p2a, q2a
        );
   }

   inline
   int CalculateUCoef(
       const SU3Irrep& irrep1, const SU3Irrep& irrep2, const SU3Irrep& irrep,
       const SU3Irrep& irrep3, const SU3Irrep& irrep12, const SU3Irrep& irrep23,
       const int& rhomaxa, const int& rhomaxb, const int& rhomaxc, const int& rhomaxd,
       double* rac
     )
   {
      // Calculate SU(3) recoupling U coefficients
      // U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
      // for given lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
      // calling Lapack subroutine dgesv to solve system of linear equations
      // Returns is 0 iff Lapack subroutine dgesv ran withou errors.
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
      int info;
      fortran::calculate_u_coef(
          irrep1, irrep2, irrep,
          irrep3, irrep12, irrep23,
          rhomaxa, rhomaxb, rhomaxc, rhomaxd,
          rac, &info
        );
      return info;
   }

   inline
   int CalculateZCoef(
       const SU3Irrep& irrep2, const SU3Irrep& irrep1, const SU3Irrep& irrep,
       const SU3Irrep& irrep3, const SU3Irrep& irrep12, const SU3Irrep& irrep13,
       const int& rhomaxa, const int& rhomaxb, const int& rhomaxc, const int& rhomaxd,
       double* Zcoeff
     )
   {
      // Calculate SU(3) recoupling Z coefficients
      // Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
      // for given lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
      // calling Lapack subroutine dgesv to solve system of linear equations
      // Returns is 0 iff Lapack subroutine dgesv ran withou errors.
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
      int info;
      fortran::calculate_z_coef(
          irrep2, irrep1, irrep,
          irrep3, irrep12, irrep13,
          rhomaxa, rhomaxb, rhomaxc, rhomaxd,
          Zcoeff, &info
        );
      return info;
   }

   inline
   int Calculate9LambdaMu(
       const SU3Irrep& irrep1, const SU3Irrep& irrep2, const SU3Irrep& irrep12,
       const SU3Irrep& irrep3, const SU3Irrep& irrep4, const SU3Irrep& irrep34,
       const SU3Irrep& irrep13, const SU3Irrep& irrep24, const SU3Irrep& irrep,
       const int& rhomax12, const int& rhomax34, const int& rhomax1234,
       const int& rhomax13, const int& rhomax24, const int& rhomax1324,
       double* ninelm
     )
   {
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
      // Returns is 0 iff Lapack subroutine dgesv ran withou errors.
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
      int info;
      fortran::calculate_9_lambda_mu(
          irrep1, irrep2, irrep12,
          irrep3, irrep4, irrep34,
          irrep13, irrep24, irrep,
          rhomax12, rhomax34, rhomax1234,
          rhomax13, rhomax24, rhomax1324,
          ninelm, &info
        );
      return info;
   }

   inline
   void CalculateCouplingSU3SO3(
       const SU3Irrep& irrep1, const int& L1,
       const SU3Irrep& irrep2, const int& L2,
       const SU3Irrep& irrep3, const int& L3,
       const int& kappa1max, const int& kappa2max, const int& kappa3max, const int& rhomax,
       double* wigner
     )
   {
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
      fortran::calculate_coupling_su3so3(
          irrep1, L1, irrep2, L2, irrep3, L3,
          kappa1max, kappa2max, kappa3max, rhomax,
          wigner
        );
   }

}  // namespace ndsu3lib

#endif  // NDSU3LIB_H_
