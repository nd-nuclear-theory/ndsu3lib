/***************************************************************

 c++wrappers.h -- C++ wrappers for Fortran subroutines and functions

 Jakub Herko
 University of Notre Dame

 SPDX-License-Identifier: MIT

***************************************************************/
namespace ndsu3lib
{
  extern "C"
  {
    typedef struct {
      int lambda, mu;
    } su3irrep;
    namespace fortran
    {
    // Fortran subroutines and functions
      extern int outer_multiplicity(const su3irrep&, const su3irrep&, const su3irrep&);
      extern int inner_multiplicity(const su3irrep&, const int&);
      extern void ndsu3lib_init(const bool&, const int&);
      extern void ndsu3lib_free(const bool&);
      extern void calculate_wigner_canonical(const su3irrep&, const su3irrep&, const su3irrep&, const int&, const int&, const int&, const int&, const int&, int&, double[], int[], int[], int[]);
      extern void u_coeff_wrapper(const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const int&, const int&, const int&, const int&, const int&, double[], int&);
      extern void z_coeff_wrapper(const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const int&, const int&, const int&, const int&, const int&, double[], int&);
      extern void nine_lambda_mu_wrapper(const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const su3irrep&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], int&);
      extern void wigner_su3so3_wrapper(const su3irrep&, const int&, const su3irrep&, const int&, const su3irrep&, const int&, const int&, const int&, const int&, const int&, const int&, double[]);
    }
  }

  int outer_multiplicity(su3irrep irrep1, su3irrep irrep2, su3irrep irrep3);
  // Multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
  // Input arguments: irrep1, irrep2, irrep3

  int inner_multiplicity(su3irrep irrep, int L);
  // Inner multiplicity of L within SU(3) irrep (lambda,mu)
  // Input arguments: irrep,L

  void ndsu3lib_init(bool wso3, int j2max);
  // ndsu3lib initialization subroutine
  // This subroutine must be called by the main program before calling ndsu3lib
  // subroutines for SU(3) Wigner or recoupling coefficients.
  //
  // Input arguments: wso3,j2max
  //
  // wso3 must be true if SU(3)-SO(3) Wigner coefficients are going to be
  // calculated.
  // If WIGXJPF is not going to be utilized, j2max is not used. Otherwise j2max
  // must be greater than or equal to two times the maximal angular momentum
  // expected in ordinary Clebsch-Gordan or SU(2) recoupling coefficients.
  // j2max should be at least the maximal expected value of lambda+mu if
  // SU(3)-SO(3) Wigner coefficients are not going to be calculated. If
  // SU(3)-SO(3) Wigner coefficients are going to be calculated, j2max should be
  // at least two times the maximal expected value of lambda+mu. If this j2max
  // is insufficient, WIGXJPF will terminate the program and display an error
  // message.

  void ndsu3lib_free(bool wso3);
  // This subroutine can be called by the main program once SU(3) Wigner or
  // recoupling coefficients are not going to be calculated anymore to free memory.
  //
  // Input argument: wso3
  //
  // wso3 should be true if ndsu3lib_init was called with the first argument
  // being true.

  void calculate_wigner_canonical(su3irrep irrep1, su3irrep irrep2, su3irrep irrep3, int epsilon3, int Lambda32, int dimpq, int dimw, int rhomax, int& numb, double wigner[], int p1a[], int p2a[], int q2a[], int& info);
  // Calculates all SU(3)-SU(2)xU(1) reduced Wigner coefficients
  // <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
  // for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3=Lambda32/2
  //
  // Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,dimpq,dimw,rhomax
  // Output arguments: numb,wigner,p1a,p2a,q2a,info
  //
  // rmomax is multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
  // dimpq is size of arrays p1a,p2a,q2a, which should be at least (max(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
  // dimw is size of array wigner, which should be at least rhomax*dimpq
  // info = 1 if epsilon3,Lambda32 are invalid, otherwise info = 0
  //
  // <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho=wigner[ind]
  //   where epsilon2=2*lambda2+mu2-3*(p2+q2)
  //         epsilon1=epsilon3-epsilon2
  //         Lambda1=(mu1+p1-q1)/2
  //         Lambda2=(mu2+p2-q2)/2
  //         p1=p1a(i) 
  //         p2=p2a(i)
  //         q2=q2a(i)
  //         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
  //         ind=i+numb*(rho-1)
  //         0<=i<=numb-1

  void calculate_u_coeff(su3irrep irrep1, su3irrep irrep2, su3irrep irrep, su3irrep irrep3, su3irrep irrep12, su3irrep irrep23, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double rac[], int& info);
  // Calculates SU(3) U-recoupling coefficients
  // U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
  //
  // Input arguments: irrep1,irrep2,irrep,irrep3,irrep12,irrep23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
  // Output argumrnts: rac,info
  //
  // rhomaxa is multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
  // rhomaxb is multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
  // rhomaxc is multiplicity of coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23)
  // rhomaxd is multiplicity of coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu)
  // dimen is size of array rac, which must be at least rhomaxa*rhomaxb*rhomaxc*rhomaxd
  // rac(ind) = U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
  // ind = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1)-1
  // info = 0 if MKL subroutine dgesv ran withou errors

  void calculate_z_coeff(su3irrep irrep2, su3irrep irrep1, su3irrep irrep, su3irrep irrep3, su3irrep irrep12, su3irrep irrep13, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double Zcoeff[], int& info);
  // Calculates SU(3) Z-recoupling coefficients
  // Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
  //
  // Input arguments: irrep2,irrep1,irrep,irrep3,irrep12,irrep13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
  // Output argumrnts: Zcoeff,info
  //
  // rhomaxa is multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
  // rhomaxb is multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
  // rhomaxc is multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
  // rhomaxd is multiplicity of coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu)
  // dimen is size of array Zcoeff, which must be at least rhomaxa*rhomaxb*rhomaxc*rhomaxd
  // Zcoeff(ind) = Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
  // ind = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1)-1
  // info = 0 if MKL subroutine dgesv ran withou errors

  void calculate_9_lambda_mu(su3irrep irrep1, su3irrep irrep2, su3irrep irrep12, su3irrep irrep3, su3irrep irrep4, su3irrep irrep34, su3irrep irrep13, su3irrep irrep24, su3irrep irrep, int rhomax12, int rhomax34, int rhomax1234, int rhomax13, int rhomax24, int rhomax1324, int dimen, double ninelm[], int& info);
  // Calculates 9-(lambda,mu) coefficients
  //
  // | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
  // | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
  // |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
  // |     rho13           rho24          rho1234            |
  //
  // Input arguments: irrep1,irrep2,irrep12,irrep3,irrep4,irrep34,irrep13,irrep24,irrep,
  //                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,dimen
  // Output arguments: ninelm,info
  //
  // rhomax12 is multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
  // rhomax34 is multiplicity of coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34)
  // rhomax1234 is multiplicity of coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu)
  // rhomax13 is multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
  // rhomax24 is multiplicity of coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24)
  // rhomax1324 is multiplicity of coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu)
  // dimen is size of array ninelm, which must be at least rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*rhomax1324
  // ninelm(ind) is the 9-(lambda,mu) coefficient for given rho12,rho34,rho1234,rho13,rho24,rho1324
  // ind = rho12+rhomax12*(rho34-1)+rhomax12*rhomax34*(rho1234-1)+rhomax12*rhomax34*rhomax1234*(rho13-1)
  //       +rhomax12*rhomax34*rhomax1234*rhomax13*(rho24-1)+rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*(rho1324-1)-1
  // info = 0 if MKL subroutine dgesv ran without errors

  void calculate_wigner_su3so3(su3irrep irrep1, int L1, su3irrep irrep2, int L2, su3irrep irrep3, int L3, int kappa1max, int kappa2max, int kappa3max, int rhomax, int dimen, double wigner[]);
  // Calculates reduced SU(3)-SO(3) Wigner coefficients
  // <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
  //
  // Input arguments: irrep1,L1,irrep2,L2,irrep3,L3,kappa1max,kappa2max,kappa3max,rhomax,dimen
  // Output argument: wigner
  //
  // kappa1max is inner multiplicity of L1 within (lambda1,mu1)
  // kappa2max is inner multiplicity of L2 within (lambda2,mu2)
  // kappa3max is inner multiplicity of L3 within (lambda3,mu3)
  // rhomax is multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
  // dimen is size of array wigner, which must be at least kappa1max*kaapa2max*kaapa3max*rhomax
  // wigner(ind) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
  // ind = kappa1+kappa1max*(kappa2-1)+kappa1max*kappa2max*(kappa3-1)+kappa1max*kappa2max*kappa3max*(rho-1)-1

}
