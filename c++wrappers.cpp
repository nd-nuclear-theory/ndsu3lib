/***************************************************************

 c++wrappers.cpp -- C++ wrappers for Fortran subroutines and functions

 Jakub Herko
 University of Notre Dame

 SPDX-License-Identifier: MIT

***************************************************************/

#include <algorithm>
#include "c++wrappers.h"

namespace ndsu3lib
{

  int outer_multiplicity(int lambda1, int mu1, int lambda2, int mu2, int lambda3, int mu3)
  {
    int rhomax;
    rhomax = fortran::outer_multiplicity_(lambda1, mu1, lambda2, mu2, lambda3, mu3);
    return rhomax;
  }

  int inner_multiplicity(int lambda, int mu, int L)
  {
    int kappamax;
    kappamax = fortran::inner_multiplicity_(lambda, mu, L);
    return kappamax;
  }

  void ndsu3lib_init(bool wso3, int j2max)
  {
    fortran::ndsu3lib_init_(wso3, j2max);
  }

  void ndsu3lib_free(bool wso3)
  {
    fortran::ndsu3lib_free_(wso3);
  }

  void calculate_wigner_canonical(int lambda1, int mu1, int lambda2, int mu2, int lambda3, int mu3, int epsilon3, int Lambda32, int dimpq, int dimw, int rhomax, int& numb, double wigner[], int p1a[], int p2a[], int q2a[], int& info)
  {
    int p3, q3;
    p3 = (2*(lambda3-mu3)-epsilon3)/3+Lambda32;
    if ((p3<0)||(p3>2*lambda3)||(2*(p3/2)!=p3))
    {
      info = 1;
      return;
    }
    q3 = (2*lambda3+mu3-epsilon3)/3+mu3-Lambda32;
    if ((q3<0)||(q3>2*mu3)||(2*(q3/2)!=q3))
    {
      info = 1;
      return;
    }
    info = 0;
    fortran::calculate_wigner_canonical_(lambda1, mu1, lambda2, mu2, lambda3, mu3, epsilon3, Lambda32, dimpq, dimw, rhomax, numb, wigner, p1a, p2a, q2a);
  }

  void calculate_U_coeff(int lambda1, int mu1, int lambda2, int mu2, int lambda, int mu, int lambda3, int mu3, int lambda12, int mu12, int lambda23, int mu23, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double rac[], int& info)
  {
    fortran::u_coeff_wrapper_(lambda1, mu1, lambda2, mu2, lambda, mu, lambda3, mu3, lambda12, mu12, lambda23, mu23, rhomaxa, rhomaxb, rhomaxc, rhomaxd, dimen, rac, info);
  }

  void calculate_Z_coeff(int lambda2, int mu2, int lambda1, int mu1, int lambda, int mu, int lambda3, int mu3, int lambda12, int mu12, int lambda13, int mu13, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double Zcoeff[], int& info)
  {
    fortran::z_coeff_wrapper_(lambda2, mu2, lambda1, mu1, lambda, mu, lambda3, mu3, lambda12, mu12, lambda13, mu13, rhomaxa, rhomaxb, rhomaxc, rhomaxd, dimen, Zcoeff, info);
  }

  void calculate_9_lambda_mu(int lambda1, int mu1, int lambda2, int mu2, int lambda12, int mu12, int lambda3, int mu3, int lambda4, int mu4, int lambda34, int mu34, int lambda13, int mu13, int lambda24, int mu24, int lambda, int mu, int rhomax12, int rhomax34, int rhomax1234, int rhomax13, int rhomax24, int rhomax1324, int dimen, double ninelm[], int& info)
  {
    fortran::nine_lambda_mu_wrapper_(lambda1, mu1, lambda2, mu2, lambda12, mu12, lambda3, mu3, lambda4, mu4, lambda34, mu34, lambda13, mu13, lambda24, mu24, lambda, mu, rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, dimen, ninelm, info);
  }

  void calculate_wigner_su3so3(int lambda1, int mu1, int L1, int lambda2, int mu2, int L2, int lambda3, int mu3, int L3, int kappa1max, int kappa2max, int kappa3max, int rhomax, int dimen, double wigner[])
  {
    fortran::wigner_su3so3_wrapper_(lambda1, mu1, L1, lambda2, mu2, L2, lambda3, mu3, L3, kappa1max, kappa2max, kappa3max, rhomax, dimen, wigner);
  }

}
