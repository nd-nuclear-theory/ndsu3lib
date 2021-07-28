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

  int outer_multiplicity(su3irrep irrep1, su3irrep irrep2, su3irrep irrep3)
  {
    int rhomax;
    rhomax = fortran::outer_multiplicity(irrep1, irrep2, irrep3);
    return rhomax;
  }

  int inner_multiplicity(su3irrep irrep, int L)
  {
    int kappamax;
    kappamax = fortran::inner_multiplicity(irrep, L);
    return kappamax;
  }

  void ndsu3lib_init(bool wso3, int j2max)
  {
    fortran::ndsu3lib_init(wso3, j2max);
  }

  void ndsu3lib_free(bool wso3)
  {
    fortran::ndsu3lib_free(wso3);
  }

  void calculate_wigner_canonical(su3irrep irrep1, su3irrep irrep2, su3irrep irrep3, int epsilon3, int Lambda32, int dimpq, int dimw, int rhomax, int& numb, double wigner[], int p1a[], int p2a[], int q2a[], int& info)
  {
    int pq3 = (2*(irrep3.lambda-irrep3.mu)-epsilon3)/3+Lambda32;
    if ((pq3<0)||(pq3>2*irrep3.lambda)||(pq3%2!=0))
    {
      info = 1;
      return;
    }
    pq3 = (2*irrep3.lambda+irrep3.mu-epsilon3)/3+irrep3.mu-Lambda32;
    if ((pq3<0)||(pq3>2*irrep3.mu)||(pq3%2!=0))
    {
      info = 1;
      return;
    }
    info = 0;
    fortran::calculate_wigner_canonical(irrep1, irrep2, irrep3, epsilon3, Lambda32, dimpq, dimw, rhomax, numb, wigner, p1a, p2a, q2a);
  }

  void calculate_u_coeff(su3irrep irrep1, su3irrep irrep2, su3irrep irrep, su3irrep irrep3, su3irrep irrep12, su3irrep irrep23, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double rac[], int& info)
  {
    fortran::u_coeff_wrapper(irrep1, irrep2, irrep, irrep3, irrep12, irrep23, rhomaxa, rhomaxb, rhomaxc, rhomaxd, dimen, rac, info);
  }

  void calculate_z_coeff(su3irrep irrep2, su3irrep irrep1, su3irrep irrep, su3irrep irrep3, su3irrep irrep12, su3irrep irrep13, int rhomaxa, int rhomaxb, int rhomaxc, int rhomaxd, int dimen, double Zcoeff[], int& info)
  {
    fortran::z_coeff_wrapper(irrep2, irrep1, irrep, irrep3, irrep12, irrep13, rhomaxa, rhomaxb, rhomaxc, rhomaxd, dimen, Zcoeff, info);
  }

  void calculate_9_lambda_mu(su3irrep irrep1, su3irrep irrep2, su3irrep irrep12, su3irrep irrep3, su3irrep irrep4, su3irrep irrep34, su3irrep irrep13, su3irrep irrep24, su3irrep irrep, int rhomax12, int rhomax34, int rhomax1234, int rhomax13, int rhomax24, int rhomax1324, int dimen, double ninelm[], int& info)
  {
    fortran::nine_lambda_mu_wrapper(irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep, rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, dimen, ninelm, info);
  }

  void calculate_wigner_su3so3(su3irrep irrep1, int L1, su3irrep irrep2, int L2, su3irrep irrep3, int L3, int kappa1max, int kappa2max, int kappa3max, int rhomax, int dimen, double wigner[])
  {
    fortran::wigner_su3so3_wrapper(irrep1, L1, irrep2, L2, irrep3, L3, kappa1max, kappa2max, kappa3max, rhomax, dimen, wigner);
  }

}
