/********************************************************************************************************************

 ndsu3lib_example_cpp.cpp
//! Simple program tabulating SU(3) reduced coupling and recoupling coefficients to demonstrate usage of C++ wrappers

 Jakub Herko
 University of Notre Dame

 SPDX-License-Identifier: MIT

********************************************************************************************************************/
#include <iostream>

#include "ndsu3lib.h"

void tabulate_coupling_canonical()
   // Tabulate SU(3)-SU(2)xU(1) reduced coupling coefficients
   // /     (6,1)            (2,1)       ||  (5,2) \
   // \epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho
   {
      ndsu3lib::SU3Irrep irrep1 = {6, 1}, irrep2 = {2, 1}, irrep3 = {5, 2};
      int epsilon3 = -9,
          Lambda3_twice = 5, // Lambda3_twice is 2*Lambda3
          rhomax = ndsu3lib::OuterMultiplicity(irrep1, irrep2, irrep3),
          dimpq = (std::max(irrep1.lambda, irrep1.mu) + 1)*(irrep2.lambda + 1)*(irrep2.mu + 1);
      int dimw = rhomax*dimpq,
          numb;
      int p1a[dimpq], p2a[dimpq], q2a[dimpq];
      double wigner[dimw];
      ndsu3lib::CalculateCouplingCanonical(
          irrep1, irrep2, irrep3, epsilon3, Lambda3_twice, dimpq, dimw, rhomax,
	  numb, wigner, p1a, p2a, q2a
        );
      std::cout << "SU(3)-SU(2)xU(1) reduced coupling coefficients" << std::endl;
      std::cout << "/     (6,1)            (2,1)       ||  (5,2) \\" << std::endl;
      std::cout << "\\epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho" << std::endl;
      std::cout << "epsilon1  2*Lambda1  epsilon2  2*Lambda2  coefficients for rho=1,...,rhomax=" << rhomax << std::endl;
      for (int i = 0; i <= numb - 1; i++)
      {
         int p1 = p1a[i],
             p2 = p2a[i],
             q2 = q2a[i];
         int q1 = (2*(irrep1.lambda + irrep2.lambda) + irrep1.mu + irrep2.mu - epsilon3)/3 - p1 - p2 - q2,
             epsilon2 = 2*irrep2.lambda + irrep2.mu - 3*(p2 + q2);
         int epsilon1 = epsilon3 - epsilon2,
             Lambda1_twice = irrep1.mu + p1 - q1, // Lambda1_twice is 2*Lambda1
             Lambda2_twice = irrep2.mu + p2 - q2; // Lambda2_twice is 2*Lambda2
         std::cout << "   " << epsilon1 << "         " << Lambda1_twice << "         "
                   << epsilon2 << "         " << Lambda2_twice << "      ";
         for (int rho = 1; rho <= rhomax; rho++)
         {
            int ind = i + numb*(rho - 1);
            std::cout << wigner[ind] << "  ";
         }
         std::cout << std::endl;
      }
   }

void tabulate_coupling_su3so3()
   // Tabulate SU(3)-SO(3) reduced coupling coefficients
   // / (6,1)    (2,1)   ||  (5,2)  \
   // \kappa1,2 kappa2,3 || kappa3,3/rho
   {
      ndsu3lib::SU3Irrep irrep1 = {6, 1}, irrep2 = {2, 1}, irrep3 = {5, 2};
      int L1 = 2, L2 = 3, L3 = 3;
      int rhomax = ndsu3lib::OuterMultiplicity(irrep1, irrep2, irrep3),
          kappa1max = ndsu3lib::InnerMultiplicity(irrep1, L1),
          kappa2max = ndsu3lib::InnerMultiplicity(irrep2, L2),
          kappa3max = ndsu3lib::InnerMultiplicity(irrep3, L3);
      int dimen = kappa1max*kappa2max*kappa3max*rhomax;
      double wigner[dimen];
      ndsu3lib::CalculateCouplingSU3SO3(
          irrep1, L1, irrep2, L2, irrep3, L3, kappa1max, kappa2max, kappa3max, rhomax,
	  wigner
        );
      std::cout << std::endl;
      std::cout << "SU(3)-SO(3) reduced coupling coefficients" << std::endl;
      std::cout << "/ (6,1)    (2,1)   ||  (5,2)  \\" << std::endl;
      std::cout << "\\kappa1,2 kappa2,3 || kappa3,3/rho" << std::endl;
      std::cout << "kappa1  kappa2  kappa3  coefficients for rho=1,...,rhomax=" << rhomax << std::endl;
      for (int kappa1 = 1; kappa1 <= kappa1max; kappa1++)
      {
         for (int kappa2 = 1; kappa2 <= kappa2max; kappa2++)
         {
            for (int kappa3 = 1; kappa3 <= kappa3max; kappa3++)
            {
               std::cout << "   " << kappa1 << "       " << kappa2 << "       " << kappa3 << "    ";
               for (int rho = 1; rho <= rhomax; rho++)
               {
                  int ind = kappa1 + kappa1max*(kappa2 - 1) + kappa1max*kappa2max*(kappa3 - 1)
                            + kappa1max*kappa2max*kappa3max*(rho - 1) - 1;
                  std::cout << wigner[ind] << "  ";
               }
               std::cout << std::endl;
            }
         }
      }
   }

void tabulate_u_coef()
   // Tabulate SU(3) recoupling U coefficients
   // U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   {
      ndsu3lib::SU3Irrep irrep1 = {9, 3}, irrep2 = {1, 1}, irrep = {6, 6},
                         irrep3 = {2, 2}, irrep12 = {9, 3}, irrep23 = {3, 3};
      int rhomaxa = ndsu3lib::OuterMultiplicity(irrep1, irrep2, irrep12),
          rhomaxb = ndsu3lib::OuterMultiplicity(irrep12, irrep3, irrep),
          rhomaxc = ndsu3lib::OuterMultiplicity(irrep2, irrep3, irrep23),
          rhomaxd = ndsu3lib::OuterMultiplicity(irrep1, irrep23, irrep);
      int dimen = rhomaxa*rhomaxb*rhomaxc*rhomaxd;
      double rac[dimen];
      int info = ndsu3lib::CalculateUCoef(
          irrep1, irrep2, irrep, irrep3, irrep12, irrep23, rhomaxa, rhomaxb, rhomaxc, rhomaxd,
	  rac
        );
      std::cout << std::endl;
      if (info != 0)
      {
         std::cout << "CalculateUCoef: Lapack subroutine dgesv ran with error: info=" << info << std::endl;
         return;
      }
      std::cout << "SU(3) recoupling coefficients U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]" << std::endl;
      std::cout << "rhoa  rhob  rhoc  rhod  coefficient" << std::endl;
      for (int rhoa = 1; rhoa <= rhomaxa; rhoa++)
      {
         for (int rhob = 1; rhob <= rhomaxb; rhob++)
         {
            for (int rhoc = 1; rhoc <= rhomaxc; rhoc++)
            {
               for (int rhod = 1; rhod <= rhomaxd; rhod++)
               {
                  int ind = rhod - 1 + rhomaxd*(rhoa - 1) + rhomaxd*rhomaxa*(rhob - 1) + rhomaxd*rhomaxa*rhomaxb*(rhoc - 1);
                  std::cout << "  " << rhoa << "     " << rhob << "     " << rhoc
                            << "     " << rhod << "   " << rac[ind] << std::endl;
               }
            }
         }
      }
   }

void tabulate_z_coef()
   // Tabulate SU(3) recoupling Z coefficients
   // Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   {
      ndsu3lib::SU3Irrep irrep2 = {9, 3}, irrep1 = {1, 1}, irrep = {6, 6},
                         irrep3 = {2, 2}, irrep12 = {9, 3}, irrep13 = {3, 3};
      int rhomaxa = ndsu3lib::OuterMultiplicity(irrep1, irrep2, irrep12),
          rhomaxb = ndsu3lib::OuterMultiplicity(irrep12, irrep3, irrep),
          rhomaxc = ndsu3lib::OuterMultiplicity(irrep1, irrep3, irrep13),
          rhomaxd = ndsu3lib::OuterMultiplicity(irrep13, irrep2, irrep);
      int dimen = rhomaxa*rhomaxb*rhomaxc*rhomaxd;
      double Zcoeff[dimen];
      int info = ndsu3lib::CalculateZCoef(
          irrep2, irrep1, irrep, irrep3, irrep12, irrep13, rhomaxa, rhomaxb, rhomaxc, rhomaxd,
	  Zcoeff
        );
      std::cout << std::endl;
      if (info != 0)
      {
         std::cout << "CalculateZCoef: Lapack subroutine dgesv ran with error: info=" << info << std::endl;
         return;
      }
      std::cout << "SU(3) recoupling coefficients Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]" << std::endl;
      std::cout << "rhoa  rhob  rhoc  rhod  coefficient" << std::endl;
      for (int rhoa = 1; rhoa <= rhomaxa; rhoa++)
      {
         for (int rhob = 1; rhob <= rhomaxb; rhob++)
         {
            for (int rhoc = 1; rhoc <= rhomaxc; rhoc++)
            {
               for (int rhod = 1; rhod <= rhomaxd; rhod++)
               {
                  int ind = rhod - 1 + rhomaxd*(rhoa - 1) + rhomaxd*rhomaxa*(rhob - 1) + rhomaxd*rhomaxa*rhomaxb*(rhoc - 1);
                  std::cout << "  " << rhoa << "     " << rhob << "     " << rhoc
                            << "     " << rhod << "   " << Zcoeff[ind] << std::endl;
               }
            }
         }
      }
   }

void tabulate_nine_lm()
   // Tabulate 9-(lambda,mu) coefficients
   // |(1,1) (0,0)  (1,1)  rho12 |
   // |(0,4) (4,2)  (2,4)  rho34 |
   // |(2,3) (4,2)  (2,4) rho1324|
   // |rho13 rho24 rho1234       |
   {
      ndsu3lib::SU3Irrep irrep1 = {1, 1}, irrep2 = {0, 0}, irrep12 = {1, 1},
                         irrep3 = {0, 4}, irrep4 = {4, 2}, irrep34 = {2, 4},
                         irrep13 = {2, 3}, irrep24 = {4, 2}, irrep = {2, 4};
      int rhomax12 = ndsu3lib::OuterMultiplicity(irrep1, irrep2, irrep12),
          rhomax34 = ndsu3lib::OuterMultiplicity(irrep3, irrep4, irrep34),
          rhomax1234 = ndsu3lib::OuterMultiplicity(irrep12, irrep34, irrep),
          rhomax13 = ndsu3lib::OuterMultiplicity(irrep1, irrep3, irrep13),
          rhomax24 = ndsu3lib::OuterMultiplicity(irrep2, irrep4, irrep24),
          rhomax1324 = ndsu3lib::OuterMultiplicity(irrep13, irrep24, irrep);
      int dimen = rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*rhomax1324;
      double ninelm[dimen];
      int info = ndsu3lib::Calculate9LambdaMu(
          irrep1, irrep2, irrep12,
          irrep3, irrep4, irrep34,
          irrep13, irrep24, irrep,
          rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324,
          ninelm
        );
      std::cout << std::endl;
      if (info != 0)
      {
         std::cout << "Calculate9LambdaMu: Lapack subroutine dgesv ran with error: info=" << info << std::endl;
         return;
      }
      std::cout << "                           |(1,1) (0,0)  (1,1)  rho12 |" << std::endl;
      std::cout << "9-(lambda,mu) coefficients |(0,4) (4,2)  (2,4)  rho34 |" << std::endl;
      std::cout << "                           |(2,3) (4,2)  (2,4) rho1324|" << std::endl;
      std::cout << "                           |rho13 rho24 rho1234       |" << std::endl;
      std::cout << "rho12  rho34  rho1234  rho13  rho24  rho1324  coefficient" << std::endl;
      for (int rho12 = 1; rho12 <= rhomax12; rho12++)
      {
         for (int rho34 = 1; rho34 <= rhomax34; rho34++)
         {
            for (int rho1234 = 1; rho1234 <= rhomax1234; rho1234++)
            {
               for (int rho13 = 1; rho13 <= rhomax13; rho13++)
               {
                  for (int rho24 = 1; rho24 <= rhomax24; rho24++)
                  {
                     for (int rho1324 = 1; rho1324 <= rhomax1324; rho1324++)
                     {
                        int ind = rho12 - 1 + rhomax12*(rho34 - 1) + rhomax12*rhomax34*(rho1234 - 1)
                                  + rhomax12*rhomax34*rhomax1234*(rho13 - 1)
                                  + rhomax12*rhomax34*rhomax1234*rhomax13*(rho24 - 1)
                                  + rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*(rho1324 - 1);
                        std::cout << "  " << rho12 << "      " << rho34 << "       " << rho1234
                                  << "       " << rho13 << "      " << rho24 << "       "
                                  << rho1324 << "     " << ninelm[ind] << std::endl;
                     }
                  }
               }
            }
         }
      }
   }

int main()
{
   ndsu3lib::InitializeNdsu3lib(true, false, 50);
   tabulate_coupling_canonical();
   tabulate_coupling_su3so3();
   tabulate_u_coef();
   tabulate_z_coef();
   tabulate_nine_lm();
   ndsu3lib::FinalizeNdsu3lib(true);
   return 0;
}
