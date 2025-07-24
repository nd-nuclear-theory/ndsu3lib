/********************************************************************************************************************

 ndsu3lib_example_cpp.cpp
//! Simple program tabulating SU(3) reduced coupling and recoupling coefficients to demonstrate usage of C wrappers

 Jakub Herko
 University of Notre Dame
 Patrick J. Fasano
 University of Notre Dame, Argonne National Laboratory

 SPDX-License-Identifier: MIT

********************************************************************************************************************/
#include <stdio.h>

#include "ndsu3lib.h"

#define max(a,b)      (((a) > (b)) ? (a) : (b))

void tabulate_coupling_canonical()
   // Tabulate SU(3)-SU(2)xU(1) reduced coupling coefficients
   // /     (6,1)            (2,1)       ||  (5,2) \
   // \epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho
   {
      su3irrep irrep1 = {6,1}, irrep2 = {2,1}, irrep3 = {5,2};
      int epsilon3 = -9,
          Lambda3_twice = 5, // Lambda3_twice is 2*Lambda3
          rhomax = outer_multiplicity(irrep1, irrep2, irrep3),
          dimpq = (max(irrep1.lambda,irrep1.mu) + 1)*(irrep2.lambda + 1)*(irrep2.mu + 1);
      int dimw = rhomax*dimpq,
          numb;
      int p1a[dimpq], p2a[dimpq], q2a[dimpq];
      double coupling[dimw];
      calculate_coupling_canonical(
          irrep1, irrep2, irrep3, epsilon3, Lambda3_twice,
          dimpq, dimw, rhomax, &numb, coupling, p1a, p2a, q2a
        );
      printf("SU(3)-SU(2)xU(1) reduced coupling coefficients\n");
      printf("/     (6,1)            (2,1)       ||  (5,2) \\\n");
      printf("\\epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho\n");
      printf("epsilon1  2*Lambda1  epsilon2  2*Lambda2  coefficients for rho=1,...,rhomax=%d\n",rhomax);
      for (int i = 0; i < numb; i++)
      {
         int p1 = p1a[i],
            p2 = p2a[i],
            q2 = q2a[i];
         int q1 = (2*(irrep1.lambda + irrep2.lambda) + irrep1.mu + irrep2.mu - epsilon3)/3 - p1 - p2 - q2,
            epsilon2 = 2*irrep2.lambda + irrep2.mu - 3*(p2 + q2);
         int epsilon1 = epsilon3 - epsilon2,
            Lambda1_twice = irrep1.mu + p1 - q1, // Lambda1_twice is 2*Lambda1
            Lambda2_twice = irrep2.mu + p2 - q2; // Lambda2_twice is 2*Lambda2
         printf(
            "   %d         %d         %d         %d      ",
            epsilon1, Lambda1_twice, epsilon2, Lambda2_twice
            );
         for (int rho = 1; rho <= rhomax; rho++)
         {
            int ind = i + numb*(rho - 1);
            printf("%f  ", coupling[ind]);
         }
         printf("\n");
      }
   }

void tabulate_coupling_su3so3()
   // Tabulate SU(3)-SO(3) reduced coupling coefficients
   // / (6,1)    (2,1)   ||  (5,2)  \
   // \kappa1,2 kappa2,3 || kappa3,3/rho
   {
      su3irrep irrep1 = {6,1}, irrep2 = {2,1}, irrep3 = {5,2};
      int L1 = 2, L2 = 3, L3 = 3;
      int rhomax = outer_multiplicity(irrep1, irrep2, irrep3),
          kappa1max = inner_multiplicity(irrep1, L1),
          kappa2max = inner_multiplicity(irrep2, L2),
          kappa3max = inner_multiplicity(irrep3, L3);
      int dimen = kappa1max*kappa2max*kappa3max*rhomax;
      double coupling[dimen];
      calculate_coupling_su3so3(
         irrep1, L1, irrep2, L2, irrep3, L3, kappa1max, kappa2max, kappa3max,
         rhomax, coupling
         );
      printf("\n");
      printf("SU(3)-SO(3) reduced coupling coefficients\n");
      printf("/ (6,1)    (2,1)   ||  (5,2)  \\\n");
      printf("\\kappa1,2 kappa2,3 || kappa3,3/rho\n");
      printf("kappa1  kappa2  kappa3  coefficients for rho=1,...,rhomax=%d\n", rhomax);
      for (int kappa1 = 1; kappa1 <= kappa1max; kappa1++)
      {
         for (int kappa2 = 1; kappa2 <= kappa2max; kappa2++)
         {
            for (int kappa3 = 1; kappa3 <= kappa3max; kappa3++)
            {
               printf("   %d       %d       %d    ", kappa1, kappa2, kappa3);
               for (int rho = 1; rho <= rhomax; rho++)
               {
                  int ind = kappa1 + kappa1max*(kappa2 - 1)
                            + kappa1max*kappa2max*(kappa3 - 1)
                            + kappa1max*kappa2max*kappa3max*(rho - 1) - 1;
                  printf("%f  ", coupling[ind]);
               }
               printf("\n");
            }
         }
      }
   }

void tabulate_u_coef()
   // Tabulate SU(3) recoupling U coefficients
   // U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   {
      su3irrep irrep1 = {9,3}, irrep2 = {1,1}, irrep = {6,6},
               irrep3 = {2,2}, irrep12 = {9,3}, irrep23 = {3,3};
      int rhomaxa = outer_multiplicity(irrep1, irrep2, irrep12),
          rhomaxb = outer_multiplicity(irrep12, irrep3, irrep),
          rhomaxc = outer_multiplicity(irrep2, irrep3, irrep23),
          rhomaxd = outer_multiplicity(irrep1, irrep23, irrep);
      int dimen = rhomaxa*rhomaxb*rhomaxc*rhomaxd;
      double rac[dimen];
      int info = calculate_u_coef(
         irrep1, irrep2, irrep, irrep3, irrep12, irrep23,
         rhomaxa, rhomaxb, rhomaxc, rhomaxd,
         rac
         );
      printf("\n");
      if (info != 0)
      {
         printf("calculate_u_coef: Lapack subroutine dgesv ran with error: info=%d\n", info);
         return;
      }
      printf("SU(3) recoupling coefficients U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]\n");
      printf("rhoa  rhob  rhoc  rhod  coefficient\n");
      for (int rhoa = 1; rhoa <= rhomaxa; rhoa++)
      {
         for (int rhob = 1; rhob <= rhomaxb; rhob++)
         {
            for (int rhoc = 1; rhoc <= rhomaxc; rhoc++)
            {
               for (int rhod = 1; rhod <= rhomaxd; rhod++)
               {
                  int ind = rhod - 1 + rhomaxd*(rhoa - 1) + rhomaxd*rhomaxa*(rhob - 1) + rhomaxd*rhomaxa*rhomaxb*(rhoc - 1);
                  printf("  %d     %d     %d     %d   %f\n", rhoa, rhob, rhoc, rhod, rac[ind]);
               }
            }
         }
      }
   }

void tabulate_z_coef()
   // Tabulate SU(3) recoupling Z coefficients
   // Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   {
      su3irrep irrep2 = {9, 3}, irrep1 = {1, 1}, irrep = {6, 6},
               irrep3 = {2, 2}, irrep12 = {9, 3}, irrep13 = {3, 3};
      int rhomaxa = outer_multiplicity(irrep1, irrep2, irrep12),
          rhomaxb = outer_multiplicity(irrep12, irrep3, irrep),
          rhomaxc = outer_multiplicity(irrep1, irrep3, irrep13),
          rhomaxd = outer_multiplicity(irrep13, irrep2, irrep);
      int dimen = rhomaxa*rhomaxb*rhomaxc*rhomaxd;
      double Zcoeff[dimen];
      int info = calculate_z_coef(
         irrep2, irrep1, irrep, irrep3, irrep12, irrep13,
         rhomaxa, rhomaxb, rhomaxc, rhomaxd,
         Zcoeff
      );
      printf("\n");
      if (info != 0)
      {
         printf("calculate_z_coef: Lapack subroutine dgesv ran with error: info=%d\n", info);
         return;
      }
      printf("SU(3) recoupling coefficients Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]\n");
      printf("rhoa  rhob  rhoc  rhod  coefficient\n");
      for (int rhoa = 1; rhoa <= rhomaxa; rhoa++)
      {
         for (int rhob = 1; rhob <= rhomaxb; rhob++)
         {
            for (int rhoc = 1; rhoc <= rhomaxc; rhoc++)
            {
               for (int rhod = 1; rhod <= rhomaxd; rhod++)
               {
                  int ind = rhod - 1 + rhomaxd*(rhoa - 1) + rhomaxd*rhomaxa*(rhob - 1) + rhomaxd*rhomaxa*rhomaxb*(rhoc - 1);
                  printf("  %d     %d     %d     %d   %f\n", rhoa, rhob, rhoc, rhod, Zcoeff[ind]);
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
      su3irrep irrep1 = {1, 1}, irrep2 = {0, 0}, irrep12 = {1, 1},
               irrep3 = {0, 4}, irrep4 = {4, 2}, irrep34 = {2, 4},
               irrep13 = {2, 3}, irrep24 = {4, 2}, irrep = {2, 4};
      int rhomax12 = outer_multiplicity(irrep1, irrep2, irrep12),
          rhomax34 = outer_multiplicity(irrep3, irrep4, irrep34),
          rhomax1234 = outer_multiplicity(irrep12, irrep34, irrep),
          rhomax13 = outer_multiplicity(irrep1, irrep3, irrep13),
          rhomax24 = outer_multiplicity(irrep2, irrep4, irrep24),
          rhomax1324 = outer_multiplicity(irrep13, irrep24, irrep);
      int dimen = rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*rhomax1324;
      double ninelm[dimen];
      int info = calculate_9_lambda_mu(
         irrep1, irrep2, irrep12,
         irrep3, irrep4, irrep34,
         irrep13, irrep24, irrep,
         rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324,
         ninelm
      );
      printf("\n");
      if(info != 0)
      {
         printf("calculate_9_lambda_mu: Lapack subroutine dgesv ran with error: info=%d\n", info);
         return;
      }
      printf("                           |(1,1) (0,0)  (1,1)  rho12 |\n");
      printf("9-(lambda,mu) coefficients |(0,4) (4,2)  (2,4)  rho34 |\n");
      printf("                           |(2,3) (4,2)  (2,4) rho1324|\n");
      printf("                           |rho13 rho24 rho1234       |\n");
      printf("rho12  rho34  rho1234  rho13  rho24  rho1324  coefficient\n");
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
                        printf(
                           "  %d      %d       %d       %d      %d       %d     %f\n",
                           rho12, rho34, rho1234, rho13, rho24, rho1324, ninelm[ind]
                        );
                     }
                  }
               }
            }
         }
      }
   }

int main()
{
   initialize_ndsu3lib(true, false, 50);
   tabulate_coupling_canonical();
   tabulate_coupling_su3so3();
   tabulate_u_coef();
   tabulate_z_coef();
   tabulate_nine_lm();
   finalize_ndsu3lib(true);
   return 0;
}
