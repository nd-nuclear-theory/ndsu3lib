!***********************************************
!
! ndsu3lib_example.f90
!! Simple program demonstrating usage of ndsu3lib
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!***********************************************
SUBROUTINE tabulate_coupling_canonical
   !---------------------------------------------------------
   !! Tabulate SU(3)-SU(2)xU(1) reduced coupling coefficients
   ! /     (6,1)            (2,1)       ||  (5,2) \
   ! \epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho
   !---------------------------------------------------------
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   IMPLICIT NONE
   TYPE(su3irrep) :: irrep1, irrep2, irrep3
   INTEGER :: epsilon3, Lambda3_twice, dimpq, dimw, rhomax, numb, i, p1, q1, p2, q2, &
              epsilon1, Lambda1_twice, epsilon2, Lambda2_twice, rho
   INTEGER, ALLOCATABLE, DIMENSION(:) :: p1a, p2a, q2a
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: wigner_block

   irrep1%lambda = 6
   irrep1%mu = 1
   irrep2%lambda = 2
   irrep2%mu = 1
   irrep3%lambda = 5
   irrep3%mu = 2
   epsilon3 = -9
   Lambda3_twice = 5 ! Lambda3_twice is 2*Lambda3

   rhomax = outer_multiplicity(irrep1, irrep2, irrep3)

   dimpq = (MAX(irrep1%lambda, irrep1%mu) + 1)*(irrep2%lambda + 1)*(irrep2%mu + 1)
   dimw = rhomax*dimpq

   ALLOCATE (wigner_block(dimw), p1a(dimpq), p2a(dimpq), q2a(dimpq))

   CALL calculate_coupling_canonical(irrep1, irrep2, irrep3, epsilon3, Lambda3_twice, &
                                     dimpq, dimw, rhomax, numb, wigner_block, p1a, p2a, q2a)

   WRITE (*, "(A)") "SU(3)-SU(2)xU(1) reduced coupling coefficients"
   WRITE (*, "(A)") "/     (6,1)            (2,1)       ||  (5,2) \ "
   WRITE (*, "(A)") "\epsilon1,Lambda1 epsilon2,Lambda2 || -9,5/2 /rho"
   WRITE (*, "(A,I1)") "epsilon1  2*Lambda1  epsilon2  2*Lambda2  coefficients for rho=1,...,rhomax=", rhomax
   DO i = 1, numb
      p1 = p1a(i)
      p2 = p2a(i)
      q2 = q2a(i)
      q1 = (2*(irrep1%lambda + irrep2%lambda) + irrep1%mu + irrep2%mu - epsilon3)/3 - p1 - p2 - q2
      epsilon2 = 2*irrep2%lambda + irrep2%mu - 3*(p2 + q2)
      epsilon1 = epsilon3 - epsilon2
      Lambda1_twice = irrep1%mu + p1 - q1 ! Lambda1_twice is 2*Lambda1
      Lambda2_twice = irrep2%mu + p2 - q2 ! Lambda2_twice is 2*Lambda2
      WRITE (*, "(I5,I10,I11,I10,2D28.16)") epsilon1, Lambda1_twice, epsilon2, Lambda2_twice, &
         (wigner_block(i + numb*(rho - 1)), rho=1, rhomax)
   END DO

   DEALLOCATE (wigner_block, p1a, p2a, q2a)

END SUBROUTINE tabulate_coupling_canonical

SUBROUTINE tabulate_coupling_su3so3
   !----------------------------------------------------
   !! Tabulate SU(3)-SO(3) reduced coupling coefficients
   ! / (6,1)    (2,1)   ||  (5,2)  \
   ! \kappa1,2 kappa2,3 || kappa3,3/rho
   !----------------------------------------------------
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   USE ndsu3lib_coupling_su3so3
   IMPLICIT NONE
   TYPE(su3irrep) :: irrep1, irrep2, irrep3
   INTEGER :: rho, rhomax, L1, L2, L3, kappa1max, kappa2max, kappa3max, kappa1, kappa2, kappa3
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :, :) :: wigner_phys

   irrep1%lambda = 6
   irrep1%mu = 1
   irrep2%lambda = 2
   irrep2%mu = 1
   irrep3%lambda = 5
   irrep3%mu = 2
   L1 = 2
   L2 = 3
   L3 = 3
   rhomax = outer_multiplicity(irrep1, irrep2, irrep3)
   kappa1max = inner_multiplicity(irrep1, L1)
   kappa2max = inner_multiplicity(irrep2, L2)
   kappa3max = inner_multiplicity(irrep3, L3)

   ALLOCATE (wigner_phys(kappa1max, kappa2max, kappa3max, rhomax))

   CALL calculate_coupling_su3so3(irrep1, L1, kappa1max, irrep2, L2, kappa2max, irrep3, L3, kappa3max, rhomax, wigner_phys)

   WRITE (*, *)
   WRITE (*, "(A)") "SU(3)-SO(3) reduced coupling coefficients"
   WRITE (*, "(A)") "/ (6,1)    (2,1)   ||  (5,2)  \ "
   WRITE (*, "(A)") "\kappa1,2 kappa2,3 || kappa3,3/rho"
   WRITE (*, "(A,I1)") "kappa1  kappa2  kappa3  coefficients for rho=1,...,rhomax=", rhomax
   DO kappa1 = 1, kappa1max
      DO kappa2 = 1, kappa2max
         DO kappa3 = 1, kappa3max
            WRITE (*, "(I4,I8,I8,2D26.16)") kappa1, kappa2, kappa3, &
               (wigner_phys(kappa1, kappa2, kappa3, rho), rho=1, rhomax)
         END DO
      END DO
   END DO

   DEALLOCATE (wigner_phys)

END SUBROUTINE tabulate_coupling_su3so3

SUBROUTINE tabulate_u_coef
   !-----------------------------------------------------
   !! Tabulate SU(3) recoupling U coefficients
   ! U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   !-----------------------------------------------------
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   USE ndsu3lib_recoupling
   IMPLICIT NONE
   TYPE(su3irrep) :: irrep1, irrep2, irrep3, irrep, irrep12, irrep23
   INTEGER :: rhomaxa, rhomaxb, rhomaxc, rhomaxd, info, rhoa, rhob, rhoc, rhod, n
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: rac

   irrep1%lambda = 9
   irrep1%mu = 3
   irrep2%lambda = 1
   irrep2%mu = 1
   irrep%lambda = 6
   irrep%mu = 6
   irrep3%lambda = 2
   irrep3%mu = 2
   irrep12%lambda = 9
   irrep12%mu = 3
   irrep23%lambda = 3
   irrep23%mu = 3

   rhomaxa = outer_multiplicity(irrep1, irrep2, irrep12)
   rhomaxb = outer_multiplicity(irrep12, irrep3, irrep)
   rhomaxc = outer_multiplicity(irrep2, irrep3, irrep23)
   rhomaxd = outer_multiplicity(irrep1, irrep23, irrep)

   ALLOCATE (rac(rhomaxd, rhomaxa*rhomaxb*rhomaxc))

   CALL calculate_u_coef(irrep1, irrep2, irrep, irrep3, irrep12, irrep23, &
                         rhomaxa, rhomaxb, rhomaxc, rhomaxd, rac, rhomaxd, info)

   WRITE (*, *)
   IF (info /= 0) THEN
      WRITE (*, "(A,I2)") "calculate_u_coef: Lapack subroutine dgesv ran with error: info=", info
      DEALLOCATE (rac)
      RETURN
   END IF
   WRITE (*, "(A)") "SU(3) recoupling coefficients U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]"
   WRITE (*, "(A)") "rhoa  rhob  rhoc  rhod  coefficient"
   DO rhoa = 1, rhomaxa
      DO rhob = 1, rhomaxb
         DO rhoc = 1, rhomaxc
            DO rhod = 1, rhomaxd
               n = rhoa + rhomaxa*(rhob - 1) + rhomaxa*rhomaxb*(rhoc - 1)
               WRITE (*, "(I3,I6,I6,I6,D25.16)") rhoa, rhob, rhoc, rhod, rac(rhod, n)
            END DO
         END DO
      END DO
   END DO

   DEALLOCATE (rac)

END SUBROUTINE tabulate_u_coef

SUBROUTINE tabulate_z_coef
   !-----------------------------------------------------
   !! Tabulate SU(3) recoupling Z coefficients
   ! Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
   !-----------------------------------------------------
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   USE ndsu3lib_recoupling
   IMPLICIT NONE
   TYPE(su3irrep) :: irrep1, irrep2, irrep3, irrep, irrep12, irrep13
   INTEGER :: rhomaxa, rhomaxb, rhomaxc, rhomaxd, info, rhoa, rhob, rhoc, rhod, n
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: Zcoeff

   irrep2%lambda = 9
   irrep2%mu = 3
   irrep1%lambda = 1
   irrep1%mu = 1
   irrep%lambda = 6
   irrep%mu = 6
   irrep3%lambda = 2
   irrep3%mu = 2
   irrep12%lambda = 9
   irrep12%mu = 3
   irrep13%lambda = 3
   irrep13%mu = 3

   rhomaxa = outer_multiplicity(irrep1, irrep2, irrep12)
   rhomaxb = outer_multiplicity(irrep12, irrep3, irrep)
   rhomaxc = outer_multiplicity(irrep1, irrep3, irrep13)
   rhomaxd = outer_multiplicity(irrep13, irrep2, irrep)

   ALLOCATE (Zcoeff(rhomaxd, rhomaxa*rhomaxb*rhomaxc))

   CALL calculate_z_coef(irrep2, irrep1, irrep, irrep3, irrep12, irrep13, &
                         rhomaxa, rhomaxb, rhomaxc, rhomaxd, Zcoeff, rhomaxd, info)

   WRITE (*, *)
   IF (info /= 0) THEN
      WRITE (*, "(A,I2)") "calculate_z_coef: Lapack subroutine dgesv ran with error: info=", info
      DEALLOCATE (Zcoeff)
      RETURN
   END IF
   WRITE (*, "(A)") "SU(3) recoupling coefficients Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]"
   WRITE (*, "(A)") "rhoa  rhob  rhoc  rhod  coefficient"
   DO rhoa = 1, rhomaxa
      DO rhob = 1, rhomaxb
         DO rhoc = 1, rhomaxc
            DO rhod = 1, rhomaxd
               n = rhoa + rhomaxa*(rhob - 1) + rhomaxa*rhomaxb*(rhoc - 1)
               WRITE (*, "(I3,I6,I6,I6,D25.16)") rhoa, rhob, rhoc, rhod, Zcoeff(rhod, n)
            END DO
         END DO
      END DO
   END DO

   DEALLOCATE (Zcoeff)

END SUBROUTINE tabulate_z_coef

SUBROUTINE tabulate_nine_lm
   !-------------------------------------
   !! Tabulate 9-(lambda,mu) coefficients
   ! |(1,1) (0,0)  (1,1)  rho12 |
   ! |(0,4) (4,2)  (2,4)  rho34 |
   ! |(2,3) (4,2)  (2,4) rho1324|
   ! |rho13 rho24 rho1234       |
   !-------------------------------------
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   USE ndsu3lib_recoupling
   IMPLICIT NONE
   TYPE(su3irrep) :: irrep1, irrep2, irrep3, irrep, irrep12, irrep13, irrep4, irrep24, irrep34
   INTEGER :: info, rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, &
              rho12, rho34, rho1234, rho13, rho24, rho1324
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :, :, :, :) :: ninelm

   irrep1%lambda = 1
   irrep1%mu = 1
   irrep2%lambda = 0
   irrep2%mu = 0
   irrep12%lambda = 1
   irrep12%mu = 1
   irrep3%lambda = 0
   irrep3%mu = 4
   irrep4%lambda = 4
   irrep4%mu = 2
   irrep34%lambda = 2
   irrep34%mu = 4
   irrep13%lambda = 2
   irrep13%mu = 3
   irrep24%lambda = 4
   irrep24%mu = 2
   irrep%lambda = 2
   irrep%mu = 4

   rhomax12 = outer_multiplicity(irrep1, irrep2, irrep12)
   rhomax34 = outer_multiplicity(irrep3, irrep4, irrep34)
   rhomax1234 = outer_multiplicity(irrep12, irrep34, irrep)
   rhomax13 = outer_multiplicity(irrep1, irrep3, irrep13)
   rhomax24 = outer_multiplicity(irrep2, irrep4, irrep24)
   rhomax1324 = outer_multiplicity(irrep13, irrep24, irrep)

   ALLOCATE (ninelm(rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324))

   CALL calculate_9_lambda_mu(irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep, &
                              rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, ninelm, info)

   WRITE (*, *)
   IF (info /= 0) THEN
      WRITE (*, "(A,I2)") "calculate_9_lambda_mu: Lapack subroutine dgesv ran with error: info=", info
      DEALLOCATE (ninelm)
      RETURN
   END IF
   WRITE (*, "(27X,A)") "|(1,1) (0,0)  (1,1)  rho12 |"
   WRITE (*, "(A)") "9-(lambda,mu) coefficients |(0,4) (4,2)  (2,4)  rho34 |"
   WRITE (*, "(27X,A)") "|(2,3) (4,2)  (2,4) rho1324|"
   WRITE (*, "(27X,A)") "|rho13 rho24 rho1234       |"
   WRITE (*, "(A)") "rho12  rho34  rho1234  rho13  rho24  rho1324  coefficient"
   DO rho12 = 1, rhomax12
      DO rho34 = 1, rhomax34
         DO rho1234 = 1, rhomax1234
            DO rho13 = 1, rhomax13
               DO rho24 = 1, rhomax24
                  DO rho1324 = 1, rhomax1324
                     WRITE (*, "(I3,I7,I8,I8,I7,I8,D27.16)") rho12, rho34, rho1234, rho13, rho24, rho1324, &
                        ninelm(rho12, rho34, rho1234, rho13, rho24, rho1324)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   DEALLOCATE (ninelm)

END SUBROUTINE tabulate_nine_lm

PROGRAM ndsu3lib_example
   !------------------------------------------------------------------------------
   !! Simple program tabulating SU(3) reduced coupling and recoupling coefficients
   !------------------------------------------------------------------------------
   USE ndsu3lib_tools
   IMPLICIT NONE
   LOGICAL(KIND=1) :: wso3 = .TRUE., openmp = .FALSE.
      !! First input argument of subroutines initialize_ndsu3lib and finalize_ndsu3lib.
      !! Must be .TRUE. if SU(3)-SO(3) reduced coupling coefficients are going to be calculated

   CALL initialize_ndsu3lib(wso3, openmp, 50)
   CALL tabulate_coupling_canonical
   CALL tabulate_coupling_su3so3
   CALL tabulate_u_coef
   CALL tabulate_z_coef
   CALL tabulate_nine_lm
   CALL finalize_ndsu3lib(wso3)

END PROGRAM ndsu3lib_example
