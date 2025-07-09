!******************************************
!
! ndsu3lib_recoupling.F90
!! Module for SU(3) recoupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!******************************************
MODULE ndsu3lib_recoupling
   !! Module for SU(3) recoupling coefficients
   USE ndsu3lib_tools
   USE ndsu3lib_coupling_canonical
   IMPLICIT NONE
CONTAINS

   FUNCTION dimen(irrep) RESULT(dm)
      !----------------------------------------------------
      !! Dimension of SU(3) irrep
      ! (lambda,mu), where lambda=irrep%lambda, mu=irrep%mu
      !----------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN) :: irrep
         !! SU(3) irrep
      INTEGER :: dm
         !! Dimension of SU(3) irrep
      dm = (irrep%lambda + 1)*(irrep%mu + 1)*(irrep%lambda + irrep%mu + 2)/2
   END FUNCTION dimen

#if defined(NDSU3LIB_RACAH_GSL)
   FUNCTION su2racah(j1, j2, j, j3, j12, j23) RESULT(w)
      !------------------------------------------------------------
      !! Calculate SU(2) Racah coefficient W(J1,J2,J,J3;J12,J23)
      ! using GSL function gsl_sf_coupling_6j calculating 6j symbol
      !------------------------------------------------------------
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTERFACE
         REAL(C_DOUBLE) FUNCTION gsl_sf_coupling_6j(l1, l2, l12, l3, l, l23) BIND(C)
            USE ISO_C_BINDING
            INTEGER(C_INT), VALUE :: l1, l2, l12, l3, l, l23
         END FUNCTION gsl_sf_coupling_6j
      END INTERFACE
      INTEGER(C_INT), INTENT(IN) :: j1
         !! Twice the value of J1
      INTEGER(C_INT), INTENT(IN) :: j2
         !! Twice the value of J2
      INTEGER(C_INT), INTENT(IN) :: j
         !! Twice the value of J
      INTEGER(C_INT), INTENT(IN) :: j3
         !! Twice the value of J3
      INTEGER(C_INT), INTENT(IN) :: j12
         !! Twice the value of J12
      INTEGER(C_INT), INTENT(IN) :: j23
         !! Twice the value of J23
      REAL(C_DOUBLE) :: w
         !! Racah coefficient
      INTEGER :: a
      w = gsl_sf_coupling_6j(j1, j2, j12, j3, j, j23)
      a = j1 + j2 + j3 + j
      IF ((a/4)*4 /= a) w = -w
      RETURN
   END FUNCTION su2racah
#endif

   SUBROUTINE calculate_u_coef(irrep1, irrep2, irrep, irrep3, irrep12, irrep23, rhomaxa, rhomaxb, rhomaxc, rhomaxd, rac, ldb, info)
      !----------------------------------------------------------------------------------------------------------------------------
      !! Calculate SU(3) recoupling U coefficients
      ! U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda23,mu23)rhoc,rhod]
      ! for given lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
      ! using Eq.(A.24) in [1] or equivalently Eq.(22),(35,1B) in [2] and
      !! calling Lapack subroutine dgesv to solve system of linear equations
      !
      ! References: [1] J.Herko et al. arXiv:2505.08993
      !             [2] J.P.Draayer, Y.Akiyama, J.Math.Phys. 14, 1904 (1973)
      !
      ! Input arguments: irrep1,irrep2,irrep,irrep3,irrep12,irrep23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
      ! Output arguments: rac,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      ! rhomaxc is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23),
      ! rhomaxd is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu),
      ! ldb is leading dimension of array rac,
      ! rac(rhod,n), where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1), is U coefficient for given rhoa,rhob,rhoc,rhod,
      ! info is 0 iff dgesv ran without errors.
      !----------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
#if defined(NDSU3LIB_RACAH_WIGXJPF)
      REAL(KIND=8), EXTERNAL :: fwig6jj
      REAL(KIND=8) :: su2racah
#endif
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! SU(3) irrep 1
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! SU(3) irrep 2
      TYPE(su3irrep), INTENT(IN) :: irrep
         !! SU(3) irrep resulting from SU(3) couplings 1x23 and 12x3
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! SU(3) irrep 3
      TYPE(su3irrep), INTENT(IN) :: irrep12
         !! SU(3) irrep 12
      TYPE(su3irrep), INTENT(IN) :: irrep23
         !! SU(3) irrep 23
      INTEGER, INTENT(IN) :: rhomaxa
         !! Multiplicity of SU(3) coupling 1x2->12
      INTEGER, INTENT(IN) :: rhomaxb
         !! Multiplicity of SU(3) coupling 12x3->irrep
      INTEGER, INTENT(IN) :: rhomaxc
         !! Multiplicity of SU(3) coupling 2x3->23
      INTEGER, INTENT(IN) :: rhomaxd
         !! Multiplicity of SU(3) coupling 1x23->irrep
      INTEGER, INTENT(IN) :: ldb
         !! Leading dimension of array rac
      INTEGER, INTENT(OUT) :: info
         !! Equals 0 iff Lapack subroutine dgesv ran without errors.
      TYPE(su3irrep) :: irrep2c
      INTEGER :: epsilon23, rhomaxabc, numba, numbb, numbc, numbd, i1, i2, inda, indd, i, &
                 Lambda122, epsilon2, Lambda22, p3, q3, n, rhoa, rhob, rhoc, I23, phase, &
                 Lambda232, Lambda32, p23, q23, p12, q12, p2, q2, m, noname3, epsilon2max, aux, p3min, pqdima, pqdimc, pqdimd
      REAL(KIND=8) :: factor1, factor2, factor3
      REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: rac
         !! Array of U coefficients.
         !! rac(rhod,n), where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1), is U coefficient for given rhoa,rhob,rhoc,rhod.
         !! Sizes are at least rhomaxd and rhomaxa*rhomaxb*rhomaxc
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: matrix ! Sizes are at least rhomaxd and rhomaxd.
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :, :) :: wignera, wignerb, wignerc, wignerd, wigner
      INTEGER, ALLOCATABLE, DIMENSION(:) :: p1aa, p2aa, q2aa, p1ac, p2ac, q2ac, p2ad, q2ad

      pqdima = (irrep12%lambda + 1)*(MAX(irrep2%mu, irrep3%lambda) + 1)*(MAX(irrep2%lambda, irrep3%mu) + 1)
      pqdimc = (MAX(irrep2%lambda, irrep2%mu) + 1)*(irrep3%lambda + 1)*(irrep3%mu + 1)
      pqdimd = (irrep1%lambda + 1)*(irrep23%lambda + 1)*(irrep23%mu + 1)
      ALLOCATE (matrix(rhomaxd, rhomaxd), wignera(0:irrep12%lambda, 0:irrep2%mu, 0:irrep2%lambda, 1:rhomaxa), &
                wignerb(0:irrep12%lambda, 0:irrep3%lambda, 0:irrep3%mu, 1:rhomaxb), &
                wignerc(0:MAX(irrep2%lambda, irrep2%mu), 0:MAX(irrep3%lambda, irrep3%mu), &
                        0:MAX(irrep3%lambda, irrep3%mu), 1:rhomaxc), &
                wignerd(0:irrep1%lambda, 0:irrep23%lambda, 0:irrep23%mu, 1:rhomaxd), &
                wigner(0:irrep2%lambda, 0:irrep3%lambda, 0:irrep3%mu, 1:rhomaxc), &
                p1aa((MAX(irrep12%lambda, irrep1%lambda) + 1)*(MAX(irrep2%mu, irrep3%lambda, irrep23%lambda) + 1) &
                     *(MAX(irrep2%lambda, irrep3%mu, irrep23%mu) + 1)), &
                p2aa(MAX(pqdima, rhomaxd)), &
                q2aa(pqdima), &
                p1ac(pqdimc), &
                p2ac(pqdimc), &
                q2ac(pqdimc), &
                p2ad(pqdimd), &
                q2ad(pqdimd))

      rhomaxabc = rhomaxa*rhomaxb*rhomaxc
      epsilon23 = -irrep%lambda - 2*irrep%mu + irrep1%lambda + 2*irrep1%mu

      IF (2*epsilon23 <= irrep23%lambda - irrep23%mu) THEN
         I23 = 1
      ELSE
         I23 = 0
      END IF

      epsilon2max = 2*irrep2%mu + irrep2%lambda
      i1 = 3*irrep1%lambda - 6*irrep12%lambda + 6*irrep1%mu - 6*irrep12%mu + 8*irrep2%lambda + 4*irrep2%mu
      factor1 = DBLE(INT(irrep1%lambda + 1, 8)*dimen(irrep12))/DBLE(dimen(irrep1))
      m = (2*(irrep12%lambda + irrep2%mu + irrep1%mu - irrep12%mu) + irrep2%lambda + irrep1%lambda)/3
      aux = 2*irrep3%lambda + irrep3%mu - epsilon23

      rac(1:rhomaxd, 1:rhomaxabc) = 0.D0

      CALL calculate_coupling_canonical_extremal(irrep1, irrep23, irrep, 1, rhomaxd, numbd, wignerd, p1aa, p2ad, q2ad)
      CALL calculate_coupling_canonical_extremal(irrep2, irrep3, irrep23, I23, rhomaxc, numbc, wignerc, p1ac, p2ac, q2ac)
      CALL calculate_coupling_canonical_extremal(irrep12, irrep3, irrep, 1, rhomaxb, numbb, wignerb, p1aa, p2aa, q2aa)
      irrep2c%lambda = irrep2%mu
      irrep2c%mu = irrep2%lambda
      CALL calculate_coupling_canonical_extremal(irrep12, irrep2c, irrep1, 1, rhomaxa, numba, wignera, p1aa, p2aa, q2aa)

      i = 0
      Lambda232 = irrep23%mu + p2ad(numbd) - q2ad(numbd)
      DO indd = numbd, numbd - rhomaxd + 1, -1 ! This is loop over Lambda23
         p23 = p2ad(indd)
         q23 = q2ad(indd)
         i = i + 1
         matrix(i, 1:rhomaxd) = wignerd(irrep1%lambda, p23, q23, 1:rhomaxd)
         factor2 = DSQRT(factor1*DBLE(Lambda232 + 1))

         CALL calculate_coupling_canonical_nonextremal(irrep2, irrep3, irrep23, epsilon23, Lambda232, I23, &
                                                       rhomaxc, numbc, wignerc, wigner, p1ac, p2ac, q2ac)

         DO inda = 1, numba ! sum over epsilon2,Lambda2,Lambda12
            p12 = p1aa(inda)
            p2 = p2aa(inda)
            q2 = q2aa(inda)
            Lambda122 = 2*p12 - m + p2 + q2
            epsilon2 = epsilon2max - 3*(p2 + q2) ! epsilon2 is -epsilon2 in the formula
            Lambda22 = irrep2%lambda + p2 - q2
            noname3 = (aux - epsilon2)/3
            i2 = i1 + epsilon2 + 3*Lambda122
            p3min = MAX(0, noname3 - irrep3%mu)
            q3 = noname3 - p3min
            Lambda32 = irrep3%mu + p3min - q3
            DO p3 = p3min, MIN(irrep3%lambda, noname3) ! sum over Lambda3
#if defined(NDSU3LIB_RACAH_GSL)
               factor3 = factor2*su2racah(irrep1%lambda, Lambda22, irrep%lambda, Lambda32, Lambda122, Lambda232)
#elif defined(NDSU3LIB_RACAH_WIGXJPF)
               su2racah = fwig6jj(irrep1%lambda, Lambda22, Lambda122, Lambda32, irrep%lambda, Lambda232)
               phase = irrep1%lambda + Lambda22 + Lambda32 + irrep%lambda
               IF ((phase/4)*4 /= phase) su2racah = -su2racah
               factor3 = factor2*su2racah
#endif
               IF (12*(i2/12) /= i2) factor3 = -factor3
               n = 0
               DO rhoc = 1, rhomaxc
                  DO rhob = 1, rhomaxb
                     DO rhoa = 1, rhomaxa
                        ! n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
                        n = n + 1
                        rac(i, n) = rac(i, n) + factor3*wignera(p12, p2, q2, rhoa) &
                                    *wignerb(p12, p3, q3, rhob)*wigner(irrep2%lambda - q2, p3, q3, rhoc)
                        ! See relation between p and \tilde{p} in Eq.(32) in [2].
                     END DO
                  END DO
               END DO
               q3 = q3 - 1
               Lambda32 = Lambda32 + 2
            END DO
         END DO
         Lambda232 = Lambda232 - 2
      END DO

      IF (rhomaxd > 1) THEN
         CALL dgesv(rhomaxd, rhomaxabc, matrix, rhomaxd, p2aa, rac, ldb, info)
      ELSE
         rac(1, 1:rhomaxabc) = rac(1, 1:rhomaxabc)/matrix(1, 1)
         info = 0
      END IF

      DEALLOCATE (matrix, wignera, wignerb, wignerc, wignerd, wigner, p1aa, p2aa, q2aa, p1ac, p2ac, q2ac, p2ad, q2ad)

   END SUBROUTINE calculate_u_coef

   SUBROUTINE calculate_z_coef(irrep2, irrep1, irrep, irrep3, irrep12, irrep13, &
                               rhomaxa, rhomaxb, rhomaxc, rhomaxd, Zcoeff, ldb, info)
      !------------------------------------------------------------------------------------------------------------------------
      !! Calculate SU(3) recoupling Z coefficients
      ! Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda13,mu13)rhoc,rhod]
      ! for given lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
      ! using Eq.(A.25) in [1] or equivalently Eq.(2) in [2] and
      !! calling Lapack subroutine dgesv to solve system of linear equations
      !
      ! References: [1] J.Herko et al. arXiv:2505.08993
      !             [2] D.J.Millener, J.Math.Phys. 19, 1513 (1978)
      !
      ! Input arguments: irrep2,irrep1,irrep,irrep3,irrep12,irrep13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
      ! Output arguments: Zcoeff,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      ! rhomaxc is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      ! rhomaxd is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu),
      ! ldb is leading dimension of array Zcoeff,
      ! Zcoeff(rhod,n), where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1), is Z coefficient for given rhoa,rhob,rhoc,rhod,
      ! info is 0 iff dgesv ran without errors.
      !------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
#if defined(NDSU3LIB_RACAH_WIGXJPF)
      REAL(KIND=8), EXTERNAL :: fwig6jj
      REAL(KIND=8) :: su2racah
#endif
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! SU(3) irrep 1
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! SU(3) irrep 2
      TYPE(su3irrep), INTENT(IN) :: irrep
         !! SU(3) irrep resulting from SU(3) couplings 12x3 and 13x2
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! SU(3) irrep 3
      TYPE(su3irrep), INTENT(IN) :: irrep12
         !! SU(3) irrep 12
      TYPE(su3irrep), INTENT(IN) :: irrep13
         !! SU(3) irrep 13
      INTEGER, INTENT(IN) :: rhomaxa
         !! Multiplicity of SU(3) coupling 1x2->12
      INTEGER, INTENT(IN) :: rhomaxb
         !! Multiplicity of SU(3) coupling 12x3->irrep
      INTEGER, INTENT(IN) :: rhomaxc
         !! Multiplicity of SU(3) coupling 1x3->13
      INTEGER, INTENT(IN) :: rhomaxd
         !! Multiplicity of SU(3) coupling 13x2->irrep
      INTEGER, INTENT(IN) :: ldb
         !! Leading dimension of array Zcoeff
      INTEGER, INTENT(OUT) :: info
         !! Equals 0 iff Lapack subroutine dgesv ran without errors.
      INTEGER :: rhomaxabc, numbahw, numbalw, numbb, numbd, i, indd, p2, q2, indb, &
                 p12, p3, q3, epsilon, epsilon12, Lambda122, Lambda32, &
                 p1, pq1, n, rhoa, rhob, rhoc, expp, Lambda22, LLambda12, p1min, &
                 aux1, mu1mpq1, aux3, aux4, lambdammu12, phase, &
                 aux5, aux6, inddmin, epsilon1lwpepsilon2, epsilonmepsilon3lw, &
                 epsilon12lw, lambdamlambda13, pqdima, pqdimb, pqdimd
      REAL(KIND=8) :: factor1, factor2, factor3
      REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: Zcoeff
         !! Array of Z coefficients.
         !! Zcoeff(rhod,n), where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1), is Z coefficient for given rhoa,rhob,rhoc,rhod.
         !! Sizes are at least rhomaxd and rhomaxa*rhomaxb*rhomaxc.
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: matrix
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :, :) :: wignerahw, wigneralw, wignerb, wignerc, wignerd, wigner
      INTEGER, ALLOCATABLE, DIMENSION(:) :: p1aa, p2aa, q2aa, p1ab, p2ab, q2ab, p2ad, q2ad

      pqdima = (irrep1%lambda + 1)*(irrep2%lambda + 1)*(irrep2%mu + 1)
      pqdimb = MAX((irrep13%lambda + 1)*(irrep2%lambda + 1)*(irrep2%mu + 1), &
                   (irrep1%lambda + 1)*(irrep3%lambda + 1)*(irrep3%mu + 1), &
                   (irrep1%lambda + 1)*(irrep2%lambda + 1)*(irrep2%mu + 1), &
                   (irrep12%lambda + 1)*(irrep3%lambda + 1)*(irrep3%mu + 1), &
                   (irrep1%mu + 1)*(irrep2%mu + 1)*(irrep2%lambda + 1))
      pqdimd = (irrep13%lambda + 1)*(irrep2%lambda + 1)*(irrep2%mu + 1)
      ALLOCATE (matrix(rhomaxd, rhomaxd), &
                wignerahw(0:irrep1%lambda, 0:irrep2%lambda, 0:irrep2%mu, rhomaxa), &
                wigneralw(0:irrep1%mu, 0:irrep2%mu, 0:irrep2%lambda, rhomaxa), &
                wignerb(0:irrep12%lambda, 0:irrep3%lambda, 0:irrep3%mu, rhomaxb), &
                wignerc(0:irrep1%lambda, 0:irrep3%lambda, 0:irrep3%mu, rhomaxc), &
                wignerd(0:irrep13%lambda, 0:irrep2%lambda, 0:irrep2%mu, rhomaxd), &
                wigner(0:irrep1%lambda, 0:irrep2%lambda, 0:irrep2%mu, rhomaxa), &
                p1aa(pqdima), &
                p2aa(MAX(pqdima, rhomaxd)), &
                q2aa(pqdima), &
                p1ab(pqdimb), &
                p2ab(pqdimb), &
                q2ab(pqdimb), &
                p2ad(pqdimd), &
                q2ad(pqdimd))

      rhomaxabc = rhomaxa*rhomaxb*rhomaxc
      epsilon = -irrep%lambda - 2*irrep%mu
      epsilon1lwpepsilon2 = 2*irrep1%lambda + irrep1%mu + epsilon + irrep13%lambda + 2*irrep13%mu
      epsilonmepsilon3lw = epsilon - 2*irrep3%lambda - irrep3%mu
      epsilon12lw = 2*irrep12%lambda + irrep12%mu
      lambdamlambda13 = irrep%lambda - irrep13%lambda
      lambdammu12 = irrep12%lambda - irrep12%mu

      CALL calculate_coupling_canonical_extremal(irrep13, irrep2, irrep, 1, rhomaxd, numbd, wignerd, p1ab, p2ad, q2ad)
      CALL calculate_coupling_canonical_extremal(irrep1, irrep3, irrep13, 1, rhomaxc, numbahw, wignerc, p1ab, p2ab, q2ab)
      CALL calculate_coupling_canonical_extremal(irrep1, irrep2, irrep12, 1, rhomaxa, numbahw, wignerahw, p1ab, p2ab, q2ab)
      CALL calculate_coupling_canonical_extremal(irrep1, irrep2, irrep12, 0, rhomaxa, numbalw, wigneralw, p1ab, p2ab, q2ab)
      CALL calculate_coupling_canonical_extremal(irrep12, irrep3, irrep, 1, rhomaxb, numbb, wignerb, p1ab, p2ab, q2ab)

      inddmin = numbd - rhomaxd + 1

      ! Construction of matrix
      i = 0
      DO indd = inddmin, numbd ! loop over Lambda2
         p2 = p2ad(indd)
         q2 = q2ad(indd)
         i = i + 1
         matrix(i, 1:rhomaxd) = wignerd(irrep13%lambda, p2, q2, 1:rhomaxd)
      END DO

      ! Construction of RHS
      Zcoeff(1:rhomaxd, 1:rhomaxabc) = 0.D0
      DO indb = 1, numbb ! loop over epsilon1,epsilon3,epsilon12,Lambda3,Lambda12
         p12 = p1ab(indb)
         p3 = p2ab(indb)
         q3 = q2ab(indb)
         epsilon12 = epsilonmepsilon3lw + 3*(p3 + q3) ! =epsilon-epsilon3
         Lambda122 = irrep12%mu - (epsilon12lw - epsilon12)/3 + 2*p12
         Lambda32 = irrep3%mu + p3 - q3
         aux1 = lambdamlambda13 - Lambda122

         IF (2*epsilon12 <= lambdammu12) THEN
            CALL calculate_coupling_canonical_nonextremal(irrep1, irrep2, irrep12, epsilon12, Lambda122, 1, &
                                                          rhomaxa, numbahw, wignerahw, wigner, p1aa, p2aa, q2aa)
         ELSE
            CALL calculate_coupling_canonical_nonextremal(irrep1, irrep2, irrep12, epsilon12, Lambda122, 0, &
                                                          rhomaxa, numbalw, wigneralw, wigner, p1aa, p2aa, q2aa)
         END IF

         factor1 = DSQRT(DBLE((Lambda122 + 1)*(irrep13%lambda + 1)))

         pq1 = (epsilon1lwpepsilon2 - epsilon12)/3 ! pq1 is p1+q1
         mu1mpq1 = irrep1%mu - pq1
         aux3 = MIN(irrep1%lambda, pq1, (Lambda32 + irrep13%lambda - mu1mpq1)/2)
         aux4 = MAX(0, -mu1mpq1, (Lambda32 - irrep13%lambda - mu1mpq1)/2, -(Lambda32 - irrep13%lambda + mu1mpq1)/2)
         aux5 = Lambda122 - mu1mpq1
         aux6 = Lambda122 + mu1mpq1

         i = 0
         DO indd = inddmin, numbd ! loop over Lambda2
            i = i + 1
            p2 = p2ad(indd)
            q2 = q2ad(indd)
            Lambda22 = irrep2%mu + p2 - q2

            p1min = MAX(aux4, (Lambda22 - aux6)/2, (aux5 - Lambda22)/2)
            LLambda12 = mu1mpq1 + 2*p1min ! LLambda12 is 2*Lambda1
            expp = LLambda12 + aux1
            IF (4*(expp/4) == expp) THEN
               factor2 = -factor1
            ELSE
               factor2 = factor1
            END IF
            DO p1 = p1min, MIN((Lambda22 + aux5)/2, aux3)
               ! Lower and upper bounds on p1 are such that:
               ! 1) 0<=q1<=mu1, where q1=pq1-p1
               ! 2) ABS(LLambda12-Lambda22)<=Lambda122<=LLambda12+Lambda22, where LLambda12=mu1+p1-q1=mu1-pq1+2*p1
               ! 3) ABS(LLambda12-Lambda32)<=lambda13<=LLambda12+Lambda32

               factor2 = -factor2
#if defined(NDSU3LIB_RACAH_GSL)
               factor3 = factor2*su2racah(Lambda22, LLambda12, irrep%lambda, Lambda32, Lambda122, irrep13%lambda)
#elif defined(NDSU3LIB_RACAH_WIGXJPF)
               su2racah = fwig6jj(Lambda22, LLambda12, Lambda122, Lambda32, irrep%lambda, irrep13%lambda)
               phase = Lambda22 + LLambda12 + Lambda32 + irrep%lambda
               IF ((phase/4)*4 /= phase) su2racah = -su2racah
               factor3 = factor2*su2racah
#endif

               n = 0
               DO rhoc = 1, rhomaxc
                  DO rhob = 1, rhomaxb
                     DO rhoa = 1, rhomaxa
                        ! n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
                        n = n + 1
                        Zcoeff(i, n) = Zcoeff(i, n) + factor3*wignerc(p1, p3, q3, rhoc) &
                                       *wigner(p1, p2, q2, rhoa)*wignerb(p12, p3, q3, rhob)
                     END DO
                  END DO
               END DO

               LLambda12 = LLambda12 + 2
            END DO
         END DO
      END DO

      ! Solution of system of linear equations
      IF (rhomaxd > 1) THEN
         CALL dgesv(rhomaxd, rhomaxabc, matrix, rhomaxd, p2aa, Zcoeff, ldb, info)
      ELSE
         Zcoeff(1, 1:rhomaxabc) = Zcoeff(1, 1:rhomaxabc)/matrix(1, 1)
         info = 0
      END IF

      DEALLOCATE (matrix, wignerahw, wigneralw, wignerb, wignerc, wignerd, wigner, p1aa, p2aa, q2aa, p1ab, p2ab, q2ab, p2ad, q2ad)

   END SUBROUTINE calculate_z_coef

   SUBROUTINE calculate_9_lambda_mu(irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep, &
                                    rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, ninelm, info)
      !----------------------------------------------------------------------------------------------------------------------------
      !! Calculate 9-(lambda,mu) coefficients
      !
      ! | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
      ! | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
      ! |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
      ! |     rho13           rho24          rho1234            |
      !
      ! for given lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,lambda24,mu24,lambda,mu
      !! from U and Z coefficients obtained by calling calculate_u_coef and calculate_z_coef
      ! using Eq.(56) in J.Herko et al. arXiv:2505.08993
      !
      ! Input arguments: irrep1,irrep2,irrep12,irrep3,irrep4,irrep34,irrep13,irrep24,irrep,
      !                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
      ! Output arguments: ninelm,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomax12 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomax34 is multiplicity of SU(3) coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34),
      ! rhomax1234 is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu),
      ! rhomax13 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      ! rhomax24 is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24),
      ! rhomax1324 is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu),
      ! ninelm(rho12,rho34,rho1234,rho13,rho24,rho1324) is 9-(lambda,mu) coefficient for given
      !   rho12,rho34,rho1234,rho13,rho24,rho1324,
      ! info is 0 iff Lapack subroutine dgesv called by calculate_u_coef and calculate_z_coef ran without errors.
      !----------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! SU(3) irrep 1
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! SU(3) irrep 2
      TYPE(su3irrep), INTENT(IN) :: irrep12
         !! SU(3) irrep 12
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! SU(3) irrep 3
      TYPE(su3irrep), INTENT(IN) :: irrep4
         !! SU(3) irrep 4
      TYPE(su3irrep), INTENT(IN) :: irrep34
         !! SU(3) irrep 34
      TYPE(su3irrep), INTENT(IN) :: irrep13
         !! SU(3) irrep 13
      TYPE(su3irrep), INTENT(IN) :: irrep24
         !! SU(3) irrep 24
      TYPE(su3irrep), INTENT(IN) :: irrep
         !! SU(3) irrep resulting from SU(3) couplings 12x34 and 13x24
      INTEGER, INTENT(IN) :: rhomax12
         !! Multiplicity of SU(3) coupling 1x2->12
      INTEGER, INTENT(IN) :: rhomax34
         !! Multiplicity of SU(3) coupling 3x4->34
      INTEGER, INTENT(IN) :: rhomax1234
         !! Multiplicity of SU(3) coupling 12x34->irrep
      INTEGER, INTENT(IN) :: rhomax13
         !! Multiplicity of SU(3) coupling 1x3->13
      INTEGER, INTENT(IN) :: rhomax24
         !! Multiplicity of SU(3) coupling 2x4->24
      INTEGER, INTENT(IN) :: rhomax1324
         !! Multiplicity of SU(3) coupling 13x24->irrep
      INTEGER, INTENT(OUT) :: info
         !! Equals 0 iff Lapack subroutine dgesv called by calculate_u_coeff and calculate_z_coeff ran without errors.
      REAL(KIND=8), DIMENSION(:, :, :, :, :, :), INTENT(OUT) :: ninelm
         !! Array of 9-(lambda,mu) coefficients.
         !! ninelm(rho12,rho34,rho1234,rho13,rho24,rho1324) is 9-(lambda,mu) coefficient for given
         !! rho12,rho34,rho1234,rho13,rho24,rho1324.
         !! Sizes are rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324.
      TYPE(su3irrep) :: irrep0
      INTEGER :: lambda0, mu0, rho132, rho04, rho123, rhomax132, rhomax04, rhomax123, &
                 rho12, rho34, rho13, rho24, rho1324, nU1, nZ, nU2, rhomax12304
      INTEGER :: i1, i2
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: U1, Z, U2

      ninelm = 0.D0

      DO lambda0 = 0, MIN(irrep13%lambda + irrep2%lambda + MIN(irrep2%mu, irrep13%lambda + irrep13%mu), &
                          irrep12%lambda + irrep3%lambda + MIN(irrep3%mu, irrep12%lambda + irrep12%mu))
         irrep0%lambda = lambda0
         DO mu0 = 0, MIN(irrep13%mu + irrep2%mu + MIN(irrep13%lambda, irrep2%lambda), &
                         irrep12%mu + irrep3%mu + MIN(irrep12%lambda, irrep3%lambda))
            irrep0%mu = mu0
            rhomax132 = outer_multiplicity(irrep13, irrep2, irrep0)
            IF (rhomax132 == 0) CYCLE
            rhomax04 = outer_multiplicity(irrep0, irrep4, irrep)
            IF (rhomax04 == 0) CYCLE
            rhomax123 = outer_multiplicity(irrep12, irrep3, irrep0)
            IF (rhomax123 == 0) CYCLE
            rhomax12304 = rhomax123*rhomax04

            ALLOCATE (U1(rhomax1324, rhomax132*rhomax04*rhomax24), Z(rhomax132, rhomax12*rhomax123*rhomax13), &
                      U2(rhomax1234, rhomax123*rhomax04*rhomax34))

            CALL calculate_u_coef(irrep13, irrep2, irrep, irrep4, irrep0, irrep24, &
                                  rhomax132, rhomax04, rhomax24, rhomax1324, U1, rhomax1324, info)
            IF (info /= 0) THEN
               DEALLOCATE (U1, Z, U2)
               RETURN
            END IF
            CALL calculate_z_coef(irrep2, irrep1, irrep0, irrep3, irrep12, irrep13, &
                                  rhomax12, rhomax123, rhomax13, rhomax132, Z, rhomax132, info)
            IF (info /= 0) THEN
               DEALLOCATE (U1, Z, U2)
               RETURN
            END IF
            CALL calculate_u_coef(irrep12, irrep3, irrep, irrep4, irrep0, irrep34, &
                                  rhomax123, rhomax04, rhomax34, rhomax1234, U2, rhomax1234, info)
            IF (info /= 0) THEN
               DEALLOCATE (U1, Z, U2)
               RETURN
            END IF

            nU1 = 0
            DO rho24 = 1, rhomax24
               i1 = -rhomax123
               DO rho04 = 1, rhomax04
                  i1 = i1 + rhomax123 ! i1=rhomax123*(rho04-1)
                  DO rho132 = 1, rhomax132
                     nU1 = nU1 + 1 ! nU1=rho132+rhomax132*(rho04-1)+rhomax132*rhomax04*(rho24-1)
                     nZ = 0
                     DO rho13 = 1, rhomax13
                        i2 = i1
                        DO rho123 = 1, rhomax123
                           i2 = i2 + 1 ! i2=rho123+i1
                           DO rho12 = 1, rhomax12
                              nZ = nZ + 1 ! nZ=rho12+rhomax12*(rho123-1)+rhomax12*rhomax123*(rho13-1)
                              nU2 = i2 - rhomax12304
                              DO rho34 = 1, rhomax34
                                 nU2 = nU2 + rhomax12304 ! nU2=rho123+rhomax123*(rho04-1)+rhomax123*rhomax04*(rho34-1)
                                 DO rho1324 = 1, rhomax1324
                                    ninelm(rho12, rho34, 1:rhomax1234, rho13, rho24, rho1324) = &
                                       ninelm(rho12, rho34, 1:rhomax1234, rho13, rho24, rho1324) &
                                       + U1(rho1324, nU1)*Z(rho132, nZ)*U2(1:rhomax1234, nU2)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DEALLOCATE (U1, Z, U2)

         END DO
      END DO

   END SUBROUTINE calculate_9_lambda_mu

   SUBROUTINE calculate_u_coef_c(irrep1, irrep2, irrep, irrep3, irrep12, irrep23, &
                                 rhomaxa, rhomaxb, rhomaxc, rhomaxd, racah_ptr, info) BIND(C, NAME="calculate_u_coef")
      !------------------------------------------------------------------------------------------------------------------------
      !! Wrapper of subroutine calculate_u_coef for SU(3) recoupling U coefficients
      ! U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
      !! to be used by C++ wrapper
      !
      ! Input arguments: irrep1,irrep2,irrep,irrep3,irrep12,irrep23,rhomaxa,rhomaxb,rhomaxc,rhomaxd
      ! Output arguments: racah_ptr,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      ! rhomaxc is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23),
      ! rhomaxd is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu),
      ! racah_ptr(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1),
      !   is U coefficient for given rhoa,rhob,rhoc,rhod,
      ! info is 0 iff Lapack subroutine dgesv called by calculate_u_coef ran withou errors.
      !-------------------------------------------------------------------------------------------------------------------------
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep1
         !! SU(3) irrep 1
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep2
         !! SU(3) irrep 2
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep
         !! SU(3) irrep resulting from SU(3) couplings 12x3 and 1x23
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep3
         !! SU(3) irrep 3
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep12
         !! SU(3) irrep 12
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep23
         !! SU(3) irrep 23
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxa
         !! Multiplicity of SU(3) coupling 1x2->12
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxb
         !! Multiplicity of SU(3) coupling 12x3->irrep
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxc
         !! Multiplicity of SU(3) coupling 2x3->23
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxd
         !! Multiplicity of SU(3) coupling 1x23->irrep
      TYPE(C_PTR), INTENT(IN), VALUE :: racah_ptr
         !! Array of U coefficients.
         !! racah_ptr(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1),
         !! is U coefficient for given rhoa,rhob,rhoc,rhod.
      INTEGER(C_INT), INTENT(OUT) :: info
         !! Equals 0 iff Lapack subroutine dgesv called by calculate_u_coef ran withou errors.
      REAL(C_DOUBLE), POINTER, DIMENSION(:, :) :: rac
      INTEGER :: rhomaxabc
      rhomaxabc = rhomaxa*rhomaxb*rhomaxc
      CALL C_F_POINTER(racah_ptr, rac, [rhomaxd, rhomaxabc])
      CALL calculate_u_coef(irrep1, irrep2, irrep, irrep3, irrep12, irrep23, &
                            rhomaxa, rhomaxb, rhomaxc, rhomaxd, rac, rhomaxd, info)
   END SUBROUTINE calculate_u_coef_c

   SUBROUTINE calculate_z_coef_c(irrep2, irrep1, irrep, irrep3, irrep12, irrep13, &
                                 rhomaxa, rhomaxb, rhomaxc, rhomaxd, Z_ptr, info) BIND(C, NAME="calculate_z_coef")
      !---------------------------------------------------------------------------------------------------------
      !! Wrapper of subroutine calculate_z_coef for SU(3) recoupling Z coefficients
      ! Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
      !! to be used by C++ wrapper
      !
      ! Input arguments: irrep2,irrep1,irrep,irrep3,irrep12,irrep13,rhomaxa,rhomaxb,rhomaxc,rhomaxd
      ! Output argumrnts: Z_ptr,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomaxa is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomaxb is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu),
      ! rhomaxc is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      ! rhomaxd is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu),
      ! Z_ptr(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1),
      !   is Z coefficient for given rhoa,rhob,rhoc,rhod,
      ! info is 0 iff Lapack subroutine dgesv called by calculate_z_coef ran withou errors.
      !---------------------------------------------------------------------------------------------------------
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep1
         !! SU(3) irrep 1
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep2
         !! SU(3) irrep 2
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep
         !! SU(3) irrep resulting from SU(3) couplings 12x3 and 13x2
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep3
         !! SU(3) irrep 3
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep12
         !! SU(3) irrep 12
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep13
         !! SU(3) irrep 13
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxa
         !! Multiplicity of SU(3) coupling 1x2->12
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxb
         !! Multiplicity of SU(3) coupling 12x3->irrep
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxc
         !! Multiplicity of SU(3) coupling 1x3->13
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomaxd
         !! Multiplicity of SU(3) coupling 13x2->irrep
      INTEGER(C_INT), INTENT(OUT) :: info
         !! Equals 0 iff Lapack subroutine dgesv called by calculate_z_coeff ran withou errors.
      TYPE(C_PTR), INTENT(IN), VALUE :: Z_ptr
         !! Array of Z coefficients.
         !! Z_ptr(ind), where ind=rhod-1+rhomaxd*(rhoa-1)+rhomaxd*rhomaxa*(rhob-1)+rhomaxd*rhomaxa*rhomaxb*(rhoc-1),
         !! is Z coefficient for given rhoa,rhob,rhoc,rhod.
      REAL(C_DOUBLE), POINTER, DIMENSION(:, :) :: Zcoeff
      INTEGER :: rhomaxabc
      rhomaxabc = rhomaxa*rhomaxb*rhomaxc
      CALL C_F_POINTER(Z_ptr, Zcoeff, [rhomaxd, rhomaxabc])
      CALL calculate_z_coef(irrep2, irrep1, irrep, irrep3, irrep12, irrep13, &
                            rhomaxa, rhomaxb, rhomaxc, rhomaxd, Zcoeff, rhomaxd, info)
   END SUBROUTINE calculate_z_coef_c

   SUBROUTINE calculate_9_lambda_mu_c(irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep, &
                                      rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, ninelm_ptr, info) &
      BIND(C, NAME="calculate_9_lambda_mu")
      !---------------------------------------------------------------------------------------------------------------------------
      !! Wrapper of subroutine calculate_9_lambda_mu for 9-(lambda,mu) coefficients
      !
      ! | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
      ! | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
      ! |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
      ! |     rho13           rho24          rho1234            |
      !
      !! to be used by C++ wrapper
      !
      ! Input arguments: irrep1,irrep2,irrep12,irrep3,irrep4,irrep34,irrep13,irrep24,irrep,
      !                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
      ! Output arguments: ninelm_ptr,info
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rhomax12 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12),
      ! rhomax34 is multiplicity of SU(3) coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34),
      ! rhomax1234 is multiplicity of SU(3) coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu),
      ! rhomax13 is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13),
      ! rhomax24 is multiplicity of SU(3) coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24),
      ! rhomax1324 is multiplicity of SU(3) coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu),
      ! ninelm_ptr(ind), where ind=rho12-1+rhomax12*(rho34-1)+rhomax12*rhomax34*(rho1234-1)+rhomax12*rhomax34*rhomax1234*(rho13-1)
      !   +rhomax12*rhomax34*rhomax1234*rhomax13*(rho24-1)+rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*(rho1324-1),
      !   is 9-(lambda,mu) coefficient for given rho12,rho34,rho1234,rho13,rho24,rho1324,
      ! info is 0 iff Lapack subroutine dgesv called by calculate_9_lambda_mu ran without errors.
      !---------------------------------------------------------------------------------------------------------------------------
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep
      INTEGER(C_INT), INTENT(IN), VALUE :: rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324
      INTEGER(C_INT), INTENT(OUT) :: info
      TYPE(C_PTR), INTENT(IN), VALUE :: ninelm_ptr
      REAL(C_DOUBLE), POINTER, DIMENSION(:, :, :, :, :, :) :: ninelm
      CALL C_F_POINTER(ninelm_ptr, ninelm, [rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324])
      CALL calculate_9_lambda_mu(irrep1, irrep2, irrep12, irrep3, irrep4, irrep34, irrep13, irrep24, irrep, &
                                 rhomax12, rhomax34, rhomax1234, rhomax13, rhomax24, rhomax1324, ninelm, info)
   END SUBROUTINE calculate_9_lambda_mu_c

END MODULE ndsu3lib_recoupling
