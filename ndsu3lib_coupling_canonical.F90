!***********************************************************
!
! ndsu3lib_coupling_canonical.F90
!! Module for SU(3)-SU(2)xU(1) reduced coupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!***********************************************************
MODULE ndsu3lib_coupling_canonical
   !! Module for SU(3)-SU(2)xU(1) reduced coupling coefficients
   USE ndsu3lib_tools
   IMPLICIT NONE

CONTAINS

   SUBROUTINE calculate_coupling_canonical_extremal(irrep1, irrep2, irrep3, I3, rhomax, i2, wigner, p1a, p2a, q2a)
      !------------------------------------------------------------------------------------------------------------------
      !! Calculate extremal-weight reduced SU(3)-SU(2)xU(1) coupling coefficients
      ! / (lambda1,mu1)    (lambda2,mu2)   ||   (lambda3,mu3)   \
      ! \epsilon1,Lambda1 epsilon2,Lambda2 || epsilon3E,Lambda3E/rho
      !! for given lambda1,mu1,lambda2,mu2,lambda3,mu3.
      ! using Eq.(...) in Ref.[1]
      !
      ! References: [1] J.Herko et al. in preparation
      !             [2] D.Goldberg, ACM Computing Surveys, Vol.23, No.1 (1991) 5
      !
      ! Input arguments: irrep1,irrep2,irrep3,I3,rhomax
      ! Output arguments: i2,wigner,p1a,p2a,q2a
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! I3=1 for E=H, I3=0 for E=L,
      ! rhomax is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0.
      !
      !! Resulting reduced coupling coefficient is wigner(p1,p2,q2,rho), where
      !!   epsilon2=2*lambda2+mu2-3*(p2+q2) for I3=1, epsilon2=-2*mu2-lambda2+3*(p2+q2) for I3=0,
      !!   epsilon1=epsilon3E-epsilon2,
      !!   Lambda1=(mu1+p1-q1)/2 for I3=1, Lambda1=(lambda1+p1-q1)/2 for I3=0,
      !!   Lambda2=(mu2+p2-q2)/2 for I3=1, Lambda2=(lambda2+p2-q2)/2 for I3=0,
      !!   p1=p1a(i),
      !!   p2=p2a(i),
      !!   q2=q2a(i),
      !!   q1=(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3-p1-p2-q2 for I3=1,
      !!   q1=(2*(mu1+mu2+lambda3)+lambda1+lambda2+mu3)/3-p1-p2-q2 for I3=0,
      !!   1<=i<=i2.
      !------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! First SU(3) irrep in SU(3) coupling, i.e., (lambda1,mu1)
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! Second SU(3) irrep in SU(3) coupling, i.e., (lambda2,mu2)
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! Resulting SU(3) irrep, i.e., (lambda3,mu3)
      INTEGER, INTENT(IN) :: I3
         !! Parameter determining weight in resulting SU(3) irrep space (I3=1(0) for highest(lowest) weight)
      INTEGER, INTENT(IN) ::  rhomax
         !! Multiplicity of SU(3) coupling which must be greater than 0
      INTEGER, INTENT(OUT) :: i2
         !! Number of reduced coupling coefficients
      TYPE(su3irrep) :: irrep1a, irrep2a, irrep3a
      INTEGER :: lambda1, mu1, lambda2, mu2, lambda3, mu3, p2, q2, Lambda22, i1, j2, j1, p2tilde, q2tilde, i, j, n, a, b, c, d, &
                 epsilon2, steps21, Lambda12, eta, p1, q1, rho, i4, Lambda22max, ABCD, phiprhomax, p2min, p2max, &
                 noname1, p1min, noname2, p1max
      INTEGER(KIND=8) :: Sq2, Rp2
      REAL(KIND=8) :: F, G, H, scalprod
      INTEGER, DIMENSION(:), INTENT(OUT) :: p1a
         !! Array of values of p1.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
      INTEGER, DIMENSION(:), INTENT(OUT) :: p2a
         !! Array of values of p2.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
      INTEGER, DIMENSION(:), INTENT(OUT) :: q2a
         !! Array of values of q2.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
      REAL(KIND=8), DIMENSION(0:, 0:, 0:, 1:), INTENT(OUT) :: wigner
         !! Array of resulting reduced coupling coefficients.
         !! Sizes are at least (0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax).

      IF (I3 == 1) THEN ! E=H
         lambda1 = irrep1%lambda
         lambda2 = irrep2%lambda
         lambda3 = irrep3%lambda
         mu1 = irrep1%mu
         mu2 = irrep2%mu
         mu3 = irrep3%mu
         irrep1a = irrep1
         irrep3a = irrep3
      ELSE ! E=L
         ! Coefficients with E=L are obtained from those with E=H and flipped lambdas and mus.
         ! See Eq.(...) in [1].
         lambda1 = irrep1%mu
         lambda2 = irrep2%mu
         lambda3 = irrep3%mu
         mu1 = irrep1%lambda
         mu2 = irrep2%lambda
         mu3 = irrep3%lambda
         irrep1a%lambda = lambda1
         irrep1a%mu = mu1
         irrep3a%lambda = lambda3
         irrep3a%mu = mu3
      END IF

      wigner(0:lambda1, 0:lambda2, 0:mu2, 1:rhomax) = 0.D0
      Lambda22max = 0

      DO WHILE (MAX(mu2, (lambda1 + lambda2 - lambda3 + 2*(mu1 + mu2 - mu3))/3 + 1 + lambda2) > upbound_binom)
         CALL reallocate_binom(50)
      END DO

      DO rho = 1, rhomax

         eta = 0
         DO
            irrep2a%lambda = lambda2 - 1
            irrep2a%mu = mu2 - 1
            IF (outer_multiplicity(irrep1a, irrep2a, irrep3a) == rho - 1) EXIT
            lambda2 = lambda2 - 1
            mu2 = mu2 - 1
            eta = eta + 1
         END DO
         ! lambda2 is \bar{lambda2}, mu2 is \bar{mu2}
         !************************************************************************
         ! Calculation of
         ! /  (lambda1,mu1)    (\bar{lambda2},\bar{mu2}) ||   (lambda3,mu3)   \
         ! \epsilon1H,Lambda1H     epsilon2,Lambda2      || epsilon3H,Lambda3H/rho
         ! using Eq.(...) in [1]
         !************************************************************************
         epsilon2 = lambda1 + 2*mu1 - lambda3 - 2*mu3
         noname1 = (2*lambda2 + mu2 - epsilon2)/3
         p2max = MIN(lambda2, noname1, (lambda1 + lambda3 - mu2 + noname1)/2)
         p2min = MAX(0, noname1 - mu2, (lambda3 - lambda1 - mu2 + noname1 + 1)/2, (lambda1 - lambda3 - mu2 + noname1 + 1)/2)
         q2 = noname1 - p2min
         IF (p2min /= p2max) THEN
            n = (lambda1 + lambda2 - lambda3 + 2*(mu1 + mu2 - mu3))/3
            a = (lambda2 + lambda3 - lambda1 - n)/2
            b = (lambda3 + lambda1 - lambda2 + n + 2)/2
            c = (lambda1 + lambda2 - lambda3 - n)/2
            d = (lambda1 + lambda2 + lambda3 + 2 - n)/2
            j1 = lambda2 - p2max
            j2 = lambda2 - p2min
            DO p2 = p2min, p2max
               ! Upper and lower bound on p2 are such that:
               ! 1) 0<=p2<=lambda2
               ! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
               ! 3) ABS(lambda1-Lambda22)<=lambda3<=lambda1+Lambda22, where Lambda22=mu2-(2*lambda2+mu2-epsilon2)/3+2*p2
               p2tilde = mu2 - q2
               q2tilde = lambda2 - p2
               IF (p2tilde == 0) THEN
                  F = 1.D0
               ELSE
                  F = DBLE((a + 1)*(b - 1))
                  DO j = 1, p2tilde - 1
                     F = F*DBLE((a + j + 1)*(b - j - 1))
                  END DO
                  i1 = p2tilde*(p2tilde + 1)/2
                  DO i = 1, p2tilde
                     i1 = i1 + 1
                     scalprod = DBLE((p2 + 1)*(mu1 + lambda2 + mu2 - n + 2))
                     DO j = 1, p2tilde - 1
                        IF (j < i) THEN
                           scalprod = scalprod*DBLE((p2 + j + 1)*(mu1 + lambda2 + mu2 - n + j + 2))
                        ELSE
                           scalprod = scalprod*DBLE((a + j + 1)*(b - j - 1))
                        END IF
                     END DO
                     F = F + binom(i1)*scalprod
                  END DO
                  IF (BTEST(p2tilde, 0)) F = -F
               END IF
               IF (j1 < q2tilde) THEN
                  G = DBLE(INT((a + n - j1)*(b - n + j1), 8)*(c + n - j1)*(d + n - j1)*(lambda2 + mu2 - j1 + 1))
               ELSE
                  G = DBLE(mu2 - n + j1 + 1)
               END IF
               DO j = j1 + 1, j2 - 1
                  IF (j < q2tilde) THEN
                     G = G*DBLE(INT((a + n - j)*(b - n + j), 8)*(c + n - j)*(d + n - j)*(lambda2 + mu2 - j + 1))
                  ELSE
                     G = G*DBLE(mu2 - n + j + 1)
                  END IF
               END DO
               G = G*DBLE(lambda2 - 2*q2tilde + n + 1)*binom(n*(n + 1)/2 + q2tilde)
               i1 = n + 1 + lambda2 - q2tilde
               H = binom(i1*(i1 + 1)/2 + lambda2 - q2tilde)
               wigner(lambda1, p2, q2, rho) = F*DSQRT(G/H)
               q2 = q2 - 1
            END DO
         ELSE
            wigner(lambda1, p2min, q2, rho) = 1.D0
         END IF
         !************************************************************************
         ! Calculation of
         ! / (lambda1,mu1)   (\bar{lambda2},\bar{mu2}) ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1    epsilon2H,Lambda2H     || epsilon3H,Lambda3H/rho
         ! from
         ! /  (lambda1,mu1)    (\bar{lambda2},\bar{mu2}) ||   (lambda3,mu3)   \
         ! \epsilon1H,Lambda1H     epsilon2,Lambda2      || epsilon3H,Lambda3H/rho
         ! using Eq.(...) in [1]
         !************************************************************************
         noname2 = lambda1 + mu1 + 1
         steps21 = (epsilon2 + lambda2 + 2*mu2)/3 ! steps21 is the number of iterations in Eq.(...)
         DO i2 = 1, steps21
            noname1 = noname1 + 1
            noname2 = noname2 - 1
            p1min = MAX(0, noname2 - mu1, (lambda1 - i2 - mu1 + noname2 + 2)/2)
            p1max = MIN(lambda1, noname2 - 1, (lambda1 - 1 + i2 - mu1 + noname2)/2)
            p2min = MAX(0, noname1 - mu2, (lambda2 - steps21 + i2 - mu2 + noname1 + 1)/2)
            q2 = noname1 - p2min
            Lambda22 = mu2 + p2min - q2 - 2
            DO p2 = p2min, MIN(lambda2, noname1, (lambda2 + steps21 - i2 - mu2 + noname1)/2)
               ! Upper and lower bound on p2 are such that:
               ! 1) 0<=p2<=lambda2
               ! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
               ! 3) lambda2-steps21+i2<=Lambda22<=lambda2+steps21-i2, where Lambda22=mu2+p2-q2
               Lambda22 = Lambda22 + 2 ! Lambda22 is 2*Lambda2 in Eq.(...)
               Sq2 = q2*(mu2 + 1 - q2)*INT(lambda2 + mu2 + 2 - q2, 8) ! Sq2 is S(q2)
               Rp2 = p2*(lambda2 + 1 - p2)*INT(mu2 + 1 + p2, 8) ! Rp2 is R(p2)
               IF (lambda1 >= i2) THEN
                  IF (Lambda22 + 1 <= lambda2 + mu2 .AND. q2 > 0) THEN
                     ABCD = (lambda1 - i2 + Lambda22 - lambda3 + 2)*(lambda1 - i2 + Lambda22 + lambda3 + 4)
                     IF (ABCD > 0) &
                        wigner(p1min - 1, p2, q2, rho) &
                        = -DSQRT(DBLE((lambda1 - i2 + 1)*Sq2*ABCD) &
                                 /DBLE(INT((lambda1 + 2 - i2)*(Lambda22 + 1)*p1min, 8) &
                                       *(lambda1 + 1 - p1min)*(mu1 + 1 + p1min)*(Lambda22 + 2))) &
                          *wigner(p1min, p2, q2 - 1, rho)
                  END IF
                  IF (Lambda22 >= 1 .AND. p2 > 0) THEN
                     ABCD = (Lambda22 + lambda3 - lambda1 + i2)*(lambda3 + lambda1 - i2 - Lambda22 + 2)
                     IF (ABCD > 0) &
                        wigner(p1min - 1, p2, q2, rho) &
                        = wigner(p1min - 1, p2, q2, rho) &
                          - DSQRT(DBLE((lambda1 - i2 + 1)*Rp2*ABCD) &
                                  /DBLE(INT((lambda1 + 2 - i2)*(Lambda22 + 1)*p1min, 8) &
                                        *(lambda1 + 1 - p1min)*(mu1 + 1 + p1min)*Lambda22)) &
                          *wigner(p1min, p2 - 1, q2, rho)
                  END IF
               END IF
               q1 = noname2 - p1min
               Lambda12 = mu1 + p1min - q1 - 2
               DO p1 = p1min, p1max
                  ! Upper and lower bound on p1 are such that:
                  ! 1) 0<=p1<=lambda1
                  ! 2) 0<q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1 cannot be 0 to avoid division by 0
                  ! 3) lambda1+1-i2<=Lambda12<=lambda1-1+i2, where Lambda12=mu1+p1-q1
                  Lambda12 = Lambda12 + 2
                  IF (Lambda22 + 1 <= lambda2 + mu2 .AND. q2 > 0) THEN
                     ABCD = (Lambda22 + lambda3 - Lambda12 + 1)*(lambda3 + Lambda12 - Lambda22 + 1)
                     IF (ABCD > 0) &
                        wigner(p1, p2, q2, rho) &
                        = -DSQRT(DBLE((Lambda12 + 2)*Sq2*ABCD) &
                                 /DBLE(INT((Lambda12 + 1)*(Lambda22 + 1)*q1, 8) &
                                       *(mu1 + 1 - q1)*(lambda1 + mu1 + 2 - q1)*(Lambda22 + 2))) &
                          *wigner(p1, p2, q2 - 1, rho)
                  END IF
                  IF (Lambda22 >= 1 .AND. p2 > 0) THEN
                     ABCD = (Lambda12 + Lambda22 - lambda3 + 1)*(Lambda12 + Lambda22 + lambda3 + 3)
                     IF (ABCD > 0) &
                        wigner(p1, p2, q2, rho) &
                        = wigner(p1, p2, q2, rho) &
                          + DSQRT(DBLE((Lambda12 + 2)*Rp2*ABCD) &
                                  /DBLE(INT((Lambda12 + 1)*(Lambda22 + 1)*q1, 8) &
                                        *(mu1 + 1 - q1)*(lambda1 + mu1 + 2 - q1)*Lambda22)) &
                          *wigner(p1, p2 - 1, q2, rho)
                  END IF
                  q1 = q1 - 1
               END DO
               q2 = q2 - 1
            END DO
         END DO
         !**********************************************************************
         ! Calculation of
         ! / (lambda1,mu1)     (lambda2,mu2)    ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1 epsilon2H,Lambda2H || epsilon3H,Lambda3H/rho
         ! from
         ! / (lambda1,mu1)   (\bar{lambda2},\bar{mu2}) ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1    epsilon2H,Lambda2H     || epsilon3H,Lambda3H/rho
         ! using Eq.(...) in [1]
         !**********************************************************************
         noname2 = noname2 - 1
         DO i1 = 1, eta
            lambda2 = lambda2 + 1
            mu2 = mu2 + 1
            noname2 = noname2 - 1
            IF (noname2 == mu1 .AND. wigner(1, lambda2 - 1, mu2 - 1, rho) /= 0.D0) &
               wigner(0, lambda2, mu2, rho) &
               = -DSQRT(DBLE(INT(lambda1*(mu1 + 2), 8) &
                             *(lambda2 + lambda3)*(lambda3 - lambda2 + 2))/2.D0) &
                 *wigner(1, lambda2 - 1, mu2 - 1, rho)
            p1min = MAX(0, noname2 - mu1, (2 - mu1 + noname2)/2)
            q1 = noname2 - p1min
            Lambda12 = mu1 + p1min - q1 - 2
            DO p1 = p1min, MIN(lambda1, noname2, (lambda1 - 1 + noname2)/2)
               ! Upper and lower bound on p1 are such that:
               ! 1) 0<=p1<=lambda1
               ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
               ! 3) 1<=Lambda12<=lambda1+mu1-1, where Lambda12=mu1+p1-q1
               Lambda12 = Lambda12 + 2
               IF (p1 < lambda1 .AND. wigner(p1 + 1, lambda2 - 1, mu2 - 1, rho) /= 0.D0) &
                  wigner(p1, lambda2, mu2, rho) &
                  = -DSQRT(DBLE(INT((p1 + 1)*(lambda1 - p1), 8) &
                                *(mu1 + 2 + p1)*(lambda2 + lambda3 - Lambda12) &
                                *(lambda3 + Lambda12 - lambda2 + 2)) &
                           /DBLE((Lambda12 + 2)*(Lambda12 + 1))) &
                    *wigner(p1 + 1, lambda2 - 1, mu2 - 1, rho)
               IF (wigner(p1, lambda2 - 1, mu2 - 1, rho) /= 0.D0) &
                  wigner(p1, lambda2, mu2, rho) &
                  = wigner(p1, lambda2, mu2, rho) &
                    + DSQRT(DBLE(INT((q1 + 1)*(mu1 - q1), 8) &
                                 *(lambda1 + mu1 + 1 - q1)*(Lambda12 + lambda2 - lambda3) &
                                 *(Lambda12 + lambda2 + lambda3 + 2)) &
                            /DBLE((Lambda12 + 1)*Lambda12)) &
                    *wigner(p1, lambda2 - 1, mu2 - 1, rho)
               q1 = q1 - 1
            END DO
            IF (noname2 == lambda1 .AND. wigner(lambda1, lambda2 - 1, mu2 - 1, rho) /= 0.D0) THEN
               Lambda12 = lambda1 + mu1
               wigner(lambda1, lambda2, mu2, rho) &
                  = DSQRT(DBLE(INT(mu1*(lambda1 + mu1 + 1), 8) &
                               *(Lambda12 + lambda2 - lambda3)*(Lambda12 + lambda2 + lambda3 + 2)) &
                          /DBLE((Lambda12 + 1)*Lambda12)) &
                    *wigner(lambda1, lambda2 - 1, mu2 - 1, rho)
            END IF
         END DO

         ! Now arrays p1a,p2a,q2a are filled with values corresponding to coefficients
         ! / (lambda1,mu1)     (lambda2,mu2)    ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1 epsilon2H,Lambda2H || epsilon3H,Lambda3H/
         IF (rho == rhomax) THEN
            i2 = 0
            p1min = MAX(0, noname2 - mu1, (lambda3 - mu1 + noname2 - lambda2 + 1)/2, (lambda2 - lambda3 - mu1 + noname2 + 1)/2)
            DO p1 = p1min, MIN(lambda1, noname2, (lambda2 + lambda3 - mu1 + noname2)/2)
               ! Upper and lower bound on p1 are such that:
               ! 1) 0<=p1<=lambda1
               ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
               ! 3) ABS(Lambda12-lambda2)<=lambda3<=Lambda12+lambda2, where Lambda12=mu1+p1-q1
               i2 = i2 + 1
               p1a(i2) = p1
               p2a(i2) = lambda2
               q2a(i2) = mu2
            END DO
            IF (-lambda1 - 2*mu1 + 3*(steps21 + eta) == -lambda1 - 2*mu1) Lambda22max = lambda2
            ! Lambda22max is greatest value of 2*Lambda2 for which epsilon1,Lambda1) are of highest weight.
            ! This is needed for setting the phase.
         END IF
         !***************************************************************
         ! Calculation of
         ! / (lambda1,mu1)    (lambda2,mu2)   ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1 epsilon2,Lambda2 || epsilon3H,Lambda3H/rho
         ! from
         ! / (lambda1,mu1)     (lambda2,mu2)    ||   (lambda3,mu3)   \
         ! \epsilon1,Lambda1 epsilon2H,Lambda2H || epsilon3H,Lambda3H/rho
         ! using Eq.(...) in [1]
         !***************************************************************
         noname1 = lambda2 + mu2
         DO i1 = 1, lambda2 + mu2 ! This is loop over epsilon2=epsilon2HW+3,...,2*lambda2+mu2; i1 is (epsilon2-epsilon2HW)/3
            noname2 = noname2 + 1
            IF (noname2 > lambda1 + mu1) EXIT
            noname1 = noname1 - 1
            p2min = MAX(0, noname1 - mu2, (lambda2 - i1 - mu2 + noname1 + 1)/2)
            q2 = noname1 - p2min
            Lambda22 = mu2 + p2min - q2 - 2 ! Lambda22 is 2*Lambda2' in Eq.(18)
            DO p2 = p2min, MIN(lambda2, noname1, (lambda2 + i1 - mu2 + noname1)/2)
               Lambda22 = Lambda22 + 2
               p1min = MAX(0, noname2 - mu1, (lambda3 - mu1 + noname2 - Lambda22 + 1)/2, (Lambda22 - lambda3 - mu1 + noname2 + 1)/2)
               q1 = noname2 - p1min
               Lambda12 = mu1 + p1min - q1 - 2
               DO p1 = p1min, MIN(lambda1, noname2, (Lambda22 + lambda3 - mu1 + noname2)/2)
                  ! Upper and lower bound on p1 are such that:
                  ! 1) 0<=p1<=lambda1
                  ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                  ! 3) ABS(Lambda12-Lambda22)<=lambda3<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
                  Lambda12 = Lambda12 + 2
                  IF (Lambda22 == lambda2 - i1 .OR. Lambda22 == 0) THEN
                     IF (Lambda12 >= 1 .AND. p1 /= 0) THEN
                        wigner(p1, p2, q2, rho) &
                           = -DSQRT(DBLE(INT((Lambda22 + 1)*p1*(lambda1 + 1 - p1), 8) &
                                         *(mu1 + 1 + p1)*(Lambda22 + lambda3 - Lambda12 + 2) &
                                         *(lambda3 + Lambda12 - Lambda22)) &
                                    /DBLE(INT((Lambda12 + 1)*(Lambda22 + 2)*(p2 + 1), 8) &
                                          *(lambda2 - p2)*(mu2 + 2 + p2)*4*Lambda12)) &
                             *wigner(p1 - 1, p2 + 1, q2, rho)
                     ELSE
                        wigner(p1, p2, q2, rho) = 0.D0
                     END IF
                     IF (Lambda12 + 1 <= lambda1 + mu1 .AND. q1 /= 0) &
                        wigner(p1, p2, q2, rho) &
                        = wigner(p1, p2, q2, rho) &
                          + DSQRT(DBLE(INT((Lambda22 + 1)*q1*(mu1 + 1 - q1), 8) &
                                       *(lambda1 + mu1 + 2 - q1)*(Lambda12 + Lambda22 - lambda3 + 2) &
                                       *(Lambda12 + Lambda22 + lambda3 + 4)) &
                                  /DBLE(INT((Lambda12 + 1)*(Lambda22 + 2)*(p2 + 1), 8) &
                                        *(lambda2 - p2)*(mu2 + 2 + p2)*4*(Lambda12 + 2))) &
                          *wigner(p1, p2 + 1, q2, rho)
                  ELSE IF (q2 < mu2) THEN
                     IF (Lambda12 >= 1 .AND. p1 /= 0) THEN
                        wigner(p1, p2, q2, rho) &
                           = -DSQRT(DBLE(INT((Lambda22 + 1)*p1*(lambda1 + 1 - p1), 8) &
                                         *(mu1 + 1 + p1)*(Lambda12 + Lambda22 - lambda3) &
                                         *(Lambda12 + Lambda22 + lambda3 + 2)) &
                                    /DBLE(INT((Lambda12 + 1)*Lambda22*(q2 + 1), 8) &
                                          *(mu2 - q2)*(lambda2 + mu2 + 1 - q2)*4*Lambda12)) &
                             *wigner(p1 - 1, p2, q2 + 1, rho)
                     ELSE
                        wigner(p1, p2, q2, rho) = 0.D0
                     END IF
                     IF (Lambda12 + 1 <= lambda1 + mu1 .and. q1 /= 0) &
                        wigner(p1, p2, q2, rho) &
                        = wigner(p1, p2, q2, rho) &
                          - DSQRT(DBLE(INT((Lambda22 + 1)*q1*(mu1 + 1 - q1), 8) &
                                       *(lambda1 + mu1 + 2 - q1)*(Lambda22 + lambda3 - Lambda12) &
                                       *(lambda3 + Lambda12 - Lambda22 + 2)) &
                                  /DBLE(INT((Lambda12 + 1)*Lambda22*(q2 + 1), 8) &
                                        *(mu2 - q2)*(lambda2 + mu2 + 1 - q2)*4*(Lambda12 + 2))) &
                          *wigner(p1, p2, q2 + 1, rho)
                  END IF

                  IF (rho == rhomax) THEN
                     i2 = i2 + 1
                     IF (p1 == lambda1 .AND. q1 == mu1) Lambda22max = MAX(Lambda22max, Lambda22)
                     ! Lambda22max is greatest value of 2*Lambda2 for which epsilon1,Lambda1 are of highest weight.
                     ! This is needed for setting the phase.
                     p1a(i2) = p1
                     p2a(i2) = p2
                     q2a(i2) = q2
                  END IF

                  q1 = q1 - 1
               END DO
               q2 = q2 - 1
            END DO
         END DO

         ! Valid values of p1, p2, and q2 are elements of arrays p1a, p2a, and q2a with indeces from 1 to i2.
      END DO ! End of loop over rho
      !*****************************************
      ! Orthonormalization according to Eq.(...)
      !*****************************************
      ! Gram-Schmidt orthonormalization is performed, where coupling coefficients for given rho represent vector,
      ! whose components are indexed by epsilon1,Lambda1,espilon2,Lambda2, and scalar product is defined in Eq.(...).
      DO rho = 1, rhomax
         DO i4 = 1, rho
            p1 = p1a(1)
            p2 = p2a(1)
            q2 = q2a(1)
            scalprod = wigner(p1, p2, q2, rho)*wigner(p1, p2, q2, i4)
            F = 0.D0
            DO i1 = 2, i2 ! Kahan summation fomula (see Ref.[2]) is used to calculate scalar product.
               p1 = p1a(i1)
               p2 = p2a(i1)
               q2 = q2a(i1)
               G = wigner(p1, p2, q2, rho)*wigner(p1, p2, q2, i4) - F
               H = scalprod + G
               F = H - scalprod
               F = F - G
               scalprod = H
            END DO
            IF (i4 /= rho) THEN
               DO i1 = 1, i2
                  p1 = p1a(i1)
                  p2 = p2a(i1)
                  q2 = q2a(i1)
                  wigner(p1, p2, q2, rho) = wigner(p1, p2, q2, rho) - scalprod*wigner(p1, p2, q2, i4)
               END DO
            ELSE
               scalprod = 1.D0/DSQRT(scalprod)
               DO i1 = 1, i2
                  p1 = p1a(i1)
                  p2 = p2a(i1)
                  q2 = q2a(i1)
                  wigner(p1, p2, q2, rho) = wigner(p1, p2, q2, rho)*scalprod
               END DO
            END IF
         END DO
      END DO
      !****************************************
      ! Setting the phase according to Eq.(...)
      !****************************************
      phiprhomax = lambda1 + lambda2 - lambda3 + mu1 + mu2 - mu3 + rhomax ! phiprhomax is phi+rhomax
      i4 = phiprhomax + (lambda1 + Lambda22max - lambda3)/2 
      ! i4 is lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax+(lambda1+2*Lambda2_max-lambda3)/2
      p2 = (noname1 - mu2 + Lambda22max)/2
      q2 = noname1 - p2
      DO rho = 1, rhomax
         IF ((BTEST(i4 - rho + 1, 0) .AND. wigner(lambda1, p2, q2, rho) < 0.D0) &
             .OR. (BTEST(i4 - rho, 0) .AND. wigner(lambda1, p2, q2, rho) > 0.D0)) THEN
            DO i1 = 1, i2
               p1 = p1a(i1)
               p2 = p2a(i1)
               q2 = q2a(i1)
               wigner(p1, p2, q2, rho) = -wigner(p1, p2, q2, rho)
            END DO
         END IF
      END DO

      IF (I3 == 0) THEN ! E=L
         DO i1 = 1, i2
            p1 = p1a(i1)
            p2 = p2a(i1)
            DO rho = 1, rhomax
               ! See Eq.(35,2B) in Draayer & Akiyama
               IF (BTEST(phiprhomax - rho + p1 + p2 &
                         + (mu1 - (2*(lambda1 + lambda2 + mu3) + mu1 + mu2 + lambda3)/3 + mu2 - lambda3)/2, 0)) THEN
                  q2 = q2a(i1)
                  wigner(p1, p2, q2, rho) = -wigner(p1, p2, q2, rho)
               END IF
            END DO
         END DO
      END IF
   END SUBROUTINE calculate_coupling_canonical_extremal

   SUBROUTINE calculate_coupling_canonical_nonextremal(irrep1, irrep2, irrep3, epsilon3, Lambda32, I3, rhomax,&
                                                       numb, wignerex, wigner, p1a, p2a, q2a)
      !------------------------------------------------------------------------------------------------------------------
      !! Calculate non-extremal-weight reduced SU(3)-SU(2)xU(1) coupling coefficients
      ! / (lambda1,mu1)    (lambda2,mu2)   ||  (lambda3,mu3)  \
      ! \epsilon1,Lambda1 epsilon2,Lambda2 || epsilon3,Lambda3/rho
      !! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3 from extremal-weight coefficients.
      ! using Eq.(...) in [1] and Table 9.1 on page 311 in [2].
      !
      ! References: [1] J.Herko et al. in preparation
      !             [2] Varshalovich, Quantum theory of angular momentum
      !
      ! Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,I3,rhomax,wignerex
      ! Input and output arguments: numb,p1a,p2a,q2a
      ! Output argument: wigner
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! Lambda32 is 2*Lambda3 (epsilon3 and Lambda32 must be valid),
      ! rhomax is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0,
      ! wignerex is array of extremal-weight coefficients calculated by subroutine calculate_coupling_canonical_extremal.
      ! If I3=1, coefficients are calculated from highest-weight (HW) coefficients.
      ! If I3=0, coefficients are calculated from lowest-weight (LW) coefficients.
      !
      !! Resulting reduced coupling coefficient is wigner(p1,p2,q2,rho),
      !!   where epsilon2=2*lambda2+mu2-3*(p2+q2),
      !!         epsilon1=epsilon3-epsilon2,
      !!         Lambda1=(mu1+p1-q1)/2,
      !!         Lambda2=(mu2+p2-q2)/2,
      !!         p1=p1a(i),
      !!         p2=p2a(i),
      !!         q2=q2a(i),
      !!         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2,
      !!         1<=i<=numb,
      !!   unless I3=0 and epsilon3=epsilon3L=2*lambda3+mu3, in this case:
      !!         q1=mu1-p1a(i),
      !!         p2=lambda2-q2a(i),
      !!         q2=mu2-p2a(i),
      !!         p1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-q1-p2-q2.
      !------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! First SU(3) irrep in SU(3) coupling, i.e., (lambda1,mu1)
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! Second SU(3) irrep in SU(3) coupling, i.e., (lambda2,mu2)
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! Resulting SU(3) irrep, i.e., (lambda3,mu3)
      INTEGER, INTENT(IN) :: epsilon3
         !! epsilon3.
         !! Must be valid
      INTEGER, INTENT(IN) :: Lambda32
         !! 2*Lambda3.
         !! Must be valid
      INTEGER, INTENT(IN) :: I3
         !! Parameter determining from what extremal-weight reduced coupling coefficients desired coefficient are calculated.
         !! If I3=1(0), coefficients are calculated from highest(lowest)-weight coefficients.
      INTEGER, INTENT(IN) :: rhomax
         !! Multiplicity of SU(3) coupling which must be greater than 0
      INTEGER :: numb
         !! Number of extremal-weight(resulting) reduced coupling coefficients for given outer multiplicity index as input(output) argument
      INTEGER :: eps3, rho, p1, q1, p2, q2, Lam32, epsilon2max, epsilon3ex, Lambda22, Lambda12, Lam32prime, p3, q3, &
                 pq1, noname1, noname2, s2, lm1, lm2, lm3, mu1p, mu2p, mu3p, q2ex, p1ex, lambda1p, lambda2p, lambda3p
      INTEGER(KIND=8) :: Rp2, Sq2
      REAL(KIND=8) :: N3
      INTEGER, DIMENSION(:) :: p1a
         !! Array of values of p1 for extremal-weight(resulting) reduced coupling coefficients as input(output) argument.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
      INTEGER, DIMENSION(:) :: p2a
         !! Array of values of p2 for extremal-weight(resulting) reduced coupling coefficients as input(output) argument.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
      INTEGER, DIMENSION(:) :: q2a
         !! Array of values of q2 for extremal-weight(resulting) reduced coupling coefficients as input(output) argument.
         !! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
      REAL(KIND=8), DIMENSION(0:, 0:, 0:, 1:), INTENT(IN) :: wignerex
         !! Array of extremal-weight reduced coupling coefficients calculated by subroutine wigner_canonical_extremal.
         !! Sizes are at least (0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax).
      REAL(KIND=8), DIMENSION(0:, 0:, 0:, 1:), INTENT(OUT) :: wigner
         !! Array of resulting reduced coupling coefficients.
         !! Sizes are at least (0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax).

      IF (I3 == 1) THEN

         DO rho = 1, rhomax
            DO q2 = 0, irrep2%mu
               DO p2 = 0, irrep2%lambda
                  wigner(0:irrep1%lambda, p2, q2, rho) = wignerex(0:irrep1%lambda, p2, q2, rho)
               END DO
            END DO
         END DO

         epsilon3ex = -irrep3%lambda - 2*irrep3%mu
         IF (epsilon3 == epsilon3ex) RETURN

         lm1 = irrep1%lambda + irrep1%mu + 1
         lm2 = irrep2%lambda + irrep2%mu + 1
         lm3 = irrep3%lambda + irrep3%mu + 1
         mu1p = irrep1%mu + 2
         mu2p = irrep2%mu + 2
         mu3p = irrep3%mu + 2

         epsilon2max = 2*irrep2%lambda + irrep2%mu
         noname1 = (epsilon2max - irrep1%lambda - 2*irrep1%mu - epsilon3ex)/3
         noname2 = (epsilon2max + 2*irrep1%lambda + irrep1%mu - epsilon3ex)/3
         epsilon3ex = epsilon3ex + 3

         numb = 0

         Lam32 = irrep3%lambda
         p3 = irrep3%lambda
         q3 = irrep3%mu

         DO eps3 = epsilon3ex, epsilon3, 3 ! eps3 is epsilon3 in Eq.(..)
            noname1 = noname1 - 1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
            noname2 = noname2 - 1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3

            Lam32prime = Lam32 ! Lam32prime is 2*Lambda3' in Eq.(...)
            IF (Lam32 < Lambda32) THEN
               Lam32 = Lam32 + 1
               q3 = q3 - 1 ! p3 and q3 correspond to eps3 and Lam32
               N3 = DSQRT(DBLE(Lam32 + 1)/DBLE(INT((q3 + 1)*(irrep3%mu - q3), 8)*(lm3 - q3)*(Lam32prime + 2)*4))
            ELSE
               Lam32 = Lam32 - 1
               p3 = p3 - 1 ! p3 and q3 correspond to eps3 and Lam32
               N3 = DSQRT(DBLE(Lam32 + 1)/DBLE(INT((p3 + 1)*(irrep3%lambda - p3), 8)*(mu3p + p3)*Lam32prime*4))
            END IF ! Lam32 is 2*Lambda3 in Eq.(...) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.

            DO p2 = MAX(0, noname1 - irrep2%mu, (Lam32 - irrep1%mu + noname2 - irrep2%mu)/2 - irrep1%lambda), &
               MIN(irrep2%lambda, noname2)
               Rp2 = INT(p2 + 1, 8)*(irrep2%lambda - p2)*(mu2p + p2)
               q2ex = MAX(0, noname1 - p2, (irrep2%mu + noname2 - Lam32 - irrep1%mu)/2 - irrep1%lambda)
               pq1 = noname2 - p2 - q2ex ! pq1 is p1+q1
               Lambda22 = irrep2%mu + p2 - q2ex ! Lambda22 is 2*Lambda2 in Eq.(...)
               DO q2 = q2ex, MIN(irrep2%mu, noname2 - p2, (irrep2%mu + noname2 + Lam32 - irrep1%mu)/2)
                  ! Lower and upper bounds on q2 are such that:
                  ! 1) 0<=q2<=mu2
                  ! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
                  !     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(...)
                  !     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(...)
                  Sq2 = INT(q2 + 1, 8)*(irrep2%mu - q2)*(lm2 - q2)
                  p1ex = MAX(0, pq1 - irrep1%mu, (Lam32 - irrep1%mu + pq1 - Lambda22)/2, (Lambda22 - Lam32 - irrep1%mu + pq1)/2)
                  q1 = pq1 - p1ex
                  Lambda12 = irrep1%mu + p1ex - q1 ! Lambda12 is 2*Lambda1 in Eq.(...)
                  DO p1 = p1ex, MIN(irrep1%lambda, pq1, (Lambda22 + Lam32 - irrep1%mu + pq1)/2)
                     ! Lower and upper bounds on p1 are such that:
                     ! 1) 0<=p1<=lambda1
                     ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                     ! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

                     IF (eps3 == epsilon3) THEN
                        numb = numb + 1
                        p1a(numb) = p1
                        p2a(numb) = p2
                        q2a(numb) = q2
                     END IF

                     s2 = Lambda12 + Lambda22 + Lam32prime + 1
                     IF (Lam32 > Lam32prime) THEN

                        IF (q1 /= irrep1%mu) THEN
                           wigner(p1, p2, q2, 1:rhomax) &
                              = DSQRT(DBLE(INT((q1 + 1)*(irrep1%mu - q1), 8) &
                                           *(lm1 - q1)*(s2 + 2)*(s2 - 2*Lambda22)) &
                                      /DBLE(Lambda12*(Lambda12 + 1))) &
                                *wigner(p1, p2, q2, 1:rhomax)
                        ELSE
                           wigner(p1, p2, q2, 1:rhomax) = 0.D0
                        END IF

                        IF (p1 /= irrep1%lambda) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(INT((p1 + 1)*(irrep1%lambda - p1), 8) &
                                          *(mu1p + p1)*(s2 - 2*Lambda12)*(s2 - 2*Lam32prime)) &
                                     /DBLE((Lambda12 + 1)*(Lambda12 + 2))) &
                             *wigner(p1 + 1, p2, q2, 1:rhomax)

                        IF (p2 /= irrep2%lambda) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             - DSQRT(DBLE(Rp2*(s2 - 2*Lam32prime)*(s2 - 2*Lambda22)) &
                                     /DBLE((Lambda22 + 1)*(Lambda22 + 2))) &
                             *wigner(p1, p2 + 1, q2, 1:rhomax)

                        IF (q2 /= irrep2%mu) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Sq2*(s2 + 2)*(s2 - 2*Lambda12)) &
                                     /DBLE(Lambda22*(Lambda22 + 1))) &
                             *wigner(p1, p2, q2 + 1, 1:rhomax)

                     ELSE

                        IF (q1 /= irrep1%mu) THEN
                           wigner(p1, p2, q2, 1:rhomax) &
                              = -DSQRT(DBLE(INT((q1 + 1)*(irrep1%mu - q1), 8) &
                                            *(lm1 - q1)*(s2 - 2*Lambda12)*(s2 - 2*Lam32prime)) &
                                       /DBLE(Lambda12*(Lambda12 + 1))) &
                                *wigner(p1, p2, q2, 1:rhomax)
                        ELSE
                           wigner(p1, p2, q2, 1:rhomax) = 0.D0
                        END IF

                        IF (p1 /= irrep1%lambda) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(INT((p1 + 1)*(irrep1%lambda - p1), 8) &
                                          *(mu1p + p1)*(s2 + 2)*(s2 - 2*Lambda22)) &
                                     /DBLE((Lambda12 + 1)*(Lambda12 + 2))) &
                             *wigner(p1 + 1, p2, q2, 1:rhomax)

                        IF (p2 /= irrep2%lambda) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Rp2*(s2 + 2)*(s2 - 2*Lambda12)) &
                                     /DBLE((Lambda22 + 1)*(Lambda22 + 2))) &
                             *wigner(p1, p2 + 1, q2, 1:rhomax)

                        IF (q2 /= irrep2%mu) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Sq2*(s2 - 2*Lam32prime)*(s2 - 2*Lambda22)) &
                                     /DBLE(Lambda22*(Lambda22 + 1))) &
                             *wigner(p1, p2, q2 + 1, 1:rhomax)

                     END IF

                     wigner(p1, p2, q2, 1:rhomax) = wigner(p1, p2, q2, 1:rhomax)*N3

                     q1 = q1 - 1
                     Lambda12 = Lambda12 + 2 ! Lambda12 is 2*Lambda1 in Eq.(...)
                  END DO
                  pq1 = pq1 - 1 ! pq1 is p1+q1
                  Lambda22 = Lambda22 - 1 ! Lambda22 is 2*Lambda2 in Eq.(...)
               END DO
            END DO
         END DO

      ELSE

         epsilon3ex = 2*irrep3%lambda + irrep3%mu
         epsilon2max = 2*irrep2%lambda + irrep2%mu
         noname1 = (epsilon2max - irrep1%lambda - 2*irrep1%mu - epsilon3ex)/3
         noname2 = (epsilon2max + 2*irrep1%lambda + irrep1%mu - epsilon3ex)/3

         DO p1 = 0, irrep1%lambda
            DO p2 = 0, irrep2%lambda
               DO q2 = 0, irrep2%mu
                  ! q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3max)/3-p1-p2-q2
                  q1 = noname2 - p1 - p2 - q2
                  DO rho = 1, rhomax
                     N3 = wignerex(irrep1%mu - q1, irrep2%mu - q2, irrep2%lambda - p2, rho)
                     IF (N3 == N3) THEN
                        wigner(p1, p2, q2, rho) = N3
                     ELSE
                        wigner(p1, p2, q2, rho) = 0.D0
                     END IF
                  END DO
               END DO
            END DO
         END DO

         IF (epsilon3 == epsilon3ex) RETURN

         lambda1p = irrep1%lambda + 1
         lambda2p = irrep2%lambda + 1
         lambda3p = irrep3%lambda + 1
         mu1p = irrep1%mu + 1
         mu2p = irrep2%mu + 1
         mu3p = irrep3%mu + 1
         lm1 = lambda1p + mu1p
         lm2 = lambda2p + mu2p
         lm3 = lambda3p + mu3p

         epsilon3ex = epsilon3ex - 3

         numb = 0

         Lam32 = irrep3%mu
         p3 = 0
         q3 = 0

         DO eps3 = epsilon3ex, epsilon3, -3 ! eps3 is epsilon3 in Eq.(...)
            noname1 = noname1 + 1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
            noname2 = noname2 + 1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3

            Lam32prime = Lam32 ! Lam32prime is 2*Lambda3' in Eq.(...)
            IF (Lam32 < Lambda32) THEN
               Lam32 = Lam32 + 1
               p3 = p3 + 1 ! p3 and q3 correspond to eps3 and Lam32
               N3 = DSQRT(DBLE(Lam32 + 1)/DBLE(INT(p3*(lambda3p - p3), 8)*(mu3p + p3)*(Lam32prime + 2)*4))
               ! N3 is sqrt((2*Lambda3+1)/(4*(2*Lambda3'+1)))/N3
            ELSE
               Lam32 = Lam32 - 1
               q3 = q3 + 1 ! p3 and q3 correspond to eps3 and Lam32
               N3 = DSQRT(DBLE(Lam32 + 1)/DBLE(INT(q3*(mu3p - q3), 8)*(lm3 - q3)*Lam32prime*4))
               ! N3 is sqrt((2*Lambda3+1)/(8*Lambda3'))/N3
            END IF ! Lam32 is 2*Lambda3 in Eq.(...)

            DO p2 = MIN(irrep2%lambda, noname2), &
               MAX(0, noname1 - irrep2%mu, (Lam32 - irrep1%mu + noname2 - irrep2%mu)/2 - irrep1%lambda), -1
               Rp2 = INT(p2, 8)*(lambda2p - p2)*(mu2p + p2) ! Rp2 is R(p2)
               q2ex = MIN(irrep2%mu, noname2 - p2, (irrep2%mu + noname2 + Lam32 - irrep1%mu)/2)
               pq1 = noname2 - p2 - q2ex ! pq1 is p1+q1
               Lambda22 = irrep2%mu + p2 - q2ex ! Lambda22 is 2*Lambda2 in Eq.(...)
               DO q2 = q2ex, MAX(0, noname1 - p2, (irrep2%mu + noname2 - Lam32 - irrep1%mu)/2 - irrep1%lambda), -1
                  ! Lower and upper bounds on q2 are such that:
                  ! 1) 0<=q2<=mu2
                  ! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
                  !     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(...)
                  !     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(...)
                  Sq2 = INT(q2, 8)*(mu2p - q2)*(lm2 - q2)
                  p1ex = MIN(irrep1%lambda, pq1, (Lambda22 + Lam32 - irrep1%mu + pq1)/2)
                  q1 = pq1 - p1ex
                  Lambda12 = irrep1%mu + p1ex - q1 ! Lambda12 is 2*Lambda1 in Eq.(...)
                  DO p1 = p1ex, &
                     MAX(0, pq1 - irrep1%mu, (Lam32 - irrep1%mu + pq1 - Lambda22)/2, (Lambda22 - Lam32 - irrep1%mu + pq1)/2), -1
                     ! Lower and upper bounds on p1 are such that:
                     ! 1) 0<=p1<=lambda1
                     ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                     ! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

                     IF (eps3 == epsilon3) THEN
                        numb = numb + 1
                        p1a(numb) = p1
                        p2a(numb) = p2
                        q2a(numb) = q2
                     END IF

                     s2 = Lambda12 + Lambda22 + Lam32prime + 1
                     IF (Lam32 > Lam32prime) THEN

                        ! term with Lambda1p=Lambda1+1/2
                        IF (q1 /= 0) THEN
                           wigner(p1, p2, q2, 1:rhomax) &
                              = -DSQRT(DBLE(INT(q1*(mu1p - q1), 8) &
                                            *(lm1 - q1)*(s2 - 2*Lambda12)*(s2 - 2*Lam32prime)) &
                                       /DBLE((Lambda12 + 1)*(Lambda12 + 2))) &
                                *wigner(p1, p2, q2, 1:rhomax)
                        ELSE
                           wigner(p1, p2, q2, 1:rhomax) = 0.D0
                        END IF

                        ! term with Lambda1p=Lambda1-1/2
                        IF (p1 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(INT(p1*(lambda1p - p1), 8) &
                                          *(mu1p + p1)*(s2 + 2)*(s2 - 2*Lambda22)) &
                                     /DBLE(Lambda12*(Lambda12 + 1))) &
                             *wigner(p1 - 1, p2, q2, 1:rhomax)

                        ! term with Lambda2p=Lambdda2+1/2
                        IF (q2 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Sq2*(s2 - 2*Lam32prime)*(s2 - 2*Lambda22)) &
                                     /DBLE((Lambda22 + 1)*(Lambda22 + 2))) &
                             *wigner(p1, p2, q2 - 1, 1:rhomax)

                        ! term with Lambda2p=Lambda2-1/2
                        IF (p2 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Rp2*(s2 + 2)*(s2 - 2*Lambda12)) &
                                     /DBLE(Lambda22*(Lambda22 + 1))) &
                             *wigner(p1, p2 - 1, q2, 1:rhomax)

                     ELSE

                        ! term with Lambda1p=Lambda1+1/2
                        IF (q1 /= 0) THEN
                           wigner(p1, p2, q2, 1:rhomax) &
                              = DSQRT(DBLE(INT(q1*(mu1p - q1), 8) &
                                           *(lm1 - q1)*(s2 + 2)*(s2 - 2*Lambda22)) &
                                      /DBLE((Lambda12 + 1)*(Lambda12 + 2))) &
                                *wigner(p1, p2, q2, 1:rhomax)
                        ELSE
                           wigner(p1, p2, q2, 1:rhomax) = 0.D0
                        END IF

                        ! term with Lambda1p=Lambda1-1/2
                        IF (p1 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(INT(p1)*(lambda1p - p1)*(mu1p + p1) &
                                          *(s2 - 2*Lambda12)*(s2 - 2*Lam32prime)) &
                                     /DBLE(Lambda12*(Lambda12 + 1))) &
                             *wigner(p1 - 1, p2, q2, 1:rhomax)

                        ! term with Lambda2p=Lambda2+1/2
                        IF (q2 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             + DSQRT(DBLE(Sq2*(s2 + 2)*(s2 - 2*Lambda12)) &
                                     /DBLE((Lambda22 + 1)*(Lambda22 + 2))) &
                             *wigner(p1, p2, q2 - 1, 1:rhomax)

                        ! term with Lambda2p=Lambda2-1/2
                        IF (p2 /= 0) &
                           wigner(p1, p2, q2, 1:rhomax) &
                           = wigner(p1, p2, q2, 1:rhomax) &
                             - DSQRT(DBLE(Rp2*(s2 - 2*Lam32prime)*(s2 - 2*Lambda22)) &
                                     /DBLE(Lambda22*(Lambda22 + 1))) &
                             *wigner(p1, p2 - 1, q2, 1:rhomax)

                     END IF

                     wigner(p1, p2, q2, 1:rhomax) = wigner(p1, p2, q2, 1:rhomax)*N3

                     q1 = q1 + 1
                     Lambda12 = Lambda12 - 2 ! Lambda12 is 2*Lambda1 in Eq.(...)
                  END DO
                  pq1 = pq1 + 1 ! pq1 is p1+q1
                  Lambda22 = Lambda22 + 1 ! Lambda22 is 2*Lambda2 in Eq.(...)
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE calculate_coupling_canonical_nonextremal

   SUBROUTINE calculate_coupling_canonical(irrep1, irrep2, irrep3, epsilon3, Lambda32, &
                                           dimpq, dimw, rhomax, numb, wigner_block, p1a, p2a, q2a) BIND(C)
      !--------------------------------------------------------------------------------------------------------
      !! Calculate SU(3)-SU(2)xU(1) reduced coupling coefficients
      ! / (lambda1,mu1)    (lambda2,mu2)   ||  (lambda3,mu3)  \
      ! \epsilon1,Lambda1 epsilon2,Lambda2 || epsilon3,Lambda3/rho
      !! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3.
      !! Quantum numbers epsilon3,Lambda3 must be valid, i.e., there must be integers p3,q3 satisfying:
      !!   0<=p3<=lambda3,
      !!   0<=q3<=mu3,
      !!   epsilon3=2*lambda3+mu3-3*(p3+q3),
      !!   Lambda3=(mu3+p3-q3)/2.
      !
      ! Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,dimpq,dimw,rhomax
      ! Output arguments: numb,wigner_block,p1a,p2a,q2a
      !
      ! irrep%lambda is lambda, irrep%mu is mu,
      ! rmomax is multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3),
      ! dimpq is size of arrays p1a,p2a,q2a, which should be at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1),
      ! dimw is size of array wigner_block, which should be at least rhomax*dimpq.
      !
      !! Resulting reduced coupling coefficient is wigner_block(ind), where
      !!   epsilon2=2*lambda2+mu2-3*(p2+q2),
      !!   epsilon1=epsilon3-epsilon2,
      !!   Lambda1=(mu1+p1-q1)/2,
      !!   Lambda2=(mu2+p2-q2)/2,
      !!   p1=p1a(i),
      !!   p2=p2a(i),
      !!   q2=q2a(i),
      !!   q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2,
      !!   ind=i+numb*(rho-1),
      !!   1<=i<=numb.
      !--------------------------------------------------------------------------------------------------------
      !USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN) :: irrep1
         !! First SU(3) irrep in SU(3) coupling, i.e., (lambda1,mu1)
      TYPE(su3irrep), INTENT(IN) :: irrep2
         !! Second SU(3) irrep in SU(3) coupling, i.e., (lambda2,mu2)
      TYPE(su3irrep), INTENT(IN) :: irrep3
         !! Resulting SU(3) irrep, i.e., (lambda3,mu3)
      INTEGER(C_INT), INTENT(IN) :: epsilon3
         !! epsilon3
      INTEGER(C_INT), INTENT(IN) :: Lambda32
         !! 2*Lambda3
      INTEGER(C_INT), INTENT(IN) :: dimpq
         !! Size of arrays p1a,p2a,q2a, which should be at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
      INTEGER(C_INT), INTENT(IN) :: dimw
         !! Size of array wigner_block, which should be at least rhomax*dimpq
      INTEGER(C_INT), INTENT(IN) :: rhomax
         !! Multiplicity of SU(3) coupling
      INTEGER(C_INT), INTENT(OUT) :: numb
         !! Number of resulting reduced coupling coefficients for given outer multiplicity index
      INTEGER(C_INT), DIMENSION(dimpq), INTENT(OUT) :: p1a
         !! Array of values of p1
      INTEGER(C_INT), DIMENSION(dimpq), INTENT(OUT) :: p2a
         !! Array of values of p2
      INTEGER(C_INT), DIMENSION(dimpq), INTENT(OUT) :: q2a
         !! Array of values of q2
      REAL(C_DOUBLE), DIMENSION(dimw), INTENT(OUT) :: wigner_block
         !! Array of resulting reduced coupling coefficients
      INTEGER :: i, rho, ind
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :, :) :: wignerex, wigner

      i = MAX(irrep2%lambda, irrep2%mu)
      ALLOCATE (wignerex(0:MAX(irrep1%lambda, irrep1%mu), 0:i, 0:i, 1:rhomax), &
                wigner(0:MAX(irrep1%lambda, irrep1%mu), 0:i, 0:i, 1:rhomax))

      IF (2*epsilon3 <= irrep3%lambda - irrep3%mu) THEN
         CALL calculate_coupling_canonical_extremal(irrep1, irrep2, irrep3, 1, rhomax, numb, wignerex, p1a, p2a, q2a)
         CALL calculate_coupling_canonical_nonextremal(irrep1, irrep2, irrep3, epsilon3, Lambda32, 1, &
                                                       rhomax, numb, wignerex, wigner, p1a, p2a, q2a)
      ELSE
         CALL calculate_coupling_canonical_extremal(irrep1, irrep2, irrep3, 0, rhomax, numb, wignerex, p1a, p2a, q2a)
         CALL calculate_coupling_canonical_nonextremal(irrep1, irrep2, irrep3, epsilon3, Lambda32, 0, &
                                                       rhomax, numb, wignerex, wigner, p1a, p2a, q2a)
         IF (epsilon3 == 2*irrep3%lambda + irrep3%mu) THEN
            rho = (2*(irrep1%lambda - irrep1%mu - irrep2%mu) - irrep2%lambda - epsilon3)/3
            DO i = 1, numb
               ind = p2a(i)
               p1a(i) = rho + p1a(i) + q2a(i) + ind
               p2a(i) = irrep2%lambda - q2a(i)
               q2a(i) = irrep2%mu - ind
            END DO
         END IF
      END IF

      ind = 0
      DO rho = 1, rhomax
         DO i = 1, numb
            ind = ind + 1 ! ind=i+numb*(rho-1)
            wigner_block(ind) = wigner(p1a(i), p2a(i), q2a(i), rho)
         END DO
      END DO

      DEALLOCATE (wignerex, wigner)

   END SUBROUTINE calculate_coupling_canonical

END MODULE ndsu3lib_coupling_canonical
