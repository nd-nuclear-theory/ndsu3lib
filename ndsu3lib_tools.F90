!***************************************************************************************************************
!
! ndsu3lib_tools.F90
!! Module with procedures for initialization and finalization of ndsu3lib and for outer and inner multiplicities
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!***************************************************************************************************************
MODULE ndsu3lib_tools
   !! Module with procedures for initialization and finalization of ndsu3lib and for outer and inner multiplicities
   USE, INTRINSIC :: ISO_C_BINDING
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
   USE, INTRINSIC :: ISO_FORTRAN_ENV, wp => REAL128
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
   USE mpmodule
#endif
   IMPLICIT NONE
   TYPE, BIND(C) :: su3irrep
      !! Derived type for SU(3) irrep labels
      INTEGER(C_INT) :: lambda, mu
   END TYPE su3irrep
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
   INTEGER :: ndig, nwds, nwdsm
   PARAMETER(ndig=37, nwds=INT(ndig/mpdpw + 2.D0), nwdsm=8) ! 37 is the minimal ndig, Medium precision set to 8 words, i.e., quad
#endif
   !----------------------------------------------------------------------------------------------------------------------------
   ! Binomial coefficients: (n choose k)=binom(n*(n+1)/2+k)=binom_quad(n*(n+1)/2+k)=binom_mp(n*(n+1)/2+k): n=0,...,upbound_binom
   !                                                                                                       k=0,...,n
   ! I(p,q,s)=\sum_{n}(-1)^{n}*(p choose s-n)*(q choose n)=Ia((p^3+(p+q)^2+q)/2+s): p=0,...,upbound_I
   !                                                                                q=0,...,p
   !                                                                                s=0,...,p+q
   ! S(p,q,s)=\sum_{n}(-1)^{n}*(p choose n)/(p+q choose s+n)
   !         =Sa((upbound_S*p*(p+upbound_S+2)+(p+q)*(p+q+1))/2+s): p=0,...,upbound_S
   !                                                               q=0,...,ubbound_S
   !                                                               s=0,...,p+q
   ! For p<q: I(p,q,s)=(-1)^{s}*I(q,p,s)
   !-----------------------------------------------------------------------------------------------------------------------------
   INTEGER :: upbound_binom
      !! Upper bound on n of stored binomial coefficients (n choose k)
   INTEGER :: upbound_I
      !! Upper bound on p of stored values of I(p,q,s)
   INTEGER :: upbound_S
      !! Upper bound on p and q of stored values of S(p,q,s)
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: binom
      !! Array of binomial coefficients in double precision.
      !! (n choose k)=binom(n*(n+1)/2+k)
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Ia
      !! Array of values of I(p,q,s) in double precision.
      !! I(p,q,s)=Ia((p^3+(p+q)^2+q)/2+s) for p>=q.
      !! For p<q: I(p,q,s)=(-1)^{s}*I(q,p,s)
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Sa
      !! Array of values of S(p,q,s) in double precision.
      !! S(p,q,s)=Sa((upbound_S*p*(p+upbound_S+2)+(p+q)*(p+q+1))/2+s)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
   REAL(wp), ALLOCATABLE, DIMENSION(:) :: binom_quad
      !! Array of binomial coefficients in quad precision
   REAL(wp), ALLOCATABLE, DIMENSION(:) :: Ia_quad
      !! Array of values of I(p,q,s) in quad precision
   REAL(wp), ALLOCATABLE, DIMENSION(:) :: Sa_quad
      !! Array of values of S(p,q,s) in quad precision
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
   TYPE(mp_realm), ALLOCATABLE, DIMENSION(:) :: binom_quad
      !! Array of binomial coefficients in medium arbitrary precision
   TYPE(mp_realm), ALLOCATABLE, DIMENSION(:) :: Ia_quad
      !! Array of values of I(p,q,s) in medium arbitrary precision
   TYPE(mp_realm), ALLOCATABLE, DIMENSION(:) :: Sa_quad
      !! Array of values of S(p,q,s) in medium arbitrary precision
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
   TYPE(mp_real), ALLOCATABLE, DIMENSION(:) :: binom_mp
      !! Array of binomial coefficients in full arbitrary precision
   TYPE(mp_real), ALLOCATABLE, DIMENSION(:) :: Ia_mp
      !! Array of values of I(p,q,s) in full arbitrary precision
   TYPE(mp_real), ALLOCATABLE, DIMENSION(:) :: Sa_mp
      !! Array of values of S(p,q,s) in fumm arbitrary precision
#endif

CONTAINS

   SUBROUTINE allocate_binom(upbound)
      !----------------------------------------------------------------------------------------------------------------
      !! Allocate arrays binom, binom_quad, and binom_mp, calculate binomial coefficients, and store them in the arrays
      !----------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: upbound
         !! Upper bound on n of stored binomial coefficients (n choose k)
      INTEGER :: n, k, ind
      upbound_binom = upbound
      ind = upbound_binom*(upbound_binom + 1)/2 + upbound_binom
#if defined(NDSU3LIB_DBL)
      ALLOCATE (binom(0:ind))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      ALLOCATE (binom(0:ind), binom_quad(0:ind), binom_mp(0:ind))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      ALLOCATE (binom(0:ind), binom_quad(0:ind))
#endif
      binom(0) = 1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      binom_quad(0) = 1.0_wp
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      binom_quad(0) = mprealm(1.D0, nwdsm)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      binom_mp(0) = mpreal(1.D0, nwds)
#endif
      ind = 1
      DO n = 1, upbound_binom
         binom(ind) = 1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
         binom_quad(ind) = 1.0_wp
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         binom_quad(ind) = mprealm(1.D0, nwdsm)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         binom_mp(ind) = mpreal(1.D0, nwds)
#endif
         ind = ind + 1
         DO k = 1, n - 1
            binom(ind) = binom(ind - n - 1) + binom(ind - n)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            binom_quad(ind) = binom_quad(ind - n - 1) + binom_quad(ind - n)
            binom_mp(ind) = binom_mp(ind - n - 1) + binom_mp(ind - n)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
            binom_quad(ind) = binom_quad(ind - n - 1) + binom_quad(ind - n)
#endif
            ind = ind + 1
         END DO
         binom(ind) = 1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
         binom_quad(ind) = 1.0_wp
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         binom_quad(ind) = mprealm(1.D0, nwdsm)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         binom_mp(ind) = mpreal(1.D0, nwds)
#endif
         ind = ind + 1
      END DO
   END SUBROUTINE allocate_binom

   SUBROUTINE deallocate_binom
      !---------------------------------------------------------------------------------
      !! Deallocate arrays binom, binom_quad, and binom_mp storing binomial coefficients
      !---------------------------------------------------------------------------------
      IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
      DEALLOCATE (binom)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      DEALLOCATE (binom, binom_quad, binom_mp)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      DEALLOCATE (binom, binom_quad)
#endif
   END SUBROUTINE deallocate_binom

   SUBROUTINE allocate_I(upboundI)
      !----------------------------------------------------------------------------------------------------------------
      !! Allocate arrays Ia, Ia_quad, and Ia_mp, calculate values of I(p,q,s) recurrently, and store them in the arrays
      ! The recursion relation is
      !
      ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
      !
      ! with starting condition
      !
      ! I(p,0,sigma)=(p choose sigma)
      !----------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: upboundI
         !! Upper bound on p of stored values of I(p,q,s)
      INTEGER :: p, q, sigma, ind1, ind2, aux
      upbound_I = upboundI
      IF (upbound_I > upbound_binom) THEN
         CALL reallocate_binom(upbound_I - upbound_binom)
      END IF
      ind1 = (upbound_I**3 + (2*upbound_I)**2 + upbound_I)/2 + 2*upbound_I
#if defined(NDSU3LIB_DBL)
      ALLOCATE (Ia(0:ind1))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      ALLOCATE (Ia(0:ind1), Ia_quad(0:ind1), Ia_mp(0:ind1))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      ALLOCATE (Ia(0:ind1), Ia_quad(0:ind1))
#endif
      aux = 0
      ind2 = 0
      DO p = 0, upbound_I
         aux = aux + p ! aux1=p*(p+1)/2
         ind1 = p*aux ! ind1=(p^3+p^2)/2
         DO sigma = 0, p
            Ia(ind1) = binom(ind2) ! I(p,0,sigma)=(p choose sigma)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Ia_quad(ind1) = binom_quad(ind2)
            Ia_mp(ind1) = binom_mp(ind2)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
            Ia_quad(ind1) = binom_quad(ind2)
#endif
            ind1 = ind1 + 1 ! ind1=(p^3+p^2)/2+sigma
            ind2 = ind2 + 1 ! ind2=p*(p+1)/2+sigma
         END DO
      END DO
      Ia(3) = Ia(0)
      Ia(4) = 0.D0
      Ia(5) = -Ia(0)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      Ia_quad(3) = Ia_quad(0)
      Ia_quad(4) = 0.0_wp
      Ia_quad(5) = -Ia_quad(0)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      Ia_quad(3) = Ia_quad(0)
      Ia_quad(4) = mprealm(0.D0, nwdsm)
      Ia_quad(5) = -Ia_quad(0)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      Ia_mp(3) = Ia_mp(0)
      Ia_mp(4) = mpreal(0.D0, nwds)
      Ia_mp(5) = -Ia_mp(0)
#endif
      ind1 = 5
      ind2 = 1
      DO p = 2, upbound_I
         ind1 = ind1 + p + 1
         DO q = 1, p
            DO sigma = 0, 1
               ind1 = ind1 + 1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
               Ia(ind1) = Ia(ind2) ! I(p,q,sigma)=I(p-1,q-1,sigma)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
               Ia_quad(ind1) = Ia_quad(ind2)
               Ia_mp(ind1) = Ia_mp(ind2)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
               Ia_quad(ind1) = Ia_quad(ind2)
#endif
               ind2 = ind2 + 1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
            END DO
            DO sigma = 2, p + q - 2
               ind1 = ind1 + 1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
               Ia(ind1) = Ia(ind2) - Ia(ind2 - 2) ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
               Ia_quad(ind1) = Ia_quad(ind2) - Ia_quad(ind2 - 2)
               Ia_mp(ind1) = Ia_mp(ind2) - Ia_mp(ind2 - 2)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
               Ia_quad(ind1) = Ia_quad(ind2) - Ia_quad(ind2 - 2)
#endif
               ind2 = ind2 + 1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
            END DO
            ind2 = ind2 - 2
            DO sigma = p + q - 1, p + q
               ind1 = ind1 + 1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
               Ia(ind1) = -Ia(ind2) ! I(p,q,sigma)=-I(p-1,q-1,sigma-2)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
               Ia_quad(ind1) = -Ia_quad(ind2)
               Ia_mp(ind1) = -Ia_mp(ind2)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
               Ia_quad(ind1) = -Ia_quad(ind2)
#endif
               ind2 = ind2 + 1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma-2
            END DO
         END DO
      END DO
   END SUBROUTINE allocate_I

   SUBROUTINE allocate_S(upboundS)
      !----------------------------------------------------------------------------------------------------------------
      !! Allocate arrays Sa, Sa_quad, and Sa_mp, calculate values of S(p,q,s) recurrently, and store them in the arrays
      ! The recursion relations are
      !
      ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1) unless q=upbound_S
      ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q) iff q=upbound_S
      !
      ! with starting condition
      !
      ! S(0,q,sigma)=1/(q choose sigma)
      !----------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: upboundS
         !! Upper bound on p and q of stored values of S(p,q,s)
      INTEGER :: p, q, sigma, ind1, ind2, aux, ind3, n, ind4, ind5, nmax, aux2, aux3
      upbound_S = upboundS
      IF (2*upbound_S > upbound_binom) THEN
         CALL reallocate_binom(2*upbound_S - upbound_binom)
      END IF
      n = (upbound_S*upbound_S*(2*upbound_S + 2) + 2*upbound_S*(2*upbound_S + 1))/2 + 2*upbound_S
#if defined(NDSU3LIB_DBL)
      ALLOCATE (Sa(0:n))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      ALLOCATE (Sa(0:n), Sa_quad(0:n), Sa_mp(0:n))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      ALLOCATE (Sa(0:n), Sa_quad(0:n))
#endif
      ind1 = 0
      DO q = 0, upbound_S
         DO sigma = 0, q
            Sa(ind1) = 1.D0/binom(ind1)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
            Sa_quad(ind1) = 1.0_wp/binom_quad(ind1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_quad(ind1) = 1.D0/binom_quad(ind1)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_mp(ind1) = 1.D0/binom_mp(ind1)
#endif
            ind1 = ind1 + 1 ! ind1=q*(q+1)/2+sigma
         END DO
      END DO
      aux = 0
      ind1 = upbound_S*(upbound_S + 3)/2 + 1
      aux2 = 1
      DO p = 1, upbound_S
         aux2 = aux2 + p ! aux2=p*(p+1)/2+1
         aux3 = (p + upbound_S)*(p + upbound_S + 1)/2 + 1
         ind2 = aux + p
         aux = ind1
         DO q = 0, upbound_S - 1
            DO sigma = 0, p + q - 1
               Sa(ind1) = Sa(ind2) - Sa(ind2 + 1) ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
               Sa_quad(ind1) = Sa_quad(ind2) - Sa_quad(ind2 + 1)
               Sa_mp(ind1) = Sa_mp(ind2) - Sa_mp(ind2 + 1)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
               Sa_quad(ind1) = Sa_quad(ind2) - Sa_quad(ind2 + 1)
#endif
               ind1 = ind1 + 1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
               ind2 = ind2 + 1 ! ind2=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q+1)*(p-1+q+1+1))/2+sigma
            END DO
            Sa(ind1) = 1.D0 ! S(p,q,p+q)=1
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
            Sa_quad(ind1) = 1.0_wp
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_quad(ind1) = mprealm(1.D0, nwdsm)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_mp(ind1) = mpreal(1.D0, nwds)
#endif
            ind1 = ind1 + 1
            ind2 = ind2 + 1
         END DO
         ! q=upbound
         ind3 = ind2 - p - upbound_S
         ind2 = ind1 - p - upbound_S
         DO sigma = 0, upbound_S - 1
            Sa(ind1) = (DBLE(upbound_S - sigma)*Sa(ind2) + DBLE(p)*Sa(ind3))/DBLE(p + upbound_S)
#if defined(NDSU3LIB_QUAD)
            Sa_quad(ind1) = (QFLOAT(upbound_S - sigma)*Sa_quad(ind2) + QFLOAT(p)*Sa_quad(ind3))/QFLOAT(p + upbound_S)
#elif defined(NDSU3LIB_QUAD_GNU)
            Sa_quad(ind1) = (REAL(upbound_S - sigma, 16)*Sa_quad(ind2) + REAL(p, 16)*Sa_quad(ind3))/REAL(p + upbound_S, 16)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_quad(ind1) = (DBLE(upbound_S - sigma)*Sa_quad(ind2) + DBLE(p)*Sa_quad(ind3))/DBLE(p + upbound_S)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_mp(ind1) = (DBLE(upbound_S - sigma)*Sa_mp(ind2) + DBLE(p)*Sa_mp(ind3))/DBLE(p + upbound_S)
#endif
            ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q)
            ind1 = ind1 + 1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
            ind2 = ind2 + 1 ! ind2=(upbound*p*(p+upbound+2)+(p+q-1)*(p+q-1+1))/2+sigma
            ind3 = ind3 + 1 ! ind3=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q)*(p-1+q+1))/2+sigma
         END DO
         nmax = p - 1
         DO sigma = upbound_S, p + upbound_S - 1
            nmax = nmax - 1
            ind4 = aux2
            ind5 = aux3 + sigma
            Sa(ind1) = 1.D0/binom(ind5 - 1)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
            Sa_quad(ind1) = 1.0_wp/binom_quad(ind5 - 1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_quad(ind1) = 1.D0/binom_quad(ind5 - 1)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
            Sa_mp(ind1) = 1.D0/binom_mp(ind5 - 1)
#endif
            DO n = 1, nmax, 2
               Sa(ind1) = Sa(ind1) - binom(ind4)/binom(ind5)
               Sa(ind1) = Sa(ind1) + binom(ind4 + 1)/binom(ind5 + 1)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
               Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
               Sa_quad(ind1) = Sa_quad(ind1) + binom_quad(ind4 + 1)/binom_quad(ind5 + 1)
               Sa_mp(ind1) = Sa_mp(ind1) - binom_mp(ind4)/binom_mp(ind5)
               Sa_mp(ind1) = Sa_mp(ind1) + binom_mp(ind4 + 1)/binom_mp(ind5 + 1)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
               Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
               Sa_quad(ind1) = Sa_quad(ind1) + binom_quad(ind4 + 1)/binom_quad(ind5 + 1)
#endif
               ind4 = ind4 + 2
               ind5 = ind5 + 2
            END DO
            IF (nmax >= -1) THEN
               IF (BTEST(nmax, 0)) THEN ! nmax is odd
                  Sa(ind1) = Sa(ind1) - binom(ind4)/binom(ind5)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
                  Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
                  Sa_mp(ind1) = Sa_mp(ind1) - binom_mp(ind4)/binom_mp(ind5)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
                  Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
#endif
               ELSE ! nmax is even
                  Sa(ind1) = Sa(ind1) - binom(ind4)/binom(ind5)
                  Sa(ind1) = Sa(ind1) + binom(ind4 + 1)/binom(ind5 + 1)
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
                  Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
                  Sa_quad(ind1) = Sa_quad(ind1) + binom_quad(ind4 + 1)/binom_quad(ind5 + 1)
                  Sa_mp(ind1) = Sa_mp(ind1) - binom_mp(ind4)/binom_mp(ind5)
                  Sa_mp(ind1) = Sa_mp(ind1) + binom_mp(ind4 + 1)/binom_mp(ind5 + 1)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
                  Sa_quad(ind1) = Sa_quad(ind1) - binom_quad(ind4)/binom_quad(ind5)
                  Sa_quad(ind1) = Sa_quad(ind1) + binom_quad(ind4 + 1)/binom_quad(ind5 + 1)
#endif
               END IF
            END IF
            ind1 = ind1 + 1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
         END DO
         Sa(ind1) = 1.D0 ! S(p,q,p+q)=1
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
         Sa_quad(ind1) = 1.0_wp
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         Sa_quad(ind1) = mprealm(1.D0, nwdsm)
#endif
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
         Sa_mp(ind1) = mpreal(1.D0, nwds)
#endif
         ind1 = ind1 + 1
      END DO
   END SUBROUTINE allocate_S

   SUBROUTINE deallocate_I
      !---------------------------------------------------------------------
      !! Deallocate arrays Ia, Ia_quad, and Ia_mp storing values of I(p,q,s)
      !---------------------------------------------------------------------
      IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
      DEALLOCATE (Ia)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      DEALLOCATE (Ia, Ia_quad, Ia_mp)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      DEALLOCATE (Ia, Ia_quad)
#endif
   END SUBROUTINE deallocate_I

   SUBROUTINE deallocate_S
      !---------------------------------------------------------------------
      !! Deallocate arrays Sa, Sa_quad, and Sa_mp storing values of S(p,q,s)
      !---------------------------------------------------------------------
      IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
      DEALLOCATE (Sa)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      DEALLOCATE (Sa, Sa_quad, Sa_mp)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      DEALLOCATE (Sa, Sa_quad)
#endif
   END SUBROUTINE deallocate_S

   SUBROUTINE initialize_ndsu3lib(wso3, openmp, lmpmu) BIND(C)
      !--------------------------------------------------------------------------------------
      !! ndsu3lib initialization subroutine
      !! Must be called by the main program before calling ndsu3lib subroutines for SU(3) coupling or recoupling coefficients.
      !! In OpenMP parallelized programs it must be called by each thread.
      !
      ! Input arguments: wso3,openmp,lmpmu
      !
      ! wso3 must be .TRUE. if SU(3)-SO(3) coupling coefficients are going to be calculated.
      ! openmp must be .TRUE. if OpenMP parallelization is used, otherwise it must be .FALSE.
      ! lmpmu should be greater than or equal to the maximal expected value of lambda+mu.
      !--------------------------------------------------------------------------------------
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
      USE fwigxjpf
#endif
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN), VALUE :: wso3
         !! Must be .TRUE. if SU(3)-SO(3) coupling coefficients are going to be calculated
      LOGICAL(C_BOOL), INTENT(IN), VALUE :: openmp
         !! Must be .TRUE. if OpenMP parallelization is used, otherwise it must be .FALSE.
      INTEGER(C_INT), INTENT(IN), VALUE :: lmpmu
         !! Should be greater than or equal to the maximal expected value of lambda+mu
!$OMP SINGLE
      CALL allocate_binom(4*lmpmu + 1)
!$OMP END SINGLE
      IF (wso3) THEN
!$OMP SINGLE
         CALL allocate_I(3*lmpmu)
         CALL allocate_S(2*lmpmu)
!$OMP END SINGLE
#if defined(NDSU3LIB_RACAH_WIGXJPF)
!$OMP SINGLE
         CALL fwig_table_init(2*lmpmu, 6)
!$OMP END SINGLE
         IF (openmp) THEN
            CALL fwig_thread_temp_init(2*lmpmu)
         ELSE
            CALL fwig_temp_init(2*lmpmu)
         END IF
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
!$OMP SINGLE
         CALL fwig_table_init(2*lmpmu, 3)
!$OMP END SINGLE
         IF (openmp) THEN
            CALL fwig_thread_temp_init(2*lmpmu)
         ELSE
            CALL fwig_temp_init(2*lmpmu)
         END IF
#endif
      ELSE
#if defined(NDSU3LIB_RACAH_WIGXJPF)
!$OMP SINGLE
         CALL fwig_table_init(lmpmu, 6)
!$OMP END SINGLE
         IF (openmp) THEN
            CALL fwig_thread_temp_init(lmpmu)
         ELSE
            CALL fwig_temp_init(lmpmu)
         END IF
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
!$OMP SINGLE
         CALL fwig_table_init(lmpmu, 3)
!$OMP END SINGLE
         IF (openmp) THEN
            CALL fwig_thread_temp_init(lmpmu)
         ELSE
            CALL fwig_temp_init(lmpmu)
         END IF
#endif
      END IF
   END SUBROUTINE initialize_ndsu3lib

   SUBROUTINE finalize_ndsu3lib(wso3) BIND(C)
      !-------------------------------------------------------------------------------------------------------------------------------------
      !! Can be called by main program once SU(3) coupling or recoupling coefficients are not going to be calculated anymore to free memory.
      !! In OpenMP parallelized programs it should be called by each thread.
      !
      ! Input argument: wso3
      !
      ! wso3 should be .TRUE. if initialize_ndsu3lib or initialize_ndsu3lib_thread was
      ! called with the first argument being .TRUE.
      !-------------------------------------------------------------------------------------------------------------------------------------
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
      USE fwigxjpf
#endif
      IMPLICIT NONE
      LOGICAL(C_BOOL), INTENT(IN), VALUE :: wso3
         !! Should be .TRUE. if initialize_ndsu3lib or initialize_ndsu3lib_thread was called with first argument being .TRUE.
      IF (wso3) THEN
!$OMP SINGLE
         CALL deallocate_I
         CALL deallocate_S
!$OMP END SINGLE
      END IF
!$OMP SINGLE
      CALL deallocate_binom
!$OMP END SINGLE
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
      CALL fwig_temp_free(); 
!$OMP SINGLE
      CALL fwig_table_free(); 
!$OMP END SINGLE
#endif
   END SUBROUTINE finalize_ndsu3lib

   SUBROUTINE reallocate_binom(incr)
      !-------------------------------------------------------------------------------------------------------------
      !! Print error message and terminate execution due to insufficient upper bound on stored binomial coefficients
      ! Originally: Deallocate arrays binom, binom_quad and binom_mp, allocate them again
      !             with upbound_binom increased by incr, and calculate entries
      !             This functionality was abandoned for tread safety.
      !-------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: incr
         !! Not used
      LOGICAL(KIND=1) :: wso3 = .TRUE.
      PRINT *, "ndsu3lib ERROR: Insufficient upper bound on lambda+mu set in initialize_ndsu3lib"
      CALL finalize_ndsu3lib(wso3)
      STOP
      CALL deallocate_binom
      CALL allocate_binom(upbound_binom + incr)
   END SUBROUTINE reallocate_binom

   SUBROUTINE reallocate_I(incr)
      !------------------------------------------------------------------------------------------------------------------------
      !! Print error message and terminate execution due to insufficient upper bound on stored values of I(p,q,s)
      ! Originally: Deallocate arrays Ia, Ia_quad and Ia_mp, reallocate them with size increased by incr, and calculate entries
      !             This functionality was abandoned for tread safety.
      !------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: incr
         !! Not used
      LOGICAL(KIND=1) :: wso3 = .TRUE.
      PRINT *, "ndsu3lib ERROR: Insufficient upper bound on lambda+mu set in initialize_ndsu3lib"
      CALL finalize_ndsu3lib(wso3)
      STOP
      CALL deallocate_I
      CALL allocate_I(upbound_I + incr)
   END SUBROUTINE reallocate_I

   SUBROUTINE reallocate_S(incr)
      !------------------------------------------------------------------------------------------------------------------------
      !! Print error message and terminate execution due to insufficient upper bound on stored values of S(p,q,s)
      ! Originally: Deallocate arrays Sa, Sa_quad and Sa_mp, reallocate them with size increased by incr, and calculate entries
      !             This functionality was abandoned for tread safety.
      !------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: incr
         !! Not used
      LOGICAL(KIND=1) :: wso3 = .TRUE.
      PRINT *, "ndsu3lib ERROR: Insufficient upper bound on lambda+mu set in initialize_ndsu3lib"
      CALL finalize_ndsu3lib(wso3)
      STOP
      CALL deallocate_S
      CALL allocate_S(upbound_S + incr)
   END SUBROUTINE reallocate_S

   FUNCTION outer_multiplicity(irrep1, irrep2, irrep3) RESULT(rhomax) BIND(C)
      !----------------------------------------------------------------------------------
      !! Calculate multiplicity of SU(3) coupling
      ! (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3),
      ! where lambda1=irrep1%lambda, mu1=irrep1%mu,
      !       lambda2=irrep2%lambda, mu2=irrep2%mu,
      !       lambda3=irrep3%lambda, mu3=irrep3%mu
      ! Reference: M.F.O'Reilly, J.Math.Phys. 23 (1982) 2022: Section 5, Proposition 7(a)
      !----------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep1
         !! First SU(3) irrep in coupling
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep2
         !! Second SU(3) irrep in coupling
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep3
         !! Resulting SU(3) irrep
      INTEGER(C_INT) :: rhomax
         !! Resulting multiplicity
      INTEGER :: C3, C, D
      C3 = irrep1%lambda - irrep1%mu + irrep2%lambda - irrep2%mu - irrep3%lambda + irrep3%mu ! C3 equals 3 times the C in the reference
      C = C3/3 ! C equals the C in the reference if 3 divides C3
      IF ((3*C == C3) .AND. (C >= -MIN(irrep2%mu, irrep1%mu)) .AND. (C <= MIN(irrep1%lambda, irrep2%lambda))) THEN
         D = C + irrep1%mu + irrep2%mu - irrep3%mu ! D equals the D in the reference
         IF ((D >= 0) .AND. (D <= MIN(irrep1%mu + irrep2%mu, irrep2%lambda + irrep2%mu, irrep1%lambda + irrep1%mu)) &
             .AND. (D + C >= 0) &
             .AND. (D + C <= MIN(irrep1%lambda + irrep1%mu, irrep2%lambda + irrep1%lambda, irrep2%lambda + irrep2%mu))) THEN
            rhomax = 1 + MIN(irrep2%mu, irrep1%lambda + irrep1%mu, D, irrep1%lambda - C) &
                     - MAX(0, D - irrep1%mu, D - irrep2%lambda, -C, D - C - irrep1%mu, D + C - irrep2%lambda)
         ELSE
            rhomax = 0
         END IF
      ELSE
         rhomax = 0
      END IF
      RETURN
   END FUNCTION outer_multiplicity

   FUNCTION inner_multiplicity(irrep, L) RESULT(kappamax) BIND(C)
      !--------------------------------------------------------
      !! Inner multiplicity of SO(3) irrep L within SU(3) irrep
      ! (lambda,mu), where lambda=irrep%lambda, mu=irrep%mu
      !--------------------------------------------------------
      IMPLICIT NONE
      TYPE(su3irrep), INTENT(IN), VALUE :: irrep
         !! SU(3) irrep
      INTEGER(C_INT), INTENT(IN), VALUE :: L
         !! SO(3) irrep label (angular momentum)
      INTEGER(C_INT) :: kappamax
         !! Inner multiplicity
      kappamax = MAX(0, (irrep%lambda + irrep%mu + 2 - L)/2) - MAX(0, (irrep%lambda + 1 - L)/2) - MAX(0, (irrep%mu + 1 - L)/2)
   END FUNCTION inner_multiplicity

END MODULE ndsu3lib_tools
