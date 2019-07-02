MODULE binomial_coeff_factorials
!-------------------------------------------------------------------
! BINOMIAL COEFFICIENTS AND FACTORIALS ***** SEE COMMENT BELOW *****
!-------------------------------------------------------------------
! UPDATE/MOD: (MTS,06-76)  H.SATO             EXPANDED RANGE        
!             (LSU,11-78)  J.P.DRAAYER        LOG FACTORIALS        
!             (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         
!             (LSU,08-81)  J.P.DRAAYER        EXPANDED RANGE        
!             (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    
!             (LSU,03-83)  J.P.DRAAYER        D,Q MIX & LOG CUTS    
!             (LSU,11-84)  J.P.DRAAYER        MODIFIED FACTORIALS   
!             (LSU,01-88)  J.P.DRAAYER        DLOGF RANGE/INDEX     
!             (LSU,10-89)  J.P.DRAAYER        BINOMIAL INVERSES     
!             (LSU,11-89)  J.P.DRAAYER        POWERS OF TWO ARRAY   
!             (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
!             (LSU,11-11)  T. DYTRYCH         MERGE WITH SU3GENBK
!             (ND,11-22)   A.E.MCCOY          ADDED FLAG SU3QUAD_GNU
!             (ND,2019)    J.HERKO            MODULE INSTEAD OF COMMON BLOCKS, EXPLICIT DECLARATIONS, REWRITTEN IN MODERN FORTRAN
!                                                                       
! DOUBLE PRECISION BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83):      
!       SCALE: BINO(I,J)=DBINO(I*(I+1)/2+J+1)                       
!       RANGE: BINO(0,0)=DBINO(1) TO BINO(128,128)=DBINO(8385)      
!       ADDED: 2**I = DTWOS(I) WHERE I=-128,128                     
! DOUBLE/QUAD SWITCHABLE BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83):      
!       SCALE: BINO(I,J)=QBINO(I*(I+1)/2+J+1)                       
!       RANGE: BINO(0,0)=QBINO(1) TO BINO(192,192)=QBINO(18721)     
!       ADDED: 2**I = QTWOS(I) WHERE I=-192,192                     
! LOG FACTORIALS (FACT) (INSERTED 11-78, MODIFIED 01-88):      
!       SCALE: LNFACT(I)=DLOGF(2*I)                                 
!       RANGE: LNFACT(0)=DLOGF(0) TO LNFACT(1000)=DLOGF(2000)       
!
! This is a new version of BLOCKS
!-------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: II,MM,N,JJ,LL,KK
REAL(KIND=8), DIMENSION(8385) :: DBINO,DBINV
REAL(KIND=8), DIMENSION(-128:128) :: DTWOS
REAL(KIND=8), DIMENSION(0:2000) :: DLOGF
#if defined(SU3DBL)
REAL(KIND=8), DIMENSION(18721) :: QBINO,QBINV
REAL(KIND=8), DIMENSION(-192:192) :: QTWOS
REAL(KIND=8) :: QGEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16), DIMENSION(18721) :: QBINO,QBINV
REAL(KIND=16), DIMENSION(-192:192) :: QTWOS
REAL(KIND=16) :: QGEN
#endif

CONTAINS
 SUBROUTINE calculate_binomial_coefficients_and_factorials
!This subroutine calculates binomial coefficients and factorials and stores them in corresponding arrays.
!This subroutine should be called at the begining of the main program (after "USE binomial_coefficients_and_factorials").
  DLOGF(0)=0.D0

  DO II=2,2000,2
   DLOGF(II)=DLOGF(II-2)+DLOG(DFLOAT(II/2))
  END DO

  DO MM=1,193
   II=MM-1
   DO N=1,MM
    JJ=N-1
    LL=MIN(JJ,II-JJ)
    QGEN=1.0
    IF (LL/=0) THEN
     DO KK=1,LL
#if defined(SU3DBL)
      QGEN=DFLOAT(II+1-KK)/DFLOAT(KK)*QGEN
#elif defined(SU3QUAD)
      QGEN=QFLOAT(II+1-KK)/QFLOAT(KK)*QGEN
#elif defined(SU3QUAD_GNU)
      QGEN=REAL(II+1-KK,16)/REAL(KK,16)*QGEN
#endif
     END DO
    END IF
    QBINO(II*(II+1)/2+JJ+1)=QGEN
   END DO
  END DO

  DO II=1,18721
   QBINV(II)=1.0/QBINO(II)
  END DO

  DO II=1,8385
   DBINO(II)=QBINO(II)
   DBINV(II)=QBINV(II)
  END DO

  QTWOS(0)=1.0
  DO II=1,192
   QTWOS(+II)=QTWOS(II-1)*2.0
   QTWOS(-II)=QTWOS(1-II)/2.0
  END DO

  DO II=-128,128
   DTWOS(II)=QTWOS(II)
  END DO

 END SUBROUTINE calculate_binomial_coefficients_and_factorials
END MODULE binomial_coeff_factorials