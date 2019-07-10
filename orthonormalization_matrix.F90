SUBROUTINE orthonormalization_matrix(I,J,su3irrepx,KAPMAX,L,DONMAT)                   
!-----------------------------------------------------------------------
!     ORTHONORMALIZATION COEFFICIENTS FOR THE U3R3 TRANSFORMATION       
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,05-83)  J.P.DRAAYER        EXTENDED PRECISION    
!                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
!                 (ND,11-22)   A.E.MCCOY          ADDED FLAG SU3QUAD_GNU
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTION "KSTART" REPLACED BY INTERNAL SUBPROGRAM,
!                                                 DERIVED DATA TYPE su3irrep INSTEAD OF LAMBDA AND MU,
!                                                 ALLOCATABLE ARRAY
!                                                                       
!     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
!                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
!                 J.D.VERGADOS, NUCL.PHYS.A111(1968)681-754             
!     PARAMETERS--(I,J) : (1,1)=GHW, (1,0)=GHW', (0,0)=GLW, (0,1)=GLW'  
!     DIMENSIONS--DONMAT(KAPMAX*KAPMAX)
!
! This is a new version of CONMAT
!
! INPUT ARGUMENTS: I,J,su3irrepx,KAPMAX,L
! OUTPUT ARGUMENT: DONMAT
!-----------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
REAL(KIND=8), EXTERNAL :: transformation_coeff
#if defined(SU3DBL)
REAL(KIND=8), EXTERNAL :: transformation_coeff_quad
REAL(KIND=8), ALLOCATABLE,DIMENSION(:) :: QONMAT
REAL(KIND=8) :: QHOLD
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16), EXTERNAL :: transformation_coeff_quad
REAL(KIND=16), ALLOCATABLE,DIMENSION(:) :: QONMAT
REAL(KIND=16) :: QHOLD
#endif
REAL(KIND=8), DIMENSION(1) :: DONMAT
INTEGER :: I,J,KAPMAX,L,N,KM,IEE,JET,MET,K,KKAQ,KKA,KKAKKA,M,KKBQ,KKB,KB,KKC,KKCKKA,&
           KKCKKB,KKBKKA,KKBKKB,IHELP,KA,KKAKKB
REAL(KIND=8) :: DHOLD
TYPE(su3irrep) :: su3irrepx

IHELP=1
ALLOCATE(QONMAT(KAPMAX*KAPMAX))

IF(I==1)THEN
 KM=KSTART(su3irrepx%lambda,su3irrepx%mu,L)-2
 IEE=-(su3irrepx%lambda+2*su3irrepx%mu)
 JET=su3irrepx%lambda
 MET=su3irrepx%lambda
ELSE
 KM=KSTART(su3irrepx%mu,su3irrepx%lambda,L)-2
 IEE=2*su3irrepx%lambda+su3irrepx%mu
 JET=su3irrepx%mu
 MET=-su3irrepx%mu
END IF
IF(I/=J)MET=-MET
K=KM
KKAQ=-KAPMAX ! N used to be here instead of KAPMAX
DO KKA=1,KAPMAX
 KKAQ=KKAQ+KAPMAX ! N used to be here instead of KAPMAX
 K=K+2
 KKAKKA=KKA+KKAQ
 M=KM
 KKBQ=-KAPMAX ! N used to be here instead of KAPMAX
 DO KKB=1,KKA
  KKBQ=KKBQ+KAPMAX ! N used to be here instead of KAPMAX
  M=M+2
  DHOLD=0.D0
  IF(KKB/=1)THEN
   KB=KKB-1
   DO KKC=1,KB
    KKCKKA=KKC+KKAQ
    KKCKKB=KKC+KKBQ
    DHOLD=DHOLD+DONMAT(KKCKKB)*DONMAT(KKCKKA)
   END DO
  END IF
  DHOLD=transformation_coeff(I,J,su3irrepx,IEE,JET,MET,K,L,M)-DHOLD
  IF(KKB==KKA)EXIT
  KKBKKA=KKB+KKAQ
  KKBKKB=KKB+KKBQ
  DONMAT(KKBKKA)=DONMAT(KKBKKB)*DHOLD
 END DO
 IF(DHOLD<=0.D0)THEN
  IHELP=0 ! IHELP is an indicator of how the end of the loop was reahed (if via EXIT then IHELP=0 else IHELP=1)
  EXIT
 END IF
 DONMAT(KKAKKA)=DSQRT(1.D0/DHOLD)
END DO
IF(IHELP/=0)THEN
 KKAQ=-KAPMAX ! N used to be here instead of KAPMAX
 DO KKA=1,KAPMAX
  KKAQ=KKAQ+KAPMAX ! N used to be here instead of KAPMAX
  IF(KKA==1)CYCLE
  KKAKKA=KKA+KKAQ
  KA=KKA-1
  KKBQ=-KAPMAX ! N used to be here instead of KAPMAX
  DO KKB=1,KA
   KKBQ=KKBQ+KAPMAX ! N used to be here instead of KAPMAX
   DHOLD=0.D0
   DO KKC=KKB,KA
    KKCKKA=KKC+KKAQ
    KKCKKB=KKC+KKBQ
    DHOLD=DHOLD-DONMAT(KKCKKB)*DONMAT(KKCKKA)
   END DO
   KKAKKB=KKA+KKBQ
   KKBKKA=KKB+KKAQ
   DONMAT(KKAKKB)=DONMAT(KKAKKA)*DHOLD
   DONMAT(KKBKKA)=-DONMAT(KKAKKA)*DONMAT(KKBKKA)
  END DO
 END DO
ELSE
 K=KM
 KKAQ=-KAPMAX ! N used to be here instead of KAPMAX
 DO KKA=1,KAPMAX
  KKAQ=KKAQ+KAPMAX ! N used to be here instead of KAPMAX
  K=K+2
  KKAKKA=KKA+KKAQ
  M=KM
  KKBQ=-KAPMAX ! N used to be here instead of KAPMAX
  DO KKB=1,KKA
   KKBQ=KKBQ+KAPMAX ! N used to be here instead of KAPMAX
   M=M+2
   QHOLD=0.D0
   IF(KKB/=1)THEN
    KB=KKB-1
    DO KKC=1,KB
     KKCKKA=KKC+KKAQ
     KKCKKB=KKC+KKBQ
     QHOLD=QHOLD+QONMAT(KKCKKB)*QONMAT(KKCKKA)
    END DO
   END IF
   QHOLD=transformation_coeff_quad(I,J,su3irrepx,IEE,JET,MET,K,L,M)-QHOLD
   IF(KKB==KKA)EXIT
   KKBKKA=KKB+KKAQ
   KKBKKB=KKB+KKBQ
   QONMAT(KKBKKA)=QONMAT(KKBKKB)*QHOLD
  END DO
#if defined(SU3DBL)
  QONMAT(KKAKKA)=DSQRT(1.D0/QHOLD)                                  
#elif defined(SU3QUAD)
  QONMAT(KKAKKA)=QSQRT(1.D0/QHOLD)                                  
#elif defined(SU3QUAD_GNU)
  QONMAT(KKAKKA)=SQRT(1.D0/QHOLD)                                  
#endif
 END DO
 KKAQ=-KAPMAX ! N used to be here instead of KAPMAX
 DO KKA=1,KAPMAX
  KKAQ=KKAQ+KAPMAX ! N used to be here instead of KAPMAX
  IF(KKA==1)CYCLE
  KKAKKA=KKA+KKAQ
  KA=KKA-1
  KKBQ=-KAPMAX ! N used to be here instead of KAPMAX
  DO KKB=1,KA
   KKBQ=KKBQ+KAPMAX ! N used to be here instead of KAPMAX
   QHOLD=0.D0
   DO KKC=KKB,KA
    KKCKKA=KKC+KKAQ
    KKCKKB=KKC+KKBQ
    QHOLD=QHOLD-QONMAT(KKCKKB)*QONMAT(KKCKKA)
   END DO
   KKAKKB=KKA+KKBQ
   KKBKKA=KKB+KKAQ
   QONMAT(KKAKKB)=QONMAT(KKAKKA)*QHOLD
   QONMAT(KKBKKA)=-QONMAT(KKAKKA)*QONMAT(KKBKKA)
   DONMAT(KKAKKB)=QONMAT(KKAKKB)
   DONMAT(KKBKKA)=QONMAT(KKBKKA)
  END DO
 END DO
END IF
DEALLOCATE(QONMAT)
RETURN
	  
CONTAINS
 FUNCTION KSTART(LAM,MU,L) RESULT(IKSTART)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: LAM,MU,L
  INTEGER :: IKSTART
  IKSTART=MOD(LAM,2)+2*(MAX(0,(L-MU)/2)+MOD(LAM+1,2)*MOD(ABS(L-MU),2))
 END FUNCTION KSTART
END SUBROUTINE orthonormalization_matrix