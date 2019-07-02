FUNCTION transformation_coeff_quad(I,J,su3irrepx,IEBX,JBT,MBTX,K,L,M) RESULT(QTU3R3)
!-----------------------------------------------------------------------
!     EXTREMAL U3R3 TRANSFORMATION--CHECKS OUT
!
! Calculates transformation coefficient for SU(3)-SO(3) reduction <G|(G_E)KLM>,
! where |G>=|(lambda,mu)epsilon,Lambda,M_Lambda>.
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER        EXPANDED RANGE        
!                 (LSU,05-80)  J.P.DRAAYER        INTEGER->REAL         
!                 (LSU,06-81)  J.P.DRAAYER        4->4.D0 IN DTS        
!                 (LSU,08-81)  J.P.DRAAYER        LOG BINOMIALS         
!                 (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    
!                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  
!                 (LSU,10-89)  J.P.DRAAYER        VS FORTRAN UPGRADE    
!                 (ND,11-11)   M.A.CAPRIO         DOUBLE/QUAD SWITCHABLE
!                 (ND,11-22)   A.E.MCCOY          ADDED FLAG SU3QUAD_GNU
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, MODULE INSTEAD OF COMMON BLOCKS, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTIONS "ID" AND "ND" REPLACED BY INTERNAL SUBPROGRAMS,
!                                                 DERIVED DATA TYPE su3irrep INSTEAD OF LAMBDA AND MU
!                                                                       
!     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
!                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
!     PARAMETERS--(I,J) : (1,1)=GHW, (1,0)=GHW', (0,0)=GLW, (0,1)=GLW'
!
! This is a new version of QTU3R3
! INPUT ARGUMENTS: (I,J)=(0,1) for lambda>=mu, then E=LW'
!                  (I,J)=(1,0) for lambda<mu, then E=HW'
!                  su3irrepx%lambda=lambda
!                  su3irrepx%mu=mu
!                  IEBX=epsilon
!                  JBT=2*Lambda
!                  MBTX=2*M_Lambda
!-----------------------------------------------------------------------
USE derived_types
USE binomial_coeff_factorials
IMPLICIT NONE
#if defined(SU3DBL)
REAL(KIND=8) :: QTU3R3,QTS,QAS,QAX,QBS,QS,QX,QCS,QCX
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16) :: QTU3R3,QTS,QAS,QAX,QBS,QS,QX,QCS,QCX
#endif
INTEGER, INTENT(IN) :: I,J,IEBX,JBT,MBTX,K,L,M
TYPE(su3irrep), INTENT(IN) :: su3irrepx
INTEGER :: LAM,MU,IEB,MBT,IPH,KBAR,MBAR,IAH,IBH,IPB,IQB,IRB,IW,IX,IY,IS,NTB,MA2,MA1,NA1,NA2,&
           IB1,MB2,MB1,NB1,NB2,IB2,MC2,MC1,NC1,NC2,NC3,NC4,IZMAX,IZMIN,IXMAX,IXMIN,MS2,MS1,&
		   NS1,NS2,IYMAX,IYMIN
!
!     SET HIGHEST WEIGHT/LOWEST WEIGHT
!
IF(I==0)THEN
 LAM=su3irrepx%mu
 MU=su3irrepx%lambda
 IEB=-IEBX
 MBT=-MBTX
ELSE
 LAM=su3irrepx%lambda
 MU=su3irrepx%mu
 IEB=IEBX
 MBT=MBTX
END IF
!
!     OVERALL PHASE
!
IF(I==J)THEN
 IPH=0
ELSE
 IPH=(LAM+K)/2
END IF
!
!     INITIALIZATIONS
!
QTU3R3=0.D0
KBAR=ABS(K)
MBAR=ABS(M)
IAH=(LAM+KBAR)/2
IBH=(JBT+MBAR)/2
IPB=(2*(LAM-MU)+3*JBT-IEB)/3
QTS=QTWOS(-IPB)
IPB=IPB/2
IQB=MU+IPB-JBT
IRB=(JBT+MBT)/2
IF(K<0)IPH=IPH+IPB+IQB+L+M
IF(M<0)IPH=IPH+IQB+IRB+L+K
IW=L+KBAR
IX=L-KBAR
IY=L-MBAR
IS=IW+IX
!
!     OVERALL FACTOR
!
QTS=QTS*DFLOAT(IS+1)*&
#if defined(SU3DBL)
    DSQRT(QBINO(ND(LAM,IPB))*QBINO(ND(MU,IQB))*&
#elif defined(SU3QUAD)
    QSQRT(QBINO(ND(LAM,IPB))*QBINO(ND(MU,IQB))*&
#elif defined(SU3QUAD_GNU)
    SQRT(QBINO(ND(LAM,IPB))*QBINO(ND(MU,IQB))*&
#endif
    QBINO(ND(LAM+MU+1,IQB))*QBINV(ND(IPB+MU+1,IQB))*&
    QBINV(ND(JBT,IRB))*QBINO(ND(IS,IX))*QBINV(ND(IS,IY)))
IF(BTEST(L-IPB+IPH,0))QTS=-QTS
!
!     SETUP CONSTANTS
!
NTB=ID(IPB)
MA2=LAM
MA1=0
NA1=1
NA2=ND(MA2,IAH)
IB1=IRB-MBT
MB2=IPB
MB1=MU-IQB
NB1=ID(MB1)
NB2=NTB+IB1
IB1=IB1+1
IB2=IPB-IB1-1
MC2=LAM+MU+L
MC1=IPB+IQB
NC1=ID(MC1)
NC2=ND(MC2,(IX+IY-L+MC2-MC1)/2)
NC3=ID(IX)
NC4=ND(IW,IY)
!
!     OUTER LOOP
!
IZMAX=IY
IZMIN=IZMAX-IW
IF(IX<IZMAX)IZMAX=IX
IF(IZMIN<0)IZMIN=0
DO IW=0,IPB
 IF(((K/=0).OR.(.NOT.BTEST(IW,0))).AND.((M/=0).OR.(.NOT.BTEST(IB1-MB2+1,0))))THEN
!
!        INNER LOOP "A"
!
  QAS=0.D0
  IXMAX=IAH
  IXMIN=IXMAX-MA2
  IF(MA1<IXMAX)IXMAX=MA1
  IF(IXMIN<0)IXMIN=0
  DO IX=IXMIN,IXMAX
   QAX=QBINO(NA1+IX)*QBINO(NA2-IX)
   IF(BTEST(IX,0))THEN
    QAS=QAS-QAX
   ELSE
    QAS=QAS+QAX
   END IF
  END DO
!
!        INNER LOOP "B"
!
  QBS=0.D0
  IXMAX=IB1-1
  IXMIN=IXMAX-MB2
  IF(MB1<IXMAX)IXMAX=MB1
  IF(IXMIN<0)IXMIN=0
  MS2=2*(IXMIN+1)+IB2-IW
  MS1=JBT-MS2
  NS1=ID(MS1)
  NS2=ND(MS2,IBH)
  DO IX=IXMIN,IXMAX
   IYMAX=IBH
   IYMIN=IYMAX-MS2
   IF(MS1<IYMAX)IYMAX=MS1
   IF(IYMIN<0)IYMIN=0
   QS=0.D0
   DO IY=IYMIN,IYMAX
    QX=QBINO(NS1+IY)*QBINO(NS2-IY)
    IF(BTEST(IY,0))THEN
     QS=QS-QX
    ELSE
     QS=QS+QX
    END IF
   END DO
   IF(IX/=IXMAX)THEN
    NS1=NS1-2*MS1+1
    MS1=MS1-2
    MS2=MS2+2
    NS2=NS2+2*MS2-1
   END IF
   QBS=QBS+QBINO(NB1+IX)*QBINO(NB2-IX)*QS
  END DO
!
!        "C" LOOP
!
  QCS=0.D0
  IYMAX=MC1
  DO IX=IZMIN,IZMAX
   QS=0.D0
   DO IY=0,IYMAX
    QX=QBINO(NC1+IY)*QBINV(NC2+IY-IX)
    IF(BTEST(IY,0))THEN
     QS=QS-QX
    ELSE
     QS=QS+QX
    END IF
   END DO
   QCX=QBINO(NC3+IX)*QBINO(NC4-IX)*QS
   IF(BTEST(IX,0))THEN
    QCS=QCS-QCX
   ELSE
    QCS=QCS+QCX
   END IF
  END DO
  QTU3R3=QTU3R3+QAS*QBS*QCS*QBINO(NTB+IW)/DFLOAT(MC2+1)

 END IF

 IF(IW/=IPB)THEN
  MA1=MA1+1
  NA1=NA1+MA1
  NA2=NA2-MA2
  MA2=MA2-1
  MB1=MB1+1
  NB1=NB1+MB1
  NB2=NB2-MB2
  MB2=MB2-1
  NC1=NC1-MC1
  MC1=MC1-1
  NC2=NC2-MC2
  MC2=MC2-1
 END IF
END DO
QTU3R3=QTS*QTU3R3
RETURN

CONTAINS
 FUNCTION ID(I) RESULT(IID)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: I
  INTEGER :: IID
  IID=I*(I+1)/2+1
 END FUNCTION ID
 
 FUNCTION ND(I,J) RESULT(IND)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: I,J
  INTEGER :: IND
  IND=I*(I+1)/2+J+1
 END FUNCTION ND
END FUNCTION transformation_coeff_quad