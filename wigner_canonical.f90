SUBROUTINE wigner_canonical(su3irrep1,su3irrep2,su3irrep3,IE,JT,NEC,DEWU3,KR0MAX,&
 INDMAX,DWU3,J1SMAX,J1TMAX,J2SMAX,J2TMAX,IESMAX,IE2MAX,INDMAT,N2,KIMAX2)                                           
!-----------------------------------------------------------------------
! CALCULATES SU(3)-SU(2)xU(1) WIGNER COEFFICIENTS            
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING OF DEWU3     
!                 (LSU,11-89)  J.P.DRAAYER        DWU3 ZERO-OUT RANGE
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTIONS "INDEX" AND "IDM" REPLACED BY INTERNAL SUBPROGRAMS,
!                                                 DERIVED DATA TYPE su3irrep INSTEAD OF LAMBDA AND MU, ALLOCATABLE ARRAYS
!                                                                       
!     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
!                 K.T.HECHT, NUCL.PHYS.62(1965)1                        
!     PARAMETERS--SEE ALSO XEWU3                                        
!       EXTERNAL--N1=MAX(KR0MAX)                                        
!                 N2=MAX(IESMAX) SAFE TO SET N2=MAX(LAM2+MU2+1)         
!                 N3=MAX(DIM(LAM2,MU2))                                 
!                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
!                 NB=MAX(J2SMAX) SAFE TO SET NB=MAX(LAM2+MU2+1)         
!*                KIMAX1=MAX(KR0MAX*INDMAX) (SEE XEWU3)                 
!*                KIMAX2=MAX(KR0MAX*DIM(LAM2,MU2))                      
!       INTERNAL--X1=ABS(N1*N3)                                         
!                 X2=ABS(N2)                                            
!                 X3=ABS(N2*NB)                                         
!     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT          
!                 ADJUST INTERNAL PARAMETERS BELOW                      
!*    DIMENSIONS--DEWU3(N1*NA),DWU3(N1*N3),J1SMAX(N2*NB),               
!                 J1TMAX(N2*NB),J2SMAX(N2),J2TMAX(N2),INDMAT(N2*NB),    
!                 DWU3P(X1),J2TMAP(X2),INDMAP(X3)                       
!*      COMMENTS--USE N1*NA->KIMAX1,N1*N3->KIMAX2                       
!                 ASSUME MAX N1=9,N2=42,N3=9030                         
!                        SET X1=27090,X2=42,X3=1764 (X1=3*N3,FIXED)
!
! This is a new version of XWU3
!
! INPUT ARGUMENTS: su3irrep1,su3irrep2,su3irrep3,IE,JT,NEC,DEWU3,KR0MAX,N1,N2,N3,KIMAX2
! OUTPUT ARGUMENTS: INDMAX,DWU3,J1SMAX,J1TMAX,J2SMAX,J2TMAX,IESMAX,IE2MAX,INDMAT
!-----------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
REAL(KIND=8), DIMENSION(1)     :: DEWU3,DWU3
INTEGER, DIMENSION(1)          :: J1SMAX,J1TMAX,J2SMAX,J2TMAX,INDMAT

!REAL(KIND=8), DIMENSION(27090) :: DWU3P
!INTEGER, DIMENSION(42)         :: J2TMAP
!INTEGER, DIMENSION(1764)       :: INDMAP

REAL(KIND=8), ALLOCATABLE,DIMENSION(:) :: DWU3P
INTEGER, ALLOCATABLE,DIMENSION(:)      :: J2TMAP,INDMAP
INTEGER                        :: IE,JT,NEC,KR0MAX,INDMAX,IESMAX,IE2MAX,N2,KIMAX2,LL1,MM1,LL2,&
                                  MM2,LL3,MM3,LM1,LM2,LLMM1,LLMM2,LLMM3,JJTD,IP,NCC,INC,IQ3,IP3,J3T,JJ3TD,N,&
								  NM,JJ2TDA,JJ2TDB,JJ2TDC,IND,IES,JJ2TD,J2TD,J1TD,IIQ2A,IIQ2B,IIQ1A,IIQ1B,J2S,&
								  IIQ2,IQ2,IP2,J2T,J23S,J23D,J23H,J1S,IIQ1,IQ1,IP1,J1T,J1TS,J2TS,INDQ,JA,JJA,&
								  JB,JJB,JC,JJC,JD,JJD,IESP,I,J2TP,J1TP,M,J12TP,IAB,ICD,IPH,J2SP,INDP,INDPQ,&
								  KR0,KI,KIP,IESJ2S,J3TP,JJ2TDP,J2SB,J2SQ,J1SB!,N1,N3,IDTEST
REAL(KIND=8)                   :: DC
TYPE(su3irrep)                 :: su3irrep1,su3irrep2,su3irrep3

!************************************************************************** BEGINNING OF A BLOCK ADDED BY J.H.
ALLOCATE(DWU3P(KIMAX2),J2TMAP(N2),INDMAP(N2*N2))
!***************************************** ZACIATOK BLOKU, KTORY ASI NIE JE NUTNY
DO I=1,KIMAX2
 DWU3(I)=0.D0
END DO
DO I=1,N2*N2
 J1SMAX(I)=0
 J1TMAX(I)=0
 INDMAT(I)=0
END DO
DO I=1,N2
 J2SMAX(I)=0
 J2TMAX(I)=0
END DO
!***************************************** KONIEC BLOKU, KTORY ASI NIE JE NUTNY
!************************************************************************** END OF THE BLOCK ADDED BY J.H.

!************************************************************************** BEGINNING OF A BLOCK COMMENTED OUT BY J.H. BECAUSE OF ITS IRRELEVANCE FOR ALLOCATABLE ARRAYS
!     DIMENSION CHECKS (LSU,6-81)-START                                 
!IF(N1>9)THEN
! WRITE(*,FMT='(35H *****XWU3 DIMENSION OVERFLOW:  N1=,I10)')N1
! STOP
!END IF
!IF(N2>42)THEN
! WRITE(*,FMT='(35H *****XWU3 DIMENSION OVERFLOW:  N2=,I10)')N2
! STOP
!END IF
!IF(N3>9030)THEN
! WRITE(*,FMT='(35H *****XWU3 DIMENSION OVERFLOW:  N3=,I10)')N3
! STOP
!END IF
!IDTEST=KR0MAX*IDM(su3irrep2)
!IF(IDTEST>KIMAX2.OR.IDTEST>27090)THEN
! WRITE(*,FMT='(39H *****XWU3 DIMENSION OVERFLOW:  IDTEST=,I10)')IDTEST
! STOP
!END IF
!     DIMENSION CHECKS (LSU,6-81)--STOP
!************************************************************************** END OF THE BLOCK COMMENTED OUT BY J.H.
LL1=su3irrep1%lambda+1
MM1=su3irrep1%mu+1
LL2=su3irrep2%lambda+1
MM2=su3irrep2%mu+1
LL3=su3irrep3%lambda+1
MM3=su3irrep3%mu+1
LM1=su3irrep1%lambda+su3irrep1%mu
LM2=su3irrep2%lambda+su3irrep2%mu
LLMM1=LL1+MM1
LLMM2=LL2+MM2
LLMM3=LL3+MM3
JJTD=(IE+su3irrep3%lambda+2*su3irrep3%mu)/3+1
IP=(JJTD+JT-LL3)/2
NCC=NEC-1
INC=1
IQ3=0
IP3=-1
J3T=su3irrep3%lambda+IP3-IQ3
DO JJ3TD=1,JJTD
 DO N=1,KIMAX2
! DO N=1,IDTEST
  DWU3(N)=0.D0
 END DO
 NCC=NCC+1
 IF(IP3==IP)INC=0
 IF(INC/=1)THEN
  IQ3=IQ3+1
  J3T=J3T-1
  NM=(LL3-IQ3)*IQ3*(LLMM3-IQ3)
 ELSE
  IP3=IP3+1
  J3T=J3T+1
  NM=(MM3-IP3)*IP3*(LL3+IP3)
 END IF
 JJ2TDA=NCC-LM1
 IF(JJ2TDA<0)JJ2TDA=0
 JJ2TDA=JJ2TDA+1
 JJ2TDB=LM2
 IF(NCC<JJ2TDB)JJ2TDB=NCC
 JJ2TDB=JJ2TDB+1
 JJ2TDC=JJ2TDA
 IND=0
 IES=0
 DO JJ2TD=JJ2TDA,JJ2TDB
  J2TD=JJ2TD-1
  J1TD=NCC-J2TD
  IES=IES+1
  IIQ2A=J2TD-su3irrep2%mu
  IF(IIQ2A<0)IIQ2A=0
  IIQ2A=IIQ2A+1
  IIQ2B=J2TD
  IF(su3irrep2%lambda<IIQ2B)IIQ2B=su3irrep2%lambda
  IIQ2B=IIQ2B+1
  IIQ1A=J1TD-su3irrep1%mu
  IF(IIQ1A<0)IIQ1A=0
  IIQ1A=IIQ1A+1
  IIQ1B=J1TD
  IF(su3irrep1%lambda<IIQ1B)IIQ1B=su3irrep1%lambda
  IIQ1B=IIQ1B+1
  J2S=0
  DO IIQ2=IIQ2A,IIQ2B
   IQ2=IIQ2-1
   IP2=J2TD-IQ2
   J2T=su3irrep2%lambda+IP2-IQ2
   J23S=J2T+J3T
   J23D=J3T-J2T
   J23H=ABS(J23D)
   J1S=0
   DO IIQ1=IIQ1A,IIQ1B
    IQ1=IIQ1-1
    IP1=J1TD-IQ1
    J1T=su3irrep1%lambda+IP1-IQ1
    IF((J1T<J23H).OR.(J1T>J23S))CYCLE
    J1TS=J1T
    J2TS=J2T
    INDQ=IND*KR0MAX
    IND=IND+1
    J1S=J1S+1
    IF(JJ3TD/=1)THEN
     JA=(J23S-J1T)/2
     JJA=JA+1
     JB=(J23D+J1T)/2
     JJB=JB+1
     JC=(J1T+J23S)/2+1
     JJC=JC+1
     JD=(J1T-J23D)/2+1
     JJD=JD-1
     IESP=J2TD-JJ2TDP
     DO I=1,4
      IF(I==1)IESP=IESP+1
      IF(I==3)IESP=IESP+1
      IF((IESP<1).OR.(IESP>IESMAX))CYCLE
	  
      IF(I==1)THEN
	  
       J2TP=J2T+1
       J1TP=J1T
       IF((J1TP<ABS(J2TP-J3TP)).OR.(J1TP>J2TP+J3TP))CYCLE
       M=IQ2
       IF(M==0)CYCLE
       N=(LLMM2-M)*(LL2-M)
       J12TP=J2T+1
       IF(INC==1)THEN
        IAB=JB
        ICD=JD
        IPH=-1
       ELSE
		IAB=JJA
        ICD=JJC
        IPH=1
       END IF
	   
	  ELSE IF(I==2)THEN
	  
       J2TP=J2T-1
       J1TP=J1T
       IF((J1TP<ABS(J2TP-J3TP)).OR.(J1TP>J2TP+J3TP))CYCLE
       M=IP2
       IF(M==0)CYCLE
       N=(LLMM2-MM2+M)*(MM2-M)
       J12TP=J2T
       IF(INC==1)THEN
        IAB=JA
        ICD=JC
        IPH=1
       ELSE
		IAB=JJB
        ICD=JJD
        IPH=1
       END IF
	  
	  ELSE IF(I==3)THEN
	  
       J2TP=J2T
       J1TP=J1T+1
       IF((J1TP<ABS(J2TP-J3TP)).OR.(J1TP>J2TP+J3TP))CYCLE
       M=IQ1
       IF(M==0)CYCLE
       N=(LLMM1-M)*(LL1-M)
       J12TP=J1T+1
       IF(INC==1)THEN
        IAB=JA
        ICD=JD
        IPH=1
       ELSE
		IAB=JJB
        ICD=JJC
        IPH=1
       END IF
	  
	  ELSE
	  
       J2TP=J2T
       J1TP=J1T-1
       IF((J1TP<ABS(J2TP-J3TP)).OR.(J1TP>J2TP+J3TP))CYCLE
       M=IP1
       IF(M==0)CYCLE
       N=(LLMM1-MM1+M)*(MM1-M)
       J12TP=J1T
       IF(INC==1)THEN
        IAB=JB
        ICD=JC
        IPH=1
       ELSE
		IAB=JJA
        ICD=JJD
        IPH=-1
	   END IF
	  
	  END IF
	  
      IF(J12TP<=0)THEN
       IF(INC==1)IAB=1
       IF(INC==0)ICD=1
       DC=1.D0
      ELSE
       DC=DFLOAT(J12TP*(J12TP+1))
      END IF
      J2SP=(J2TMAP(IESP)-J2TP+2)/2
      INDP=(INDMAP(IESP+(J2SP-1)*N2)-J1TP)/2
      DC=DSQRT(DFLOAT(IAB*ICD*M*N)/(DFLOAT(NM)*DC))
      IF(IPH<0)DC=-DC
      INDPQ=(INDP-1)*KR0MAX
      DO KR0=1,KR0MAX
       KI=KR0+INDQ
       KIP=KR0+INDPQ
       DWU3(KI)=DWU3(KI)+DC*DWU3P(KIP)
      END DO
     END DO
    ELSE
     INDPQ=(INDEX(J1TD,su3irrep1%lambda,J1T,J2TD,su3irrep2%lambda,J2T)-1)*KR0MAX
     DO KR0=1,KR0MAX
      KI=KR0+INDQ
      KIP=KR0+INDPQ
      DWU3(KI)=DEWU3(KIP)
     END DO
    END IF
   END DO
   IF(J1S/=0)THEN
    IESJ2S=IES+J2S*N2
    J2S=J2S+1
    J1SMAX(IESJ2S)=J1S
    J1TMAX(IESJ2S)=J1TS+2*(J1S-1)
    INDMAT(IESJ2S)=2*IND+J1TS
   END IF
  END DO
  IF(J2S==0)THEN
   IES=IES-1
   IF(IES==0)JJ2TDC=JJ2TDC+1
  ELSE
   J2SMAX(IES)=J2S
   J2TMAX(IES)=J2TS+2*(J2S-1)
  END IF
 END DO
 IESMAX=IES
 IF(JJ3TD==JJTD)CYCLE
 J3TP=J3T
 JJ2TDP=JJ2TDC
 IND=0
 DO IES=1,IESMAX
  J2TMAP(IES)=J2TMAX(IES)
  J2SB=J2SMAX(IES)
  J2SQ=-N2
  DO J2S=1,J2SB
   J2SQ=J2SQ+N2
   IESJ2S=IES+J2SQ
   INDMAP(IESJ2S)=INDMAT(IESJ2S)
   J1SB=J1SMAX(IESJ2S)
   DO J1S=1,J1SB
    INDQ=IND*KR0MAX
    IND=IND+1
    DO KR0=1,KR0MAX
     KI=KR0+INDQ
     DWU3P(KI)=DWU3(KI)
	END DO
   END DO
  END DO
 END DO
END DO
INDMAX=IND
IE2MAX=-(su3irrep2%lambda+2*su3irrep2%mu)+3*(JJ2TDC-1)+3*(IES-1)
DEALLOCATE(DWU3P,J2TMAP,INDMAP)
RETURN

CONTAINS
 FUNCTION INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T) RESULT(IINDEX)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: J1TD,LAM1,J1T,J2TD,LAM2,J2T
  INTEGER :: IINDEX
  IINDEX=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2
 END FUNCTION INDEX
 
 FUNCTION IDM(su3irrepx) RESULT(IIDM)
  IMPLICIT NONE
  TYPE(su3irrep),INTENT (IN) :: su3irrepx
  INTEGER :: IIDM
  IIDM=(su3irrepx%lambda+1)*(su3irrepx%mu+1)*(su3irrepx%lambda+su3irrepx%mu+2)/2
 END FUNCTION IDM
END SUBROUTINE wigner_canonical