SUBROUTINE wigner_canonical_extremal&
 (su3irrep1x,su3irrep2x,su3irrep3x,I3,NEC,KR0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,N1,N2,KIMAX1)                          
!-----------------------------------------------------------------------
! CALCULATES EXTREMAL SU(3)-SU(2)xU(1) WIGNER COEFFICIENTS
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         
!                 (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        
!                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  
!                 (LSU,02-87)  J.P.DRAAYER        OVERFLOW CORRECTION   
!                 (LSU,10-89)  J.P.DRAAYER        ZERO OUT RELOCATED    
!                 (UoT,04-97)  C.BAHRI            NX=82
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, MODULE INSTEAD OF COMMON BLOCKS, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTION "INDEX" REPLACED BY INTERNAL SUBPROGRAM, DERIVED DATA TYPE su3irrep
!                                                 INSTEAD OF LAMBDA AND MU, ALLOCATABLE ARRAYS
!                                                                       
!     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
!                 K.T.HECHT, NUCL.PHYS.62(1965)1                        
!     PARAMETERS--(I3) : (1)=GHW, (0)=GLW                               
!       EXTERNAL--N1=MAX(KR0MAX)                                        
!                 N2=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
!*                KIMAX1=MAX(KR0MAX*INDMAX)                             
!       INTERNAL--X1=ABS(N1*NX)                                         
!                 X2=ABS(NX)                                            
!     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL                    
!                 ADJUST INTERNAL PARAMETERS BELOW                      
!*    DIMENSIONS--DEWU3(N1*N2->KIMAX1),J1TA(N2),J2TA(N2),IEA(N2),       
!                 DEWU3P(X1),DZ(X2),J1TAP(X2),IAB(X2),ICD(X2)           
!B      COMMENTS--ASSUME MAX N1=9,NX=42,N2=13244                        
!B                       SET X1=378,X2=42                               
!-                           N1=9,NX=82,N2=95284
!-                           X1=738,X2=82
!                 DZ ARRAY ADDED FOR CORRECTING THE OVERFLOW PROBLEM    
!
! This is a new version of XEWU3
!
! INPUT ARGUMENTS: su3irrep1x,su3irrep2x,su3irrep3x,I3,N1,N2,KIMAX1
! OUTPUT ARGUMENTS: NEC,KR0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA
!-----------------------------------------------------------------------
USE derived_types
USE binomial_coeff_factorials
IMPLICIT NONE
INTEGER, EXTERNAL :: outer_multiplicity,outer_multiplicity_theory
REAL(KIND=8), DIMENSION(1)      :: DEWU3
INTEGER, DIMENSION(1)           :: J1TA,J2TA,IEA

!!REAL(KIND=8), DIMENSION(738)    :: DEWU3P
!!REAL(KIND=8), DIMENSION(82)     :: DZ
!!INTEGER, DIMENSION(82)          :: J1TAP,IAB,ICD

REAL(KIND=8), ALLOCATABLE,DIMENSION(:)    :: DEWU3P,DZ
INTEGER, ALLOCATABLE,DIMENSION(:)          :: J1TAP,IAB,ICD

!B REAL(KIND=8), DIMENSION(378) :: DEWU3P
!B REAL(KIND=8), DIMENSION(42)  :: DZ
!B INTEGER, DIMENSION(42)       :: J1TAP,IAB,ICD

INTEGER                         :: I3,NEC,KR0MAX,INDMAX,N1,N2,KIMAX1,J1TD,J1T,J2TD,J2T,NX,IAH,IBH,ICH,&
                                   IDH,I,NCDMAX,NCDMIN,LL1,MM1,LL2,MM2,IA1,IB1,IC1,IA2,&
                                   IB2,IC2,IS1,IS2,ISS,IE3,IEH,KR0CNT,NCD,NNC,LN1,LN2,INN,IND,IE2,IIE,&
                                   IE1,JJ2TA,JJ2TB,JJ1TA,JJ1TB,J,JJ2T,L,M,JJ1T,INDQ,IQ1,J1TP,INDP,INDPQ,&
                                   KR0,KI,KIP,IP1,IIQ2,IIQ2B,IIQ2A,IQ2B,IZ,IHELP,IX,IY,IN,ID,IP2,IIP2,&
                                   KR0A,KR0B,INC,KR0PA,KR0P,IPH!,NNCMAX,KITEST
REAL(KIND=8)                    :: DC,DN,DD,DS,DMIN
TYPE(su3irrep)                  :: su3irrep1x,su3irrep2x,su3irrep3x,su3irrep1,su3irrep2,su3irrep3

!************************************************************************** BEGINNING OF A BLOCK ADDED BY J.H.
DO I=1,N2
 J1TA(I)=0
 J2TA(I)=0
 IEA(I)=0
END DO
NX=su3irrep2x%lambda+su3irrep2x%mu+1
ALLOCATE(DEWU3P(N1*NX),DZ(NX),J1TAP(NX),IAB(NX),ICD(NX))
!************************************************************************** END OF THE BLOCK ADDED BY J.H.

!************************************************************************** BEGINNING OF A BLOCK COMMENTED OUT BY J.H. BECAUSE OF ITS IRRELEVANCE FOR ALLOCATABLE ARRAYS
! DIMENSION CHECKS (LSU,6-81)-START                                 
!IF(N1>9)THEN
! WRITE(*,FMT='(36H ***** XEWU3 DIMENSION OVERFLOW: N1=,I10)')N1                                                    
! STOP
!END IF
!NX=su3irrep2x%lambda+su3irrep2x%mu+1
!B IF(NX>42)THEN
!IF(NX>82)THEN
! WRITE(*,FMT='(36H ***** XEWU3 DIMENSION OVERFLOW: NX=,I10)')NX                                                    
! STOP
!END IF
! DIMENSION CHECKS (LSU,6-81)-START
!************************************************************************** END OF THE BLOCK COMMENTED OUT BY J.H.
KR0MAX=outer_multiplicity(su3irrep1x,su3irrep2x,su3irrep3x)                   
IF(KR0MAX==0)RETURN
IF(I3.EQ.1)THEN
 su3irrep1=su3irrep1x
 su3irrep2=su3irrep2x
 su3irrep3=su3irrep3x
ELSE	  
 su3irrep1%lambda=su3irrep1x%mu
 su3irrep2%lambda=su3irrep2x%mu
 su3irrep3%lambda=su3irrep3x%mu
 su3irrep1%mu=su3irrep1x%lambda
 su3irrep2%mu=su3irrep2x%lambda
 su3irrep3%mu=su3irrep3x%lambda
END IF
NEC=(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda+2*(su3irrep1%mu+su3irrep2%mu-su3irrep3%mu))/3
IAH=(su3irrep2%lambda+su3irrep3%lambda-su3irrep1%lambda-NEC)/2
IBH=(su3irrep3%lambda+su3irrep1%lambda-su3irrep2%lambda+NEC+2)/2
ICH=(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda-NEC)/2
IDH=(su3irrep1%lambda+su3irrep2%lambda+su3irrep3%lambda-NEC+2)/2
DO I=1,NEC
 IAB(I)=(IAH+I)*(IBH-I)
 ICD(I)=(ICH+I)*(IDH+I)
END DO
NCDMAX=outer_multiplicity_theory(su3irrep1,su3irrep2,su3irrep3) ! outer_multiplicity_theory >= outer_multiplicity
NEC=NEC-NCDMAX
su3irrep2=su3irrep2+(-NCDMAX) ! Operator + (defined in module derived_types) adds an integer to lambda and mu
NCDMIN=1
DO WHILE ((NCDMIN/=NCDMAX).AND.(outer_multiplicity(su3irrep1,su3irrep2+1,su3irrep3)<=0))
 NEC=NEC+1
 su3irrep2=su3irrep2+1
 NCDMIN=NCDMIN+1
END DO
!!************************************************************************** BEGINNING OF A BLOCK COMMENTED OUT BY J.H. BECAUSE OF ITS IRRELEVANCE FOR ALLOCATABLE ARRAYS
! DIMENSION MODIFICATION (LSU,6-81)-START
!NNCMAX=NEC+NCDMAX-NCDMIN+2
!KITEST=KR0MAX*(NNCMAX)*(NNCMAX+1)*(NNCMAX+2)/6
!IF(KITEST>KIMAX1)THEN !DIMENSION CHECKS (LSU,6-81)
! WRITE(*,FMT='(40H ***** XEWU3 DIMENSION OVERFLOW: KITEST=,I10,5X,7HKIMAX1=,I10)')KITEST,KIMAX1
! STOP
!END IF
! DIMENSION MODIFICATION (LSU,6-81)--STOP
!*************************************************************************** END OF THE BLOCK COMMENTED OUT BY J.H.
!DO I=1,KITEST !************************************************************ THIS HAS BEEN COMMENTED OUT BY J.H.
DO I=1,KIMAX1 !************************************************************* THIS HAS BEEN ADDED BY J.H.
 DEWU3(I)=0.D0
ENDDO
LL1=su3irrep1%lambda+1
MM1=su3irrep1%mu+1
LL2=su3irrep2%lambda+1
MM2=su3irrep2%mu+1
IA1=2*su3irrep1%lambda+4*su3irrep1%mu
IB1=4*su3irrep1%lambda+2*su3irrep1%mu
IC1=IB1-IA1
IA2=2*su3irrep2%lambda+4*su3irrep2%mu
IB2=4*su3irrep2%lambda+2*su3irrep2%mu
IC2=IB2-IA2
IS1=LL1+MM1
IS2=LL2+MM2
ISS=MM1+su3irrep2%lambda+su3irrep2%mu-NEC
IE3=-(su3irrep3%lambda+2*su3irrep3%mu)
IEH=-(su3irrep2%lambda+2*su3irrep2%mu+3)
KR0CNT=0
DO NCD=NCDMIN,NCDMAX
 NEC=NEC+1
 su3irrep2=su3irrep2+1
 NNC=NEC+1
 INDMAX=NNC*(NNC+1)*(NNC+2)/6
 IA2=IA2+6
 IB2=IB2+6
 IS2=IS2+2
 ISS=ISS+1
 IEH=IEH-3
 LL2=su3irrep2%lambda+1
 MM2=su3irrep2%mu+1
 LN1=su3irrep1%lambda+NEC
 LN2=su3irrep2%lambda+NEC
 INN=NEC*NNC/2
 IF(NCD/=NCDMIN)THEN
!!  DO I=1,KITEST !********************************************************** THIS HAS BEEN COMMENTED OUT BY J.H.
  DO I=1,KIMAX1 !************************************************************ THIS HAS BEEN ADDED BY J.H.
   DEWU3(I)=0.D0
  END DO
 END IF
 DO IND=1,INDMAX
  IEA(IND)=-1000
  J2TA(IND)=-1000
  J1TA(IND)=-1000
 END DO
 IE2=IEH
 I=1000
 DO IIE=1,NNC
  IE2=IE2+3
  IE1=IE3-IE2
  J2TD=IIE-1
  J1TD=NNC-IIE
  JJ2TA=IA2-IE2
  JJ2TB=IB2+IE2
  IF(JJ2TB<JJ2TA)JJ2TA=JJ2TB
  JJ2TA=JJ2TA/3+1
  JJ2TB=JJ2TA-ABS(IC2-IE2)/3
  JJ1TA=IA1-IE1
  JJ1TB=IB1+IE1
  IF(JJ1TB<JJ1TA)JJ1TA=JJ1TB
  JJ1TA=JJ1TA/3+1
  JJ1TB=JJ1TA-ABS(IC1-IE1)/3
  J=0
  DO JJ2T=1,JJ2TB,2
   J2T=JJ2TA-JJ2T
   L=ABS(J2T-su3irrep3%lambda)
   M=J2T+su3irrep3%lambda
   DO JJ1T=1,JJ1TB,2                                              
    J1T=JJ1TA-JJ1T
    IF(J1T<L)CYCLE
    IF(J1T>M)CYCLE
    IND=INDEX(J1TD,su3irrep1%lambda,J1T,J2TD,su3irrep2%lambda,J2T)
    IEA(IND)=IE2
    J2TA(IND)=J2T
    J1TA(IND)=J1T
    J=J+1
   END DO
  END DO
  IF(J<I)I=J
 END DO
 IF(I==0)CYCLE
 IF(KR0CNT/=0)THEN
 
!     GENERATE <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>
!     FROM <(LAM1,MU1)????;(LAM2-1,MU2-1)HIGH::KR0(LAM3,MU3)HIGH>

  INDQ=-KR0MAX
  DO IND=1,NNC
   INDQ=INDQ+KR0MAX
   J1T=J1TA(IND)                                                     
   IF(J1T<0)CYCLE
   IQ1=(LN1-J1T)/2
   IF(IQ1>0)THEN
    J1TP=J1T+1
    INDP=(LN1-J1TP-1)/2+1
    IF(J1TAP(INDP)>=0)THEN
     I=IAB(IQ1)*IQ1*(LL1-IQ1)*(IS1-IQ1)
     DC=-DSQRT(DFLOAT(I)/DFLOAT((J1T+2)*J1TP))
     INDPQ=(INDP-1)*KR0MAX
     DO KR0=1,KR0CNT
      KI=KR0+INDQ
      KIP=KR0+INDPQ
      DEWU3(KI)=DC*DEWU3P(KIP)
     END DO
    END IF
   END IF
   IP1=NEC-IQ1
   IF(IP1<=0)CYCLE
   J1TP=J1T-1
   INDP=(LN1-J1TP-1)/2+1
   IF(J1TAP(INDP)<0)CYCLE
   I=ICD(NNC-IND)*IP1*(MM1-IP1)*(LL1+IP1)
   DC=DSQRT(DFLOAT(I)/DFLOAT((J1TP+2)*J1T))
   INDPQ=(INDP-1)*KR0MAX
   DO KR0=1,KR0CNT
    KI=KR0+INDQ
    KIP=KR0+INDPQ
    DEWU3(KI)=DEWU3(KI)+DC*DEWU3P(KIP)
   END DO
  END DO
  INC=0
 END IF
 IF(KR0CNT/=KR0MAX)THEN
 
!     EVALUATE <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!     WITH (LAM2,MU2) A MINIMUM FOR KR0=KR0CNT

  KR0CNT=KR0CNT+1
  I=0
  IND=INDEX(0,su3irrep1%lambda,su3irrep1%lambda,NEC,su3irrep2%lambda,LN2)-1
  INDQ=-KR0MAX
  DO IIQ2=1,NNC                                                  
   INDQ=INDQ+KR0MAX
   IND=IND+1
   KI=KR0CNT+INDQ
   DEWU3P(KI)=0.D0
   IF(J1TA(IND)<0)CYCLE
   I=I+1
   IIQ2B=IIQ2
  END DO
!                                                                       
!     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)-START                  
!      ****DLOGB CHANGED TO SAVE SPACE (LSU,3-83)****                   
!      ****FURTHER OVERFLOW CORRECTION (LSU,2-87)****                   
!                                                                       
  IIQ2A=IIQ2B-I+1
  IQ2B=IIQ2B-1
  INDQ=(IIQ2A-2)*KR0MAX
  IZ=0
  DO IIQ2=IIQ2A,IIQ2B
   IZ=IZ+1
   DZ(IZ)=1.D0
   INDQ=INDQ+KR0MAX
   L=LL2-IIQ2
!                                                                       
! --> START NUMERATOR PRODUCT LOOP                                      
!                                                                       
   IX=L-IIQ2+NNC+1
   IF(IX==0)THEN !"GO TO 120" was originally here
	IHELP=0 !IHELP is an indicator of how the end of the loop was reahed (if via EXIT then IHELP=0 else IHELP=1)
	EXIT
   ELSE
	IHELP=1
   END IF
   IY=ABS(IX)
   IN=IX/IY
   DN=DLOG(DFLOAT(IY))
   IF(IIQ2A/=IIQ2B)THEN
    DO I=IIQ2A,IQ2B
     J=NNC-I
     IF(I<IIQ2)THEN
      K=IAB(J)*ICD(J)*(IS2-I)
     ELSE
      K=MM2-J
     ENDIF
     IF(K==0)EXIT !"GO TO 120" was originally here
     IF(K<0)IN=-IN
     DN=DLOG(DFLOAT(ABS(K)))+DN
    END DO
    IF(K==0)THEN
	 IHELP=0
	 EXIT
!	ELSE !This is probably redundant
!	 IHELP=1 !This is probably redundant
	END IF
   END IF
   DN=DN+DLOG(DBINO(INN+IIQ2))
!                                                                       
! --> END NUMERATOR PRODUCT LOOP & START DENOMINATOR PRODUCT LOOP       
!                                                                       
   ID=1
   DD=0.D0
   DO I=1,NNC
    IX=I+L
    IF(IX<0)ID=-ID
    DD=DLOG(DFLOAT(I+L))+DD
   END DO
!                                                                       
! --> END DENOMINATOR PRODUCT LOOP & START INNER PRODUCT/SUM LOOP       
!                                                                       
   IP2=NNC-IIQ2
!                                                                       
!     MULTIPLY BY SMALL NUMBER --> DEXP(-172) LIMIT FOR IBM SYSTEMS     
!                                                                       
   DZ(IZ)=DEXP(-DMIN1(DLOGF(2*IP2),172.D0))
   IIP2=IP2+1
   M=IP2*IIP2/2
   DS=0.D0
   DO I=1,IIP2
    DC=DZ(IZ)*DBINO(I+M)
    IF(IIP2/=1)THEN
     DO J=1,IP2
      IF(J<I)THEN
       K=(J+L)*(ISS+J)
      ELSE
       K=IAB(J)
      ENDIF
      DC=DFLOAT(K)*DC
     END DO
    END IF
    DS=DS+DC
   END DO
!                                                                       
! --> END INNER PRODUCT/SUM LOOP & ASSIGN unnormalized DEWU3P VALUE     
!                                                                       
   IF(2*(IP2/2)/=IP2)DS=-DS
   KI=KR0CNT+INDQ
   DEWU3P(KI)=DFLOAT(IN*ID)*DS*DEXP((DN-DD)/2.D0)
  END DO
  
  IF(IHELP/=0)THEN ! IHELP is an indicator of how the end of the loop was reahed (if via EXIT then IHELP=0)
  
!                                                                       
! --> START renormalization PROCEDURE                                   
!                                                                       
   DMIN=1.D0
   IZ=0
   DO IIQ2=IIQ2A,IIQ2B
    IZ=IZ+1
    DMIN=DMIN1(DMIN,DZ(IZ))
   END DO
   INDQ=(IIQ2A-2)*KR0MAX
   IZ=0
   DO IIQ2=IIQ2A,IIQ2B
    IZ=IZ+1
    INDQ=INDQ+KR0MAX
    KI=KR0CNT+INDQ
    DEWU3P(KI)=(DMIN/DZ(IZ))*DEWU3P(KI)
   END DO
!                                                                       
!     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)--STOP                  
!                                                                       
   KR0A=KR0CNT
   KR0B=KR0CNT
   
!     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!     FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>

   CALL support_wigner_canonical&
    (1,su3irrep1,su3irrep2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,IAB,ICD,INDMAX,DEWU3,KR0MAX)
   INC=1

  ELSE

   KR0CNT=KR0CNT-1
   WRITE(1,FMT='(28H *****U3 COUPLING ERROR*****,3X,3(4X,2I3),3X,15HRO(ABSOLUTE) = ,&
    I2,3X,15HRO(RELATIVE) = ,I2,4X,26H*****REPORT TO AUTHOR*****)')&
    su3irrep1x%lambda,su3irrep1x%mu,su3irrep2x%lambda,su3irrep2x%mu,su3irrep3x%lambda,su3irrep3x%mu,KR0MAX,KR0CNT

  END IF
	  
 END IF
	  
 IF(KR0CNT==0)CYCLE
 INDQ=-KR0MAX
 DO IND=1,NNC
  INDQ=INDQ+KR0MAX
  J1TAP(IND)=J1TA(IND)
  DO KR0=1,KR0CNT
   KI=KR0+INDQ
   DEWU3P(KI)=DEWU3(KI)
  END DO
 END DO
END DO
IF(KR0CNT==0)RETURN
KR0A=1
KR0B=KR0CNT-INC
IF(KR0B/=0)THEN

!     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!     FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>

 CALL support_wigner_canonical&
  (0,su3irrep1,su3irrep2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,IAB,ICD,INDMAX,DEWU3,KR0MAX)
END IF

!     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)-START

DO KR0=1,KR0CNT
 DC=1.D0
 INDQ=-KR0MAX
 DO IND=1,INDMAX
  INDQ=INDQ+KR0MAX
  KI=KR0+INDQ
  DC=DMAX1(DC,DABS(DEWU3(KI)))
 END DO
 INDQ=-KR0MAX
 DO IND=1,INDMAX
  INDQ=INDQ+KR0MAX
  KI=KR0+INDQ
  DEWU3(KI)=DEWU3(KI)/DC
 END DO
END DO

!     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)--STOP
!     ORTHONORMALIZATION OF SOLUTIONS

DO KR0=1,KR0CNT
 KR0PA=1
 IF(INC==1.AND.KR0==KR0CNT)KR0PA=KR0CNT
 DO KR0P=KR0PA,KR0
  DN=0.D0
  INDQ=-KR0MAX
  DO IND=1,INDMAX
   INDQ=INDQ+KR0MAX
   KI=KR0+INDQ
   KIP=KR0P+INDQ
   DN=DN+DEWU3(KI)*DEWU3(KIP)
  END DO
  IF(KR0P/=KR0)THEN
   IF(DABS(DN)<1.D-12)CYCLE
   INDQ=-KR0MAX
   DO IND=1,INDMAX
    INDQ=INDQ+KR0MAX
    KI=KR0+INDQ
    KIP=KR0P+INDQ
    DEWU3(KI)=DEWU3(KI)-DN*DEWU3(KIP)
   END DO
  ELSE
   DN=1.D0/DSQRT(DN)
   INDQ=-KR0MAX
   DO IND=1,INDMAX
    INDQ=INDQ+KR0MAX
    KI=KR0+INDQ
    DEWU3(KI)=DN*DEWU3(KI)
   END DO
  END IF
 END DO
END DO

!     SET PHASE CONVENTION (K.T.HECHT, NUCL.PHYS.62(1965)1)

IE2=IE3+(su3irrep1%lambda+2*su3irrep1%mu)
IPH=2*(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda+su3irrep1%mu+su3irrep2%mu-su3irrep3%mu+KR0MAX)                         
INDQ=-KR0MAX
DO IND=1,INDMAX
 INDQ=INDQ+KR0MAX
 IF(IEA(IND)/=IE2)CYCLE
 I=IPH+J1TA(IND)+J2TA(IND)-su3irrep3%lambda
 DO KR0=1,KR0MAX
  I=I-2
  KI=KR0+INDQ
  J=I
  IF(DEWU3(KI)<0.D0)J=J-2
  IF(4*(J/4)==J)CYCLE
  INDPQ=-KR0MAX
  DO INDP=1,INDMAX
   INDPQ=INDPQ+KR0MAX
   KIP=KR0+INDPQ
   DEWU3(KIP)=-DEWU3(KIP)
  END DO
 END DO
 EXIT
END DO
IF(I3==1)RETURN
INDQ=-KR0MAX
DO IND=1,INDMAX
 INDQ=INDQ+KR0MAX
 IEA(IND)=-IEA(IND)
 I=IPH+J1TA(IND)+J2TA(IND)-su3irrep3%lambda
 DO KR0=1,KR0MAX
  I=I-2
  KI=KR0+INDQ
  IF(4*(I/4)==I)DEWU3(KI)=-DEWU3(KI)
 END DO
END DO
DEALLOCATE(DEWU3P,DZ,J1TAP,IAB,ICD)
RETURN

CONTAINS
 FUNCTION INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T) RESULT(IINDEX)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: J1TD,LAM1,J1T,J2TD,LAM2,J2T
  INTEGER :: IINDEX
  IINDEX=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2
 END FUNCTION INDEX
END SUBROUTINE wigner_canonical_extremal