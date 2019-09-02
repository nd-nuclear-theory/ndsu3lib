SUBROUTINE support_wigner_canonical&
 (INC,su3irrep1,su3irrep2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,IAB,ICD,INDMAX,DEWU3,KR0MAX)
!---------------------------------------------------------------------------------
! SUPPORT ROUTINE FOR SUBROUTINE wigner_canonical_extremal
!---------------------------------------------------------------------------------
!     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3
!                 (LSU,03-97)  J.P.DRAAYER        INTEGER OVERFLOW FIX
!                              C.BAHRI
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTION "INDEX" REPLACED BY INTERNAL SUBPROGRAM,
!                                                 DERIVED DATA TYPE su3irrep INSTEAD OF LAMBDA AND MU
!
!     PARAMETERS--
!       INC=0:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!             FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>
!       INC=1:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!             FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!
! This is a new version of XEWU3S
!
! INPUT ARGUMENTS: INC,su3irrep1,su3irrep2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,IAB,ICD,INDMAX,KR0MAX
! OUTPUT ARGUMENT: DEWU3
!
! INC=0 corresponds to formula (18) and INC=1 corresponds to formula (21)
! In the comments in expression like a(b) a corresponds to formula (18) and b corresponds to formula (21)
! In formulae (18) and (21) "the first term" corresponds to Lambda'1(2)=Lambda1(2)+1/2 and "the second term" corresponds to Lambda'1(2)=Lambda1(2)-1/2
!---------------------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
REAL(KIND=8), DIMENSION(1) :: DEWU3,DEWU3P
INTEGER, DIMENSION(1)      :: J1TA,IAB,ICD
INTEGER                    :: INC,NEC,NNC,KR0A,KR0B,INDMAX,KR0MAX,INDQ,INDPQ,IIQ,KR0,KI,KIP,L1,M1,L2,M2,&
                              LL1,MM1,LL2,MM2,LM1,LM2,J2TD,IIQ2,IIQ1,J1TD,J2D,J1D,IIQ2A,IIQ2B,IIQ1A,IIQ1B,&
							  IQ2,IP2,J2T,JJ2T,IQ,IP,IQ2P,IP2P,IQ2D,IP2D,NM,J2TP,JTA,JTB,NQD,IQ1,IP1,J1T,&
							  IND,IQ1P,IP1P,J1TP,INDP,IHELP,J123,I
REAL(KIND=8)               :: DM,DN,DC,DX
TYPE(su3irrep)             :: su3irrep1,su3irrep2
!W    write(6,*) INC,su3irrep1%lambda,su3irrep1%mu,su3irrep2%lambda,su3irrep2%mu,&
!W    NEC,NNC,KR0A,KR0B,DEWU3P(1),J1TA(1),IAB(1),ICD(1),INDMAX,DEWU3(1),KR0MAX
INDQ=-KR0MAX
IF(INC==1)INDQ=INDQ+(INDMAX-NNC)*KR0MAX
INDPQ=-KR0MAX
DO IIQ=1,NNC
 INDQ=INDQ+KR0MAX
 INDPQ=INDPQ+KR0MAX
 DO KR0=KR0A,KR0B    ! For INC=1: KR0A=KR0B=rho
  KI=KR0+INDQ
  KIP=KR0+INDPQ
  DEWU3(KI)=DEWU3P(KIP)
 END DO
END DO
IF(NEC==0)RETURN
IF(INC==1)THEN
 L1=su3irrep2%lambda ! L1=lam2
 M1=su3irrep2%MU     ! M1=mu2
 L2=su3irrep1%lambda ! L2=lam1
 M2=su3irrep1%mu     ! M2=mu1
ELSE
 L1=su3irrep1%lambda ! L1=lam1
 M1=su3irrep1%mu     ! M1=mu1
 L2=su3irrep2%lambda ! L2=lam2
 M2=su3irrep2%mu     ! M2=mu2
END IF
LL1=L1+1             ! LL1=lam1(2)+1
MM1=M1+1             ! MM1=mu1(2)+1
LL2=L2+1             ! LL2=lam2(1)+1
MM2=M2+1             ! MM2=mu2(1)+1
LM1=LL1+MM1          ! LM1=lam1(2)+mu1(2)+2
LM2=LL2+MM2          ! LM2=lam2(1)+mu2(1)+2
DO J2TD=1,NEC        ! J2TD=(eps2-eps2_HW)/3+1 ?
 J1TD=NEC-J2TD       ! J1TD=(eps1-eps1_HW)/3-1 ?
 J2D=J2TD-1          ! J2D=(eps2-eps2_HW)/3 ?
 J1D=J1TD+1          ! J1D=(eps1-eps1_HW)/3 ?
 IIQ2A=J2TD-M2
 IF(IIQ2A<0)IIQ2A=0
 IIQ2A=IIQ2A+1
 IIQ2B=J2TD
 IF(L2<IIQ2B)IIQ2B=L2
 IIQ2B=IIQ2B+1
 IIQ1A=J1TD-M1
 IF(IIQ1A<0)IIQ1A=0
 IIQ1A=IIQ1A+1
 IIQ1B=J1TD
 IF(L1<IIQ1B)IIQ1B=L1
 IIQ1B=IIQ1B+1
 DO IIQ2=IIQ2A,IIQ2B
  IQ2=IIQ2-1
  IP2=J2TD-IQ2
  J2T=L2+IP2-IQ2                                  ! J2T=2*Lam'2(1) - see Eq.(1) for Lambda: 2*Lambda=lambda+\tilde{p}-\tilde{q}
  JJ2T=J2T+1                                      ! JJ2T=2*Lam'2(1)+1
  IQ=-1
  IP=-1
  IF(IP2/=0)THEN
   IQ2P=IQ2
   IP2P=IP2-1
   IQ2D=0
   IP2D=1
   IF(INC==1)IP=1
!B    NM=IP2*(M2-IP2P)*(LL2+IP2)                                        
   DM=DFLOAT(IP2)*DFLOAT(M2-IP2P)*DFLOAT(LL2+IP2) ! Eq.(17): DM=S(q2(1)) where IP2P=\tilde{p2(1)}=mu2(1)-q2(1) and IP2=IP2P+1
   DN=DFLOAT(J2T)                                 ! DN=2*Lam'2(1)=2*Lam2(1)+1
  ELSE
   IQ2P=IQ2-1
   IP2P=IP2
   IQ2D=1
   IP2D=0
   IF(INC==0)IQ=1
!B    NM=IQ2*(L2-IQ2P)*(LM2-IQ2)                                        
   DM=DFLOAT(IQ2)*DFLOAT(L2-IQ2P)*DFLOAT(LM2-IQ2) ! Eq.(17): NM=R(p2(1)) where IQ2P=\tilde{q2(1)}=lam2(1)-p2(1) and IQ2=IQ2P+1
   DN=DFLOAT(J2T+2)                               ! DN=2*Lam'2(1)+2=2*Lam2(1)+1
  END IF
  J2TP=L2+IP2P-IQ2P
  JTA=J2TD-IQ2P
  JTB=NNC-JTA
  NQD=NEC-IQ2P
  DO IIQ1=IIQ1A,IIQ1B
   IQ1=IIQ1-1                                     ! Eq.(18)((21)): IQ1=\tilde{q1(2)}=lam1(2)-p1(2)
   IP1=J1TD-IQ1                                   ! Eq.(18)((21)): IP1=\tilde{p1(2)}=mu1(2)-q1(2)
   J1T=L1+IP1-IQ1                                 ! J1T=2*Lam1(2) - see formula (1) for Lambda
   IF(INC==0)IND=INDEX(J1TD,L1,J1T,J2TD,L2,J2T)
   IF(INC==1)IND=INDEX(J2TD,L2,J2T,J1TD,L1,J1T)
   IF(J1TA(IND)<0)CYCLE
   IF(IP1/=M1)THEN
    IQ1P=IQ1
    IP1P=IP1+1                                    ! Eq.(18)((21)): IP1P=\tilde{p1(2)}+1=mu1(2)-q1(2)+1
    J1TP=J1T+1                                    ! J1TP=2*Lam1(2)+1
    IF(INC/=0)THEN
     INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)
     IF(J1TA(INDP)>=0)THEN
	  IHELP=1
      J123=JTB-IQ1
	 ELSE
	  IHELP=0
	 END IF
    ELSE
     INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)
     IF(J1TA(INDP)>=0)THEN
	  IHELP=1
      J123=JTA+IQ1
	 ELSE
	  IHELP=0
	 END IF
	END IF
	IF(IHELP==1)THEN
     IF(IP2D==1)I=IAB(J123)                       ! Eq.(18)/(21): I=[A(...)+1/2][B(...)+1/2]
     IF(IQ2D==1)I=ICD(NQD-IQ1)                    ! Eq.(18)/(21): I=[C(...)+1/2][D(...)+1/2]
!     I=JJ2T*IP1P*(MM1-IP1P)*(LL1+IP1P)*I
!     DC=DSQRT(DFLOAT(I)/(DFLOAT((J1T+2)*J1TP*NM)*DN))
!B    DX=DFLOAT(JJ2T*IP1P*(MM1-IP1P))*DFLOAT((LL1+IP1P)*I)
!B    DC=DSQRT(DX/(DFLOAT((J1T+2)*J1TP*NM)*DN))                  
     DX=DFLOAT(JJ2T)*DFLOAT(IP1P)*DFLOAT(MM1-IP1P)*DFLOAT((LL1+IP1P)*I) ! Eq.(18)((21)): DX=(2*Lam'2(1)+1)*(X(Y)(Lam'1,Lam'2))^2 in the first term on the right-hand side; where JJ2T=2*Lam'2(1)+1 and the rest is X(Y)^2, where IP1P*(MM1-IP1P)*(LL1+IP1P)=S(q1(2)) - in the first term there is always S(q1(2))
     DC=DSQRT(DX/(DFLOAT((J1T+2)*J1TP)*DM*DN))                          ! Eq.(18)((21)): DC=factor of the Wigner coefficient in the first term on the right-hand side; where DM=N(Lam'1(2)) and DN=2*Lam2(1)+1
!W    write(6,*) ' line', 40, dx, dc
     IF(IQ<0)DC=-DC
     INDQ=(IND-1)*KR0MAX
     INDPQ=(INDP-1)*KR0MAX
     DO KR0=KR0A,KR0B
      KI=KR0+INDQ
      KIP=KR0+INDPQ
      DEWU3(KI)=DC*DEWU3(KIP)                     ! Eq.(18)/(21): DEWU3(KI)=the first term on the right-hand side
     END DO
    END IF
   END IF
   IF(IQ1==L1)CYCLE
   IQ1P=IQ1+1                                     ! Eq.(18)((21)): IQ1P=\tilde{q1(2)}+1=lam1(2)-p1(2)+1
   IP1P=IP1
   J1TP=J1T-1                                     ! J1TP=2*Lam1(2)-1
   IF(INC/=0)THEN
    INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)
    IF(J1TA(INDP)<0)CYCLE
    J123=JTB-IQ1
   ELSE
    INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)
    IF(J1TA(INDP)<0)CYCLE
    J123=JTA+IQ1
   END IF
   IF(IP2D==1)I=ICD(NQD-IQ1)                      ! Eq.(18)/(21): I=[C(...)+1/2][D(...)+1/2]
   IF(IQ2D==1)I=IAB(J123)                         ! Eq.(18)/(21): I=[A(...)+1/2][B(...)+1/2]
!     I=JJ2T*IQ1P*(LL1-IQ1P)*(LM1-IQ1P)*I
!     DC=DSQRT(DFLOAT(I)/(DFLOAT((J1TP+2)*J1T*NM)*DN))
!B    DX=DFLOAT(JJ2T*IQ1P*(LL1-IQ1P))*DFLOAT((LM1-IQ1P)*I)
!B    DC=DSQRT(DX/(DFLOAT((J1TP+2)*J1T*NM)*DN))                  
   DX=DFLOAT(JJ2T)*DFLOAT(IQ1P)*DFLOAT(LL1-IQ1P)*DFLOAT((LM1-IQ1P)*I) ! Eq.(18)((21)): DX=(2*Lam'2(1)+1)*(X(Y)(Lam'1,Lam'2))^2 in the second term on the right-hand side; where JJ2T=2*Lam'2(1)+1 and the rest is X(Y)^2, where IQ1P*(LL1-IQ1P)*(LM1-IQ1P)=R(p1(2)) - in the second term there is always R(p1(2))
   DC=DSQRT(DX/(DFLOAT((J1TP+2)*J1T)*DM*DN))                          ! Eq.(18)((21)): DC=factor of the Wigner coefficient in the second term on the right-hand side; where DM=N(Lam'1(2)) and DN=2*Lam2(1)+1
!W    write(6,*) ' line', 60, dx, dc
   IF(IP<0)DC=-DC
   INDQ=(IND-1)*KR0MAX
   INDPQ=(INDP-1)*KR0MAX
   DO KR0=KR0A,KR0B
    KI=KR0+INDQ
    KIP=KR0+INDPQ
    DEWU3(KI)=DEWU3(KI)+DC*DEWU3(KIP)             ! Eq.(18)/(21): DEWU3(KI)=Wigner coefficient on the left-hand side
   END DO
  END DO
 END DO
END DO
RETURN
	  
CONTAINS
 FUNCTION INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T) RESULT(IINDEX)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: J1TD,LAM1,J1T,J2TD,LAM2,J2T
  INTEGER :: IINDEX
  IINDEX=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2
 END FUNCTION INDEX
END SUBROUTINE support_wigner_canonical