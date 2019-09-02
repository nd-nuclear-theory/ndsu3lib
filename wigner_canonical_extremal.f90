SUBROUTINE wigner_canonical_extremal&
 (su3irrep1x,su3irrep2x,su3irrep3x,I3,NEC,KR0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,NX)                          
!-----------------------------------------------------------------------
! CALCULATES EXTREMAL SU(3)-SU(2)xU(1) REDUCED WIGNER COEFFICIENTS
! <(lam1,mu1)eps1,Lam1;(lam2,mu2)eps2,Lam2||(lam3,mu3)E>_rho
! FOR GIVEN lam1,mu1,lam2,mu2,lam3,mu3
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         
!                 (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        
!                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  
!                 (LSU,02-87)  J.P.DRAAYER        OVERFLOW CORRECTION   
!                 (LSU,10-89)  J.P.DRAAYER        ZERO OUT RELOCATED    
!                 (UoT,04-97)  C.BAHRI            NX=82
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, MODULE INSTEAD OF COMMON BLOCKS, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTION "INDEX" REPLACED BY INTERNAL SUBPROGRAM, DERIVED DATA TYPE su3irrep
!                                                 INSTEAD OF LAMBDA AND MU, AUTOMATIC ARRAYS, MINOR CHANGES
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
!-----------------------------------------------------------------------
! This is a new version of XEWU3
!
! INPUT ARGUMENTS: su3irrep1x,su3irrep2x,su3irrep3x,I3,KR0MAX,NX
! OUTPUT ARGUMENTS: NEC,INDMAX,DEWU3,J1TA,J2TA,IEA
!
! su3irrep1x%lambda=lam1
! su3irrep1x%mu=mu1
! su3irrep2x%lambda=lam2
! su3irrep2x%mu=mu2
! su3irrep3x%lambda=lam3
! su3irrep3x%mu=mu3
! I3=1 FOR E=HW, I3=0 FOR E=LW
! NX=lam2+mu2+1
! NEC=n=(lam1+lam2-lam3+2*(mu1+mu2-mu3))/3
! KR0MAX=MULTIPLICITY OF COUPLING (lam1,mu1)x(lam2,mu2)->(lam3,mu3)
! J1TA(IND)=2*Lam1 WHERE 1<=IND<=INDMAX=(n+1)(n+2)(n+3)/6
! J2TA(IND)=2*Lam2 WHERE 1<=IND<=INDMAX
! IEA(IND)=eps2 WHERE 1<=IND<=INDMAX
! eps1=eps3E-eps2
! FOR GIVEN eps1,Lam1,eps2,Lam2 THE EXTREMAL REDUCED WIGNER COEFFICIENTS ARE THE ELEMENTS OF ARRAY
!  DEWU3 WITH INDECES rho+KR0MAX*(IND-1)
!
! The comments perain to the I3=1 case only!
!-----------------------------------------------------------------------
USE derived_types
USE binomial_coeff_factorials
IMPLICIT NONE
INTEGER, EXTERNAL :: outer_multiplicity,outer_multiplicity_theory
REAL(KIND=8), DIMENSION(1)      :: DEWU3
INTEGER, DIMENSION(1)           :: J1TA,J2TA,IEA

!REAL(KIND=8), DIMENSION(738)    :: DEWU3P
!REAL(KIND=8), DIMENSION(82)     :: DZ
!INTEGER, DIMENSION(82)          :: J1TAP,IAB,ICD

REAL(KIND=8), DIMENSION(KR0MAX*NX)    :: DEWU3P
REAL(KIND=8), DIMENSION(NX)     :: DZ
INTEGER, DIMENSION(NX)          :: J1TAP,IAB,ICD

!B REAL(KIND=8), DIMENSION(378) :: DEWU3P
!B REAL(KIND=8), DIMENSION(42)  :: DZ
!B INTEGER, DIMENSION(42)       :: J1TAP,IAB,ICD

INTEGER                         :: I3,NEC,KR0MAX,INDMAX,J1TD,J1T,J2TD,J2T,NX,IAH,IBH,ICH,&
                                   IDH,I,NCDMAX,NCDMIN,LL1,MM1,LL2,MM2,IA1,IB1,IC1,IA2,&
                                   IB2,IC2,IS1,IS2,ISS,IE3,IEH,KR0CNT,NCD,NNC,LN1,LN2,INN,IND,IE2,IIE,&
                                   IE1,JJ2TA,JJ2TB,JJ1TA,JJ1TB,J,JJ2T,L,M,JJ1T,INDQ,IQ1,J1TP,INDP,INDPQ,&
                                   KR0,KI,KIP,IP1,IIQ2,IIQ2B,IIQ2A,IQ2B,IZ,IHELP,IX,IY,IN,ID,IP2,IIP2,&
                                   KR0A,KR0B,INC,KR0PA,KR0P,IPH,K,NNCMAX,KITEST!,N2,N1,KIMAX1
REAL(KIND=8)                    :: DC,DN,DD,DS,DMIN
TYPE(su3irrep)                  :: su3irrep1x,su3irrep2x,su3irrep3x,su3irrep1,su3irrep2,su3irrep3

!KR0MAX=outer_multiplicity(su3irrep1x,su3irrep2x,su3irrep3x)                   
IF(KR0MAX==0)RETURN

!************************************************************************** BEGINNING OF A BLOCK IRRELEVANT FOR AUTOMATIC ARRAYS
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
!************************************************************************** END OF THE BLOCK IRRELEVANT FOR AUTOMATIC ARRAYS
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
NEC=(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda+2*(su3irrep1%mu+su3irrep2%mu-su3irrep3%mu))/3 ! NEC=n=(lam1+lam2-lam3+2*(mu1+mu2-mu3))/3
IAH=(su3irrep2%lambda+su3irrep3%lambda-su3irrep1%lambda-NEC)/2   ! IAH=a in Eq.(20) ?
IBH=(su3irrep3%lambda+su3irrep1%lambda-su3irrep2%lambda+NEC+2)/2 ! IBH=b in Eq.(20) ?
ICH=(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda-NEC)/2   ! ICH=c in Eq.(20) ?
IDH=(su3irrep1%lambda+su3irrep2%lambda+su3irrep3%lambda-NEC+2)/2 ! IDH=d in Eq.(20) ?
DO I=1,NEC
 IAB(I)=(IAH+I)*(IBH-I) ! IAB(I)=(a+I)*(b-I)
 ICD(I)=(ICH+I)*(IDH+I) ! ICD(I)=(c+I)*(d+I)
END DO
NCDMAX=outer_multiplicity_theory(su3irrep1,su3irrep2,su3irrep3) ! NCDMAX>=eta_max
! outer_multiplicity_theory>=outer multiplicity
NEC=NEC-NCDMAX ! NEC<=n-eta_max
su3irrep2=su3irrep2+(-NCDMAX) ! su3irrep2%lambda<=lam2-eta_max, su3irrep2%mu<=mu2-eta_max
! Operator + (defined in module derived_types) adds an integer to lambda and mu
NCDMIN=1
DO WHILE ((NCDMIN/=NCDMAX).AND.(outer_multiplicity(su3irrep1,su3irrep2+1,su3irrep3)<=0))
 NEC=NEC+1
 su3irrep2=su3irrep2+1
 NCDMIN=NCDMIN+1
END DO
! NEC=n-eta_max
! su3irrep2%lambda=lam2-eta_max, su3irrep2%mu=mu2-eta_max
! NCDMAX-NCDMIN+1=eta_max

! DIMENSION MODIFICATION (LSU,6-81)-START
NNCMAX=NEC+NCDMAX-NCDMIN+2 ! NNCMAX=n+1
KITEST=KR0MAX*(NNCMAX)*(NNCMAX+1)*(NNCMAX+2)/6 ! KITEST=rho_max*(n+1)*(n+2)*(n+3)/6=rho_max*INDMAX=maximal needed index of DEWU3=minimal needed size of DEWU3
!*************************************************************************** BEGINNING OF A BLOCK IRRELEVANT FOR AUTOMATIC DEWU3
!IF(KITEST>KIMAX1)THEN !DIMENSION CHECKS (LSU,6-81)
! WRITE(*,FMT='(40H ***** XEWU3 DIMENSION OVERFLOW: KITEST=,I10,5X,7HKIMAX1=,I10)')KITEST,KIMAX1
! STOP
!END IF
! DIMENSION MODIFICATION (LSU,6-81)--STOP
!*************************************************************************** END OF THE BLOCK IRRELEVANT FOR AUTOMATIC DEWU3

DO I=1,KITEST
 DEWU3(I)=0.D0
ENDDO
LL1=su3irrep1%lambda+1                    ! LL1=lam1+1
MM1=su3irrep1%mu+1                        ! MM1=mu1+1
LL2=su3irrep2%lambda+1                    ! LL2=lam2-eta_max+1
MM2=su3irrep2%mu+1                        ! MM2=mu2-eta_max+1
IA1=2*su3irrep1%lambda+4*su3irrep1%mu     ! IA1=2*lam1+4*mu1
IB1=4*su3irrep1%lambda+2*su3irrep1%mu     ! IB1=4*lam1+2*mu1
IC1=IB1-IA1                               ! IC1=2*lam1-2*mu1
IA2=2*su3irrep2%lambda+4*su3irrep2%mu     ! IA2=2*(lam2-eta_max)+4*(mu2-eta_max)=2*lam2+4*mu2-6*eta_max
IB2=4*su3irrep2%lambda+2*su3irrep2%mu     ! IB2=4*(lam2-eta_max)+2*(mu2-eta_max)=4*lam2+2*mu2-6*eta_max
IC2=IB2-IA2                               ! IC2=2*(lam2-eta_max)-2*(mu2-eta_max)=2*lam2-2*mu2
IS1=LL1+MM1                               ! IS1=lam1+mu1+2
IS2=LL2+MM2                               ! IS2=lam2-eta_max+mu2-eta_max+2=lam2+mu2-2*eta_max
ISS=MM1+su3irrep2%lambda+su3irrep2%mu-NEC ! ISS=mu1+1+lam2-eta_max+mu2-eta_max-n+eta_max=mu1+1+lam2+mu2-n-eta_max
IE3=-(su3irrep3%lambda+2*su3irrep3%mu)    ! IE3=-lam3-2*mu3=eps3_min=eps3_HW
IEH=-(su3irrep2%lambda+2*su3irrep2%mu+3)  ! IEH=-lam2-eta_max-2*(mu2-eta_max)-3=-lam2-2*mu2+eta_max-3=eps2_HW+eta_max-3
KR0CNT=0
DO NCD=NCDMIN,NCDMAX ! This loop is a loop over eta from eta_max-1 to 0 (which is related to rho=eta_max-eta).
 NEC=NEC+1                    ! NEC=n-eta=\bar{n}=(lam1+\bar{lam2}-lam3+2*(mu1+\bar{mu2}-mu3))/3 (see the bottom of Eq.(20))
 su3irrep2=su3irrep2+1        ! su3irrep2%lambda=lam2-eta=\bar{lam2}, su3irrep2%mu=mu2-eta=\bar{mu2}
 NNC=NEC+1                    ! NNC=n-eta+1=\bar{n}+1
 INDMAX=NNC*(NNC+1)*(NNC+2)/6 ! INDMAX=INDMAX for \bar{lam2} and \bar{mu2}
 IA2=IA2+6                    ! IA2=2*lam2+4*mu2-6*eta=2*\bar{lam2}+4*\bar{mu2}
 IB2=IB2+6                    ! IB2=4*lam2+2*mu2-6*eta=4*\bar{lam2}+2*\bar{mu2}
 IS2=IS2+2                    ! IS2=lam2+mu2-2*eta=\bar{lam2}+\bar{mu2}
 ISS=ISS+1                    ! ISS=mu1+1+lam2+mu2-n-eta=mu1+1+\bar{lam2}+\bar{mu2}-\bar{n}
 IEH=IEH-3                    ! IEH=\bar{eps2}_min-3=\bar{eps2}_HW-3=eps2_HW+eta-3
 LL2=su3irrep2%lambda+1       ! LL2=\bar{lam2}+1
 MM2=su3irrep2%mu+1           ! \bar{mu2}+1
 LN1=su3irrep1%lambda+NEC     ! LN1=lam1+\bar{n}
 LN2=su3irrep2%lambda+NEC     ! LN2=\bar{lam2}+\bar{n}
 INN=NEC*NNC/2                ! INN=\bar{n}*(\bar{n}+1)/2
 IF(NCD/=NCDMIN)THEN ! If NCD=NCDMIN, DEWU3 does not need to be zeroized because it already is zeroized - see above
  DO I=1,KR0MAX*INDMAX ! KITEST used to be here instead of KR0MAX*INDMAX (<=KITEST)
   DEWU3(I)=0.D0
  END DO
 END IF
 DO IND=1,INDMAX
  IEA(IND)=-1000
  J2TA(IND)=-1000
  J1TA(IND)=-1000
 END DO
 IE2=IEH ! IE2=\bar{eps2}_HW-3
 I=1000
 DO IIE=1,NNC
  IE2=IE2+3 ! IE2=\bar{eps2}
  IE1=IE3-IE2 ! IE1=eps3_HW-\bar{eps2}=eps1
  J2TD=IIE-1 ! J2TD=(\bar{eps2}-\bar{eps2}_HW)/3=0,...,\bar{n}
  J1TD=NNC-IIE ! J1TD=(eps1-eps1_HW)/3=\bar{n},...,0 , It holds: J1TD+J2TD=\bar{n}
  JJ2TA=IA2-IE2 ! JJ2TA=2*\bar{lam2}+4*\bar{mu2}-\bar{eps2}
  JJ2TB=IB2+IE2 ! JJ2TB=4*\bar{lam2}+2*\bar{mu2}+\bar{eps2}
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
   J2T=JJ2TA-JJ2T ! J2T=2*\bar{Lam2}
   L=ABS(J2T-su3irrep3%lambda) ! L=the lowest possible 2*Lam1 according to the triangular rule: 2*Lam3_HW=lam3
   M=J2T+su3irrep3%lambda ! M=the greatest possible 2*Lam1 according to the triangular rule
   DO JJ1T=1,JJ1TB,2                                              
    J1T=JJ1TA-JJ1T ! J1T=2*Lam1
    IF(J1T<L)CYCLE
    IF(J1T>M)CYCLE
    IND=INDEX(J1TD,su3irrep1%lambda,J1T,J2TD,su3irrep2%lambda,J2T) ! IND is a reference to <(lam1,mu1)eps1,Lam1;(bar{lam2},bar{mu2})\bar{eps2},\bar{Lam2}||(lam1,mu1)HW>_rho
    IEA(IND)=IE2 ! IEA(IND)=\bar{eps2}
    J2TA(IND)=J2T ! J2TA(IND)=2*\bar{Lam2}
    J1TA(IND)=J1T ! J1TA(IND)=2*Lam1
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
   J1T=J1TA(IND)                               ! Eq.(17): J1T=2*Lam1
   IF(J1T<0)CYCLE
   IQ1=(LN1-J1T)/2                             ! Eq.(17): IQ1=\tilde{q_1}=lam1-p_1
   IF(IQ1>0)THEN
    J1TP=J1T+1                                 ! Eq.(17): J1TP=2*Lam1+1
    INDP=(LN1-J1TP-1)/2+1
    IF(J1TAP(INDP)>=0)THEN
     I=IAB(IQ1)*IQ1*(LL1-IQ1)*(IS1-IQ1)        ! Eq.(17): I=numerator in the first square root: IAB(IQ1)=A(...)[B(...)+1] (because Lam1=-lam1/2+n/2+p_1), (LL1-IQ1)*IQ1*(IS1-IQ1)=R(p_1+1)
     DC=-DSQRT(DFLOAT(I)/DFLOAT((J1T+2)*J1TP)) ! Eq.(17): DC=-the first square root
     INDPQ=(INDP-1)*KR0MAX
     DO KR0=1,KR0CNT
      KI=KR0+INDQ
      KIP=KR0+INDPQ
      DEWU3(KI)=DC*DEWU3P(KIP)                 ! Eq.(17): DEWU3(KI)=the first term on the right-hand side
     END DO
    END IF
   END IF
   IP1=NEC-IQ1                                 ! Eq.(17): IP1=\tilde{p_1}=mu1-q_1
   IF(IP1<=0)CYCLE
   J1TP=J1T-1                                  ! Eq.(17): J1TP=2*Lam1-1
   INDP=(LN1-J1TP-1)/2+1
   IF(J1TAP(INDP)<0)CYCLE
   I=ICD(NNC-IND)*IP1*(MM1-IP1)*(LL1+IP1)      ! Eq.(17): I=numerator in the second square root: ICD(NNC-IND)=C(...)D(...), (MM1-IP1)*IP1*(LL1+IP1)=S(q_1+1)
   DC=DSQRT(DFLOAT(I)/DFLOAT((J1TP+2)*J1T))    ! Eq.(17): DC=the second square root
   INDPQ=(INDP-1)*KR0MAX
   DO KR0=1,KR0CNT
    KI=KR0+INDQ
    KIP=KR0+INDPQ
    DEWU3(KI)=DEWU3(KI)+DC*DEWU3P(KIP)         ! Eq.(17): DEWU3(KI)=Wigner coefficient on the left-hand side
   END DO
  END DO
  INC=0
 END IF
 IF(KR0CNT/=KR0MAX)THEN
 
!     EVALUATE <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>
!     WITH (LAM2,MU2) A MINIMUM FOR KR0=KR0CNT

  KR0CNT=KR0CNT+1              ! KR0CNT=rho=1,...,rho_max
  I=0
  IND=INDEX(0,su3irrep1%lambda,su3irrep1%lambda,NEC,su3irrep2%lambda,LN2)-1 ! IND is a reference to <(lam1,mu)HW;(\bar{lam2},\bar{mu2})3*\bar{n}+\bar{eps2}_HW,LN2/2||(lam3,mu3)HW>_rho
  INDQ=-KR0MAX
  DO IIQ2=1,NNC                                                  
   INDQ=INDQ+KR0MAX
   IND=IND+1
   KI=KR0CNT+INDQ
   DEWU3P(KI)=0.D0
   IF(J1TA(IND)<0)CYCLE
   I=I+1
   IIQ2B=IIQ2                  ! Eq.(20): IIQ2B=j_2+1 ?
  END DO
!                                                                       
!     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)-START                  
!      ****DLOGB CHANGED TO SAVE SPACE (LSU,3-83)****                   
!      ****FURTHER OVERFLOW CORRECTION (LSU,2-87)****                   
!                                                                       
  IIQ2A=IIQ2B-I+1              ! Eq.(20): IIQ2A=j_1+1 ?
  IQ2B=IIQ2B-1                 ! Eq.(20): IQ2B=j_2 ?
  INDQ=(IIQ2A-2)*KR0MAX
  IZ=0
  DO IIQ2=IIQ2A,IIQ2B          ! Eq.(20): IIQ2=\tilde{\bar{q_2}}+1
   IZ=IZ+1
   DZ(IZ)=1.D0
   INDQ=INDQ+KR0MAX
   L=LL2-IIQ2                  ! Eq.(20): L=(\bar{lam2}+1)-(\tilde{\bar{q_2}}+1)=\bar{lam2}-\tilde{\bar{q_2}} ?
!                                                                       
! --> START NUMERATOR PRODUCT LOOP                                      
!                                                                       
   IX=L-IIQ2+NNC+1             ! Eq.(20): IX=the expression in front of the binomial coefficient in G
   IF(IX==0)THEN !"GO TO 120" was originally here
	IHELP=0 !IHELP is an indicator of how the end of the loop was reahed (if via EXIT then IHELP=0 else IHELP=1)
	EXIT
   ELSE
	IHELP=1
   END IF
   IY=ABS(IX)
   IN=IX/IY                    ! IN=the sign of IX
   DN=DLOG(DFLOAT(IY))         ! DN=log(|IX|)
   IF(IIQ2A/=IIQ2B)THEN        ! If j_2/=j_1, G is calculated. If j_2=j_1, G=1 i.e. there is no G.
    DO I=IIQ2A,IQ2B            ! Eq.(20): I=j+1 in g(j) ; This loop calculates the logarithm of G in Eq.(20) without the binomial coefficient.
     J=NNC-I                   ! Eq.(20): J=(\bar{n}+1)-(j+1)=\bar{n}-j
     IF(I<IIQ2)THEN            ! if j<\tilde{\bar{q_2}}
      K=IAB(J)*ICD(J)*(IS2-I)  ! Eq.(20): K=g(j)=... for j<\tilde{\bar{q_2}}
     ELSE                      ! if j>=\tilde{\bar{q_2}}
      K=MM2-J                  ! Eq.(20): K=g(j)=... for j>=\tilde{\bar{q_2}}
     ENDIF
     IF(K==0)EXIT !"GO TO 120" was originally here
     IF(K<0)IN=-IN
     DN=DLOG(DFLOAT(ABS(K)))+DN
    END DO                     ! DN=the logarithm of G in Eq.(20) without the binomial coefficient
    IF(K==0)THEN
	 IHELP=0
	 EXIT
!	ELSE !This is probably redundant
!	 IHELP=1 !This is probably redundant
	END IF
   END IF
   DN=DN+DLOG(DBINO(INN+IIQ2)) ! DN=the logarithm of G in Eq.(20)
!                                                                       
! --> END NUMERATOR PRODUCT LOOP & START DENOMINATOR PRODUCT LOOP       
!                                                                       
   ID=1
   DD=0.D0
   DO I=1,NNC                  ! I=1,...,\bar{n}+1
    IX=I+L
    IF(IX<0)ID=-ID
    DD=DLOG(DFLOAT(I+L))+DD    ! DD=log(L+1)+log(L+2)+...+log(L+\bar{n}+1)=log((L+\bar{n}+1)*(L+\bar{n})*...*(L+1))=log((L+\bar{n}+1)!/L!)
   END DO                      ! Eq.(20): DD=log(H*(\bar{n}+1)!) ?????
!                                                                       
! --> END DENOMINATOR PRODUCT LOOP & START INNER PRODUCT/SUM LOOP       
!                                                                       
   IP2=NNC-IIQ2                ! Eq.(20): IP2=\tilde{\bar{p_2}} (=\bar{n}-\tilde{\bar{q_2}})
!                                                                       
!     MULTIPLY BY SMALL NUMBER --> DEXP(-172) LIMIT FOR IBM SYSTEMS     
!                                                                       
   DZ(IZ)=DEXP(-DMIN1(DLOGF(2*IP2),172.D0))
   IIP2=IP2+1                  ! Eq.(20): IIP2=\tilde{\bar{p_2}}+1
   M=IP2*IIP2/2                ! Eq.(20): M=\tilde{\bar{p_2}}*(\tilde{\bar{p_2}}+1)/2
   DS=0.D0
   DO I=1,IIP2                 ! Eq.(20): I=i+1 ; This loop calculates the sum in F in Eq.(20)
    DC=DZ(IZ)*DBINO(I+M)       ! Eq.(20): DC=binomial coefficient in F multiplied by DZ(IZ) whose effect is renormalized below
    IF(IIP2/=1)THEN
     DO J=1,IP2                ! Eq.(20): J=j+1 (j is the index in the product in F) ; This loop calculates the binomial coefficient times the product in F in Eq.(20)
      IF(J<I)THEN              ! if j<i
       K=(J+L)*(ISS+J)         ! Eq.(20): K=f(j)=... for j<i (L=\bar{lam2}-\tilde{\bar{q_2}}=\bar{p_2})
      ELSE                     ! if j>=i
       K=IAB(J)                ! Eq.(20): K=f(j)=... for j>=i
      ENDIF
      DC=DFLOAT(K)*DC
     END DO                    ! Eq.(20): DC=the binomial coefficient times the product in F
    END IF
    DS=DS+DC
   END DO                      ! Eq.(20): DS=the sum in F
!                                                                       
! --> END INNER PRODUCT/SUM LOOP & ASSIGN unnormalized DEWU3P VALUE     
!                                                                       
   IF(2*(IP2/2)/=IP2)DS=-DS    ! Eq.(20): DS=F
   KI=KR0CNT+INDQ
   DEWU3P(KI)=DFLOAT(IN*ID)*DS*DEXP((DN-DD)/2.D0) ! Eq.(20): DEWU3P(KI)=the unnormalized Wigner coefficient (IN and ID take care of the sign)
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
    DEWU3P(KI)=(DMIN/DZ(IZ))*DEWU3P(KI) ! Eq.(20): DEWU3P(KI)=the Wigner coefficient
   END DO
!                                                                       
!     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)--STOP                  
!                                                                       
   KR0A=KR0CNT ! KR0A=rho
   KR0B=KR0CNT ! KR0B-rho
   
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
  IF(4*(I/4)/=I)DEWU3(KI)=-DEWU3(KI)
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
END SUBROUTINE wigner_canonical_extremal
