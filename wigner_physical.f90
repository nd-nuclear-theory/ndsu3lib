SUBROUTINE wigner_physical(I1,J1,su3irrep1,KA1MAX,L1,DON1,I2,J2,su3irrep2,KA2MAX,L2,&
 DON2,I3,J3,su3irrep3,KA3MAX,L3,DON3,KR0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,DWU3R3,N01,N012,N0123,N12,N123)                                    
!-----------------------------------------------------------------------
! Calculates reduced SU(3)-SO(3) Wigner coefficients <(lam1,mu1)kap1,L1;(lam2,mu2)kap2,L2||(lam3,mu3)kap3,L3>_rho
! for fixed lam1,mu1,L1,lam2,mu2,L2,lam3,mu3,L3
!-----------------------------------------------------------------------
!     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3
!                 (ND,2019)    J.HERKO            REWRITTEN IN MODERN FORTRAN, EXPLICIT DECLARATIONS,
!                                                 STATEMENT FUNCTION "KSTART" REPLACED BY INTERNAL SUBPROGRAM,
!                                                 DERIVED DATA TYPE su3irrep INSTEAD OF LAMBDA AND MU,
!                                                 AUTOMATIC ARRAYS
!                                                                       
!     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   
!                 J.P.DRAAYER, NUCL.PHYS.A129(1969)647-665              
!     PARAMETERS--SEE ALSO XEWU3 AND DTU3R3 AND/OR CONMAT               
!       EXTERNAL--N0=MAX(KR0MAX)                                        
!                 N1=MAX(KA1MAX)                                        
!                 N2=MAX(KA2MAX)                                        
!                 N3=MAX(KA3MAX)                                        
!                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) 
!       INTERNAL--X1=ABS(N1)                                            
!                 X2=ABS(N2)                                            
!                 X3=ABS(N1*N2*N3)                                      
!                 X4=ABS(N1*N2)                                         
!     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT AND      
!                 ADJUST INTERNAL PARAMETERS IN STATEMENT 28            
!*    DIMENSIONS--DON1(N1*N1),DON2(N2*N2),DON3(N3*N3),DEWU3(N0*NA),     
!                 J1TA(NA),J2TA(NA),IEA(NA),DWU3R3(N0*N1*N2*N3),        
!                 DT1A(X1),DT2A(X2),DS1A(X3),DS2A(X4)                   
!       COMMENTS--ASSUME MAX N0=9,N1=9,N2=9,N3=9                        
!                        SET X1=9,X2=9,X3=729,X4=81
!
! This is a new version of XWU3R3
!
! INPUT ARGUMENTS: I1,J1,su3irrep1,KA1MAX,L1,DON1,I2,J2,su3irrep2,KA2MAX,L2,DON2,
!                  I3,J3,su3irrep3,KA3MAX,L3,DON3,KR0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,N12,N123
!
! OUTPUT ARGUMENT: DWU3R3
!
! su3irrep1%lambda=lam1
! su3irrep1%mu=mu1
! su3irrep2%lambda=lam2
! su3irrep2%mu=mu2
! su3irrep3%lambda=lam3
! su3irrep3%mu=mu3
! (I,J)=(0,1) FOR lam>=mu, (I,J)=(1,0) FOR lam<mu
! KA1MAX IS THE NUMBER OF OCCURENCES OF L1 IN SU(3) IRREP (lam1,mu1), i.e. THE MAXIMAL VALUE OF kap1
!  (SIMILARLY FOR K2MAX AND K3MAX)
! DON IS AN ARRAY PROVIDED BY SUBROUTINE orthonormalization_matrix FOR GIVEN lam,mu,L
! KR0MAX IS THE MULTIPLICITY OF COUPLING (lam1,mu1)x(lam2,mu2)->(lam3,mu3)
! INDMAX IS PROVIDED BY SUBROUTINE wigner_canonical_extremal
! DEWU3 IS AN ARRAY PROVIDED BY SUBROUTINE wigner_canonical_extremal CONTAINING THE EXTREMAL SU(3)-SU(2)xU(1) REDUCED WIGNER COEFFICIENTS
! J1TA,J2TA,IEA ARE ARRAYS PROVIDED BY SUBROUTINE wigner_canonical_extremal - SEE THE COMMENTS THEREIN
! DWU3R3(N)=<(lam1,mu1)kap1,L1;(lam2,mu2)kap2,L2||(lam3,mu3)kap3,L3>_rho WHERE:
!  N=rho+KR0MAX*(kap1-1)+KR0MAX*KA1MAX*(kap2-1)+KR0MAX*KA1MAX*KA2MAX*(kap3-1)
! N01=KR0MAX*KA1MAX
! N012=N01*KA2MAX
! N0123=N012*KA3MAX
! N12=KA1MAX*KA2MAX
! N123=N12*KA3MAX
!-----------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
REAL(KIND=8), EXTERNAL :: transformation_coeff,clebsch_gordan
REAL(KIND=8), DIMENSION(1) :: DON1,DON2,DON3,DEWU3,DWU3R3
INTEGER, DIMENSION(1) :: J1TA,J2TA,IEA

REAL(KIND=8), DIMENSION(KA1MAX) :: DT1A
REAL(KIND=8), DIMENSION(KA2MAX) :: DT2A
REAL(KIND=8), DIMENSION(N123) :: DS1A
REAL(KIND=8), DIMENSION(N12) :: DS2A

!REAL(KIND=8), DIMENSION(9) :: DT1A,DT2A
!REAL(KIND=8), DIMENSION(729) :: DS1A
!REAL(KIND=8), DIMENSION(81) :: DS2A

INTEGER :: I1,J1,KA1MAX,L1,I2,J2,KA2MAX,L2,I3,J3,KA3MAX,L3,KR0MAX,INDMAX,N01,N012,N0123,&
           N12,N123,N,K3S,IE3,J3T,M3T,I2X,K2S,I2Y,I1X,K1S,I1Y,IS,INDP,IND,IE2,IE1,I2Z,I1Z,J2T,J1T,&
		   MM2TA,MM2TB,K3,KA3Q,KA3,MM2A,MM2B,MM2,M2,I2S,M1,I1S,MM2T,M2T,M1T,K1,KA1Q,KA1,KAS,KASKA1,&
		   K2,KA2Q,KA2,KASKA2,KA1KA2,KA23Q,KA123,KA3P,KA2P,KA23P,KA1P,KA123P,KR0,KR0123,KR0IND,&
		   KA12P,KR012,KASP,KASKA3,KR012S
REAL(KIND=8) ::	DS,DC
TYPE(su3irrep) :: su3irrep1,su3irrep2,su3irrep3

!************************************************************************** BEGINNING OF A BLOCK COMMENTED OUT BY J.H. BECAUSE OF ITS IRRELEVANCE FOR AUTOMATIC ARRAYS
! DIMENSION CHECKS (LSU,6-81)-START
! icount = 0
!IF(N1>9)THEN
! WRITE(*,FMT='(37H *****CWU3R3 DIMENSION OVERFLOW:  N1=,I10)')N1
! STOP
!END IF
!IF(N2>9)THEN
! WRITE(*,FMT='(37H *****CWU3R3 DIMENSION OVERFLOW:  N2=,I10)')N2
! STOP
!END IF
!IF(N3>9)THEN
! WRITE(*,FMT='(37H *****CWU3R3 DIMENSION OVERFLOW:  N3=,I10)')N3
! STOP
!END IF
! DIMENSION CHECKS (LSU,6-81)--STOP
!************************************************************************** END OF THE BLOCK COMMENTED OUT BY J.H.
!N01=KR0MAX*KA1MAX ! N0,N1 used to be here instead of KR0MAX,KA1MAX
!N012=N01*KA2MAX ! N2 used to be here instead of KA2MAX
!N0123=N012*KA3MAX ! N3 used to be here instead of KA3MAX
!N12=KA1MAX*KA2MAX ! N2 used to be here instead of KA2MAX
!N123=N12*KA3MAX ! N3 used to be here instead of KA3MAX

DO N=1,N0123
 DWU3R3(N)=0.D0
END DO
IF(KR0MAX==0)RETURN
IF(I3==1)THEN
 K3S=KSTART(su3irrep3%lambda,su3irrep3%mu,L3)-2
 IE3=-(su3irrep3%lambda+2*su3irrep3%mu)
 J3T=su3irrep3%lambda
 M3T=su3irrep3%lambda
ELSE
 K3S=KSTART(su3irrep3%mu,su3irrep3%lambda,L3)-2
 IE3=2*su3irrep3%lambda+su3irrep3%mu
 J3T=su3irrep3%mu
 M3T=-su3irrep3%mu
END IF
IF(I3/=J3)M3T=-M3T
I2X=su3irrep2%lambda+su3irrep2%mu+L2
IF(I2==1)THEN
 K2S=KSTART(su3irrep2%lambda,su3irrep2%mu,L2)-2
 IF(I2==J2)I2X=I2X-su3irrep2%lambda
 I2Y=0
ELSE
 K2S=KSTART(su3irrep2%mu,su3irrep2%lambda,L2)-2
 IF(I2==J2)I2X=I2X-su3irrep2%mu
 I2Y=6*(su3irrep2%lambda+su3irrep2%mu)
END IF
I1X=su3irrep1%lambda+su3irrep1%mu+L1
IF(I1==1)THEN
 K1S=KSTART(su3irrep1%lambda,su3irrep1%mu,L1)-2
 IF(I1==J1)I1X=I1X-su3irrep1%lambda
 I1Y=0
ELSE
 K1S=KSTART(su3irrep1%mu,su3irrep1%lambda,L1)-2
 IF(I1==J1)I1X=I1X-su3irrep1%mu
 I1Y=6*(su3irrep1%lambda+su3irrep1%mu)
END IF
IS=I2X
I2X=1
IF(2*(IS/2)==IS)I2X=0
IS=I1X
I1X=1
IF(2*(IS/2)==IS)I1X=0
I2Y=2*su3irrep2%lambda+4*su3irrep2%mu+6*(L2+K2S)-I2Y
I1Y=2*su3irrep1%lambda+4*su3irrep1%mu+6*(L1+K1S)-I1Y
INDP=-KR0MAX
DO IND=1,INDMAX
 INDP=INDP+KR0MAX
 IF(J1TA(IND)<0)CYCLE
 IE2=IEA(IND)
 IE1=IE3-IE2
 I2Z=(I2Y-IE2)/3
 I1Z=(I1Y-IE1)/3
 J2T=J2TA(IND)
 J1T=J1TA(IND)
 MM2TA=J1T+M3T
 IF(J2T<MM2TA)MM2TA=J2T
 MM2TA=MM2TA+1
 MM2TB=J1T-M3T
 IF(J2T<MM2TB)MM2TB=J2T
 MM2TB=MM2TA+MM2TB
 DO N=1,N123
  DS1A(N)=0.D0
 END DO
 K3=K3S
 KA3Q=-N12
 DO KA3=1,KA3MAX
  KA3Q=KA3Q+N12
  K3=K3+2
  MM2A=L1+K3
  IF(L2<MM2A)MM2A=L2
  MM2A=MM2A+1
  MM2B=L1-K3
  IF(L2<MM2B)MM2B=L2
  MM2B=MM2A+MM2B
  DO MM2=1,MM2B
   M2=MM2A-MM2
   IF(ABS(M2)>J2T)CYCLE
   I2S=J2T+M2
   IF(2*(I2S/2)/=I2S)CYCLE
   M1=K3-M2
   IF(ABS(M1)>J1T)CYCLE
   I1S=J1T+M1
   IF(2*(I1S/2)/=I1S)CYCLE
   DO N=1,N12
    DS2A(N)=0.D0
   END DO
   DO MM2T=1,MM2TB,2
    M2T=MM2TA-MM2T
    IF((M2T==0).AND.(4*(I2S/4)/=I2S))CYCLE
    IF(M2==0)THEN
     IS=I2Z+M2T
     IF(4*(IS/4)/=IS)CYCLE
     IS=(J2T+M2T)/2
     IF(2*(IS/2)/=IS)CYCLE
    END IF
    M1T=M3T-M2T
    IF((M1T==0).AND.(4*(I1S/4)/=I1S))CYCLE
    IF(M1==0)THEN
     IS=I1Z+M1T
     IF(4*(IS/4)/=IS)CYCLE
     IS=(J1T+M1T)/2
     IF(2*(IS/2)/=IS)CYCLE
	END IF
    K1=K1S
    KA1Q=-KA1MAX ! N1 used to be here instead of KA1MAX
    DO KA1=1,KA1MAX
     KA1Q=KA1Q+KA1MAX ! N1 used to be here instead of KA1MAX
     K1=K1+2
     IF(K1==0)THEN
      DT1A(KA1)=0.D0
      IF(I1X/=0)CYCLE
	 END IF
     DT1A(KA1)=transformation_coeff(I1,J1,su3irrep1,IE1,J1T,M1T,K1,L1,M1)
!     icount = icount + 1
!W    write(6,9000) ie1,j1t,m1t,k1,l1,m1,icount, 1, dt1a(ka1)
!9000 format(' Calling for DTU3R3 <',3i4,'|',3i4,'>', i10,'(',i1,')',&
!        f12.5 )
     DS=0.D0
     DO KAS=1,KA1
      KASKA1=KAS+KA1Q
      DS=DS+DON1(KASKA1)*DT1A(KAS)
     END DO
     DT1A(KA1)=DS
    END DO
    K2=K2S
    KA2Q=-KA2MAX ! N2 used to be here instead of KA2MAX
    DO KA2=1,KA2MAX
     KA2Q=KA2Q+KA2MAX ! N2 used to be here instead of KA2MAX
     K2=K2+2
     IF(K2==0)THEN
      DT2A(KA2)=0.D0
      IF(I2X/=0)CYCLE
	 END IF
     DT2A(KA2)=transformation_coeff(I2,J2,su3irrep2,IE2,J2T,M2T,K2,L2,M2)
!     icount = icount + 1
!     write(6,9000) ie2,j2t,m2t,k2,l2,m2,icount, 2, dt2a(ka2)
     DS=0.D0
     DO KAS=1,KA2
      KASKA2=KAS+KA2Q
      DS=DS+DON2(KASKA2)*DT2A(KAS)
     END DO
     DT2A(KA2)=DS
    END DO
    DC=clebsch_gordan(J1T,J2T,J3T,-M1T,-M2T,-M3T)
    KA2Q=-KA1MAX ! N1 used to be here instead of KA1MAX
    DO KA2=1,KA2MAX
     KA2Q=KA2Q+KA1MAX ! N1 used to be here instead of KA1MAX
     DO KA1=1,KA1MAX
      KA1KA2=KA1+KA2Q
      DS2A(KA1KA2)=DS2A(KA1KA2)+DC*DT1A(KA1)*DT2A(KA2)
     END DO
	END DO
   END DO
   DC=clebsch_gordan(2*L1,2*L2,2*L3,2*M1,2*M2,2*K3)
   KA2Q=-KA1MAX ! N1 used to be here instead of KA1MAX
   DO KA2=1,KA2MAX
    KA2Q=KA2Q+KA1MAX ! N1 used to be here instead of KA1MAX
    KA23Q=KA2Q+KA3Q
    DO KA1=1,KA1MAX
     KA1KA2=KA1+KA2Q
     KA123=KA1+KA23Q
     DS1A(KA123)=DS1A(KA123)+DC*DS2A(KA1KA2)
    END DO
   END DO
  END DO
 END DO
 KA3P=-N012
 KA3Q=-N12
 DO KA3=1,KA3MAX
  KA3Q=KA3Q+N12
  KA3P=KA3P+N012
  KA2P=-N01
  KA2Q=-KA1MAX ! N1 used to be here instead of KA1MAX
  DO KA2=1,KA2MAX
   KA2Q=KA2Q+KA1MAX ! N1 used to be here instead of KA1MAX
   KA2P=KA2P+N01
   KA23P=KA2P+KA3P
   KA23Q=KA2Q+KA3Q
   KA1P=-KR0MAX ! N0 used to be here instead of KR0MAX
   DO KA1=1,KA1MAX
    KA1P=KA1P+KR0MAX ! N0 used to be here instead of KR0MAX
    KA123P=KA1P+KA23P
    KA123=KA1+KA23Q
    DO KR0=1,KR0MAX
     KR0123=KR0+KA123P
     KR0IND=KR0+INDP
     DWU3R3(KR0123)=DWU3R3(KR0123)+DEWU3(KR0IND)*DS1A(KA123)
    END DO
   END DO
  END DO
 END DO
END DO
KA3P=-N012
KA3Q=-KA3MAX ! N3 used to be here instead of KA3MAX
DO KA3=1,KA3MAX
 KA3Q=KA3Q+KA3MAX ! N3 used to be here instead of KA3MAX
 KA3P=KA3P+N012
 KA2P=-N01
 DO KA2=1,KA2MAX
  KA2P=KA2P+N01
  KA23P=KA2P+KA3P
  KA1P=-KR0MAX ! N0 used to be here instead of KR0MAX
  DO KA1=1,KA1MAX
   KA1P=KA1P+KR0MAX ! N0 used to be here instead of KR0MAX
   KA123P=KA1P+KA23P
   KA12P=KA1P+KA2P
   DO KR0=1,KR0MAX
    KR0123=KR0+KA123P
    KR012=KR0+KA12P
    DS=0.D0
    KASP=-N012
    DO KAS=1,KA3
     KASP=KASP+N012
     KASKA3=KAS+KA3Q
     KR012S=KR012+KASP
     DS=DS+DON3(KASKA3)*DWU3R3(KR012S)
    END DO
    DWU3R3(KR0123)=DS
   END DO
  END DO
 END DO
END DO
!     write(6,*) ' icount ', icount

RETURN

CONTAINS
 FUNCTION KSTART(LAM,MU,L) RESULT(RKSTART)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: LAM,MU,L
  INTEGER :: RKSTART
  RKSTART=MOD(LAM,2)+2*(MAX(0,(L-MU)/2)+MOD(LAM+1,2)*MOD(ABS(L-MU),2))
 END FUNCTION KSTART
END SUBROUTINE wigner_physical