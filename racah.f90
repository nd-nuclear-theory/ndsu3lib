SUBROUTINE racah(su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,&
                 su3irrep23,KR0A,KR0B,KR0C,KR0D,RHS,INFO,KIMAX2,&
                 X1,X2,X3,NCE2A,NCE2B,NCE2D,NABC)
!--------------------------------------------------------------------------------------
! SUBROUTINE TO CALCULATE SU(3) RACAH COEFFICIENTS
! U((LAM1,MU1)(LAM2,MU2)(LAM,MU)(LAM3,MU3);(LAM12,MU12)RHOA,RHOB(LAM23,MU23)RHOC,RHOD)
! BY SOLVING THE SET OF EQUATIONS (22) USING MKL LAPACK SUBROUTINE dgesv
!
! REFERENCE: DRAAYER, AKIYAMA, J. MATH. PHYS., VOL. 14, NO. 12, 1973
!
! INPUT ARGUMENTS: su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,su3irrep23,KR0A,
!                  KR0B,KR0C,KR0D,KIMAX2,X1,X2,X3,NCE2A,NCE2B,NCE2D,NABC
! OUTPUT ARGUMENTS: RHS,INFO
!
! su3irrep1%lambda=LAM1
! su3irrep1%mu=MU1
! su3irrep2%lambda=LAM2
! su3irrep2%mu=MU2
! su3irrepx%lambda=LAM
! su3irrepx%mu=MU
! su3irrep3%lambda=LAM3
! su3irrep3%mu=MU3
! su3irrep12%lambda=LAM12
! su3irrep12%mu=MU12
! su3irrep23%lambda=LAM23
! su3irrep23%mu=MU23
! KR0A=MULTIPLICITY OF COUPLING (LAM1,MU1)x(LAM2,MU2)->(LAM12,MU12)
! KR0B=MULTIPLICITY OF COUPLING (LAM12,MU12)x(LAM3,MU3)->(LAM,MU)
! KR0C=MULTIPLICITY OF COUPLING (LAM2,MU2)x(LAM3,MU3)->(LAM23,MU23)
! KR0D=MULTIPLICITY OF COUPLING (LAM1,MU1)x(LAM23,MU23)->(LAM,MU)
! KIMAX2=KR0C*DIM(LAM3,MU3)
! X1=LAM3+MU3+1
! X2=LAM2+MU2+1
! X3=LAM23+MU23+1
! NCE2A=X2*(X2+1)*(X2+2)/6
! NCE2B=X1*(X1+1)*(X1+2)/6
! NCE2D=X3*(X3+1)*(X3+2)/6
! NABC=KR0A*KR0B*KR0C
! RHS(RHOD,N)=U((LAM1,MU1)(LAM2,MU2)(LAM,MU)(LAM3,MU3);(LAM12,MU12)RHOA,RHOB(LAM23,MU23)RHOC,RHOD)
!  WHERE N=RHOA+KR0A*(RHOB-1)+KR0A*KR0B*(RHOC-1)
! INFO=0 IFF dgesv RAN WITHOUT ERRORS
!--------------------------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
REAL(KIND=8), EXTERNAL :: su2racah
INTEGER :: KR0A,KR0B,KR0C,KR0D,INFO,NABC,NECA,NECB,NECC,NECD,INDMAX,I1,I2,I3,I4,I5,I6,&
           JTDMIN,JTDMAX,IE12,IE23,KA,KB,KC,KD,IND,IES,J,N,KR0IND,MINNUM,MAXNUM,J23T,&
           IESMAX,IE3MAX,J12TD,JJ12TA,JJ12TB,IS,J2TD,J2T,J3TD,J3T,J3S,J3SB,J3SQ,IESJ3S,&
           INDA,INDB,INDC,J2S,J2SB,IAQ,IBQ,ICQ,J12T,JJ12T,KAIA,KBIB,KCIC,KIMAX2,&
           X1,X2,X3,NCE2A,NCE2B,NCE2D
REAL(KIND=8) :: D1,D2,DC
TYPE(su3irrep) :: su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,su3irrep23,su3irrep2s

!REAL(KIND=8), DIMENSION(27090) :: DWU3
!REAL(KIND=8), DIMENSION(39732) :: DEWU3A,DEWU3B,DEWU3C,DEWU3D
!REAL(KIND=8), DIMENSION(9,9) :: MATRIX
!REAL(KIND=8), DIMENSION(9,729) :: RHS
!INTEGER, DIMENSION(13244) :: JXTA,JYTD,IE,JXTB,JXTD
!INTEGER, DIMENSION(1764) :: J2SMAX,J2TMAX,INDMAT
!INTEGER, DIMENSION(42) :: J3SMAX,J3TMAX
!INTEGER, DIMENSION(9) :: IPIV

REAL(KIND=8), DIMENSION(KIMAX2) :: DWU3
REAL(KIND=8), DIMENSION(KR0A*NCE2A) :: DEWU3A
REAL(KIND=8), DIMENSION(KR0B*NCE2B) :: DEWU3B
REAL(KIND=8), DIMENSION(KR0C*NCE2B) :: DEWU3C
REAL(KIND=8), DIMENSION(KR0D*NCE2D) :: DEWU3D
REAL(KIND=8), DIMENSION(KR0D,KR0D) :: MATRIX
REAL(KIND=8), DIMENSION(KR0D,NABC) :: RHS
INTEGER, DIMENSION(NCE2A) :: JXTA
INTEGER, DIMENSION(MAX(NCE2A,NCE2B,NCE2D)) :: JYTD,IE
INTEGER, DIMENSION(NCE2B) :: JXTB
INTEGER, DIMENSION(NCE2D) :: JXTD
INTEGER, DIMENSION(X1*X1) :: J2SMAX,J2TMAX,INDMAT
INTEGER, DIMENSION(X1) :: J3SMAX,J3TMAX
INTEGER, DIMENSION(KR0D) :: IPIV

!      EQUIVALENCE(JXTB(1),JXTC(1),JXTD(1))                              
!      EQUIVALENCE(JYTB(1),JYTC(1),JYTD(1))                              
!      EQUIVALENCE(IEB(1),IEC(1),IED(1))                                 
                                       
!      COMMON/BKSAVE/DEWU3A,DEWU3B,DEWU3C,DEWU3D
!!$omp threadprivate(/BKSAVE/)
!$omp threadprivate(DEWU3A,DEWU3B,DEWU3C,DEWU3D)                       
                                                
su3irrep2s%lambda=su3irrep2%mu
su3irrep2s%mu=su3irrep2%lambda

CALL wigner_canonical_extremal(su3irrep2,su3irrep3,su3irrep23,1,NECC,KR0C,INDMAX,DEWU3C,JXTB,JYTD,IE,X1)
CALL wigner_canonical_extremal(su3irrep12,su3irrep3,su3irrepx,1,NECB,KR0B,INDMAX,DEWU3B,JXTB,JYTD,IE,X1)
CALL wigner_canonical_extremal(su3irrep12,su3irrep2s,su3irrep1,1,NECA,KR0A,INDMAX,DEWU3A,JXTA,JYTD,IE,X2)
CALL wigner_canonical_extremal(su3irrep1,su3irrep23,su3irrepx,1,NECD,KR0D,INDMAX,DEWU3D,JXTD,JYTD,IE,X3)

I1=su3irrepx%lambda+2*su3irrepx%mu                                                       
I2=su3irrep12%lambda+2*su3irrep12%mu                                                   
I3=4*su3irrep12%lambda+2*su3irrep12%mu                                                 
I4=2*I2                                                           
I5=2*(su3irrep12%lambda-su3irrep12%mu)                                                 
JTDMIN=MAX(0,NECA-su3irrep2%lambda-su3irrep2%mu,NECB-su3irrep3%lambda-su3irrep3%mu)                        
JTDMAX=MIN(NECA,NECB,su3irrep12%lambda+su3irrep12%mu)                                 
D1=DFLOAT((su3irrep1%lambda+1)*IDM(su3irrep12))/DFLOAT(IDM(su3irrep1))         
IE23=su3irrep1%lambda+2*su3irrep1%mu-I1                                                

KD=0
DO IND=INDMAX,1,-1
 IF(JXTD(IND)<0)CYCLE
 MINNUM=1+KR0D*(IND-1)
 MAXNUM=KR0D+MINNUM-1
 KD=KD+1
 J=0
 DO KR0IND=MINNUM,MAXNUM
  J=J+1
  MATRIX(KD,J)=DEWU3D(KR0IND)
 END DO

 J23T=JYTD(IND)
 D2=DSQRT(D1*DFLOAT(J23T+1))                                       
 CALL wigner_canonical(su3irrep2,su3irrep3,su3irrep23,IE23,J23T,NECC,DEWU3C,KR0C,INDMAX,DWU3,&
                       J2SMAX,J2TMAX,J3SMAX,J3TMAX,IESMAX,IE3MAX,INDMAT,X1,KIMAX2)
 IE12=3*IESMAX-IE3MAX-I1                                           
 J12TD=(IE12+I2)/3                                                 
 DO IES=1,IESMAX                                                
  IE12=IE12-3                                                       
  J12TD=J12TD-1                                                     
  IF(J12TD<JTDMIN)CYCLE                                       
  IF(J12TD>JTDMAX)CYCLE                                     
  JJ12TA=I3+IE12                                                    
  IS=I4-IE12                                                        
  IF(IS<JJ12TA)JJ12TA=IS                                         
  JJ12TA=JJ12TA/3+1                                                 
  IS=(I5-IE12)/3                                                    
  JJ12TB=JJ12TA-ABS(IS)                                            
  I6=2*su3irrep1%lambda+IS                                                      
  J2TD=NECA-J12TD                                                   
  J3TD=NECB-J12TD                                                   
  J3T=J3TMAX(IES)+2                                                 
  J3SB=J3SMAX(IES)                                                  
  J3SQ=-X1                                                         
  DO J3S=1,J3SB                                                  
   J3SQ=J3SQ+X1                                                    
   IESJ3S=IES+J3SQ                                                   
   J3T=J3T-2                                                         
   J2T=J2TMAX(IESJ3S)+2                                              
   INDC=(INDMAT(IESJ3S)-J2T)/2                                       
   J2SB=J2SMAX(IESJ3S)                                               
   DO J2S=1,J2SB                                                  
    J2T=J2T-2                                                         
    ICQ=INDC*KR0C                                                     
    INDC=INDC+1 
    DO JJ12T=1,JJ12TB,2                                            
     J12T=JJ12TA-JJ12T                                                 
     INDA=INDEX(J12TD,su3irrep12%lambda,J12T,J2TD,su3irrep2%mu,J2T)                         
     IF(JXTA(INDA)<0)CYCLE                                      
     INDB=INDEX(J12TD,su3irrep12%lambda,J12T,J3TD,su3irrep3%lambda,J3T)                        
     IF(JXTB(INDB)<0)CYCLE                                        
     DC=D2*su2racah(su3irrep1%lambda,J2T,su3irrepx%lambda,J3T,J12T,J23T)                            
     IS=J12T+I6                                                        
     IF(4*(IS/4)/=IS)DC=-DC                                          
     IAQ=(INDA-1)*KR0A                                                 
     IBQ=(INDB-1)*KR0B                                                 

     N=0
     DO KC=1,KR0C                                                                  
      KCIC=KC+ICQ                                                       
      DO KB=1,KR0B                                                                 
       KBIB=KB+IBQ                                                       
       DO KA=1,KR0A                                                    
        KAIA=KA+IAQ                                                       
        N=N+1
        RHS(KD,N)=RHS(KD,N)+DC*DEWU3A(KAIA)*DEWU3B(KBIB)*DWU3(KCIC)
       END DO
      END DO
     END DO

    END DO
   END DO
  END DO
 END DO
 IF(KD==KR0D)EXIT      
END DO

IF(KR0D>1)THEN
 CALL dgesv(KR0D,NABC,MATRIX,KR0D,IPIV,RHS,KR0D,INFO)
ELSE
 DO N=1,NABC
  RHS(1,N)=RHS(1,N)/MATRIX(1,1)
 END DO
END IF

RETURN

CONTAINS
 FUNCTION IDM(su3irrepx) RESULT(IIDM)
  IMPLICIT NONE
  TYPE(su3irrep),INTENT (IN) :: su3irrepx
  INTEGER :: IIDM
  IIDM=(su3irrepx%lambda+1)*(su3irrepx%mu+1)*(su3irrepx%lambda+su3irrepx%mu+2)/2
 END FUNCTION IDM

 FUNCTION INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T) RESULT(IINDEX)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: J1TD,LAM1,J1T,J2TD,LAM2,J2T
  INTEGER :: IINDEX
  IINDEX=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2
 END FUNCTION INDEX     
END SUBROUTINE racah
