SUBROUTINE wigner_physical_mid(su3irrep1,su3irrep2,su3irrep3,L1,L2,L3,K0MAX,K1MAX,K2MAX,K3MAX,DVALU,NX,NCE2,N01,N012,N0123)
!-----------------------------------------------------------------------
! Calculates reduced SU(3)-SO(3) Wigner coefficients <(lam1,mu1)kap1,L1;(lam2,mu2)kap2,L2||(lam3,mu3)kap3,L3>_rho
!
! INPUT ARGUMENTS: su3irrep1: su3irrep1%lambda=lam1, su3irrep1%mu=mu1
!                  su3irrep2: su3irrep2%lambda=lam2, su3irrep2%mu=mu2
!                  su3irrep3: su3irrep3%lambda=lam3, su3irrep3%mu=mu3
!                  L1,L2,L3
!                  K0MAX = multiplicity of SU(3) coupling (lam1,mu1)x(lam2,mu2)->(lam3,mu3), i.e. the maximal value of rho
!                  K1MAX = the number of occurences of L1 within SU(3) irrep (lam1,mu1), i.e. the maximal value of kap1
!                  K2MAX = the number of occurences of L2 within SU(3) irrep (lam2,mu2), i.e. the maximal value of kap2
!                  K3MAX = the number of occurences of L3 within SU(3) irrep (lam3,mu3), i.e. the maximal value of kap3
!                  NX=lam2+mu2+1
!                  NCE2=NX*(NX+1)*(NX+2)/6
!                  N01=K0MAX*K1MAX
!                  N012=N01*K2MAX
!                  N0123=N012*K3MAX
!
! OUTPUT ARGUMENT: DVALU: DVALU(rho,kap1,kap2,kap3)=<(lam1,mu1)kap1,L1;(lam2,mu2)kap2,L2||(lam3,mu3)kap3,L3>_rho
!-----------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
!INTEGER, EXTERNAL :: outer_multiplicity
INTEGER :: L1,L2,L3,K0MAX,K1MAX,K2MAX,K3MAX,I1,J1,I2,J2,I3,J3,NEC,INDMAX,NX,NCE2,N01,N012,N0123,N12,N123,K0,K1,K2,K3,KOVER!,KIMAX1
INTEGER, DIMENSION(NCE2) :: J1TA,J2TA,IEA
REAL(KIND=8), DIMENSION(K0MAX*NCE2) :: DEWU3
REAL(KIND=8), DIMENSION(K1MAX*K1MAX) :: DON1
REAL(KIND=8), DIMENSION(K2MAX*K2MAX) :: DON2
REAL(KIND=8), DIMENSION(K3MAX*K3MAX) :: DON3
REAL(KIND=8), DIMENSION(N0123) :: DWU3R3
REAL(KIND=8), DIMENSION(K0MAX,K1MAX,K2MAX,K3MAX) :: DVALU
TYPE(su3irrep) :: su3irrep1,su3irrep2,su3irrep3

! CHECK IF THIS IS AN ALLOWED CASE
!K0MAX=outer_multiplicity(su3irrep1,su3irrep2,su3irrep3)
!IF(K0MAX==0)THEN
! PRINT*,"wigner_so3: THE COUPLING IS NOT ALLOWED!"
! RETURN
!END IF

! CHECK ANGULAR MOMENTUM CONTENT
!K1MAX=MULT(su3irrep1,L1)
!IF(K1MAX==0)THEN
! PRINT*,"wigner_so3: THE FIRST SU(3) IRREP DOES NOT CONTAIN THE VALUE OF L1!"
! RETURN
!END IF
!K2MAX=MULT(su3irrep2,L2)
!IF(K2MAX==0)THEN
! PRINT*,"wigner_so3: THE SECOND SU(3) IRREP DOES NOT CONTAIN THE VALUE OF L2!"
! RETURN
!END IF
!K3MAX=MULT(su3irrep3,L3)
!IF(K3MAX==0)THEN
! PRINT*,"wigner_so3: THE THIRD SU(3) IRREP DOES NOT CONTAIN THE VALUE OF L3!"
! RETURN
!END IF

! CHECK ANGULAR MOMENTUM COUPLING
!IF(L1+L2<L3.OR.ABS(L1-L2)>L3)THEN
! PRINT*,"wigner_so3: THE ANGULAR MOMENTUM COUPLING IS NOT ALLOWED!"
! RETURN
!END IF

IF(su3irrep1%lambda>=su3irrep1%mu)THEN
 I1=0
 J1=1
ELSE
 I1=1
 J1=0
END IF
IF(su3irrep2%lambda>=su3irrep2%mu)THEN
 I2=0
 J2=1
ELSE
 I2=1
 J2=0
END IF
IF(su3irrep3%lambda>=su3irrep3%mu)THEN
 I3=0
 J3=1
ELSE
 I3=1
 J3=0
END IF

N12=K1MAX*K2MAX
N123=N12*K3MAX

CALL wigner_canonical_extremal(su3irrep1,su3irrep2,su3irrep3,I3,NEC,K0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,NX)

CALL orthonormalization_matrix(I1,J1,su3irrep1,K1MAX,L1,DON1)
CALL orthonormalization_matrix(I2,J2,su3irrep2,K2MAX,L2,DON2)
CALL orthonormalization_matrix(I3,J3,su3irrep3,K3MAX,L3,DON3)

CALL wigner_physical(I1,J1,su3irrep1,K1MAX,L1,DON1,I2,J2,su3irrep2,K2MAX,L2,DON2,&
                     I3,J3,su3irrep3,K3MAX,L3,DON3,K0MAX,INDMAX,DEWU3,J1TA,J2TA,IEA,DWU3R3,N01,N012,N0123,N12,N123)

KOVER=0
DO K3=1,K3MAX
 DO K2=1,K2MAX
  DO K1=1,K1MAX
!   MINNUM=K0MAX*(K1-1)+K0MAX*K1MAX*(K2-1)+K0MAX*K1MAX*K2MAX*(K3-1)+1
!   MAXNUM=MINNUM+K0MAX-1
!   K0=0
!   DO I=MINNUM,MAXNUM
   DO K0=1,K0MAX
    KOVER=KOVER+1
!	K0=K0+1
!    DVALU(K0,K1,K2,K3)=DWU3R3(I)
    DVALU(K0,K1,K2,K3)=DWU3R3(KOVER)
   END DO
  END DO
 END DO
END DO

RETURN

!CONTAINS
! FUNCTION MULT(su3irrepx,L) RESULT(RMULT)
!  IMPLICIT NONE
!  INTEGER, INTENT(IN) :: L
!  TYPE(su3irrep), INTENT(IN) :: su3irrepx
!  INTEGER :: RMULT
!  RMULT=MAX(0,(su3irrepx%lambda+su3irrepx%mu+2-L)/2)-MAX(0,(su3irrepx%lambda+1-L)/2)-MAX(0,(su3irrepx%mu+1-L)/2)
! END FUNCTION MULT
END SUBROUTINE wigner_physical_mid