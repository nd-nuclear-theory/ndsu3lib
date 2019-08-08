SUBROUTINE wigner_physical_wrap(su3irrep1,su3irrep2,su3irrep3,L1,L2,L3,K0MAX,K1MAX,K2MAX,K3MAX,DVALU)
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
!
! OUTPUT ARGUMENTS: DVALU: DVALU(rho,kap1,kap2,kap3)=<(lam1,mu1)kap1,L1;(lam2,mu2)kap2,L2||(lam3,mu3)kap3,L3>_rho
!-----------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
!INTEGER, EXTERNAL :: outer_multiplicity
INTEGER :: L1,L2,L3,K0MAX,K1MAX,K2MAX,K3MAX,NX,NCE2,N01,N012,N0123
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
IF(L1+L2<L3.OR.ABS(L1-L2)>L3)THEN
! PRINT*,"wigner_so3: THE ANGULAR MOMENTUM COUPLING IS NOT ALLOWED!"
 RETURN
END IF

NX=su3irrep2%lambda+su3irrep2%mu+1
NCE2=NX*(NX+1)*(NX+2)/6
N01=K0MAX*K1MAX
N012=N01*K2MAX
N0123=N012*K3MAX

CALL wigner_physical_mid(su3irrep1,su3irrep2,su3irrep3,L1,L2,L3,K0MAX,K1MAX,K2MAX,K3MAX,DVALU,NX,NCE2,N01,N012,N0123)

RETURN
END SUBROUTINE wigner_physical_wrap