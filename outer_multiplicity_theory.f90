FUNCTION outer_multiplicity_theory(su3irrep1,su3irrep2,su3irrep3) RESULT(rhomax)
!-----------------------------------------------------------------------------------------
! Calculates the multiplicity (theory) of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! lambda1=su3irrep1%lambda, mu1=su3irrep1%mu etc.
! This is a new version of MULTHY
!-----------------------------------------------------------------------------------------
USE derived_types ! This module contains definitions of derived data types and operators
IMPLICIT NONE
TYPE(su3irrep), INTENT(IN) :: su3irrep1,su3irrep2,su3irrep3
INTEGER                    :: rhomax,IXMIN,I,IXDB3
INTEGER, DIMENSION(6)      :: IX
rhomax=0
IX(1)=su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda+2*(su3irrep1%mu+su3irrep2%mu-su3irrep3%mu)
IX(2)=su3irrep1%mu+su3irrep2%mu-su3irrep3%mu+2*(su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda)
IX(3)=2*su3irrep2%lambda+su3irrep2%mu+su3irrep1%mu-su3irrep1%lambda-su3irrep3%mu+su3irrep3%lambda
IX(4)=2*su3irrep2%mu+su3irrep2%lambda+su3irrep1%lambda-su3irrep1%mu-su3irrep3%lambda+su3irrep3%mu
IX(5)=su3irrep3%lambda+su3irrep2%mu-su3irrep1%lambda+2*(su3irrep3%mu+su3irrep2%lambda-su3irrep1%mu)                                       
IX(6)=su3irrep3%mu+su3irrep2%lambda-su3irrep1%mu+2*(su3irrep3%lambda+su3irrep2%mu-su3irrep1%lambda)                                       
IXMIN=1000                                                        
DO I=1,6
 IXDB3=IX(I)/3
 IF(3*IXDB3<IX(I))RETURN
 IF(IXDB3<IXMIN)IXMIN=IXDB3
END DO
IF(IXMIN<0)RETURN
rhomax=MIN(IXMIN,su3irrep2%lambda,su3irrep2%mu)+1
RETURN                                                            
END FUNCTION outer_multiplicity_theory
