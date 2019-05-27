FUNCTION outer_multiplicity(su3irrep1,su3irrep2,su3irrep3) RESULT(rhomax)
!-----------------------------------------------------------------------------------------
! Calculates the multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! lambda1=su3irrep1%lambda, mu1=su3irrep1%mu etc.
! This is a new version of MULTU3
!-----------------------------------------------------------------------------------------
USE derived_types  ! This module contains definitions of derived data types and operators
IMPLICIT NONE
TYPE(su3irrep), INTENT(IN) :: su3irrep1,su3irrep2,su3irrep3
INTEGER                    :: rhomax,NX,MX,L1,L2,L3,M1,M2,M3,MU,NU,MY,NY
rhomax=0                                                          
NX=su3irrep1%lambda+su3irrep2%lambda-su3irrep3%lambda-su3irrep1%mu-su3irrep2%mu+su3irrep3%mu
MX=NX/3            ! If NX is not a multiple of 3, MX is not equal to NX/3
IF(3*MX/=NX)RETURN ! If NX is not a multiple of 3, rhomax=0
IF(MX>=0)THEN
 L1=su3irrep1%lambda
 L2=su3irrep2%lambda
 L3=su3irrep3%lambda
 M1=su3irrep1%mu
 M2=su3irrep2%mu
 M3=su3irrep3%mu
ELSE
 L1=su3irrep1%mu
 L2=su3irrep2%mu
 L3=su3irrep3%mu
 M1=su3irrep1%lambda
 M2=su3irrep2%lambda
 M3=su3irrep3%lambda
 MX=-MX
END IF
NX=MX+M1+M2-M3                                                    
MU=MIN(L1-MX,M2)                                                          
IF(MU<0)RETURN                                                 
NU=MIN(L2-MX,M1)
IF(NU<0)RETURN                                                 
rhomax=MAX(MIN(NX,NU)-MAX(NX-MU,0)+1,0)
RETURN                                                            
END FUNCTION outer_multiplicity