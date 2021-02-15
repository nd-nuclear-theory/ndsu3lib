!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! outer_multiplicity.f90 -- multiplicity of SU(3) coupling
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda3,mu3) RESULT(rhomax)
!-------------------------------------------------------------------------------------
! Calculates multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! Reference: M.F.O'Reilly, J.Math.Phys. 23 (1982) 2022: Section 5, Proposition 7(a)
!-------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda3,mu3
INTEGER :: rhomax,C3,C,D
C3=lambda1-mu1+lambda2-mu2-lambda3+mu3 ! C3 equals 3 times the C in the reference
C=C3/3 ! C equals the C in the reference if 3 divides C3
IF((3*C==C3).AND.(C>=-MIN(mu2,mu1)).AND.(C<=MIN(lambda1,lambda2)))THEN
  D=C+mu1+mu2-mu3 ! D equals the D in the reference
  IF((D>=0).AND.(D<=MIN(mu1+mu2,lambda2+mu2,lambda1+mu1))&
   .AND.(D+C>=0).AND.(D+C<=MIN(lambda1+mu1,lambda2+lambda1,lambda2+mu2)))THEN
    rhomax=1+MIN(mu2,lambda1+mu1,D,lambda1-C)-MAX(0,D-mu1,D-lambda2,-C,D-C-mu1,D+C-lambda2)
  ELSE
    rhomax=0
  END IF
ELSE
  rhomax=0
END IF
RETURN                                                            
END FUNCTION outer_multiplicity
