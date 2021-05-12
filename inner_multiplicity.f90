!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! inner_multiplicity.f90 -- number of occurences of L within SU(3) irrep
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION inner_multiplicity(lambda,mu,L) RESULT(kappamax)
!-------------------------------------------------------
! Inner multiplicity of L within SU(3) irrep (lambda,mu)
!-------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda,mu,L
INTEGER :: kappamax
kappamax=MAX(0,(lambda+mu+2-L)/2)-MAX(0,(lambda+1-L)/2)-MAX(0,(mu+1-L)/2)
END FUNCTION inner_multiplicity
