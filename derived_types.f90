MODULE derived_types
!---------------------------------------------------------
! Contains definitions of derived data types and operators
!---------------------------------------------------------
IMPLICIT NONE
TYPE su3irrep
 INTEGER :: lambda,mu
END TYPE su3irrep

!! TYPE(su3irrep) :: su3irrep11=su3irrep(lambda=1,mu=1)

INTERFACE OPERATOR (+)
 PROCEDURE add_integer_to_lambda_and_mu
END INTERFACE OPERATOR (+)

CONTAINS
FUNCTION add_integer_to_lambda_and_mu(old_irrep,offset) RESULT(new_irrep)
IMPLICIT NONE
TYPE(su3irrep), INTENT(IN) :: old_irrep
INTEGER, INTENT(IN) :: offset
TYPE(su3irrep) :: new_irrep

new_irrep%lambda=old_irrep%lambda+offset
new_irrep%mu=old_irrep%mu+offset
END FUNCTION add_integer_to_lambda_and_mu

END MODULE derived_types