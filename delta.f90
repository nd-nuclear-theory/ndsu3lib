FUNCTION delta(J1T,J2T,J3T) RESULT(DELT)
!------------------------------------------------------
! DELTA FOR R3 ROUTINES--TRIANGLE RELATIONS CHECKED                 
!------------------------------------------------------
USE binomial_coeff_factorials
IMPLICIT NONE
INTEGER, INTENT(IN) :: J1T,J2T,J3T
REAL(KIND=8) :: DELT
INTEGER :: I1,I2,I3
DELT=12345.D0
I1=J1T+J2T-J3T
IF(BTEST(I1,0).OR.I1<0)RETURN
I2=J2T+J3T-J1T
IF(I2<0)RETURN
I3=J3T+J1T-J2T
IF(I3<0)RETURN
DELT=DLOGF(I1)+DLOGF(I2)+DLOGF(I3)-DLOGF(J1T+J2T+J3T+2)
RETURN
END FUNCTION delta