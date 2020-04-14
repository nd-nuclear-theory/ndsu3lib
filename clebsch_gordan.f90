FUNCTION clebsch_gordan(j1,m1,j2,m2,j3,m3) RESULT(cg)
!----------------------------------------------------------------------------
! Calculates SU(2) Clebsch-Gordan coefficient <j1/2,m1/2,j2/2,m2/2|j3/2,m3/2>
! using GSL function gsl_sf_coupling_3j calculating 3j symbol.
!----------------------------------------------------------------------------
USE iso_c_binding
IMPLICIT NONE

INTERFACE
REAL(C_DOUBLE) FUNCTION gsl_sf_coupling_3j(l1,n1,l2,n2,l3,n3) BIND(C)
USE iso_c_binding
INTEGER(C_INT), VALUE :: l1,n1,l2,n2,l3,n3
END FUNCTION gsl_sf_coupling_3j
END INTERFACE

INTEGER(C_INT),INTENT(IN) :: j1,m1,j2,m2,j3,m3
REAL(C_DOUBLE) :: cg
INTEGER :: a

cg=gsl_sf_coupling_3j(j1,j2,j3,m1,m2,-m3)*DSQRT(DFLOAT(j3+1))
a=j1-j2+m3
IF((a/4)*4/=a)cg=-cg

RETURN
END FUNCTION clebsch_gordan
