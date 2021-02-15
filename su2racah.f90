!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! su2racah.f90 -- bindings for GSL 6j coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION su2racah(j1,j2,j,j3,j12,j23) RESULT(w)
!---------------------------------------------------------------------
! Calculates SU(2) Racah coefficient W(j1/2,j2/2,j/2,j3/2,j12/2,j23/2)
! using GSL function gsl_sf_coupling_6j calculating 6j symbol.
!---------------------------------------------------------------------
USE iso_c_binding
IMPLICIT NONE

INTERFACE
REAL(C_DOUBLE) FUNCTION gsl_sf_coupling_6j(l1,l2,l12,l3,l,l23) BIND(C)
USE iso_c_binding
INTEGER(C_INT), VALUE :: l1,l2,l12,l3,l,l23
END FUNCTION gsl_sf_coupling_6j
END INTERFACE

INTEGER(C_INT),INTENT(IN) :: j1,j2,j,j3,j12,j23
REAL(C_DOUBLE) :: w
INTEGER :: a

w=gsl_sf_coupling_6j(j1,j2,j12,j3,j,j23)
a=j1+j2+j3+j
IF((a/4)*4/=a)w=-w

RETURN
END FUNCTION su2racah
