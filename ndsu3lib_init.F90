!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ndsu3lib_init.F90 -- ndsu3lib initialization
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ndsu3lib_init(wso3,j2max)
!----------------------------------------------------------------------------
! ndsu3lib initialization subroutine
! This subroutine must be called by the main program before calling ndsu3lib
! subroutines for SU(3) Wigner or recoupling coefficients.
!
! Input arguments: wso3,j2max
!
! wso3 must be .TRUE. if SU(3)-SO(3) Wigner coefficients are going to be
! calculated.
! If WIGXJPF is not going to be utilized, j2max is not used. Otherwise j2max
! must be greater than or equal to two times the maximal angular momentum
! expected in ordinary Clebsch-Gordan or SU(2) recoupling coefficients.
! j2max should be at least the maximal expected value of lambda+mu if
! SU(3)-SO(3) Wigner coefficients are not going to be calculated. If
! SU(3)-SO(3) Wigner coefficients are going to be calculated, j2max should be
! at least two times the maximal expected value of lambda+mu. If this j2max
! is insufficient, WIGXJPF will terminate the program and display an error
! message.
!----------------------------------------------------------------------------
USE binomial_coeff
USE I_S_module
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
USE fwigxjpf
#endif
IMPLICIT NONE
LOGICAL,INTENT(IN) :: wso3
INTEGER,INTENT(IN) :: j2max

CALL allocate_binom(100)
IF(wso3)THEN
  CALL allocate_I(50)
  CALL allocate_S(50)
END IF
#if defined(NDSU3LIB_RACAH_WIGXJPF)
CALL fwig_table_init(j2max,6)
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
CALL fwig_table_init(j2max,3)
#endif
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
CALL fwig_temp_init(j2max)
#endif

END SUBROUTINE ndsu3lib_init
