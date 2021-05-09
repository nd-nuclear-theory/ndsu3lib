!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ndsu3lib_free.F90 -- deallocation of global arrays
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ndsu3lib_free(wso3)
!-------------------------------------------------------------------------------
! This subroutine can be called by the main program once SU(3) Wigner or
! recoupling coefficients are not going to be calculated anymore to free memory.
!
! Input argument: wso3
!
! wso3 should be .TRUE. if ndsu3lib_init was called with the first argument
! being .TRUE.
!-------------------------------------------------------------------------------
USE binomial_coeff
USE I_S_module
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
USE fwigxjpf
#endif
IMPLICIT NONE
LOGICAL,INTENT(IN) :: wso3

IF(wso3)THEN
  CALL deallocate_I
  CALL deallocate_S
END IF
CALL deallocate_binom
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
CALL fwig_temp_free();
CALL fwig_table_free();
#endif

END SUBROUTINE ndsu3lib_free
