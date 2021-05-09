!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! binomial_coeff.F90 -- binomial coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE binomial_coeff
!---------------------------------------------------------
! Binomial coefficients: (n choose k) = binom(n*(n+1)/2+k)
! n=0,...,upbound_binom
! k=0,...,n
!---------------------------------------------------------
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
USE mpmodule
#endif
IMPLICIT NONE
INTEGER :: upbound_binom
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: binom
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: binom_quad
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: binom_quad
TYPE(mp_real),ALLOCATABLE,DIMENSION(:) :: binom_mp
#endif
CONTAINS
  SUBROUTINE allocate_binom(upbound)
  ! This subroutine allocates arrays binom, binom_quad and binom_mp, calculates binomial coefficients and stores them in the arrays.
  ! This subroutine should be called at the begining of the main program (after "USE binomial_coeff").
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  USE mpmodule
  USE precision_level
#endif
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: upbound
  INTEGER :: n,k,ind
  upbound_binom=upbound
  ind=upbound_binom*(upbound_binom+1)/2+upbound_binom
#if defined(NDSU3LIB_DBL)
  ALLOCATE(binom(0:ind))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
  ALLOCATE(binom(0:ind),binom_quad(0:ind))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  ALLOCATE(binom(0:ind),binom_quad(0:ind),binom_mp(0:ind))
#endif
  binom(0)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
  binom_quad(0)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  binom_quad(0)=1.Q0
  binom_mp(0)=mpreal(1.D0,nwds)
#endif
  ind=1
  DO n=1,upbound_binom
    binom(ind)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    binom_quad(ind)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    binom_quad(ind)=1.Q0
    binom_mp(ind)=mpreal(1.D0,nwds)
#endif
    ind=ind+1
    DO k=1,n-1
      binom(ind)=binom(ind-n-1)+binom(ind-n)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
      binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
      binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
      binom_mp(ind)=binom_mp(ind-n-1)+binom_mp(ind-n)
#endif
      ind=ind+1
    END DO
    binom(ind)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    binom_quad(ind)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    binom_quad(ind)=1.Q0
    binom_mp(ind)=mpreal(1.D0,nwds)
#endif
    ind=ind+1
  END DO
  END SUBROUTINE allocate_binom
  SUBROUTINE deallocate_binom
  ! This subroutine dallocates arrays binom, binom_quad and binom_mp.
  IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
  DEALLOCATE(binom)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
  DEALLOCATE(binom,binom_quad)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  DEALLOCATE(binom,binom_quad,binom_mp)
#endif
  END SUBROUTINE deallocate_binom
  SUBROUTINE reallocate_binom(incr)
  ! This subroutine dallocates arrays binom, binom_quad and binom_mp, allocates them
  ! again with upbound_binom increased by incr and calculates the entries.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: incr
  CALL deallocate_binom
  CALL allocate_binom(upbound_binom+incr)
  END SUBROUTINE reallocate_binom
END MODULE binomial_coeff
