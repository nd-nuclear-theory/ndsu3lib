MODULE binomial_coeff
!---------------------------------------------------------
! Binomial coefficients: (n choose k) = binom(n*(n+1)/2+k)
! n=0,...,upbound_binom
! k=0,...,n
!---------------------------------------------------------
IMPLICIT NONE
INTEGER :: upbound_binom,size_binom
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: binom
#if defined(SU3DBL)
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: binom_quad
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: binom_quad
#endif
CONTAINS
  SUBROUTINE allocate_binom(upbound)
  ! This subroutine allocates arrays binom and binom_quad, calculates binomial coefficients and stores them in the arrays.
  ! This subroutine should be called at the begining of the main program (after "USE binomial_coeff").
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: upbound
  INTEGER :: n,k,ind
  upbound_binom=upbound
  size_binom=upbound_binom*(upbound_binom+1)/2+upbound_binom
  ALLOCATE(binom(0:size_binom),binom_quad(0:size_binom))
  binom(0)=1.D0
#if defined(SU3DBL)
  binom_quad(0)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  binom_quad(0)=1.Q0
#endif
  ind=1
  DO n=1,upbound_binom
    binom(ind)=1.D0
#if defined(SU3DBL)
    binom_quad(ind)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
    binom_quad(ind)=1.Q0
#endif
    ind=ind+1
    DO k=1,n-1
      binom(ind)=binom(ind-n-1)+binom(ind-n)
      binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
      ind=ind+1
    END DO
    binom(ind)=1.D0
#if defined(SU3DBL)
    binom_quad(ind)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
    binom_quad(ind)=1.Q0
#endif
    ind=ind+1
  END DO
  END SUBROUTINE allocate_binom
  SUBROUTINE deallocate_binom
  IMPLICIT NONE
  DEALLOCATE(binom,binom_quad)
  END SUBROUTINE deallocate_binom
  SUBROUTINE reallocate_binom(incr)
  ! This subroutine dallocates arrays binom and binom_quad, allocates them
  ! again with upbound_binom increased by incr and calculates the entries.
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: incr
  DEALLOCATE(binom,binom_quad)
  CALL allocate_binom(upbound_binom+incr)
  END SUBROUTINE reallocate_binom
END MODULE binomial_coeff
