!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! binomial_coeff.F90 -- precalculation of binomial coefficients
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
! n=0,...,400
! k=0,...,n
!---------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),DIMENSION(0:80600) :: binom!,binom_inv
#if defined(SU3DBL)
REAL(KIND=8),DIMENSION(0:80600) :: binom_quad!,binom_inv_quad
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16),DIMENSION(0:80600) :: binom_quad!,binom_inv_quad
#endif
!REAL(KIND=8),DIMENSION(0:100) :: inv_pow_4 ! inv_pow_4(p)=1/4^p
CONTAINS
  SUBROUTINE calculate_binomial_coeff
  ! This subroutine calculates binomial coefficients and stores them in array binom.
  ! This subroutine should be called at the begining of the main program (after "USE binomial_coeff").
  IMPLICIT NONE
  INTEGER :: n,k,ind
  binom(0)=1.D0
!  binom_inv(0)=1.D0
#if defined(SU3DBL)
  binom_quad(0)=1.D0
!  binom_inv_quad(0)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
   binom_quad(0)=1.Q0
!   binom_inv_quad(0)=1.Q0
#endif
  ind=1
  DO n=1,400
    binom(ind)=1.D0
!    binom_inv(ind)=1.D0
#if defined(SU3DBL)
    binom_quad(ind)=1.D0
!    binom_inv_quad(ind)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
    binom_quad(ind)=1.Q0
!    binom_inv_quad(ind)=1.Q0
#endif
    ind=ind+1
    DO k=1,n-1
      binom(ind)=binom(ind-n-1)+binom(ind-n)
!      binom_inv(ind)=1.D0/binom(ind)
      binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
      ind=ind+1
    END DO
    binom(ind)=1.D0
!    binom_inv(ind)=1.D0
#if defined(SU3DBL)
    binom_quad(ind)=1.D0
!    binom_inv_quad(ind)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
    binom_quad(ind)=1.Q0
!    binom_inv_quad(ind)=1.Q0
#endif
    ind=ind+1
  END DO
!  inv_pow_4(0)=1.D0
!  DO nn=1,100
!    inv_pow_4(nn)=inv_pow_4(nn-1)/4.D0
!  END DO
  END SUBROUTINE calculate_binomial_coeff
END MODULE binomial_coeff
