MODULE binomial_coeff
!---------------------------------------------------------
! Binomial coefficients: (n choose k) = binom(n*(n+1)/2+k)
! n=0,...,200
! k=0,...,n
!---------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),DIMENSION(0:20300) :: binom!,binom_inv
!REAL(KIND=8),DIMENSION(0:100) :: inv_pow_4 ! inv_pow_4(p)=1/4^p
CONTAINS
  SUBROUTINE calculate_binomial_coeff
  ! This subroutine calculates binomial coefficients and stores them in array binom.
  ! This subroutine should be called at the begining of the main program (after "USE binomial_coeff").
  IMPLICIT NONE
  INTEGER :: n,k,ind
  binom(0)=1.D0
!  binom_inv(0)=1.D0
  ind=1
  DO n=1,200
    binom(ind)=1.D0
!    binom_inv(ind)=1.D0
    ind=ind+1
    DO k=1,n-1
      binom(ind)=binom(ind-n-1)+binom(ind-n)
!      binom_inv(ind)=1.D0/binom(ind)
      ind=ind+1
    END DO
    binom(ind)=1.D0
!    binom_inv(ind)=1.D0
    ind=ind+1
  END DO
!  inv_pow_4(0)=1.D0
!  DO nn=1,100
!    inv_pow_4(nn)=inv_pow_4(nn-1)/4.D0
!  END DO
  END SUBROUTINE calculate_binomial_coeff
END MODULE binomial_coeff
