MODULE binomial_coeff
!-------------------------------------------------
! Binomial coefficients: (n choose k) = binom(n,k)
!-------------------------------------------------
IMPLICIT NONE
INTEGER :: nn,kk
REAL(KIND=8),DIMENSION(0:200,0:200) :: binom,binom_inv
!REAL(KIND=8),DIMENSION(0:100) :: inv_pow_4 ! inv_pow_4(p)=1/4^p
CONTAINS
  SUBROUTINE calculate_binomial_coeff
  ! This subroutine calculates binomial coefficients from Pascal's triangle and stores them in array binom.
  ! This subroutine should be called at the begining of the main program (after "USE binomial_coeff").
    binom(0:200,0)=1.D0
    binom_inv(0:200,0)=1.D0
    DO nn=1,200
      binom(nn,nn)=1.D0
      binom_inv(nn,nn)=1.D0
    END DO
    DO nn=2,200
      DO kk=1,nn-1
        binom(nn,kk)=binom(nn-1,kk-1)+binom(nn-1,kk)
        binom_inv(nn,kk)=1.D0/binom(nn,kk)
      END DO
    END DO
!    inv_pow_4(0)=1.D0
!    DO nn=1,100
!      inv_pow_4(nn)=inv_pow_4(nn-1)/4.D0
!    END DO
  END SUBROUTINE calculate_binomial_coeff
END MODULE binomial_coeff
