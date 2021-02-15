!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! orthonormalization_matrix.f90 -- orthonormalization matrix for SO(3) reduced basis states
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE orthonormalization_matrix(I,J,lambda,mu,L,kappamax,matrix)
!-------------------------------------------------------------------------------------------
! Calculates the bottom triangle, including the diagonal, of the orthonormalization matrix O
! for given lambda,mu,L using equations (6a),(6b),(6c),(27) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: I,J,lambda,mu,L,kappamax
! Output argument: matrix
!
! (I,J)=(1,1) => E=HW, (I,J)=(1,0) => E=HW', (I,J)=(0,0) => E=LW, (I,J)=(0,1) => E=LW'
! kappamax = the number of occurences of L in SU(3) irrep (lambda,mu)
! matrix(i,j) = O_ij
!
! Note: There is a typo in Eq.(6b): the square root should not be there!
!-------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: transformation_coeff
REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: matrix
REAL(KIND=8) :: sum
INTEGER,INTENT(IN) :: I,J,lambda,mu,L,kappamax
INTEGER :: ii,jj,k,epsilon,Ki,Kj,Kminm2,Lambda2,MLambda2

IF(I==1)THEN ! E=HW or HW'
  epsilon=-lambda-2*mu ! the HW epsilon
  Kminm2=Kmin(lambda,mu,L)-2 ! Kmin(lambda,mu,L) is the lowest K
  Lambda2=lambda ! Lambda2 is 2*Lambda
  MLambda2=lambda ! MLambda2 is 2*M_Lambda
ELSE ! E=LW or LW'
  epsilon=2*lambda+mu ! the LW epsilon
  Kminm2=Kmin(mu,lambda,L)-2 ! Kmin(mu,lambda,L) is the lowest K
  Lambda2=mu ! Lambda2 is 2*Lambda
  MLambda2=-mu ! MLambda2 is 2*M_Lambda
END IF
IF(I/=J)MLambda2=-MLambda2

Ki=Kminm2

! First, the upper triangle is calculated including the diagonal.
DO ii=1,kappamax ! loop over columns
  Ki=Ki+2 ! Ki is K_i in Eq.(6a),(6b)
  Kj=Kminm2
  DO jj=1,ii-1 ! loop over rows
    Kj=Kj+2 ! Kj is K_j in Eq.(6b)
    sum=0.D0
    DO k=1,jj-1
      sum=sum+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6b)
    END DO
!    matrix(jj,ii)=matrix(jj,jj)*DSQRT(transformation_coeff(I,J,lambda,mu,epsilon,Lambda2,MLambda2,Ki,L,Kj)-sum) ! Eq.(6b)
    matrix(jj,ii)=matrix(jj,jj)*(transformation_coeff(I,J,lambda,mu,epsilon,Lambda2,MLambda2,Ki,L,Kj)-sum) ! Eq.(6b)
  END DO
  sum=0.D0
  DO jj=1,ii-1
    sum=sum+matrix(jj,ii)*matrix(jj,ii) ! sum is the sum in Eq.(6a)
  END DO
  matrix(ii,ii)=1.D0/DSQRT(transformation_coeff(I,J,lambda,mu,epsilon,Lambda2,MLambda2,Ki,L,Ki)-sum) ! Eq.(6a)
END DO

! Now, the bottom triangle is calculated.
DO jj=1,kappamax-1 ! loop over columns
  DO ii=jj+1,kappamax ! loop over rows
    sum=0.D0
    DO k=jj,ii-1
      sum=sum+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6c)
    END DO
    matrix(ii,jj)=-matrix(ii,ii)*sum ! Eq.(6c)
  END DO
END DO

CONTAINS
  FUNCTION Kmin(lambda,mu,L) RESULT(res) ! This function calculates the lowest K for given lambda,mu,L as given by Eq.(4a) in the reference.
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda,mu,L
    INTEGER :: res
    res=MAX(0,L-mu)
    res=res+MOD(res+lambda,2) ! K and lambda must have the same parity.
    IF(res==0)res=res+MOD(L+mu,2)*2 ! If K=0, L and mu must have the same parity.
  END FUNCTION Kmin
END SUBROUTINE orthonormalization_matrix
