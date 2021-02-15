!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wigner_physical_wrap.f90 -- inner wrapper for SU(3)-SO(3) coupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wigner_physical_wrap(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
!------------------------------------------------------------------------------------------------
! Calculates physical Wigner coefficients for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3.
!
! Input arguments: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3
!                  rhomax = multipicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
!                  kappa1max = number of occurences of L1 in (lambda1,mu1)
!                  kappa2max,kappa3max
!
! Output argument: wigner_phys(kappa1,kappa2,kappa3,rho) = the Wigner coefficient
!------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax
INTEGER :: I1,J1,I2,J2,I3,J3,numb
REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys ! The indeces are kappa1,kappa2,kappa3,rho, respectively.
!REAL(KIND=8),DIMENSION(0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax) :: wigner_can
!REAL(KIND=8),DIMENSION(kappa1max,kappa1max) :: matrix1
!REAL(KIND=8),DIMENSION(kappa2max,kappa2max) :: matrix2
!REAL(KIND=8),DIMENSION(kappa3max,kappa3max) :: matrix3
!INTEGER,DIMENSION((MAX(lambda1,mu1)+1)*((MAX(lambda2,mu2)+1)**2)) :: p1a,p2a,q2a
INTEGER :: pqdim
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_can
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix1,matrix2,matrix3
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1a,p2a,q2a

INTERFACE
  SUBROUTINE wigner_canonical_extremal(lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax,i2,wigner,p1a,p2a,q2a)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
    INTEGER,INTENT(OUT) :: i2
    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
  END SUBROUTINE wigner_canonical_extremal
  SUBROUTINE wigner_physical(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
                             I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: I1,J1,lambda1,mu1,L1,kappa1max,I2,J2,lambda2,mu2,L2,kappa2max,I3,lambda3,mu3,L3,kappa3max,rhomax,numb
    INTEGER,DIMENSION(:),INTENT(IN) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2,matrix3
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wigner_can
    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
  END SUBROUTINE wigner_physical
  SUBROUTINE orthonormalization_matrix(I,J,lambda,mu,L,kappamax,matrix)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: I,J,lambda,mu,L,kappamax
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: matrix
  END SUBROUTINE orthonormalization_matrix
END INTERFACE

pqdim=(MAX(lambda1,mu1)+1)*((MAX(lambda2,mu2)+1)**2)
ALLOCATE(wigner_can(0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax),&
         matrix1(kappa1max,kappa1max),matrix2(kappa2max,kappa2max),matrix3(kappa3max,kappa3max),&
         p1a(pqdim),p2a(pqdim),q2a(pqdim))

IF(lambda1<mu1)THEN
  I1=1
  J1=0
ELSE
  I1=0
  J1=1
END IF
IF(lambda2<mu2)THEN
  I2=1
  J2=0
ELSE
  I2=0
  J2=1
END IF
IF(lambda3<mu3)THEN
  I3=1
  J3=0
ELSE
  I3=0
  J3=1
END IF

CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda3,mu3,I3,rhomax,numb,wigner_can,p1a,p2a,q2a)

CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1)
CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2)
CALL orthonormalization_matrix(I3,J3,lambda3,mu3,L3,kappa3max,matrix3)

CALL wigner_physical(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
        I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)

DEALLOCATE(wigner_can,matrix1,matrix2,matrix3,p1a,p2a,q2a)

END SUBROUTINE wigner_physical_wrap
