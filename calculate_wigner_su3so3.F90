!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculate_wigner_su3so3.F90 -- SU(3)-SO(3) reduced Wigner coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_wigner_su3so3(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
!-----------------------------------------------------------------------------------------------------------
! Calculates SU(3)-SO(3) reduced Wigner coefficients for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3.
!
! Input arguments: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3
!                  rhomax = multipicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
!                  kappa1max = number of occurences of L1 in (lambda1,mu1)
!                  kappa2max,kappa3max
!
! Output argument: wigner_phys(kappa1,kappa2,kappa3,rho) = the Wigner coefficient
!-----------------------------------------------------------------------------------------------------------
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
USE mpmodule
#endif
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax
INTEGER :: I1,J1,I2,J2,I3,J3,numb,indicator1,indicator2,sum1,sum2,pqdim
REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys ! The indeces are kappa1,kappa2,kappa3,rho, respectively.
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_can
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix1,matrix2,matrix3
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
TYPE(mp_real),ALLOCATABLE,DIMENSION(:,:) :: matrix1mp,matrix2mp
#endif
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1a,p2a,q2a

INTERFACE
  SUBROUTINE wigner_canonical_extremal(lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax,i2,wigner,p1a,p2a,q2a)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
    INTEGER,INTENT(OUT) :: i2
    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
  END SUBROUTINE wigner_canonical_extremal
  SUBROUTINE wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
                           I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
    IMPLICIT NONE
#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    CLASS(*),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
#endif
    REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix3
    INTEGER,INTENT(IN) :: I1,J1,lambda1,mu1,L1,kappa1max,I2,J2,lambda2,mu2,L2,kappa2max,I3,lambda3,mu3,L3,kappa3max,rhomax,numb
    INTEGER,DIMENSION(:),INTENT(IN) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wigner_can
    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
  END SUBROUTINE wigner_su3so3
  SUBROUTINE orthonormalization_matrix(I,J,lambda,mu,L,kappamax,matrix)
    IMPLICIT NONE
#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: matrix
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    CLASS(*),DIMENSION(:,:),INTENT(OUT) :: matrix
#endif
    INTEGER,INTENT(IN) :: I,J,lambda,mu,L,kappamax
  END SUBROUTINE orthonormalization_matrix
END INTERFACE

pqdim=(MAX(lambda1,mu1)+1)*((MAX(lambda2,mu2)+1)**2)
ALLOCATE(wigner_can(0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax),&
         matrix3(kappa3max,kappa3max),p1a(pqdim),p2a(pqdim),q2a(pqdim))

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
CALL orthonormalization_matrix(I3,J3,lambda3,mu3,L3,kappa3max,matrix3)

#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))

ALLOCATE(matrix1(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1)
CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2)
CALL wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
       I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
DEALLOCATE(wigner_can,matrix1,matrix2,matrix3,p1a,p2a,q2a)

#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))

sum1=lambda1+mu1+L1+2
IF(sum1<62)THEN
  indicator1=1
ELSE IF(sum1<=77)THEN
  IF(lambda1<=sum1/3.AND.mu1<=(sum1-1)/2-3*lambda1/4+1)THEN
    indicator1=2
  ELSE
    indicator1=1
  END IF
ELSE
  IF(lambda1<=sum1/3)THEN
    IF(mu1<=(sum1-1)/2-lambda1/2+MAX(0,(sum1-82)/3))THEN
      indicator1=2
    ELSE
      indicator1=1
    END IF
  ELSE IF(mu1<=(sum1-1)/2-sum1/6-2*(lambda1-sum1/3)+MAX(0,2*(sum1-86)/3))THEN
    indicator1=2
  ELSE
    indicator1=1
  END IF
END IF

sum2=lambda2+mu2+L2+2
IF(sum2<62)THEN
  indicator2=1
ELSE IF(sum2<=77)THEN
  IF(lambda2<=sum2/3.AND.mu2<=(sum2-1)/2-3*lambda2/4+1)THEN
    indicator2=2
  ELSE
    indicator2=1
  END IF
ELSE
  IF(lambda2<=sum2/3)THEN
    IF(mu2<=(sum2-1)/2-lambda2/2+MAX(0,(sum2-82)/3))THEN
      indicator2=2
    ELSE
      indicator2=1
    END IF
  ELSE IF(mu2<=(sum2-1)/2-sum2/6-2*(lambda2-sum2/3)+MAX(0,2*(sum2-86)/3))THEN
    indicator2=2
  ELSE
    indicator2=1
  END IF
END IF

IF(indicator1==1.AND.indicator2==1)THEN
  ALLOCATE(matrix1(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1)
  CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2)
  CALL wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
         I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1,matrix2,matrix3,p1a,p2a,q2a)
ELSE IF(indicator1==2.AND.indicator2==1)THEN
  ALLOCATE(matrix1mp(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1mp)
  CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2)
  CALL wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1mp,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
         I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1mp,matrix2,matrix3,p1a,p2a,q2a)
ELSE IF(indicator1==1.AND.indicator2==2)THEN
  ALLOCATE(matrix1(kappa1max,kappa1max),matrix2mp(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1)
  CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2mp)
  CALL wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2mp,&
         I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1,matrix2mp,matrix3,p1a,p2a,q2a)
ELSE
  ALLOCATE(matrix1mp(kappa1max,kappa1max),matrix2mp(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,lambda1,mu1,L1,kappa1max,matrix1mp)
  CALL orthonormalization_matrix(I2,J2,lambda2,mu2,L2,kappa2max,matrix2mp)
  CALL wigner_su3so3(I1,J1,lambda1,mu1,L1,kappa1max,matrix1mp,I2,J2,lambda2,mu2,L2,kappa2max,matrix2mp,&
         I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1mp,matrix2mp,matrix3,p1a,p2a,q2a)
END IF

#endif

END SUBROUTINE calculate_wigner_su3so3
