!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculate_wigner_canonical.f90 -- SU(3)-SU(2)xU(1) reduced Wigner coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,dimpq,dimw,&
                                      rhomax,numb,wigner_block,p1a,p2a,q2a)
!-------------------------------------------------------------------------------------------------------------------
! Calculates all SU(3)-SU(2)xU(1) reduced Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3=Lambda32/2. The quantum numbers epsilon3,Lambda32
! must be valid, i.e. there must be integers p3,q3 satisfying:
!   0<=p3<=lambda3
!   0<=q3<=mu3
!   epsilon3=2*lambda3+mu3-3*(p3+q3)
!   Lambda3=(mu3+p3-q3)/2
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,dimpq,dimw,rhomax
! Output arguments: numb,wigner_block,p1a,p2a,q2a
!
! rmomax is the multiplicity of the coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3).
! dimpq is the size of the arrays p1a,p2a,q2a, which should be at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
! dimw is the size of the array wigner_block, which should be at least rhomax*dimpq.
!
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho=wigner_block(ind)
!   where epsilon2=2*lambda2+mu2-3*(p2+q2)
!         epsilon1=epsilon3-epsilon2
!         Lambda1=(mu1+p1-q1)/2
!         Lambda2=(mu2+p2-q2)/2
!         p1=p1a(i) 
!         p2=p2a(i)
!         q2=q2a(i)
!         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
!         ind=i+numb*(rho-1)
!         1<=i<=numb
!-------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,dimpq,dimw,rhomax
INTEGER,INTENT(OUT) :: numb
INTEGER,DIMENSION(dimpq),INTENT(OUT) :: p1a,p2a,q2a
REAL(KIND=8),DIMENSION(dimw),INTENT(OUT) :: wigner_block
INTEGER :: i,rho,ind
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wignerex,wigner

INTERFACE
  SUBROUTINE wigner_canonical_extremal(lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax,i2,wigner,p1a,p2a,q2a)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
    INTEGER,INTENT(OUT) :: i2
    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
  END SUBROUTINE wigner_canonical_extremal
  SUBROUTINE wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,I3,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,I3,rhomax
    INTEGER :: numb
    INTEGER,DIMENSION(:) :: p1a,p2a,q2a
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wignerex
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
  END SUBROUTINE wigner_canonical
END INTERFACE

ALLOCATE(wignerex(0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax),&
         wigner(0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax))

IF(2*epsilon3<=lambda3-mu3)THEN
  CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda3,mu3,1,rhomax,numb,wignerex,p1a,p2a,q2a)
  CALL wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,1,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
ELSE
  CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda3,mu3,0,rhomax,numb,wignerex,p1a,p2a,q2a)
  CALL wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,0,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
  IF(epsilon3==2*lambda3+mu3)THEN
    DO i=1,numb
      ind=p2a(i)
      p1a(i)=(2*(lambda1-mu1-mu2)-lambda2-epsilon3)/3+p1a(i)+q2a(i)+ind
      p2a(i)=lambda2-q2a(i)
      q2a(i)=mu2-ind
    END DO
  END IF
END IF

ind=0
DO rho=1,rhomax
  DO i=1,numb
    ind=ind+1 ! ind=i+numb*(rho-1)
    wigner_block(ind)=wigner(p1a(i),p2a(i),q2a(i),rho)
  END DO
END DO

DEALLOCATE(wignerex,wigner)

END SUBROUTINE calculate_wigner_canonical
