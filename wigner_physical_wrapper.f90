!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wigner_physical_wrapper.f90 -- wrapper for SU(3)-SO(3) coupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wigner_physical_wrapper(lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3,&
                                   kappa1max,kappa2max,kappa3max,rhomax,dimen,wigner_phys_block)
!---------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine wigner_physical_wrap calculating reduced SU(3)-SO(3) coupling coefficients
! <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
!
! Input arguments: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3
! Output argumrnts: kappa1max,kappa2max,kappa3max,rhomax,dimen,wigner_phys_block
!
! kappa1max = inner multiplicity od L1 within (lambda1,mu1)
! kappa2max = inner multiplicity od L2 within (lambda2,mu2)
! kappa3max = inner multiplicity od L3 within (lambda3,mu3)
! rhomax = outer multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! wigner_phys_block(n) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! n = kappa1+kappa1max*(kappa2-1)+kappa1max*kappa2max*(kappa3-1)+kappa1max*kappa2max*kappa3max*(rho-1) = 1,...,dimen
!
! If at least one of kappa1max,kappa2max,kappa3max,rhomax equals 0, dimen=1 with wigner_phys_block(1)=0.D0 is returned.
!
! Note: This wrapper accepts unallocated allocatable array wigner_phys_block, allocates it and returns it.
!----------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER,INTENT(IN) :: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3
INTEGER,INTENT(OUT) :: kappa1max,kappa2max,kappa3max,rhomax,dimen
REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: wigner_phys_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_phys
INTEGER :: ind,rho,kappa1,kappa2,kappa3

INTERFACE
  SUBROUTINE wigner_physical_wrap(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax
    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
  END SUBROUTINE wigner_physical_wrap
END INTERFACE

rhomax=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda3,mu3)
kappa1max=inner_multiplicity(lambda1,mu1,L1)
kappa2max=inner_multiplicity(lambda2,mu2,L2)
kappa3max=inner_multiplicity(lambda3,mu3,L3)
IF(rhomax==0.OR.kappa1max==0.OR.kappa2max==0.OR.kappa3max==0)THEN
  dimen=1
  ALLOCATE(wigner_phys_block(1))
  wigner_phys_block(1)=0.D0
ELSE
  dimen=kappa1max*kappa2max*kappa3max*rhomax
  ALLOCATE(wigner_phys(kappa1max,kappa2max,kappa3max,rhomax),wigner_phys_block(dimen))
  CALL wigner_physical_wrap(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
  ind=0
  DO rho=1,rhomax
    DO kappa3=1,kappa3max
      DO kappa2=1,kappa2max
        DO kappa1=1,kappa1max
          ind=ind+1
          wigner_phys_block(ind)=wigner_phys(kappa1,kappa2,kappa3,rho)
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(wigner_phys)
END IF

CONTAINS
  FUNCTION inner_multiplicity(lambda,mu,L) RESULT(kappamax)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda,mu,L
    INTEGER :: kappamax
    kappamax=MAX(0,(lambda+mu+2-L)/2)-MAX(0,(lambda+1-L)/2)-MAX(0,(mu+1-L)/2)
  END FUNCTION inner_multiplicity
END SUBROUTINE wigner_physical_wrapper
