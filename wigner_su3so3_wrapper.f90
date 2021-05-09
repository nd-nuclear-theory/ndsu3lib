!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wigner_su3so3_wrapper.f90 -- wrapper for SU(3)-SO(3) reduced Wigner coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wigner_su3so3_wrapper(lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3,&
                                 kappa1max,kappa2max,kappa3max,rhomax,dimen,wigner_phys_block)
!------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine calculating reduced SU(3)-SO(3) Wigner coefficients
! <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
!
! Input arguments: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3,kappa1max,kappa2max,kappa3max,rhomax,dimen
! Output argument: wigner_phys_block
!
! kappa1max is the inner multiplicity of L1 within (lambda1,mu1).
! kappa2max is the inner multiplicity of L2 within (lambda2,mu2).
! kappa3max is the inner multiplicity of L3 within (lambda3,mu3).
! rhomax is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3).
! dimen is the size of the array wigner_phys_block. It must be at least kappa1max*kaapa2max*kaapa3max*rhomax.
! wigner_phys_block(ind) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! ind = kappa1+kappa1max*(kappa2-1)+kappa1max*kappa2max*(kappa3-1)+kappa1max*kappa2max*kappa3max*(rho-1)
!-------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3,kappa1max,kappa2max,kappa3max,rhomax,dimen
REAL(KIND=8),DIMENSION(dimen),INTENT(OUT) :: wigner_phys_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_phys
INTEGER :: ind,rho,kappa1,kappa2,kappa3

INTERFACE
  SUBROUTINE calculate_wigner_su3so3(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax
    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
  END SUBROUTINE calculate_wigner_su3so3
END INTERFACE

ALLOCATE(wigner_phys(kappa1max,kappa2max,kappa3max,rhomax))
CALL calculate_wigner_su3so3(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax,wigner_phys)
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

END SUBROUTINE wigner_su3so3_wrapper
