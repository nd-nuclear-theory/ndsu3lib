!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! u_coeff_wrapper.f90 -- wrapper for SU(3) U recoupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE u_coeff_wrapper(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                           rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,racah_block,info)
!------------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine calculating SU(3) recoupling coefficients
! U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,
!                  rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
! Output argumrnts: racah_block,info
!
! rhomaxa is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12).
! rhomaxb is the multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu).
! rhomaxc is the multiplicity of coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23).
! rhomaxd is the multiplicity of coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu).
! dimen is the size of the array racah_block. It must be at least rhomaxa*rhomaxb*rhomaxc*rhomaxd.
! racah_block(ind) = U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
! ind = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1)
! info = 0 if MKL subroutine dgesv ran withou errors.
!-------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                      rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
INTEGER,INTENT(OUT) :: info
REAL(KIND=8),DIMENSION(dimen),INTENT(OUT) :: racah_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: rac
INTEGER :: rhomaxabc,rhoabc,rhod,ind

INTERFACE
  SUBROUTINE calculate_u_coeff(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                               rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                          rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: rac
  END SUBROUTINE calculate_u_coeff
END INTERFACE

rhomaxabc=rhomaxa*rhomaxb*rhomaxc
ALLOCATE(rac(rhomaxd,rhomaxabc))
CALL calculate_u_coeff(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,&
                       lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,rhomaxd,info)
ind=0
DO rhod=1,rhomaxd
  DO rhoabc=1,rhomaxabc
    ind=ind+1
    racah_block(ind)=rac(rhod,rhoabc)
  END DO
END DO
DEALLOCATE(rac)

END SUBROUTINE u_coeff_wrapper
