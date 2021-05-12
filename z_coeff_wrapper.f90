!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_coeff_wrapper.f90 -- wrapper for SU(3) Z recoupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE z_coeff_wrapper(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                           rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,Z_block,info)
!---------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine calculating SU(3) recoupling coefficients
! Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
!
! Input arguments: lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,
!                  rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
! Output argumrnts: Z_block,info
!
! rhomaxa is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12).
! rhomaxb is the multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu).
! rhomaxc is the multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13).
! rhomaxd is the multiplicity of coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu).
! dimen is the size of the array Z_block. It must be at least rhomaxa*rhomaxb*rhomaxc*rhomaxd.
! Z_block(ind) = Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
! ind = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1)
! info = 0 if MKL subroutine dgesv called by Z_coeff ran withou errors.
!----------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                      rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen
INTEGER,INTENT(OUT) :: info
REAL(KIND=8),DIMENSION(dimen),INTENT(OUT) :: Z_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: Zcoeff
INTEGER :: rhomaxabc,rhoabc,rhod,ind

INTERFACE
  SUBROUTINE calculate_z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                               rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,&
                          lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff
  END SUBROUTINE calculate_z_coeff
END INTERFACE

rhomaxabc=rhomaxa*rhomaxb*rhomaxc
ALLOCATE(Zcoeff(rhomaxd,rhomaxabc))
CALL calculate_z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                       rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,rhomaxd,info)
ind=0
DO rhod=1,rhomaxd
  DO rhoabc=1,rhomaxabc
    ind=ind+1
    Z_block(ind)=Zcoeff(rhod,rhoabc)
  END DO
END DO
DEALLOCATE(Zcoeff)

END SUBROUTINE z_coeff_wrapper
