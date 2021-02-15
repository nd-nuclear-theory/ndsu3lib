!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Z_coeff_wrapper.f90 -- wrapper for SU(3) Z-recoupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Z_coeff_wrapper(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                           rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,Z_block,info)
!----------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine Z_coeff calculating SU(3) Z-recoupling coefficients
! Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
!
! Input arguments: lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
! Output argumrnts: rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,Z_block,info
!
! rhomaxa = multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomaxb = multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
! rhomaxc = multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
! rhomaxd = multiplicity of coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu)
! Z_block(n) = Z[(lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda13,mu13)rhoc,rhod]
! n = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1) = 1,...,dimen
! info = 0 if MKL subroutine dgesv called by Z_coeff ran withou errors
!
! If at least one of rhomaxa,rhomaxb,rhomaxc,rhomaxd equals 0, dimen=1 with Z_block(1)=0.D0 is returned.
!
! Note: This wrapper accepts unallocated allocatabla array Z_block, allocates it and returns it.
!-----------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
INTEGER,INTENT(OUT) :: rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,info
REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: Z_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: Zcoeff
INTEGER :: rhomaxabc,rhoabc,rhod,ind

INTERFACE
  SUBROUTINE Z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                     rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff
  END SUBROUTINE Z_coeff
END INTERFACE

info=0
rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
rhomaxc=outer_multiplicity(lambda1,mu1,lambda3,mu3,lambda13,mu13)
rhomaxd=outer_multiplicity(lambda13,mu13,lambda2,mu2,lambda,mu)
IF(rhomaxa==0.OR.rhomaxb==0.OR.rhomaxc==0.OR.rhomaxd==0)THEN
  dimen=1
  ALLOCATE(Z_block(1))
  Z_block(1)=0.D0
  RETURN
END IF
rhomaxabc=rhomaxa*rhomaxb*rhomaxc
dimen=rhomaxabc*rhomaxd
ALLOCATE(Zcoeff(rhomaxd,rhomaxabc),Z_block(dimen))
CALL Z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
             rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,rhomaxd,info)
ind=0
DO rhod=1,rhomaxd
  DO rhoabc=1,rhomaxabc
    ind=ind+1
    Z_block(ind)=Zcoeff(rhod,rhoabc)
  END DO
END DO
DEALLOCATE(Zcoeff)
END SUBROUTINE Z_coeff_wrapper
