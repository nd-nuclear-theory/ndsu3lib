SUBROUTINE racah_wrapper(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                         rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,racah_block,info)
!----------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine racah calculating SU(3) U-recoupling coefficients
! U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
! Output argumrnts: rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,racah_block,info
!
! rhomaxa = multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomaxb = multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
! rhomaxc = multiplicity of coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23)
! rhomaxd = multiplicity of coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu)
! racah_block(n) = U[(lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3)rhoa,rhob(lambda12,mu12)(lambda23,mu23)rhoc,rhod]
! n = rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)+rhomaxa*rhomaxb*rhomaxc*(rhod-1) = 1,...,dimen
! info = 0 if MKL subroutine dgesv called by racah ran withou errors
!
! If at least one of rhomaxa,rhomaxb,rhomaxc,rhomaxd equals 0, dimen=1 with racah_block(1)=0.D0 is returned.
!
! Note: This wrapper accepts unallocated allocatabla array racah_block, allocates it and returns it.
!-----------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
INTEGER,INTENT(OUT) :: rhomaxa,rhomaxb,rhomaxc,rhomaxd,dimen,info
REAL(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: racah_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: rac
INTEGER :: rhomaxabc,rhoabc,rhod,ind

INTERFACE
  SUBROUTINE racah_allocatable_gsl(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                                   rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: rac
  END SUBROUTINE racah_allocatable_gsl
END INTERFACE

info=0
rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
rhomaxc=outer_multiplicity(lambda2,mu2,lambda3,mu3,lambda23,mu23)
rhomaxd=outer_multiplicity(lambda1,mu1,lambda23,mu23,lambda,mu)
IF(rhomaxa==0.OR.rhomaxb==0.OR.rhomaxc==0.OR.rhomaxd==0)THEN
  dimen=1
  ALLOCATE(racah_block(1))
  racah_block(1)=0.D0
ELSE
  rhomaxabc=rhomaxa*rhomaxb*rhomaxc
  dimen=rhomaxabc*rhomaxd
  ALLOCATE(rac(rhomaxd,rhomaxabc),racah_block(dimen))
  CALL racah_allocatable_gsl(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,&
                             lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,rhomaxd,info)
  ind=0
  DO rhod=1,rhomaxd
    DO rhoabc=1,rhomaxabc
      ind=ind+1
      racah_block(ind)=rac(rhod,rhoabc)
    END DO
  END DO
  DEALLOCATE(rac)
END IF
END SUBROUTINE racah_wrapper
