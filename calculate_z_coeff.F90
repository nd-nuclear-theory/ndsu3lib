!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculate_z_coeff.F90 -- SU(3) Z recoupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                             rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,ldb,info)
!---------------------------------------------------------------------------------------------------------------------------
! Calsulates SU(3) recoupling coefficients
! Z((lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda13,mu13)rhoc,rhod)
! for given lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
! using Eq.(2) in the reference and MKL subroutine dgesv solving a system of linear equations.
!
! Reference: D.J.Millener, J.Math.Phys., Vol.19, No.7 (1978) 1513
!
! Input arguments: lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
! Output arguments: Zcoeff,info
!
! rhomaxa = multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomaxb = multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
! rhomaxc = multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
! rhomaxd = multiplicity of coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu)
! ldb = the leading dimension of the array Zcoeff
!
! Zcoeff(rhod,n)=Z((lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda13,mu13)rhoc,rhod)
!   where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
! info=0 if dgesv ran without errors
!---------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
#if defined(NDSU3LIB_RACAH_GSL)
REAL(KIND=8),EXTERNAL :: su2racah
#elif defined(NDSU3LIB_RACAH_WIGXJPF)
REAL(KIND=8),EXTERNAL :: fwig6jj
REAL(KIND=8) :: su2racah
#endif
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
INTEGER,INTENT(OUT) :: info
INTEGER :: rhomaxabc,numbahw,numbalw,numbb,numbd,i,indd,p2,q2,indb,p12,p3,q3,epsilon,epsilon12,Lambda122,Lambda32,&
           p1,pq1,n,rhoa,rhob,rhoc,expp,Lambda22,LLambda12,p1min,aux1,mu1mpq1,aux3,aux4,lambdammu12,phase,&
           aux5,aux6,inddmin,epsilon1lwpepsilon2,epsilonmepsilon3lw,epsilon12lw,lambdamlambda13,pqdima,pqdimb,pqdimd
REAL(KIND=8) :: factor1,factor2,factor3
REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff ! Sizes are at least rhomaxd and rhomaxa*rhomaxb*rhomaxc
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wignerahw,wigneralw,wignerb,wignerc,wignerd,wigner
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1aa,p2aa,q2aa,p1ab,p2ab,q2ab,p2ad,q2ad

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

pqdima=(lambda1+1)*(lambda2+1)*(mu2+1)
pqdimb=MAX((lambda13+1)*(lambda2+1)*(mu2+1),(lambda1+1)*(lambda3+1)*(mu3+1),&
           (lambda1+1)*(lambda2+1)*(mu2+1),(lambda12+1)*(lambda3+1)*(mu3+1),(mu1+1)*(mu2+1)*(lambda2+1))
pqdimd=(lambda13+1)*(lambda2+1)*(mu2+1)
ALLOCATE(matrix(rhomaxd,rhomaxd),&
         wignerahw(0:lambda1,0:lambda2,0:mu2,rhomaxa),&
         wigneralw(0:mu1,0:mu2,0:lambda2,rhomaxa),&
         wignerb(0:lambda12,0:lambda3,0:mu3,rhomaxb),&
         wignerc(0:lambda1,0:lambda3,0:mu3,rhomaxc),&
         wignerd(0:lambda13,0:lambda2,0:mu2,rhomaxd),&
         wigner(0:lambda1,0:lambda2,0:mu2,rhomaxa),&
         p1aa(pqdima),&
         p2aa(MAX(pqdima,rhomaxd)),&
         q2aa(pqdima),&
         p1ab(pqdimb),&
         p2ab(pqdimb),&
         q2ab(pqdimb),&
         p2ad(pqdimd),&
         q2ad(pqdimd))

rhomaxabc=rhomaxa*rhomaxb*rhomaxc
epsilon=-lambda-2*mu
epsilon1lwpepsilon2=2*lambda1+mu1+epsilon+lambda13+2*mu13
epsilonmepsilon3lw=epsilon-2*lambda3-mu3
epsilon12lw=2*lambda12+mu12
lambdamlambda13=lambda-lambda13
lambdammu12=lambda12-mu12

CALL wigner_canonical_extremal(lambda13,mu13,lambda2,mu2,lambda,mu,1,rhomaxd,numbd,wignerd,p1ab,p2ad,q2ad)
CALL wigner_canonical_extremal(lambda1,mu1,lambda3,mu3,lambda13,mu13,1,rhomaxc,numbahw,wignerc,p1ab,p2ab,q2ab)
CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda12,mu12,1,rhomaxa,numbahw,wignerahw,p1ab,p2ab,q2ab)
CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda12,mu12,0,rhomaxa,numbalw,wigneralw,p1ab,p2ab,q2ab)
CALL wigner_canonical_extremal(lambda12,mu12,lambda3,mu3,lambda,mu,1,rhomaxb,numbb,wignerb,p1ab,p2ab,q2ab)

inddmin=numbd-rhomaxd+1

! Construction of matrix
i=0
DO indd=inddmin,numbd ! loop over Lambda2
  p2=p2ad(indd)
  q2=q2ad(indd)
  i=i+1
  matrix(i,1:rhomaxd)=wignerd(lambda13,p2,q2,1:rhomaxd)
END DO

! Construction of RHS
Zcoeff(1:rhomaxd,1:rhomaxabc)=0.D0
DO indb=1,numbb ! loop over epsilon1,epsilon3,epsilon12,Lambda3,Lambda12
  p12=p1ab(indb)
  p3=p2ab(indb)
  q3=q2ab(indb)
  epsilon12=epsilonmepsilon3lw+3*(p3+q3)!=epsilon-epsilon3
  Lambda122=mu12-(epsilon12lw-epsilon12)/3+2*p12
  Lambda32=mu3+p3-q3
  aux1=lambdamlambda13-Lambda122

  IF(2*epsilon12<=lambdammu12)THEN
    CALL wigner_canonical(lambda1,mu1,lambda2,mu2,lambda12,mu12,epsilon12,Lambda122,1,&
                          rhomaxa,numbahw,wignerahw,wigner,p1aa,p2aa,q2aa)
  ELSE
    CALL wigner_canonical(lambda1,mu1,lambda2,mu2,lambda12,mu12,epsilon12,Lambda122,0,&
                          rhomaxa,numbalw,wigneralw,wigner,p1aa,p2aa,q2aa)
  END IF

  factor1=DSQRT(DFLOAT((Lambda122+1)*(lambda13+1)))

  pq1=(epsilon1lwpepsilon2-epsilon12)/3 ! pq1 is p1+q1
  mu1mpq1=mu1-pq1
  aux3=MIN(lambda1,pq1,(Lambda32+lambda13-mu1mpq1)/2)
  aux4=MAX(0,-mu1mpq1,(Lambda32-lambda13-mu1mpq1)/2,-(Lambda32-lambda13+mu1mpq1)/2)
  aux5=Lambda122-mu1mpq1
  aux6=Lambda122+mu1mpq1

  i=0
  DO indd=inddmin,numbd ! loop over Lambda2
    i=i+1
    p2=p2ad(indd)
    q2=q2ad(indd)
    Lambda22=mu2+p2-q2

    p1min=MAX(aux4,(Lambda22-aux6)/2,(aux5-Lambda22)/2)
    LLambda12=mu1mpq1+2*p1min ! LLambda12 is 2*Lambda1
    expp=LLambda12+aux1
    IF(4*(expp/4)==expp)THEN
      factor2=-factor1
    ELSE
      factor2=factor1
    END IF
    DO p1=p1min,MIN((Lambda22+aux5)/2,aux3)
    ! Lower and upper bounds on p1 are such that:
    ! 1) 0<=q1<=mu1, where q1=pq1-p1
    ! 2) ABS(LLambda12-Lambda22)<=Lambda122<=LLambda12+Lambda22, where LLambda12=mu1+p1-q1=mu1-pq1+2*p1
    ! 3) ABS(LLambda12-Lambda32)<=lambda13<=LLambda12+Lambda32
      
      factor2=-factor2
#if defined(NDSU3LIB_RACAH_GSL)
      factor3=factor2*su2racah(Lambda22,LLambda12,lambda,Lambda32,Lambda122,lambda13)
#elif defined(NDSU3LIB_RACAH_WIGXJPF)
      su2racah=fwig6jj(Lambda22,LLambda12,Lambda122,Lambda32,lambda,lambda13)
      phase=Lambda22+LLambda12+Lambda32+lambda
      IF((phase/4)*4/=phase)su2racah=-su2racah
      factor3=factor2*su2racah
#endif

      n=0
      DO rhoc=1,rhomaxc
        DO rhob=1,rhomaxb
          DO rhoa=1,rhomaxa
            ! n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
            n=n+1
            Zcoeff(i,n)=Zcoeff(i,n)+factor3*wignerc(p1,p3,q3,rhoc)&
                        *wigner(p1,p2,q2,rhoa)*wignerb(p12,p3,q3,rhob)
          END DO
        END DO
      END DO

      LLambda12=LLambda12+2
    END DO
  END DO
END DO

! Solution of system of linear equations
IF(rhomaxd>1)THEN
  CALL dgesv(rhomaxd,rhomaxabc,matrix,rhomaxd,p2aa,Zcoeff,ldb,info)
ELSE
  Zcoeff(1,1:rhomaxabc)=Zcoeff(1,1:rhomaxabc)/matrix(1,1)
END IF

DEALLOCATE(matrix,wignerahw,wigneralw,wignerb,wignerc,wignerd,wigner,p1aa,p2aa,q2aa,p1ab,p2ab,q2ab,p2ad,q2ad)

END SUBROUTINE calculate_z_coeff
