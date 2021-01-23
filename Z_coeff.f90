SUBROUTINE Z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                   rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,info)
!---------------------------------------------------------------------------------------------------------------------------
! Calsulates SU(3) recoupling coefficients
! Z((lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda13,mu13)rhoc,rhod)
! for given lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13
! using Eq.(2) in the reference and MKL subroutine dgesv solving a system of linear equations.
!
! Reference: D.J.Millener, J.Math.Phys., Vol.19, No.7 (1978) 1513
!
! Input arguments: lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd
! Output arguments: Zcoeff,info
!
! rhomaxa = multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomaxb = multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
! rhomaxc = multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
! rhomaxd = multiplicity of coupling (lambda13,mu13)x(lambda2,mu2)->(lambda,mu)
! 
! Zcoeff(rhod,n)=Z((lambda2,mu2)(lambda1,mu1)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda13,mu13)rhoc,rhod)
!   where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
! info=0 if dgesv ran without errors
!---------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: su2racah ! TO DO: Replace DRR3
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd
INTEGER,INTENT(OUT) :: info
INTEGER :: rhomaxabc,numba,numbb,numbc,numbd,i,indd,p2,q2,indb,p12,p3,q3,epsilon,epsilon3,epsilon12,q12,Lambda122,Lambda32,&
           epsilon13,epsilon2,p1,q1,pq1,epsilon1,n,rhoa,rhob,rhoc,expp,Lambda22,LLambda12
REAL(KIND=8) :: factor1,factor2
REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff ! Dimensions are at least rhomaxd and rhomaxa*rhomaxb*rhomaxc
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wignera,wignerb,wignerc,wignerd,wigner
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1aa,p2aa,q2aa,p1ab,p2ab,q2ab,p1ac,p2ac,q2ac,p1ad,p2ad,q2ad

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

ALLOCATE(matrix(rhomaxd,rhomaxd),wignera(0:10,0:10,0:10,1:9),&
                     wignerb(0:10,0:10,0:10,1:9),&
                     wignerc(0:10,0:10,0:10,1:9),&
                     wignerd(0:10,0:10,0:10,1:9),&
                     wigner(0:10,0:10,0:10,1:9),&
         p1aa(1331),&
         p2aa(1331),&
         q2aa(1331),&
         p1ab(1331),&
         p2ab(1331),&
         q2ab(1331),&
         p1ac(1331),&
         p2ac(1331),&
         q2ac(1331),&
         p1ad(1331),&
         p2ad(1331),&
         q2ad(1331))

rhomaxabc=rhomaxa*rhomaxb*rhomaxc
epsilon=-lambda-2*mu
epsilon13=-lambda13-2*mu13
epsilon2=epsilon-epsilon13

CALL wigner_canonical_extremal(lambda13,mu13,lambda2,mu2,lambda,mu,1,rhomaxd,numbd,wignerd,p1ad,p2ad,q2ad)
CALL wigner_canonical_extremal(lambda1,mu1,lambda3,mu3,lambda13,mu13,1,rhomaxc,numbc,wignerc,p1ac,p2ac,q2ac)
CALL wigner_canonical_extremal(lambda12,mu12,lambda3,mu3,lambda,mu,1,rhomaxb,numbb,wignerb,p1ab,p2ab,q2ab)
CALL wigner_canonical_extremal(lambda1,mu1,lambda2,mu2,lambda12,mu12,1,rhomaxa,numba,wignera,p1aa,p2aa,q2aa)

! Construction of matrix
i=0
DO indd=numbd,numbd-rhomaxd+1,-1 ! This is a loop over Lambda2
  p2=p2ad(indd)
  q2=q2ad(indd)
  i=i+1
  matrix(i,1:rhomaxd)=wignerd(lambda13,p2,q2,1:rhomaxd)
END DO

! Construction of RHS
Zcoeff(1:rhomaxd,1:rhomaxabc)=0.D0
DO indb=1,numbb ! This is a loop over epsilon1,epsilon3,epsilon12,Lambda3,Lambda12
  p12=p1ab(indb)
  p3=p2ab(indb)
  q3=q2ab(indb)
  epsilon3=2*lambda3+mu3-3*(p3+q3)
  epsilon12=epsilon-epsilon3
  q12=(2*lambda12+mu12-epsilon12)/3-p12
  Lambda122=mu12+p12-q12
  Lambda32=mu3+p3-q3

  CALL wigner_canonical(lambda1,mu1,lambda2,mu2,lambda12,mu12,epsilon12,Lambda122,1,rhomaxa,numba,wignera,wigner,p1aa,p2aa,q2aa)

  factor1=DSQRT(DFLOAT((Lambda122+1)*(lambda13+1)))

  epsilon1=epsilon12-epsilon2!=epsilon13-epsilon3
  pq1=(2*lambda1+mu1-epsilon1)/3 ! pq1 is p1+q1

  i=0
  DO indd=numbd,numbd-rhomaxd+1,-1 ! This is a loop over Lambda2
    i=i+1
    p2=p2ad(indd)
    q2=q2ad(indd)
    Lambda22=mu2+p2-q2

    DO p1=0,lambda1 ! Medze treba urcit, aby neboli potrebne kontroly nizsie.
      q1=pq1-p1
      IF(q1<0.OR.q1>mu1)CYCLE
      LLambda12=mu1+p1-q1
      IF(ABS(LLambda12-Lambda22)>Lambda122.OR.Lambda122>LLambda12+Lambda22)CYCLE
      IF(ABS(LLambda12-Lambda32)>lambda13.OR.lambda13>LLambda12+Lambda32)CYCLE

      expp=LLambda12+lambda-Lambda122-lambda13
      IF(4*(expp/4)==expp)THEN
        factor2=factor1*su2racah(Lambda22,LLambda12,lambda,Lambda32,Lambda122,lambda13)
      ELSE
        factor2=-factor1*su2racah(Lambda22,LLambda12,lambda,Lambda32,Lambda122,lambda13)
      END IF

      n=0
      DO rhoc=1,rhomaxc
        DO rhob=1,rhomaxb
          DO rhoa=1,rhomaxa
!            n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
            n=n+1
            Zcoeff(i,n)=Zcoeff(i,n)+factor2*wignerc(p1,p3,q3,rhoc)&
                        *wigner(p1,p2,q2,rhoa)*wignerb(p12,p3,q3,rhob)
          END DO
        END DO
      END DO

    END DO

  END DO

END DO

! Solution of system of linear equations
IF(rhomaxd>1)THEN
  CALL dgesv(rhomaxd,rhomaxabc,matrix,rhomaxd,p2aa,Zcoeff,9,info)
ELSE
  Zcoeff(1,1:rhomaxabc)=Zcoeff(1,1:rhomaxabc)/matrix(1,1)
END IF

DEALLOCATE(matrix,wignera,wignerb,wignerc,wignerd,wigner,p1aa,p2aa,q2aa,p1ab,p2ab,q2ab,p1ac,p2ac,q2ac,p1ad,p2ad,q2ad)

RETURN
CONTAINS
  FUNCTION dimen(lambdax,mux) RESULT(dm) ! Function dimen(lambdax,mux) calculates the dimension of the SU(3) irrep (lambdax,mux).
    IMPLICIT NONE
    INTEGER,INTENT (IN) :: lambdax,mux
    INTEGER :: dm
    dm=(lambdax+1)*(mux+1)*(lambdax+mux+2)/2
  END FUNCTION dimen
END SUBROUTINE Z_coeff
