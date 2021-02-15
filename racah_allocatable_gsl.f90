SUBROUTINE racah_allocatable_gsl(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                 rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,ldb,info)
!---------------------------------------------------------------------------------------------------------------------------
! Calsulates SU(3) recoupling coefficients
! U((lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda23,mu23)rhoc,rhod)
! for given lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23
! using Eq.(22),(35,1B) in the reference and MKL subroutine dgesv solving a system of linear equations.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
! Output arguments: rac,info
!
! rhomaxa = multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomaxb = multiplicity of coupling (lambda12,mu12)x(lambda3,mu3)->(lambda,mu)
! rhomaxc = multiplicity of coupling (lambda2,mu2)x(lambda3,mu3)->(lambda23,mu23)
! rhomaxd = multiplicity of coupling (lambda1,mu1)x(lambda23,mu23)->(lambda,mu)
! ldb = the leading dimension of the array rac
!
! rac(rhod,n)=U((lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda23,mu23)rhoc,rhod)
!   where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
! info=0 if dgesv ran without errors
!---------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!INTEGER,EXTERNAL :: outer_multiplicity
REAL(KIND=8),EXTERNAL :: su2racah ! TO DO: Replace DRR3
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
INTEGER,INTENT(OUT) :: info
INTEGER :: epsilon23,rhomaxabc,numba,numbb,numbc,numbd,i1,i2,inda,indd,i,&
           Lambda122,epsilon2,Lambda22,p3,q3,n,rhoa,rhob,rhoc,I23,&
           Lambda232,Lambda32,p23,q23,p12,q12,p2,q2,m,noname3,epsilon2max,aux,p3min!,j,epsilon2ind,epsilon23ind,epsilon3ind,p3min,p3max
REAL(KIND=8) :: factor1,factor2,factor3
REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: rac ! Dimensions are at least rhomaxd and rhomaxa*rhomaxb*rhomaxc
!REAL(KIND=8),DIMENSION(9,9) :: matrix ! Dimensions are at least rhomaxd and rhomaxd
!REAL(KIND=8),DIMENSION(0:20,0:20,0:20,1:9) :: wignera,wignerb,wignerc,wignerd,wigner
!INTEGER,DIMENSION(9261) :: p1aa,p2aa,q2aa,p1ac,p2ac,q2ac,p2ad,q2ad
INTEGER :: pqdima,pqdimc,pqdimd
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wignera,wignerb,wignerc,wignerd,wigner
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1aa,p2aa,q2aa,p1ac,p2ac,q2ac,p2ad,q2ad

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

pqdima=(lambda12+1)*(MAX(mu2,lambda3)+1)*(MAX(lambda2,mu3)+1)
pqdimc=(MAX(lambda2,mu2)+1)*(lambda3+1)*(mu3+1)
pqdimd=(lambda1+1)*(lambda23+1)*(mu23+1)
! How about automatic arrays?
ALLOCATE(matrix(rhomaxd,rhomaxd),wignera(0:lambda12,0:mu2,0:lambda2,1:rhomaxa),&
                                 wignerb(0:lambda12,0:lambda3,0:mu3,1:rhomaxb),&
                                 wignerc(0:MAX(lambda2,mu2),0:MAX(lambda3,mu3),0:MAX(lambda3,mu3),1:rhomaxc),&
                                 wignerd(0:lambda1,0:lambda23,0:mu23,1:rhomaxd),&
                                 wigner(0:lambda2,0:lambda3,0:mu3,1:rhomaxc),&
         p1aa((MAX(lambda12,lambda1)+1)*(MAX(mu2,lambda3,lambda23)+1)*(MAX(lambda2,mu3,mu23)+1)),&
         p2aa(MAX(pqdima,rhomaxd)),&
         q2aa(pqdima),&
         p1ac(pqdimc),&
         p2ac(pqdimc),&
         q2ac(pqdimc),&
         p2ad(pqdimd),&
         q2ad(pqdimd))
!ALLOCATE(matrix(9,9),wignera(0:20,0:20,0:20,1:9),&
!                     wignerb(0:20,0:20,0:20,1:9),&
!                     wignerc(0:20,0:20,0:20,1:9),&
!                     wignerd(0:20,0:20,0:20,1:9),&
!                     wigner(0:20,0:20,0:20,1:9),&
!         p1aa(9261),&
!         p2aa(9261),&
!         q2aa(9261),&
!         p1ac(9261),&
!         p2ac(9261),&
!         q2ac(9261),&
!         p2ad(9261),&
!         q2ad(9261))

!print*,"rhomaxa=",rhomaxa
!print*,"rhomaxb=",rhomaxb
!print*,"rhomaxc=",rhomaxc
!print*,"rhomaxd=",rhomaxd

!rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
!rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
!rhomaxc=outer_multiplicity(lambda2,mu2,lambda3,mu3,lambda23,mu23)
!rhomaxd=outer_multiplicity(lambda1,mu1,lambda23,mu23,lambda,mu)
rhomaxabc=rhomaxa*rhomaxb*rhomaxc
epsilon23=-lambda-2*mu+lambda1+2*mu1

IF(2*epsilon23<=lambda23-mu23)THEN
  I23=1
ELSE
  I23=0
END IF

!epsilon23ind=(epsilon23+lambda23+2*mu23)/3!mozno netreba
epsilon2max=2*mu2+lambda2
i1=3*lambda1-6*lambda12+6*mu1-6*mu12+8*lambda2+4*mu2
factor1=DFLOAT((lambda1+1)*dimen(lambda12,mu12))/DFLOAT(dimen(lambda1,mu1))
m=(2*(lambda12+mu2+mu1-mu12)+lambda2+lambda1)/3
aux=2*lambda3+mu3-epsilon23

rac(1:rhomaxd,1:rhomaxabc)=0.D0

CALL wigner_canonical_extremal(lambda1,mu1,lambda23,mu23,lambda,mu,1,rhomaxd,numbd,wignerd,p1aa,p2ad,q2ad)
CALL wigner_canonical_extremal(lambda2,mu2,lambda3,mu3,lambda23,mu23,I23,rhomaxc,numbc,wignerc,p1ac,p2ac,q2ac)

!print*,"-----"
!do inda=1,numbc
!print*,"wignerc=",wignerc(p1ac(inda),p2ac(inda),q2ac(inda),1)
!end do

CALL wigner_canonical_extremal(lambda12,mu12,lambda3,mu3,lambda,mu,1,rhomaxb,numbb,wignerb,p1aa,p2aa,q2aa)
CALL wigner_canonical_extremal(lambda12,mu12,mu2,lambda2,lambda1,mu1,1,rhomaxa,numba,wignera,p1aa,p2aa,q2aa)

i=0
!j=rhomaxd+1
Lambda232=mu23+p2ad(numbd)-q2ad(numbd)
DO indd=numbd,numbd-rhomaxd+1,-1 ! This is a loop over Lambda23
!  Lambda232=Lambda22ad(indd)
  p23=p2ad(indd)
  q23=q2ad(indd)
!  Lambda232=mu23+p23-q23
!  j=j-1
  i=i+1
!  matrix(i,1:rhomaxd)=wignerd(lambda1,epsilon23ind,Lambda232,1:rhomaxd)
  matrix(i,1:rhomaxd)=wignerd(lambda1,p23,q23,1:rhomaxd)
  factor2=DSQRT(factor1*DFLOAT(Lambda232+1))

  CALL wigner_canonical(lambda2,mu2,lambda3,mu3,lambda23,mu23,epsilon23,Lambda232,I23,&
                        rhomaxc,numbc,wignerc,wigner,p1ac,p2ac,q2ac)

!print*,"-----"
!do inda=1,numbc
!print*,"wigner=",wigner(p1ac(inda),p2ac(inda),q2ac(inda),1)
!end do
!print*,"epsilon23,Lambda232=",epsilon23,Lambda232

  DO inda=1,numba ! sum over epsilon2,Lambda2,Lambda12
!    Lambda122=Lambda12aa(inda)
    p12=p1aa(inda)
    p2=p2aa(inda)
    q2=q2aa(inda)
!    q12=n-p12-p2-q2
!    Lambda122=mu12+p12-q12
    Lambda122=2*p12-m+p2+q2
!    epsilon2=epsilon2aa(inda) ! WARNING: epsilon2 is -epsilon2 in the formula!
    epsilon2=epsilon2max-3*(p2+q2)!mozno netreba
!    epsilon2ind=(epsilon2+mu2+2*lambda2)/3!mozno netreba
!    Lambda22=Lambda22aa(inda)
    Lambda22=lambda2+p2-q2
!    epsilon3ind=(epsilon23+epsilon2+lambda3+2*mu3)/3!mozno netreba
    noname3=(aux-epsilon2)/3
    i2=i1+epsilon2+3*Lambda122
!    p3min=MAX(0,lambda3-epsilon3ind)
!    p3max=MIN(lambda3,lambda3+mu3-epsilon3ind)
!    IF(p3min<=p3max)THEN
!      DO p3=p3min,p3max ! sum over Lambda3
      p3min=MAX(0,noname3-mu3)
      q3=noname3-p3min
      Lambda32=mu3+p3min-q3
      DO p3=p3min,MIN(lambda3,noname3)
!        q3=lambda3+mu3-epsilon3ind-p3
!        q3=noname3-p3
!        Lambda32=mu3+p3-q3
!        Lambda32=2*p3-lambda3+epsilon3ind
        factor3=factor2*su2racah(lambda1,Lambda22,lambda,Lambda32,Lambda122,Lambda232) ! TO DO: Replace DRR3

!print*,"factor2=",factor2
!print*,"su2racah=",su2racah(lambda1,Lambda22,lambda,Lambda32,Lambda122,Lambda232)

        IF(12*(i2/12)/=i2)factor3=-factor3

        n=0
        DO rhoc=1,rhomaxc
          DO rhob=1,rhomaxb
            DO rhoa=1,rhomaxa
!              n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
              n=n+1
              rac(i,n)=rac(i,n)+factor3*wignera(p12,p2,q2,rhoa)&
                       *wignerb(p12,p3,q3,rhob)*wigner(lambda2-q2,p3,q3,rhoc) ! See the relation between p and \tilde{p} in Eq.(32)

!print*,"factor3=",factor3
!print*,"wignera(p12,p2,q2,rhoa)=",wignera(p12,p2,q2,rhoa)
!print*,"wignerb(p12,p3,q3,rhob)=",wignerb(p12,p3,q3,rhob)
!print*,"wigner(lambda2-q2,p3,q3,rhoc)=",wigner(lambda2-q2,p3,q3,rhoc)

            END DO
          END DO
        END DO
        
        q3=q3-1
        Lambda32=Lambda32+2
      END DO
!    else
!      print*,"p3min>p3max"
!    END IF
  END DO
  Lambda232=Lambda232-2
END DO

!print*,matrix(1,1),matrix(1,2),rac(1,1)
!print*,matrix(2,1),matrix(2,2),rac(2,1)

IF(rhomaxd>1)THEN
  CALL dgesv(rhomaxd,rhomaxabc,matrix,rhomaxd,p2aa,rac,ldb,info)
ELSE
  rac(1,1:rhomaxabc)=rac(1,1:rhomaxabc)/matrix(1,1)

!print*,"rac(1,1)=",rac(1,1),"matrix(1,1)=",matrix(1,1)

END IF

DEALLOCATE(matrix,wignera,wignerb,wignerc,wignerd,wigner,p1aa,p2aa,q2aa,p1ac,p2ac,q2ac,p2ad,q2ad)

RETURN
CONTAINS
  FUNCTION dimen(lambdax,mux) RESULT(dm) ! Function dimen(lambdax,mux) calculates the dimension of the SU(3) irrep (lambdax,mux).
    IMPLICIT NONE
    INTEGER,INTENT (IN) :: lambdax,mux
    INTEGER :: dm
    dm=(lambdax+1)*(mux+1)*(lambdax+mux+2)/2
  END FUNCTION dimen
END SUBROUTINE racah_allocatable_gsl
