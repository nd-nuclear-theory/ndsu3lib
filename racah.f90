SUBROUTINE racah(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                 rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,info)
!---------------------------------------------------------------------------------------------------------------------
! Input arguments: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd
! Output arguments: rac,info
! 
! rac(rhod,n)=U((lambda1,mu1)(lambda2,mu2)(lambda,mu)(lambda3,mu3);(lambda12,mu12)rhoa,rhob(lambda23,mu23)rhoc,rhod)
!   where n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
! info=0 if dgesv ran without errors
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
REAL(KIND=8),EXTERNAL :: su2racah
INTEGER :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,info,epsilon23,&
           rhomaxa,rhomaxb,rhomaxc,rhomaxd,rhomaxabc,numba,numbb,numbc,numbd,i1,i2,inda,indd,i,&
           epsilon23ind,Lambda122,epsilon2,Lambda22,epsilon3ind,p3,q3,n,rhoa,rhob,rhoc,epsilon2ind,&
           Lambda232,Lambda32,p3min,p3max!,j
REAL(KIND=8) :: factor1,factor2,factor3
REAL(KIND=8),DIMENSION(9,729) :: rac
REAL(KIND=8),DIMENSION(9,9) :: matrix
REAL(KIND=8),DIMENSION(0:20,0:20,0:20,1:9) :: wignera,wignerb,wignerc,wignerd,wigner
INTEGER,DIMENSION(9261) :: Lambda12aa,Lambda22aa,epsilon2aa,Lambda12ac,Lambda22ac,epsilon2ac,Lambda22ad!,Lambda12ad,Lambda12ac,Lambda22ac,epsilon2ac

!rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
!rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
!rhomaxc=outer_multiplicity(lambda2,mu2,lambda3,mu3,lambda23,mu23)
!rhomaxd=outer_multiplicity(lambda1,mu1,lambda23,mu23,lambda,mu)
rhomaxabc=rhomaxa*rhomaxb*rhomaxc
epsilon23=-lambda-2*mu+lambda1+2*mu1
epsilon23ind=(epsilon23+lambda23+2*mu23)/3
i1=3*lambda1+8*lambda2-6*lambda12+6*mu1+4*mu2-6*mu12
factor1=DFLOAT((lambda1+1)*dimen(lambda12,mu12))/DFLOAT(dimen(lambda1,mu1))

rac(1:rhomaxd,1:rhomaxabc)=0.D0

CALL wigner_canonical_extremal(lambda1,mu1,lambda23,mu23,lambda,mu,1,rhomaxd,numbd,wignerd,Lambda12aa,Lambda22ad,epsilon2aa)
CALL wigner_canonical_extremal(lambda2,mu2,lambda3,mu3,lambda23,mu23,1,rhomaxc,numbc,wignerc,Lambda12aa,Lambda22aa,epsilon2aa)
CALL wigner_canonical_extremal(lambda12,mu12,lambda3,mu3,lambda,mu,1,rhomaxb,numbb,wignerb,Lambda12aa,Lambda22aa,epsilon2aa)
CALL wigner_canonical_extremal(lambda12,mu12,mu2,lambda2,lambda1,mu1,1,rhomaxa,numba,wignera,Lambda12aa,Lambda22aa,epsilon2aa)

i=0
!j=rhomaxd+1
DO indd=numbd,numbd-rhomaxd+1,-1 ! This is a loop over Lambda23
  Lambda232=Lambda22ad(indd)
!  j=j-1
  i=i+1
  matrix(i,1:rhomaxd)=wignerd(lambda1,epsilon23ind,Lambda232,1:rhomaxd)
  factor2=DSQRT(factor1*DFLOAT(Lambda232+1))

  CALL wigner_canonical(lambda2,mu2,lambda3,mu3,lambda23,mu23,epsilon23,Lambda232,&
                        rhomaxc,numbc,wignerc,wigner,Lambda12ac,Lambda22ac,epsilon2ac)

  DO inda=1,numba ! sum over epsilon2,Lambda2,Lambda12
    Lambda122=Lambda12aa(inda)
    epsilon2=epsilon2aa(inda) ! WARNING: epsilon2 is -epsilon2 in the formula!
    epsilon2ind=(epsilon2+mu2+2*lambda2)/3
    Lambda22=Lambda22aa(inda)
    epsilon3ind=(-lambda-2*mu+lambda1+2*mu1+epsilon2+lambda3+2*mu3)/3
    i2=i1+epsilon2+3*Lambda122
    p3min=MAX(0,lambda3-epsilon3ind)
    p3max=MIN(lambda3,lambda3+mu3-epsilon3ind)
    IF(p3min<=p3max)THEN
      DO p3=p3min,p3max ! sum over Lambda3
!        q3=lambda3+mu3-epsilon3ind-p3
!        Lambda32=mu3+p3-q3
        Lambda32=2*p3-lambda3+epsilon3ind
        factor3=factor2*su2racah(lambda1,Lambda22,lambda,Lambda32,Lambda122,Lambda232)
        IF(12*(i2/12)/=i2)factor3=-factor3

        n=0
        DO rhoc=1,rhomaxc
          DO rhob=1,rhomaxb
            DO rhoa=1,rhomaxa
!              n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
              n=n+1
              rac(i,n)=rac(i,n)+factor3*wignera(Lambda122,epsilon2ind,Lambda22,rhoa)&
                       *wignerb(Lambda122,epsilon3ind,Lambda32,rhob)*wigner(Lambda22,epsilon3ind,Lambda32,rhoc)

            END DO
          END DO
        END DO

      END DO
    END IF
  END DO
END DO

!print*,matrix(1,1),matrix(1,2),rac(1,1)
!print*,matrix(2,1),matrix(2,2),rac(2,1)

IF(rhomaxd>1)THEN
  CALL dgesv(rhomaxd,rhomaxabc,matrix,9,epsilon2aa,rac,9,info)
ELSE
  rac(1,1:rhomaxabc)=rac(1,1:rhomaxabc)/matrix(1,1)
END IF

RETURN
CONTAINS
  FUNCTION dimen(lambdax,mux) RESULT(dm) ! Function dimen(lambdax,mux) calculates the dimension of the SU(3) irrep (lambdax,mux).
    IMPLICIT NONE
    INTEGER,INTENT (IN) :: lambdax,mux
    INTEGER :: dm
    dm=(lambdax+1)*(mux+1)*(lambdax+mux+2)/2
  END FUNCTION dimen
END SUBROUTINE racah
