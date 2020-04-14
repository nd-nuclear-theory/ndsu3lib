FUNCTION transformation_coeff(I,J,lambdax,mux,epsilonx,Lambda2p,MLambda2px,M,L,Mp) RESULT(coeff)
!--------------------------------------------------------------------------------------
! Calculates transformation coefficient for SU(3)-SO(3) reduction <G|(G_E)MLM'>, where
! |G>=|(lambda,mu)epsilon,Lambda',M'_Lambda>, using equations (26),(32) and (33,6B) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input parameters: (I,J)=(1,1) => E=HW, (I,J)=(1,0) => E=HW', (I,J)=(0,0) => E=LW, (I,J)=(0,1) => E=LW'
!                   lambdax=lambda, mux=mu, epsilonx=epsilon, Lambda2p=2*Lambda', MLambda2px=2*M'_Lambda, Mp=M'
!
! Note: There are 2 typos in equation (26):
!       1) In the overall factor C, the factor of (2L+1) should not be squared.
!       2) The second S_1 should be S_1(M_Lambda=Lambda,Lambda,N_Lambda,M) instead of
!          S_1(N_Lambda,Lambda,M_Lambda=Lambda,M).
!--------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: I,J,lambdax,mux,epsilonx,Lambda2p,MLambda2px,M,L,Mp
REAL(KIND=8) :: coeff,C,S2,S2b!,r
REAL(KIND=8) :: S11,S12,S1b ! These are actually integers but might be very large
INTEGER :: lambda,mu,epsilon,MLambda2p,p,q,pq,a1,a2,a3,a4,LmM,LpM,LmMp,MmLambda,gamma,NLambda2p,k2L,k2Lo,kkappa,alpha,beta,&
           alphamin,alphamax,betamin,betamax,alphaminS2,alphamaxS2,LambdaM,MNLambda2p,LambdapMp,x,y,xm,xn!,NLambda2

IF(I==1)THEN
  lambda=lambdax
  mu=mux
  epsilon=epsilonx
  MLambda2p=MLambda2px
ELSE ! See Eq.(32) and the text below it.
  lambda=mux
  mu=lambdax
  epsilon=-epsilonx
  MLambda2p=-MLambda2px
ENDIF

p=(2*(lambda-mu)-epsilon)/3+Lambda2p
IF((p>=0).AND.(p<=2*lambda).AND.(2*(p/2)==p))THEN
  p=p/2
ELSE
  PRINT*,"transformation_coeff: Invalid values of epsilon,Lambda2"
  RETURN
END IF

!q=(2*lambda+mu-epsilon)/3+mu-Lambda2p
!IF((q>=0).AND.(q<=2*mu).AND.(2*(q/2)==q))THEN
!  q=q/2
!ELSE
!  PRINT*,"transformation_coeff: Invalid values of epsilon,Lambda2"
!  RETURN
!END IF

coeff=0.D0
! Nasledovny test parnosti 2*Lambda'+M' asi nie je nutny, lebo tato funkcia sa pre pripady neparneho 2*Lambda'+M' volat nebude.
LambdapMp=Lambda2p+Mp
IF(2*(LambdapMp/2)==LambdapMp)THEN
  LambdapMp=LambdapMp/2
ELSE
  RETURN ! If (Lambda'+M'/2) is not integer, the coefficient is zero.
END IF

q=p+mu-Lambda2p
LmM=L-M
LmMp=L-Mp

alphaminS2=MAX(0,-Mp-M) ! alphaminS2 in the lower bound for alpha in S_2 in Eq.(26)
alphamaxS2=MIN(LmM,LmMp) ! alphamaxS2 in the upper bound for alpha in S_2 in Eq.(26)
!print*,"S2: alphamin,alphamax=",alphaminS2,alphamaxS2
IF(alphaminS2>alphamaxS2)RETURN

LpM=L+M
LambdaM=(lambda+M)/2
MmLambda=(M-lambda)/2
xm=(Lambda2p-MLambda2p)/2
k2Lo=lambda+mu+L
pq=p+q
a1=q+(lambda+M+Lambda2p+Mp)/2
!print*,"q+Lambda+M/2+Lambda'+M'/2=",a1

DO gamma=0,p
!  IF(gamma>Lambda2)EXIT !????????????????????????????

  !kappa2= k2-2*(lambda-p+mu-q) ! kappa2 is 2*kappa'

  NLambda2p=2*(p-gamma)-Lambda2p ! NLambda2p is 2*N'_Lambda
  MNLambda2p=(MLambda2p+NLambda2p)/2
  xn=(Lambda2p-NLambda2p)/2
  alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
  alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
  IF(alphamin>alphamax)CYCLE

  betamin=MAX(0,MmLambda+gamma) ! betamin is the lower bound on beta in the second S_1 in Eq.(26)
  betamax=MIN(LambdaM,gamma) ! betamax is the upper bound on beta in the second S_1 in Eq.(26)
!print*,"S12: gamma,betamin,betamax=",gamma,betamin,betamax
  IF(betamin<=betamax)THEN
!    NLambda2=Lambda2-2*gamma !NLambda2 is 2*N_Lambda
    S12=0.D0 ! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION: There is a typo: it should be
!              S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)
    DO beta=betamin,betamax
      IF((beta/2)*2==beta)THEN
        S12=S12+binom(gamma,beta)*binom(lambda-gamma,LambdaM-beta)
      ELSE
        S12=S12-binom(gamma,beta)*binom(lambda-gamma,LambdaM-beta)
      END IF
    END DO
    !S12=S12*binom(lambda,gamma)
  ELSE
    CYCLE
  END IF
!print*,"S12,gamma=",S12,gamma

  S11=0.D0 ! S11 is the first S_1 in Eq.(26)
!print*,"S11: gamma,alphamin,alphamax=",gamma,alphamin,alphamax
  DO alpha=alphamin,alphamax
    x=2*alpha+MNLambda2p
    y=Lambda2p-x
    betamin=MAX(0,LambdapMp-x)
    betamax=MIN(y,LambdapMp)
!print*,"S11b: gamma,apha,betamin,betamax=",gamma,alpha,betamin,betamax
    IF(betamin<=betamax)THEN
      S1b=0.D0 !S1b is the sum over beta in S11
      DO beta=betamin,betamax
        IF((beta/2)*2==beta)THEN
          S1b=S1b+binom(y,beta)*binom(x,LambdapMp-beta)
        ELSE
          S1b=S1b-binom(y,beta)*binom(x,LambdapMp-beta)
        END IF
      END DO
!print*,"S11b,gamma=",S1b,gamma
    ELSE
      CYCLE
    END IF
    S11=S11+binom(xn,alpha)*binom(p-gamma,xm-alpha)*S1b
  END DO
!print*,"S11,gamma=",S11,gamma

  k2L=k2Lo-gamma ! k2L is 2*k+L
!print*,"2k+L,gamma=",k2L,gamma
  kkappa=pq-gamma ! kkappa is k+kappa'
!print*,"k+kappa',gamma=",kkappa,gamma
  a2=a1-gamma
  a3=a2-k2L
  S2=0.D0 ! S2 is S_2 is Eq.(26)
  DO alpha=alphaminS2,alphamaxS2
    a4=a2+alpha
!print*,"q+Lambda+M/2+Lambda'+M'/2+alpha-gamma,gamma,alpha=",a4,gamma,alpha
    betamin=MAX(0,a3+alpha)
    betamax=MIN(a4,kkappa)
!print*,"S2: gamma,alpha,betamin,betamax",gamma,alpha,betamin,betamax
    IF(betamin<=betamax)THEN
      S2b=0.D0 ! S2b is the sum over beta in S2
      DO beta=betamin,betamax
!print*,"binom(k+kappa',beta),binom(2k+L,...),2k+L,a4-beta=",DFLOAT(binom(kkappa,beta)),DFLOAT(binom(k2L,a4-beta)),k2L,a4-beta
        IF((beta/2)*2==beta)THEN
          S2b=S2b+binom(kkappa,beta)/binom(k2L,a4-beta)
        ELSE
          S2b=S2b-binom(kkappa,beta)/binom(k2L,a4-beta)
        END IF
      END DO
    ELSE
      CYCLE
    END IF
!print*,"S2b,gamma,alpha=",S2b,gamma,alpha
    IF((alpha/2)*2==alpha)THEN
      S2=S2+binom(LmM,alpha)*binom(LpM,LmMp-alpha)*S2b
    ELSE
      S2=S2-binom(LmM,alpha)*binom(LpM,LmMp-alpha)*S2b
    END IF
  END DO
  S2=S2/DFLOAT(k2L+1)
!print*,"S2,gamma=",S2,gamma

  coeff=coeff+binom(p,gamma)*S11*S12*S2
END DO

IF(coeff/=0.D0)THEN
!  r=DFLOAT(2*L+1)/DFLOAT(2**p) ! This factor appears in C in Eq.(26) as squared but that is probably a typo. The factor of (2L+1) should not be squared!
  C=DFLOAT(2*L+1)/DFLOAT(2**(2*p))*DSQRT((binom(lambda,p)*binom(mu,q)*binom(lambda+mu+1,q)*binom(2*L,L-M))&
    /(binom(2*L,L-Mp)*binom(Lambda2p,(Lambda2p+MLambda2p)/2)*binom(p+mu+1,q)))
  IF(2*((L-p)/2)/=L-p)C=-C
  coeff=coeff*C
END IF

!print*,"C=",C

! Now coeff is the coefficient for E=HW. For E=HW' there is additinal phase factor of (-1)^((lambda+M)/2) according to Eq.(33,6B), therefore:
IF(I/=J.AND.BTEST((lambda+M)/2,0))coeff=-coeff

RETURN

CONTAINS

 function factorial (n) result (res)

    implicit none
    integer, intent (in) :: n
    real(KIND=8) :: res
    integer :: i

!    res = product ((/(i, i = 1, n)/))

    res=1.d0
    do i=2,n
      res=res*dfloat(i)
    end do

  end function factorial

  function binom (n, k) result (res)

    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    real(KIND=8) :: res

    res = factorial (n) / (factorial (k) * factorial (n - k))

  end function binom
                                              
END FUNCTION transformation_coeff
