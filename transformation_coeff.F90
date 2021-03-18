!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transformation_coeff.F90 -- coefficient of transformation between U(2) and SO(3) reduced bases
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION transformation_coeff(I,J,lambdax,mux,epsilonx,Lambda2p,MLambda2px,M,L,Mp) RESULT(coeff)
!--------------------------------------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------------------------------------
USE binomial_coeff
USE I_S_module
IMPLICIT NONE
INTEGER,INTENT(IN) :: I,J,lambdax,mux,epsilonx,Lambda2p,MLambda2px,M,L,Mp
REAL(KIND=8) :: coeff,S11,S12,S2
#if defined(SU3DBL)
REAL(KIND=8) :: coeffq,S11q,S12q,S2q
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16) :: coeffq,S11q,S12q,S2q
#endif
INTEGER :: lambda,mu,epsilon,MLambda2p,p,q,a4,LmM,LpM,LmMp,gamma,NLambda2p,k2L,kkappa,alpha,beta,&
           alphamin,alphamax,alphaminS2,alphamaxS2,LambdaM,MNLambda2p,LambdapMp,x,y,xm,xn,&
           aux1,ind1,aux2,ind2,aux3,ind3,aux4,xpys,lambdas,aux5
INTEGER(KIND=8) :: p8

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

!!p=(2*(lambda-mu)-epsilon)/3+Lambda2p
!!IF((p>=0).AND.(p<=2*lambda).AND.(2*(p/2)==p))THEN
!!  p=p/2
!!ELSE
!!  PRINT*,"transformation_coeff: Invalid values of epsilon,Lambda2"
!!  RETURN
!!END IF

p=((2*(lambda-mu)-epsilon)/3+Lambda2p)/2
p8=p

!q=(2*lambda+mu-epsilon)/3+mu-Lambda2p
!IF((q>=0).AND.(q<=2*mu).AND.(2*(q/2)==q))THEN
!  q=q/2
!ELSE
!  PRINT*,"transformation_coeff: Invalid values of epsilon,Lambda2"
!  RETURN
!END IF

!coeff=0.D0
! The following test is not necessary because this function will not be called for odd 2*Lambda'+M'.
!!LambdapMp=Lambda2p+Mp
!!IF(2*(LambdapMp/2)==LambdapMp)THEN
!!  LambdapMp=LambdapMp/2
!!ELSE
!!  RETURN ! If (Lambda'+M'/2) is not integer, the coefficient is zero.
!!END IF

LambdapMp=(Lambda2p+Mp)/2
q=p+mu-Lambda2p
LmM=L-M
LmMp=L-Mp
alphaminS2=MAX(0,-Mp-M) ! alphaminS2 in the lower bound for alpha in S_2 in Eq.(26)
alphamaxS2=MIN(LmM,LmMp)-2 ! alphamaxS2 in the upper bound for alpha in S_2 in Eq.(26) minus 2
LpM=L+M
LambdaM=(lambda+M)/2
xm=(Lambda2p-MLambda2p)/2
k2L=lambda+mu+L
kkappa=p+q
lambdas=lambda**2
NLambda2p=2*(p+1)-Lambda2p
MNLambda2p=(MLambda2p+NLambda2p)/2
xn=(Lambda2p-NLambda2p)/2
aux1=xn*(xn+1)/2
ind3=p*(p+1)/2
aux2=ind3+xm
aux3=LmM*(LmM+1)/2
aux4=LpM*(LpM+1)/2+LmMp
aux5=k2L*(k2L+1)/2+k2L-q-(lambda+M+Lambda2p+Mp)/2-alphaminS2

IF(lambda+mu+L<=17)THEN
coeff=0.D0

DO gamma=0,p
  MNLambda2p=MNLambda2p-1
  xn=xn+1
  aux1=aux1+xn

  IF(lambda-gamma>=gamma)THEN ! S12=I(lambda-gamma,gamma,LambdaM)
    S12=Ia(((lambda-gamma)**3+lambdas+gamma)/2+LambdaM)
  ELSE IF(BTEST(LambdaM,0))THEN ! S12=(-1)^(LambdaM)*I(lambda-gamma,gamma,LambdaM)
    S12=-Ia((gamma**3+lambdas+lambda-gamma)/2+LambdaM)
  ELSE
    S12=Ia((gamma**3+lambdas+lambda-gamma)/2+LambdaM)
  END IF
! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION: There is a typo: it should be
! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

  IF(S12/=0.D0)THEN

  alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
  alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
  S11=0.D0 ! S11 is the first S_1 in Eq.(26)
  x=2*alphamin+MNLambda2p
  y=Lambda2p-x
  xpys=(x+y)**2
  ind1=aux1+alphamin
  ind2=aux2-alphamin
  DO alpha=alphamin,alphamax
    IF(MAX(0,LambdapMp-x)<=MIN(y,LambdapMp))THEN
      IF(x>=y)THEN
        S11=S11+binom(ind1)*binom(ind2)*Ia((x**3+xpys+y)/2+LambdapMp)
      ELSE IF(BTEST(LambdapMp,0))THEN
        S11=S11-binom(ind1)*binom(ind2)*Ia((y**3+xpys+x)/2+LambdapMp)
      ELSE
        S11=S11+binom(ind1)*binom(ind2)*Ia((y**3+xpys+x)/2+LambdapMp)
      END IF
    ELSE
      EXIT
    END IF
    x=x+2
    y=y-2
    ind1=ind1+1
    ind2=ind2-1
  END DO
  IF(S11/=0.D0)THEN

  S2=0.D0 ! S2 is S_2 is Eq.(26)
  a4=upbound*kkappa*(kkappa+upbound+2)/2+aux5
  ind1=aux3+alphaminS2
  ind2=aux4-alphaminS2
  IF(BTEST(alphaminS2,0))THEN ! if alphaminS2 is odd
    DO alpha=alphaminS2,alphamaxS2,2
      S2=S2-binom(ind1)*binom(ind2)*Sa(a4)
      S2=S2+binom(ind1+1)*binom(ind2-1)*Sa(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2=S2-binom(ind1)*binom(ind2)*Sa(a4)
    ELSE
      S2=S2-binom(ind1)*binom(ind2)*Sa(a4)
      S2=S2+binom(ind1+1)*binom(ind2-1)*Sa(a4-1)
    END IF
  ELSE
    DO alpha=alphaminS2,alphamaxS2,2
      S2=S2+binom(ind1)*binom(ind2)*Sa(a4)
      S2=S2-binom(ind1+1)*binom(ind2-1)*Sa(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2=S2+binom(ind1)*binom(ind2)*Sa(a4)
      S2=S2-binom(ind1+1)*binom(ind2-1)*Sa(a4-1)
    ELSE
      S2=S2+binom(ind1)*binom(ind2)*Sa(a4)
    END IF
  END IF
  IF(S2/=0.D0)THEN
    coeff=coeff+binom(ind3)*S11*S12*S2/DFLOAT(k2L+1)
  END IF

  END IF
  END IF

  aux5=aux5-k2L
  k2L=k2L-1
  kkappa=kkappa-1
  aux2=aux2-p+gamma
  ind3=ind3+1
END DO

IF(coeff/=0.D0)THEN
! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo: (2L+1) should not be squared!
  aux1=lambda+mu+1
  aux2=p+mu+1
  aux3=2*L*(L+1)
  S2=DFLOAT(2*L+1)*DSQRT((binom((lambdas+lambda)/2+p)/binom(aux3-Mp))*(binom(mu*(mu+1)/2+q)&
    /binom((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom(aux1*(aux1+1)/2+q)/binom(aux2*(aux2+1)/2+q))*binom(aux3-M))/DFLOAT(4**p8)
! S2 is C
  IF(BTEST(L-p,0))S2=-S2
  coeff=coeff*S2
ELSE
  RETURN
END IF

ELSE ! lambda+mu+L>17
coeffq=0.Q0

DO gamma=0,p
  MNLambda2p=MNLambda2p-1
  xn=xn+1
  aux1=aux1+xn

  IF(lambda-gamma>=gamma)THEN ! S12=I(lambda-gamma,gamma,LambdaM)
    S12q=Ia(((lambda-gamma)**3+lambdas+gamma)/2+LambdaM)
  ELSE IF(BTEST(LambdaM,0))THEN ! S12=(-1)^(LambdaM)*I(lambda-gamma,gamma,LambdaM)
    S12q=-Ia((gamma**3+lambdas+lambda-gamma)/2+LambdaM)
  ELSE
    S12q=Ia((gamma**3+lambdas+lambda-gamma)/2+LambdaM)
  END IF
! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION:
! There is a typo: it should be
! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not
! S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

#if defined(SU3DBL)
  IF(S12q/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  IF(S12q/=0.Q0)THEN
#endif

  alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
  alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
#if defined(SU3DBL)
  S11q=0.D0 ! S11 is the first S_1 in Eq.(26)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  S11q=0.Q0
#endif
  x=2*alphamin+MNLambda2p
  y=Lambda2p-x
  xpys=(x+y)**2
  ind1=aux1+alphamin
  ind2=aux2-alphamin
  DO alpha=alphamin,alphamax
    IF(MAX(0,LambdapMp-x)<=MIN(y,LambdapMp))THEN
      IF(x>=y)THEN
        S11q=S11q+binom_quad(ind1)*binom_quad(ind2)*Ia((x**3+xpys+y)/2+LambdapMp)
      ELSE IF(BTEST(LambdapMp,0))THEN
        S11q=S11q-binom_quad(ind1)*binom_quad(ind2)*Ia((y**3+xpys+x)/2+LambdapMp)
      ELSE
        S11q=S11q+binom_quad(ind1)*binom_quad(ind2)*Ia((y**3+xpys+x)/2+LambdapMp)
      END IF
    ELSE
      EXIT
    END IF
    x=x+2
    y=y-2
    ind1=ind1+1
    ind2=ind2-1
  END DO

#if defined(SU3DBL)
  IF(S11q/=0.D0)THEN
  S2q=0.D0 ! S2 is S_2 is Eq.(26)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  IF(S11q/=0.Q0)THEN
  S2q=0.Q0
#endif
  a4=upbound*kkappa*(kkappa+upbound+2)/2+aux5
  ind1=aux3+alphaminS2
  ind2=aux4-alphaminS2
  IF(BTEST(alphaminS2,0))THEN ! if alphaminS2 is odd
    DO alpha=alphaminS2,alphamaxS2,2
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
      S2q=S2q+binom_quad(ind1+1)*binom_quad(ind2-1)*Sa(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
    ELSE
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
      S2q=S2q+binom_quad(ind1+1)*binom_quad(ind2-1)*Sa(a4-1)
    END IF
  ELSE
    DO alpha=alphaminS2,alphamaxS2,2
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
      S2q=S2q-binom_quad(ind1+1)*binom_quad(ind2-1)*Sa(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
      S2q=S2q-binom_quad(ind1+1)*binom_quad(ind2-1)*Sa(a4-1)
    ELSE
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa(a4)
    END IF
  END IF
#if defined(SU3DBL)
  IF(S2q/=0.D0)THEN
    coeffq=coeffq+binom(ind3)*S11q*S12q*S2q/DFLOAT(k2L+1)
#elif defined(SU3QUAD)
  IF(S2q/=0.Q0)THEN
    coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/QFLOAT(k2L+1)
#elif defined(SU3QUAD_GNU)
  IF(S2q/=0.Q0)THEN
    coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/REAL(k2L+1,16)
#endif
  END IF

  END IF
  END IF

  aux5=aux5-k2L
  k2L=k2L-1
  kkappa=kkappa-1
  aux2=aux2-p+gamma
  ind3=ind3+1
END DO

#if defined(SU3DBL)
IF(coeffq/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
IF(coeffq/=0.Q0)THEN
#endif
! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo:
! (2L+1) should not be squared!
  aux1=lambda+mu+1
  aux2=p+mu+1
  aux3=2*L*(L+1)
#if defined(SU3DBL)
  S2q=DFLOAT(2*L+1)*DSQRT((binom((lambdas+lambda)/2+p)/binom(aux3-Mp))*(binom(mu*(mu+1)/2+q)&
    /binom((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom(aux1*(aux1+1)/2+q)/binom(aux2*(aux2+1)/2+q))*binom(aux3-M))/(4.D0**DFLOAT(p))
#elif defined(SU3QUAD)
  S2q=QFLOAT(2*L+1)*QSQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
    /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.Q0**QFLOAT(p))
#elif defined(SU3QUAD_GNU)
  S2q=REAL(2*L+1,16)*SQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
    /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.Q0**REAL(p,16))
#endif
! S2q is C
  IF(BTEST(L-p,0))THEN
    coeff=-coeffq*S2q
  ELSE
    coeff=coeffq*S2q
  END IF
ELSE
  coeff=0.D0
  RETURN
END IF

END IF

! Now coeff is the coefficient for E=HW. For E=HW' there is additinal phase factor of (-1)^((lambda+M)/2) according to Eq.(33,6B), therefore:
IF(I/=J.AND.BTEST(lambdaM,0))coeff=-coeff

RETURN
END FUNCTION transformation_coeff
