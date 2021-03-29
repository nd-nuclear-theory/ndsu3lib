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
USE mpmodule
USE precision_level
IMPLICIT NONE
INTEGER,INTENT(IN) :: I,J,lambdax,mux,epsilonx,Lambda2p,MLambda2px,M,L,Mp
REAL(KIND=8) :: coeff,S11,S12,S2
#if defined(SU3DBL)
REAL(KIND=8) :: coeffq,S11q,S12q,S2q,coeffmp,S11mp,S12mp,S2mp
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=16) :: coeffq,S11q,S12q,S2q,coeffmp,S11mp,S12mp,S2mp
#elif (defined(SU3MP) || defined(SU3MP_GNU))
REAL(KIND=16) :: coeffq,S11q,S12q,S2q
TYPE(mp_real) :: coeffmp,S11mp,S12mp,S2mp
#endif
INTEGER :: lambda,mu,epsilon,MLambda2p,p,q,a4,LmM,LpM,LmMp,gama,NLambda2p,k2L,kkappa,alpha,&
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

DO
  IF(MAX(Lambda2p,2*L,lambda+mu+1)>upbound_binom)THEN
    CALL reallocate_binom(50)
  ELSE
    EXIT
  END IF
END DO
DO
  IF(MAX(lambda,2*MIN(xm,xn+p+1)-xm+p)>upbound_I)THEN
    CALL reallocate_I(50)
  ELSE
    EXIT
  END IF
END DO
DO
  IF(kkappa>upbound_S.OR.k2L-kkappa>upbound_S)THEN
    CALL reallocate_S(50)
  ELSE
    EXIT
  END IF
END DO

IF(lambda+mu+L<=17)THEN ! Double precision is used.
coeff=0.D0

DO gama=0,p
  MNLambda2p=MNLambda2p-1
  xn=xn+1
  aux1=aux1+xn

  IF(lambda-gama>=gama)THEN ! S12=I(lambda-gama,gama,LambdaM)
    S12=Ia(((lambda-gama)**3+lambdas+gama)/2+LambdaM)
  ELSE IF(BTEST(LambdaM,0))THEN ! S12=(-1)^(LambdaM)*I(lambda-gama,gama,LambdaM)
    S12=-Ia((gama**3+lambdas+lambda-gama)/2+LambdaM)
  ELSE
    S12=Ia((gama**3+lambdas+lambda-gama)/2+LambdaM)
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
  a4=upbound_S*kkappa*(kkappa+upbound_S+2)/2+aux5
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
  aux2=aux2-p+gama
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

ELSE IF(lambda+m+L<600)THEN ! If 17<lambda+mu+L<60, quad precision is used.

#if defined(SU3DBL)
  coeffq=0.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU) || defined(SU3MP) || defined(SU3MP_GNU))
  coeffq=0.Q0
#endif

DO gama=0,p
  MNLambda2p=MNLambda2p-1
  xn=xn+1
  aux1=aux1+xn

  IF(lambda-gama>=gama)THEN ! S12=I(lambda-gama,gama,LambdaM)
    S12q=Ia_quad(((lambda-gama)**3+lambdas+gama)/2+LambdaM)
  ELSE IF(BTEST(LambdaM,0))THEN ! S12=(-1)^(LambdaM)*I(lambda-gama,gama,LambdaM)
    S12q=-Ia_quad((gama**3+lambdas+lambda-gama)/2+LambdaM)
  ELSE
    S12q=Ia_quad((gama**3+lambdas+lambda-gama)/2+LambdaM)
  END IF
! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION:
! There is a typo: it should be
! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not
! S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

#if defined(SU3DBL)
  IF(S12q/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU) || defined(SU3MP) || defined(SU3MP_GNU))
  IF(S12q/=0.Q0)THEN
#endif

  alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
  alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
#if defined(SU3DBL)
  S11q=0.D0 ! S11 is the first S_1 in Eq.(26)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU) || defined(SU3MP) || defined(SU3MP_GNU))
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
        S11q=S11q+binom_quad(ind1)*binom_quad(ind2)*Ia_quad((x**3+xpys+y)/2+LambdapMp)
      ELSE IF(BTEST(LambdapMp,0))THEN
        S11q=S11q-binom_quad(ind1)*binom_quad(ind2)*Ia_quad((y**3+xpys+x)/2+LambdapMp)
      ELSE
        S11q=S11q+binom_quad(ind1)*binom_quad(ind2)*Ia_quad((y**3+xpys+x)/2+LambdapMp)
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
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU) || defined(SU3MP) || defined(SU3MP_GNU))
  IF(S11q/=0.Q0)THEN
  S2q=0.Q0
#endif

  a4=upbound_S*kkappa*(kkappa+upbound_S+2)/2+aux5
  ind1=aux3+alphaminS2
  ind2=aux4-alphaminS2
  IF(BTEST(alphaminS2,0))THEN ! if alphaminS2 is odd
    DO alpha=alphaminS2,alphamaxS2,2
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
      S2q=S2q+binom_quad(ind1+1)*binom_quad(ind2-1)*Sa_quad(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
    ELSE
      S2q=S2q-binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
      S2q=S2q+binom_quad(ind1+1)*binom_quad(ind2-1)*Sa_quad(a4-1)
    END IF
  ELSE
    DO alpha=alphaminS2,alphamaxS2,2
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
      S2q=S2q-binom_quad(ind1+1)*binom_quad(ind2-1)*Sa_quad(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
      S2q=S2q-binom_quad(ind1+1)*binom_quad(ind2-1)*Sa_quad(a4-1)
    ELSE
      S2q=S2q+binom_quad(ind1)*binom_quad(ind2)*Sa_quad(a4)
    END IF
  END IF

!do alpha=alphaminS2,alphamaxS2+2
!  if(2*(alpha/2)==alpha)then
!    S2q=S2q+binom_quad((L-M)*(L-M+1)/2+alpha)*binom_quad((L+M)*(L+M+1)/2+L-Mp-alpha)&
!!        *Sa((upbound_S*(p+q-gama)*(p+q-gama+upbound_S+2)+(p+q-gama+lambda+mu+L-p-q)*(p+q-gama+lambda+mu+L-p-q+1))/2&
!!            +(lambda-M)/2+mu+L-q-(Lambda2p+Mp)/2-alpha)
!        *Sfu(p+q-gama,lambda+mu+L-p-q,(lambda-M)/2+mu+L-q-(Lambda2p+Mp)/2-alpha)
!  else
!    S2q=S2q-binom_quad((L-M)*(L-M+1)/2+alpha)*binom_quad((L+M)*(L+M+1)/2+L-Mp-alpha)&
!!        *Sa((upbound_S*(p+q-gama)*(p+q-gama+upbound_S+2)+(p+q-gama+lambda+mu+L-p-q)*(p+q-gama+lambda+mu+L-p-q+1))/2&
!!            +(lambda-M)/2+mu+L-q-(Lambda2p+Mp)/2-alpha)
!        *Sfu(p+q-gama,lambda+mu+L-p-q,(lambda-M)/2+mu+L-q-(Lambda2p+Mp)/2-alpha)
!  end if
!end do

#if defined(SU3DBL)
  IF(S2q/=0.D0)THEN
    coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/DFLOAT(k2L+1)
#elif (defined(SU3QUAD) || defined(SU3MP))
  IF(S2q/=0.Q0)THEN
    coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/QFLOAT(k2L+1)
#elif (defined(SU3QUAD_GNU) || defined(SU3MP_GNU))
  IF(S2q/=0.Q0)THEN
    coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/REAL(k2L+1,16)
#endif
  END IF

  END IF
  END IF

  aux5=aux5-k2L
  k2L=k2L-1
  kkappa=kkappa-1
  aux2=aux2-p+gama
  ind3=ind3+1
END DO

#if defined(SU3DBL)
IF(coeffq/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU) || defined(SU3MP) || defined(SU3MP_GNU))
IF(coeffq/=0.Q0)THEN
#endif
! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo:
! (2L+1) should not be squared!
  aux1=lambda+mu+1
  aux2=p+mu+1
  aux3=2*L*(L+1)
#if defined(SU3DBL)
  S2q=DFLOAT(2*L+1)*DSQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
    /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.D0**DFLOAT(p))
#elif (defined(SU3QUAD) || defined(SU3MP))
  S2q=QFLOAT(2*L+1)*QSQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
    /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.Q0**QFLOAT(p))
#elif (defined(SU3QUAD_GNU) || defined(SU3MP_GNU))
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

ELSE ! If lambda+mu+L>..., multi precision is used.

#if defined(SU3DBL)
  coeffmp=0.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  coeffmp=0.Q0
#elif(defined(SU3MP) || defined(SU3MP_GNU))
  coeffmp=mpreal(0.D0,nwds)
#endif

DO gama=0,p
  MNLambda2p=MNLambda2p-1
  xn=xn+1
  aux1=aux1+xn

  IF(lambda-gama>=gama)THEN ! S12=I(lambda-gama,gama,LambdaM)
    S12mp=Ia_mp(((lambda-gama)**3+lambdas+gama)/2+LambdaM)
  ELSE IF(BTEST(LambdaM,0))THEN ! S12=(-1)^(LambdaM)*I(lambda-gama,gama,LambdaM)
    S12mp=-Ia_mp((gama**3+lambdas+lambda-gama)/2+LambdaM)
  ELSE
    S12mp=Ia_mp((gama**3+lambdas+lambda-gama)/2+LambdaM)
  END IF
! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION:
! There is a typo: it should be
! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not
! S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

#if (defined(SU3DBL) || defined(SU3MP) || defined(SU3MP_GNU))
  IF(S12mp/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  IF(S12mp/=0.Q0)THEN
#endif

  alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
  alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
#if defined(SU3DBL)
  S11mp=0.D0 ! S11 is the first S_1 in Eq.(26)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  S11mp=0.Q0
#elif(defined(SU3MP) || defined(SU3MP_GNU))
  S11mp=mpreal(0.D0,nwds)
#endif
  x=2*alphamin+MNLambda2p
  y=Lambda2p-x
  xpys=(x+y)**2
  ind1=aux1+alphamin
  ind2=aux2-alphamin
  DO alpha=alphamin,alphamax
    IF(MAX(0,LambdapMp-x)<=MIN(y,LambdapMp))THEN
      IF(x>=y)THEN
        S11mp=S11mp+binom_mp(ind1)*binom_mp(ind2)*Ia_mp((x**3+xpys+y)/2+LambdapMp)
      ELSE IF(BTEST(LambdapMp,0))THEN
        S11mp=S11mp-binom_mp(ind1)*binom_mp(ind2)*Ia_mp((y**3+xpys+x)/2+LambdapMp)
      ELSE
        S11mp=S11mp+binom_mp(ind1)*binom_mp(ind2)*Ia_mp((y**3+xpys+x)/2+LambdapMp)
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
  IF(S11mp/=0.D0)THEN
  S2mp=0.D0 ! S2 is S_2 is Eq.(26)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  IF(S11mp/=0.Q0)THEN
  S2mp=0.Q0
#elif(defined(SU3MP) || defined(SU3MP_GNU))
  IF(S11mp/=0.D0)THEN
  S2mp=mpreal(0.D0,nwds)
#endif
  a4=upbound_S*kkappa*(kkappa+upbound_S+2)/2+aux5
  ind1=aux3+alphaminS2
  ind2=aux4-alphaminS2
  IF(BTEST(alphaminS2,0))THEN ! if alphaminS2 is odd
    DO alpha=alphaminS2,alphamaxS2,2
      S2mp=S2mp-binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
      S2mp=S2mp+binom_mp(ind1+1)*binom_mp(ind2-1)*Sa_mp(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2mp=S2mp-binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
    ELSE
      S2mp=S2mp-binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
      S2mp=S2mp+binom_mp(ind1+1)*binom_mp(ind2-1)*Sa_mp(a4-1)
    END IF
  ELSE
    DO alpha=alphaminS2,alphamaxS2,2
      S2mp=S2mp+binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
      S2mp=S2mp-binom_mp(ind1+1)*binom_mp(ind2-1)*Sa_mp(a4-1)
      a4=a4-2
      ind1=ind1+2
      ind2=ind2-2
    END DO
    IF(BTEST(alphamaxS2,0))THEN
      S2mp=S2mp+binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
      S2mp=S2mp-binom_mp(ind1+1)*binom_mp(ind2-1)*Sa_mp(a4-1)
    ELSE
      S2mp=S2mp+binom_mp(ind1)*binom_mp(ind2)*Sa_mp(a4)
    END IF
  END IF

#if (defined(SU3DBL) || defined(SU3MP) || defined(SU3MP_GNU))
  IF(S2mp/=0.D0)THEN
    coeffmp=coeffmp+binom_mp(ind3)*S11mp*S12mp*S2mp/DFLOAT(k2L+1)
#elif defined(SU3QUAD)
  IF(S2mp/=0.Q0)THEN
    coeffmp=coeffmp+binom_mp(ind3)*S11mp*S12mp*S2mp/QFLOAT(k2L+1)
#elif defined(SU3QUAD_GNU)
  IF(S2mp/=0.Q0)THEN
    coeffmp=coeffmp+binom_mp(ind3)*S11mp*S12mp*S2mp/REAL(k2L+1,16)
#endif
  END IF

  END IF
  END IF

  aux5=aux5-k2L
  k2L=k2L-1
  kkappa=kkappa-1
  aux2=aux2-p+gama
  ind3=ind3+1
END DO

#if (defined(SU3DBL) || defined(SU3MP) || defined(SU3MP_GNU))
IF(coeffmp/=0.D0)THEN
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
IF(coeffmp/=0.Q0)THEN
#endif
! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo:
! (2L+1) should not be squared!
  aux1=lambda+mu+1
  aux2=p+mu+1
  aux3=2*L*(L+1)
#if defined(SU3DBL)
  S2mp=DFLOAT(2*L+1)*DSQRT((binom_mp((lambdas+lambda)/2+p)/binom_mp(aux3-Mp))*(binom_mp(mu*(mu+1)/2+q)&
    /binom_mp((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_mp(aux1*(aux1+1)/2+q)/binom_mp(aux2*(aux2+1)/2+q))*binom_mp(aux3-M))/(4.D0**DFLOAT(p))
#elif defined(SU3QUAD)
  S2mp=QFLOAT(2*L+1)*QSQRT((binom_mp((lambdas+lambda)/2+p)/binom_mp(aux3-Mp))*(binom_mp(mu*(mu+1)/2+q)&
    /binom_mp((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_mp(aux1*(aux1+1)/2+q)/binom_mp(aux2*(aux2+1)/2+q))*binom_mp(aux3-M))/(4.Q0**QFLOAT(p))
#elif defined(SU3QUAD_GNU)
  S2mp=REAL(2*L+1,16)*SQRT((binom_mp((lambdas+lambda)/2+p)/binom_mp(aux3-Mp))*(binom_mp(mu*(mu+1)/2+q)&
    /binom_mp((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_mp(aux1*(aux1+1)/2+q)/binom_mp(aux2*(aux2+1)/2+q))*binom_mp(aux3-M))/(4.Q0**REAL(p,16))
#elif (defined(SU3MP) || defined(SU3MP_GNU))
  S2mp=DFLOAT(2*L+1)*SQRT((binom_mp((lambdas+lambda)/2+p)/binom_mp(aux3-Mp))*(binom_mp(mu*(mu+1)/2+q)&
    /binom_mp((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_mp(aux1*(aux1+1)/2+q)/binom_mp(aux2*(aux2+1)/2+q))*binom_mp(aux3-M))/(4.D0**DFLOAT(p))
#endif
! S2q is C
  IF(BTEST(L-p,0))THEN
    coeff=-coeffmp*S2mp
  ELSE
    coeff=coeffmp*S2mp
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
