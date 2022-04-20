!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ndsu3lib_wigner_su3so3.F90 -- module for SU(3)-SO(3) reduced Wigner coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE ndsu3lib_wigner_su3so3
USE ndsu3lib_wigner_canonical
IMPLICIT NONE
CONTAINS

FUNCTION inner_multiplicity(irrep,L) RESULT(kappamax) BIND(C)
!--------------------------------------------------------
! Inner multiplicity of L within SU(3) irrep (lambda,mu),
! where lambda=irrep%lambda, mu=irrep%mu
!--------------------------------------------------------
USE iso_c_binding
IMPLICIT NONE
TYPE(su3irrep),INTENT(IN) :: irrep
INTEGER(C_INT),INTENT(IN) :: L
INTEGER(C_INT) :: kappamax
kappamax=MAX(0,(irrep%lambda+irrep%mu+2-L)/2)-MAX(0,(irrep%lambda+1-L)/2)-MAX(0,(irrep%mu+1-L)/2)
END FUNCTION inner_multiplicity

FUNCTION Kmin(lambda,mu,L) RESULT(res)
!----------------------------------------------------------------
! The lowest K for given lambda,mu,L as given by Eq.(4a) in
! J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!----------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda,mu,L
INTEGER :: res
res=MAX(0,L-mu)
res=res+MOD(res+lambda,2) ! K and lambda must have the same parity.
IF(res==0)res=MOD(L+mu,2)*2 ! If K=0, L and mu must have the same parity.
END FUNCTION Kmin

SUBROUTINE transformation_coeff(I,J,irrep,epsilonx,Lambda2p,MLambda2px,M,L,Mp,coeff)
!--------------------------------------------------------------------------------------------------------------
! Calculates transformation coefficient for SU(3)-SO(3) reduction <G|(G_E)MLM'>, where
! |G>=|(lambda,mu)epsilon,Lambda',M'_Lambda>, using equations (26),(32) and (33,6B) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input parameters: (I,J)=(1,1) => E=HW, (I,J)=(1,0) => E=HW', (I,J)=(0,0) => E=LW, (I,J)=(0,1) => E=LW'
!                   irrep, epsilonx=epsilon, Lambda2p=2*Lambda', MLambda2px=2*M'_Lambda, Mp=M'
!
! lambda=irrep%lambda, mu=irrep%mu
!
! Note: There are 2 typos in equation (26):
!       1) In the overall factor C, the factor of (2L+1) should not be squared.
!       2) The second S_1 should be S_1(M_Lambda=Lambda,Lambda,N_Lambda,M) instead of
!          S_1(N_Lambda,Lambda,M_Lambda=Lambda,M).
!--------------------------------------------------------------------------------------------------------------
!#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!USE mpmodule
!#endif
IMPLICIT NONE
#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
REAL(KIND=8),INTENT(OUT) :: coeff
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
CLASS(*),TARGET,INTENT(OUT) :: coeff
TYPE(mp_real),POINTER :: point
#endif
TYPE(su3irrep),INTENT(IN) :: irrep
INTEGER,INTENT(IN) :: I,J,epsilonx,Lambda2p,MLambda2px,M,L,Mp
REAL(KIND=8) :: S11,S12,S2
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
REAL(KIND=16) :: coeffq,S11q,S12q,S2q
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
REAL(KIND=16) :: coeffq,S11q,S12q,S2q
TYPE(mp_real) :: S11mp,S12mp,S2mp
#endif
INTEGER :: lambda,mu,epsilon,MLambda2p,p,q,a4,LmM,LpM,LmMp,gama,NLambda2p,k2L,kkappa,alpha,&
           alphamin,alphamax,alphaminS2,alphamaxS2,LambdaM,MNLambda2p,LambdapMp,x,y,xm,xn,&
           aux1,ind1,aux2,ind2,aux3,ind3,aux4,xpys,lambdas,aux5,indicator

IF(I==1)THEN
  lambda=irrep%lambda
  mu=irrep%mu
  epsilon=epsilonx
  MLambda2p=MLambda2px
ELSE ! See Eq.(32) and the text below it.
  lambda=irrep%mu
  mu=irrep%lambda
  epsilon=-epsilonx
  MLambda2p=-MLambda2px
ENDIF

p=((2*(lambda-mu)-epsilon)/3+Lambda2p)/2
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

#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
SELECT TYPE(coeff)
TYPE IS (REAL(KIND=8))
#endif
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
IF(k2L<=17)THEN
  indicator=1
ELSE IF(k2L<60)THEN
  IF(mu<=(k2L+1)/2-lambda+k2L/3)THEN
    indicator=2
  ELSE
    indicator=1
  END IF
ELSE
  indicator=2
END IF

IF(indicator==1)THEN ! double precision used
#endif
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

        IF(S2/=0.D0)coeff=coeff+binom(ind3)*S11*S12*S2/DFLOAT(k2L+1)

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
      *(binom(aux1*(aux1+1)/2+q)/binom(aux2*(aux2+1)/2+q))*binom(aux3-M))/(4.D0**DFLOAT(p))
    ! S2 is C
    IF(BTEST(L-p,0))THEN
      coeff=-coeff*S2
    ELSE
      coeff=coeff*S2
    END IF
  ELSE
    RETURN
  END IF

#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
ELSE ! quad precision used
  coeffq=0.Q0
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
    ! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION: There is a typo: it should be
    ! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

    IF(S12q/=0.Q0)THEN
      alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
      alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
      S11q=0.Q0
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

      IF(S11q/=0.Q0)THEN
        S2q=0.Q0
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

        IF(S2q/=0.Q0)THEN
#endif
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_MP))
          coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/QFLOAT(k2L+1)
#elif (defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP_GNU))
          coeffq=coeffq+binom_quad(ind3)*S11q*S12q*S2q/REAL(k2L+1,16)
#endif
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
        END IF

      END IF
    END IF

    aux5=aux5-k2L
    k2L=k2L-1
    kkappa=kkappa-1
    aux2=aux2-p+gama
    ind3=ind3+1
  END DO

  IF(coeffq/=0.Q0)THEN
  ! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo: (2L+1) should not be squared!
    aux1=lambda+mu+1
    aux2=p+mu+1
    aux3=2*L*(L+1)
#endif
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_MP))
    S2q=QFLOAT(2*L+1)*QSQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
      /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
      *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.Q0**QFLOAT(p))
#elif (defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP_GNU))
    S2q=REAL(2*L+1,16)*SQRT((binom_quad((lambdas+lambda)/2+p)/binom_quad(aux3-Mp))*(binom_quad(mu*(mu+1)/2+q)&
      /binom_quad((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
      *(binom_quad(aux1*(aux1+1)/2+q)/binom_quad(aux2*(aux2+1)/2+q))*binom_quad(aux3-M))/(4.Q0**REAL(p,16))
#endif
    ! S2q is C
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU) || defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
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

#endif

! Now coeff is the coefficient for E=HW. For E=HW' there is additinal phase
! factor of (-1)^((lambda+M)/2) according to Eq.(33,6B), therefore:
IF(I/=J.AND.BTEST(LambdaM,0))coeff=-coeff

#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
CLASS DEFAULT ! CLASS IS (mp_real)
point=>coeff

point=mpreal(0.D0,nwds) ! coeff=mpreal(0.D0,nwds)
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
  ! S12 is the second S_1 in Eq.(26), where only alpha=0 contributes. CAUTION: There is a typo: it should be
  ! S_1(M_Lambda=Lambda,Lambda,N_Lambda,M), not S_1(N_Lambda,Lambda,M_Lambda=Lambda,M)

  IF(S12mp/=0.D0)THEN
    alphamin=MAX(0,-MNLambda2p) ! alphamin is the lower bound for alpha in the first S_1 in Eq.(26)
    alphamax=MIN(xm,xn) ! alphamin is the upper bound for alpha in the first S_1 in Eq.(26)
    S11mp=mpreal(0.D0,nwds)
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

    IF(S11mp/=0.D0)THEN
      S2mp=mpreal(0.D0,nwds)
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

      IF(S2mp/=0.D0)point=point+binom_mp(ind3)*S11mp*S12mp*S2mp/DFLOAT(k2L+1) ! coeff=coeff+binom_mp(ind3)*S11mp*S12mp*S2mp/DFLOAT(k2L+1)

    END IF
  END IF

  aux5=aux5-k2L
  k2L=k2L-1
  kkappa=kkappa-1
  aux2=aux2-p+gama
  ind3=ind3+1
END DO

IF(point/=0.D0)THEN ! IF(coeff/=0.D0)THEN
! Factor (2*L+1)/(2**p) appears in C in Eq.(26) as squared but that is a typo: (2L+1) should not be squared!
  aux1=lambda+mu+1
  aux2=p+mu+1
  aux3=2*L*(L+1)
  S2mp=DFLOAT(2*L+1)*SQRT((binom_mp((lambdas+lambda)/2+p)/binom_mp(aux3-Mp))*(binom_mp(mu*(mu+1)/2+q)&
    /binom_mp((Lambda2p*(Lambda2p+2)+MLambda2p)/2))&
    *(binom_mp(aux1*(aux1+1)/2+q)/binom_mp(aux2*(aux2+1)/2+q))*binom_mp(aux3-M))/(4.D0**DFLOAT(p))
  ! S2mp is C
  IF(BTEST(L-p,0))THEN
    p=-p*S2mp ! coeff=-coeff*S2mp
  ELSE
    p=p*S2mp ! coeff=coeff*S2mp
  END IF
ELSE
  RETURN
END IF
! Now coeff (or point) is the coefficient for E=HW. For E=HW' there is additinal phase
! factor of (-1)^((lambda+M)/2) according to Eq.(33,6B), therefore:
IF(I/=J.AND.BTEST(LambdaM,0))p=-p ! coeff=-coeff

END SELECT
#endif

END SUBROUTINE transformation_coeff

SUBROUTINE orthonormalization_matrix(I,J,irrep,L,kappamax,matrix)
!-------------------------------------------------------------------------------------------
! Calculates the bottom triangle, including the diagonal, of the orthonormalization matrix O
! for given lambda,mu,L using equations (6a),(6b),(6c),(27) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: I,J,irrep,L,kappamax
! Output argument: matrix
!
! (I,J)=(1,1) => E=HW, (I,J)=(1,0) => E=HW', (I,J)=(0,0) => E=LW, (I,J)=(0,1) => E=LW'
! lambda=irrep%lambda, mu=irrep%mu
! kappamax = the number of occurences of L in SU(3) irrep (lambda,mu)
! matrix(i,j) = O_ij
!
! Note: There is a typo in Eq.(6b): the square root should not be there!
!-------------------------------------------------------------------------------------------
!#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!USE mpmodule
!#endif
IMPLICIT NONE
#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: matrix
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
CLASS(*),TARGET,DIMENSION(:,:),INTENT(OUT) :: matrix
TYPE(mp_real),POINTER,DIMENSION(:,:) :: point
#endif
TYPE(su3irrep),INTENT(IN) :: irrep
INTEGER,INTENT(IN) :: I,J,L,kappamax
INTEGER :: ii,jj,k,epsilon,Ki,Kj,Kminm2,Lambda2,MLambda2
REAL(KIND=8) :: sum,coeff
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
TYPE(mp_real) :: summp,coeffmp
#endif

!INTERFACE
!  SUBROUTINE transformation_coeff(I,J,irrep,epsilonx,Lambda2p,MLambda2px,M,L,Mp,coeff)
!    IMPLICIT NONE
!#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
!    REAL(KIND=8),INTENT(OUT) :: coeff
!#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!    CLASS(*),INTENT(OUT) :: coeff
!#endif
!    TYPE(su3irrep),INTENT(IN) :: irrep
!    INTEGER,INTENT(IN) :: I,J,epsilonx,Lambda2p,MLambda2px,M,L,Mp
!  END SUBROUTINE transformation_coeff
!END INTERFACE

IF(I==1)THEN ! E=HW or HW'
  epsilon=-irrep%lambda-2*irrep%mu ! the HW epsilon
  Kminm2=Kmin(irrep%lambda,irrep%mu,L)-2 ! Kmin(lambda,mu,L) is the lowest K
  Lambda2=irrep%lambda ! Lambda2 is 2*Lambda
  MLambda2=irrep%lambda ! MLambda2 is 2*M_Lambda
ELSE ! E=LW or LW'
  epsilon=2*irrep%lambda+irrep%mu ! the LW epsilon
  Kminm2=Kmin(irrep%mu,irrep%lambda,L)-2 ! Kmin(mu,lambda,L) is the lowest K
  Lambda2=irrep%mu ! Lambda2 is 2*Lambda
  MLambda2=-irrep%mu ! MLambda2 is 2*M_Lambda
END IF
IF(I/=J)MLambda2=-MLambda2

Ki=Kminm2

#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
SELECT TYPE(matrix)
TYPE IS (REAL(KIND=8))
#endif

! First, the upper triangle is calculated including the diagonal.
DO ii=1,kappamax ! loop over columns
  Ki=Ki+2 ! Ki is K_i in Eq.(6a),(6b)
  Kj=Kminm2
  DO jj=1,ii-1 ! loop over rows
    Kj=Kj+2 ! Kj is K_j in Eq.(6b)
    sum=0.D0
    DO k=1,jj-1
      sum=sum+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6b)
    END DO
    CALL transformation_coeff(I,J,irrep,epsilon,Lambda2,MLambda2,Ki,L,Kj,coeff)
    matrix(jj,ii)=matrix(jj,jj)*(coeff-sum) ! Eq.(6b)
  END DO
  sum=0.D0
  DO jj=1,ii-1
    sum=sum+matrix(jj,ii)*matrix(jj,ii) ! sum is the sum in Eq.(6a)
  END DO
  CALL transformation_coeff(I,J,irrep,epsilon,Lambda2,MLambda2,Ki,L,Ki,coeff)
  matrix(ii,ii)=1.D0/DSQRT(coeff-sum) ! Eq.(6a)
END DO

! Now, the bottom triangle is calculated.
DO jj=1,kappamax-1 ! loop over columns
  DO ii=jj+1,kappamax ! loop over rows
    sum=0.D0
    DO k=jj,ii-1
      sum=sum+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6c)
    END DO
    matrix(ii,jj)=-matrix(ii,ii)*sum ! Eq.(6c)
  END DO
END DO

#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
CLASS DEFAULT ! CLASS IS (mp_real)
point=>matrix

! First, the upper triangle is calculated including the diagonal.
DO ii=1,kappamax ! loop over columns
  Ki=Ki+2 ! Ki is K_i in Eq.(6a),(6b)
  Kj=Kminm2
  DO jj=1,ii-1 ! loop over rows
    Kj=Kj+2 ! Kj is K_j in Eq.(6b)
    summp=mpreal(0.D0,nwds)
    DO k=1,jj-1
      summp=summp+point(k,jj)*point(k,ii) ! summp=summp+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6b)
    END DO
    CALL transformation_coeff(I,J,irrep,epsilon,Lambda2,MLambda2,Ki,L,Kj,coeffmp)
    point(jj,ii)=point(jj,jj)*(coeffmp-summp) ! matrix(jj,ii)=matrix(jj,jj)*(coeffmp-summp) ! Eq.(6b)
  END DO
  summp=mpreal(0.D0,nwds)
  DO jj=1,ii-1
    summp=summp+point(jj,ii)*point(jj,ii) ! summp=summp+matrix(jj,ii)*matrix(jj,ii) ! sum is the sum in Eq.(6a)
  END DO
  CALL transformation_coeff(I,J,irrep,epsilon,Lambda2,MLambda2,Ki,L,Ki,coeffmp)
  point(ii,ii)=1.D0/SQRT(coeffmp-summp) ! matrix(ii,ii)=1.D0/SQRT(coeffmp-summp) ! Eq.(6a)
END DO

! Now, the bottom triangle is calculated.
DO jj=1,kappamax-1 ! loop over columns
  DO ii=jj+1,kappamax ! loop over rows
    summp=mpreal(0.D0,nwds)
    DO k=jj,ii-1
      summp=summp+point(k,jj)*point(k,ii) ! summp=summp+matrix(k,jj)*matrix(k,ii) ! sum is the sum in Eq.(6c)
    END DO
    point(ii,jj)=-point(ii,ii)*summp ! matrix(ii,jj)=-matrix(ii,ii)*summp ! Eq.(6c)
  END DO
END DO

END SELECT
#endif

END SUBROUTINE orthonormalization_matrix

#if defined(NDSU3LIB_WSO3_GSL)
FUNCTION clebsch_gordan(j1,m1,j2,m2,j3,m3) RESULT(cg)
!----------------------------------------------------------------------------
! Calculates SU(2) Clebsch-Gordan coefficient <j1/2,m1/2,j2/2,m2/2|j3/2,m3/2>
! using GSL function gsl_sf_coupling_3j calculating 3j symbol.
!----------------------------------------------------------------------------
USE iso_c_binding
IMPLICIT NONE
INTERFACE
REAL(C_DOUBLE) FUNCTION gsl_sf_coupling_3j(l1,n1,l2,n2,l3,n3) BIND(C)
USE iso_c_binding
INTEGER(C_INT), VALUE :: l1,n1,l2,n2,l3,n3
END FUNCTION gsl_sf_coupling_3j
END INTERFACE
INTEGER(C_INT),INTENT(IN) :: j1,m1,j2,m2,j3,m3
REAL(C_DOUBLE) :: cg
INTEGER :: a
cg=gsl_sf_coupling_3j(j1,j2,j3,m1,m2,-m3)*DSQRT(DFLOAT(j3+1))
a=j1-j2+m3
IF((a/4)*4/=a)cg=-cg
RETURN
END FUNCTION clebsch_gordan
#endif

SUBROUTINE wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2,&
                         I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
!-------------------------------------------------------------------------------------------------------------------------------
! Calculates reduced SU(3)-SO(3) Wigner coefficients <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3 using equations (31),(25),(5) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2,I3,irrep3,L3,kappa3max,
!                  matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a
! Output argument: wigner_phys
!
! (I,J)=(1,0) for lambda<mu => E=HW', (I,J)=(0,1) for lambda>=mu => E=LW'
! lambda1=irrep1%lambda, mu1=irrep1%mu, lambda2=irrep2%lambda, mu2=irrep2%mu, lambda3=irrep3%lambda, mu3=irrep3%mu
! kappa1max = the number of occurences of L1 in SU(3) irrep (lambda1,mu1) (similarly kappa2max,kappa3max)
! matrix1(i,j) = element O_ij of the orthonormalization matrix for lambda1,mu1,L1 (similarly matrix2,matrix3)
! rhomax = the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! numb = the number of extremal canonical reduced Wigner coefficients
!        <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)HW> excluding the outer multiplicity
! wigner_can = array containing the extremal canonical reduced Wigner coefficients (see wigner_canonical_extremal.f90 for details)
! p1a = array containing the values of p1
! p2a = array containing the values of p2
! q2a = array containing the values of q2
! wigner_phys(kappa1,kappa2,kappa3,rho) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
!
! Note: triangular inequality for L1,L2,L3 is not checked.
!-------------------------------------------------------------------------------------------------------------------------------
!#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!USE mpmodule
!#endif
IMPLICIT NONE
#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
CLASS(*),TARGET,DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
TYPE(mp_real),POINTER,DIMENSION(:,:) :: point1,point2
#endif
REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix3
#if defined(NDSU3LIB_WSO3_WIGXJPF)
REAL(KIND=8),EXTERNAL :: fwig3jj
#endif
TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
INTEGER,INTENT(IN) :: I1,J1,L1,kappa1max,I2,J2,L2,kappa2max,I3,L3,kappa3max,rhomax,numb
INTEGER :: kappa1,kappa2,kappa3,rho,K1,K2,K3,K32,M1p,M2p,L12,L22,L32,i,M1pmin,epsilon1,epsilon2,epsilon3,&
           Lambda12,Lambda22,MLambda12,MLambda22,j,K1min,K2min,MLambda12min,MLambda32,Lambda32,K3minm2,MLambda12max,&
           p1,p2,q1,q2,hovadina,Lambda12pM1p,Lambda22pM2p,phase331a,phase332a,phase331b,phase332b,phase
INTEGER,DIMENSION(:),INTENT(IN) :: p1a,p2a,q2a
REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wigner_can
REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: transcoeff1,transcoeff2
REAL(KIND=8) :: cg1,cg2,fact,cg,coeff
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
TYPE(mp_real) :: coeffmp
TYPE(mp_real),ALLOCATABLE,DIMENSION(:) :: transcoeff1mp,transcoeff2mp
#endif

!INTERFACE
!  SUBROUTINE transformation_coeff(I,J,irrep,epsilonx,Lambda2p,MLambda2px,M,L,Mp,coeff)
!    IMPLICIT NONE
!#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
!    REAL(KIND=8),INTENT(OUT) :: coeff
!#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!    CLASS(*),INTENT(OUT) :: coeff
!#endif
!    TYPE(su3irrep),INTENT(IN) :: irrep
!    INTEGER,INTENT(IN) :: I,J,epsilonx,Lambda2p,MLambda2px,M,L,Mp
!  END SUBROUTINE transformation_coeff
!END INTERFACE

#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
ALLOCATE(transcoeff1(kappa1max),transcoeff2(kappa2max))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
ALLOCATE(transcoeff1(kappa1max),transcoeff2(kappa2max),transcoeff1mp(kappa1max),transcoeff2mp(kappa2max))
#endif

wigner_phys(1:kappa1max,1:kappa2max,1:kappa3max,1:rhomax)=0.D0

IF(I3==1)THEN ! G_3E=G_3HW'
  epsilon3=-irrep3%lambda-2*irrep3%mu
  K3minm2=Kmin(irrep3%lambda,irrep3%mu,L3)-2
  MLambda32=-irrep3%lambda ! MLambda32 is 2*M_Lambda_3
  Lambda32=irrep3%lambda ! Lambda32 is 2*Lambda_3
  hovadina=(2*(irrep1%lambda+irrep2%lambda+irrep3%mu)+irrep1%mu+irrep2%mu+irrep3%lambda)/3
ELSE ! G_3E=G_3LW'
  epsilon3=2*irrep3%lambda+irrep3%mu
  K3minm2=Kmin(irrep3%mu,irrep3%lambda,L3)-2
  MLambda32=irrep3%mu ! MLambda32 is 2*M_Lambda_3
  Lambda32=irrep3%mu ! Lambda32 is 2*Lambda_3
  hovadina=(2*(irrep1%mu+irrep2%mu+irrep3%lambda)+irrep1%lambda+irrep2%lambda+irrep3%mu)/3
END IF

L12=2*L1
L22=2*L2
L32=2*L3

IF(I1==1)THEN ! G_1E=G_1HW'
  K1min=Kmin(irrep1%lambda,irrep1%mu,L1)
  phase331a=2*(irrep1%lambda+2*irrep1%mu)+6*(L1+K1min)
ELSE ! G_1E=G_1LW'
  K1min=Kmin(irrep1%mu,irrep1%lambda,L1)
  phase331a=-2*(2*irrep1%lambda+irrep1%mu)+6*(L1+K1min)
END IF
IF(I2==1)THEN ! G_2E=G_2HW'
  K2min=Kmin(irrep2%lambda,irrep2%mu,L2)
  phase332a=2*(irrep2%lambda+2*irrep2%mu)+6*(L2+K2min)
ELSE ! G_2E=G_2LW'
  K2min=Kmin(irrep2%mu,irrep2%lambda,L2)
  phase332a=-2*(2*irrep2%lambda+irrep2%mu)+6*(L2+K2min)
END IF

DO i=1,numb ! sum over epsilon1,Lambda_1,epsilon2,Lambda_2
  p1=p1a(i)
  p2=p2a(i)
  q2=q2a(i)
  IF(I3==1)THEN
    q1=hovadina-p1-p2-q2
    epsilon2=2*irrep2%lambda+irrep2%mu-3*(p2+q2)
    Lambda12=irrep1%mu+p1-q1
    Lambda22=irrep2%mu+p2-q2
  ELSE
    q1=hovadina-p1-p2-q2
    epsilon2=-2*irrep2%mu-irrep2%lambda+3*(p2+q2)
    Lambda12=irrep1%lambda+p1-q1
    Lambda22=irrep2%lambda+p2-q2
  END IF
  epsilon1=epsilon3-epsilon2

  phase331b=(phase331a-epsilon1)/3
  phase332b=(phase332a-epsilon2)/3

  K3=K3minm2
  K32=2*K3
  MLambda12min=MAX(-Lambda12,MLambda32-Lambda22)
  MLambda12max=MIN(Lambda12,MLambda32+Lambda22)

  DO kappa3=1,kappa3max
    K3=K3+2 ! K3 is K_3
    K32=K32+4 ! K32 is 2*K_3
    M1pmin=MAX(-L1,K3-L2,-Lambda12,K3-Lambda22) ! M1pmin is the minimal M'_1
    IF(BTEST(Lambda12+M1pmin,0))M1pmin=M1pmin+1 ! M'_1 has to have the same parity as Lambda12 for <G_1|(G_1E)K1,L1,M'1> to be nonzero.

    M2p=K3-M1pmin ! M2p is M'_2
    Lambda12pM1p=Lambda12+M1pmin
    Lambda22pM2p=Lambda22+M2p
    DO M1p=M1pmin,MIN(L1,K3+L2,Lambda12,K3+Lambda22),2 ! M1p is M'_1.
    ! For <G_1|(G_1E)K1,L1,M'1> to be nonzero, M1p has to have the same parity
    ! as Lambda12 and abs(M1p) has to be less than or equal to Lambda12.
    ! This is taken care of in the lower and upper bounds on M1p.

    ! For <G_2|(G_2E)K2,L2,M'2> to be nonzero, M2p has to have the same parity
    ! as Lambda22 and abs(M2p) has to be less than or equal to Lambda22.
    ! This is taken care of in the lower and upper bounds on M1p.

    ! Lambda12pM1p=Lambda12+M1p
    ! Lambda22pM2p=Lambda22+M2p

#if defined(NDSU3LIB_WSO3_GSL)
      cg1=clebsch_gordan(L12,2*M1p,L22,2*M2p,L32,K32) ! cg1 is <L1 M'_1;L2 M'_2|L3 K3>
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
      cg1=fwig3jj(L12,L22,L32,2*M1p,2*M2p,-K32)*DSQRT(DFLOAT(L32+1))
      phase=L12-L22+K32
      IF((phase/4)*4/=phase)cg1=-cg1
#endif

      DO MLambda12=MLambda12min,MLambda12max,2 ! MLambda12 is 2*M_Lambda_1

        IF(MLambda12==0.AND.BTEST(Lambda12pM1p/2,0))CYCLE
        ! If M_Lambda1=0 and Lambda1+M1'/2 is odd, <G_1|(G_1E)K1,L1,M'1>=0. See Eq.(33,6A)

        IF(M1p==0)THEN
          IF(BTEST((phase331b+MLambda12)/2,0))CYCLE
          ! If M1'=0 and (phase331b+MLambda12)/2 is odd, <G_1|(G_1E)K1,L1,M'1>=0. Can be shown using Eq.(33).
          IF(BTEST((Lambda12+MLambda12)/2,0))CYCLE
          ! If M1'=0 and Lambda1+M_Lambda1 is odd, <G_1|(G_1E)K1,L1,M'1>=0. Can be shown using Eq.(33).
        END IF

        MLambda22=MLambda32-MLambda12 ! MLambda22 is 2*M_Lambda_2

        IF(MLambda22==0.AND.BTEST(Lambda22pM2p/2,0))CYCLE
        ! If M_Lambda2=0 and Lambda2+M2'/2 is odd, <G_2|(G_2E)K2,L2,M'2>=0. See Eq.(33,6A).

        IF(M2p==0)THEN
          IF(BTEST((phase332b+MLambda22)/2,0))CYCLE
          ! If M2'=0 and (phase332b+MLambda22)/2 is odd, <G_2|(G_2E)K2,L2,M'2>=0. Can be shown using Eq.(33).
          IF(BTEST((Lambda22+MLambda22)/2,0))CYCLE
          ! If M2'=0 and Lambda2+M_Lambda2 is odd, <G_2|(G_2E)K2,L2,M'2>=0. Can be shown using Eq.(33).
        END IF

        ! Calculation of <G_1|(G_1E)K1,L1,M'1> = transcoeff1(kappa1)
        K1=K1min
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
        SELECT TYPE(matrix1)
        TYPE IS (REAL(KIND=8))
#endif
        DO kappa1=1,kappa1max
          CALL transformation_coeff(I1,J1,irrep1,epsilon1,Lambda12,MLambda12,K1,L1,M1p,coeff)
          transcoeff1(kappa1)=coeff
          K1=K1+2
        END DO
        DO kappa1=kappa1max,1,-1
          transcoeff1(kappa1)=matrix1(kappa1,kappa1)*transcoeff1(kappa1)
          DO j=1,kappa1-1
            transcoeff1(kappa1)=transcoeff1(kappa1)+matrix1(kappa1,j)*transcoeff1(j)
          END DO
        END DO
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
        CLASS DEFAULT ! CLASS IS (mp_real)
        point1=>matrix1
        DO kappa1=1,kappa1max
          CALL transformation_coeff(I1,J1,irrep1,epsilon1,Lambda12,MLambda12,K1,L1,M1p,coeffmp)
          transcoeff1mp(kappa1)=coeffmp
          K1=K1+2
        END DO
        DO kappa1=kappa1max,1,-1
          !transcoeff1mp(kappa1)=matrix1(kappa1,kappa1)*transcoeff1mp(kappa1)
          transcoeff1mp(kappa1)=point1(kappa1,kappa1)*transcoeff1mp(kappa1)
          DO j=1,kappa1-1
            !transcoeff1mp(kappa1)=transcoeff1mp(kappa1)+matrix1(kappa1,j)*transcoeff1mp(j)
            transcoeff1mp(kappa1)=transcoeff1mp(kappa1)+point1(kappa1,j)*transcoeff1mp(j)
          END DO
          transcoeff1(kappa1)=transcoeff1mp(kappa1)
        END DO
        END SELECT
#endif

        ! Calculation of <G_2|(G_2E)K2,L2,M'2> = transcoeff2(kappa2)
        K2=K2min
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
        SELECT TYPE(matrix2)
        TYPE IS (REAL(KIND=8))
#endif
        DO kappa2=1,kappa2max
          CALL transformation_coeff(I2,J2,irrep2,epsilon2,Lambda22,MLambda22,K2,L2,M2p,coeff)
          transcoeff2(kappa2)=coeff
          K2=K2+2
        END DO
        DO kappa2=kappa2max,1,-1
          transcoeff2(kappa2)=matrix2(kappa2,kappa2)*transcoeff2(kappa2)
          DO j=1,kappa2-1
            transcoeff2(kappa2)=transcoeff2(kappa2)+matrix2(kappa2,j)*transcoeff2(j)
          END DO
        END DO
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
        CLASS DEFAULT ! CLASS IS (mp_real)
        point2=>matrix2
        DO kappa2=1,kappa2max
          CALL transformation_coeff(I2,J2,irrep2,epsilon2,Lambda22,MLambda22,K2,L2,M2p,coeffmp)
          transcoeff2mp(kappa2)=coeffmp
          K2=K2+2
        END DO
        DO kappa2=kappa2max,1,-1
          !transcoeff2mp(kappa2)=matrix2(kappa2,kappa2)*transcoeff2mp(kappa2)
          transcoeff2mp(kappa2)=point2(kappa2,kappa2)*transcoeff2mp(kappa2)
          DO j=1,kappa2-1
            !transcoeff2mp(kappa2)=transcoeff2mp(kappa2)+matrix2(kappa2,j)*transcoeff2mp(j)
            transcoeff2mp(kappa2)=transcoeff2mp(kappa2)+point2(kappa2,j)*transcoeff2mp(j)
          END DO
          transcoeff2(kappa2)=transcoeff2mp(kappa2)
        END DO
        END SELECT
#endif

#if defined(NDSU3LIB_WSO3_GSL)
        cg2=cg1*clebsch_gordan(Lambda12,-MLambda12,Lambda22,-MLambda22,Lambda32,-MLambda32)
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
        cg=fwig3jj(Lambda12,Lambda22,Lambda32,-MLambda12,-MLambda22,MLambda32)*DSQRT(DFLOAT(Lambda32+1))
        phase=Lambda12-Lambda22-MLambda32
        IF((phase/4)*4/=phase)cg=-cg
        cg2=cg1*cg
#endif

        DO kappa1=1,kappa1max
          fact=transcoeff1(kappa1)*cg2
          DO kappa2=1,kappa2max
            wigner_phys(kappa1,kappa2,kappa3,1:rhomax)=wigner_phys(kappa1,kappa2,kappa3,1:rhomax)&
                                           +transcoeff2(kappa2)*fact*wigner_can(p1,p2,q2,1:rhomax)
          END DO
        END DO

      END DO
      M2p=M2p-2
      Lambda12pM1p=Lambda12pM1p+2
      Lambda22pM2p=Lambda22pM2p-2
    END DO
  END DO
END DO

! Orthonormalization (w.r.t. K3)
DO kappa1=1,kappa1max
  DO kappa2=1,kappa2max
    DO kappa3=kappa3max,1,-1
      wigner_phys(kappa1,kappa2,kappa3,1:rhomax)=matrix3(kappa3,kappa3)*wigner_phys(kappa1,kappa2,kappa3,1:rhomax)
      DO j=kappa3-1,1,-1
        wigner_phys(kappa1,kappa2,kappa3,1:rhomax)=wigner_phys(kappa1,kappa2,kappa3,1:rhomax)&
                                              +matrix3(kappa3,j)*wigner_phys(kappa1,kappa2,j,1:rhomax)
      END DO
    END DO
  END DO
END DO

#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
DEALLOCATE(transcoeff1,transcoeff2)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
DEALLOCATE(transcoeff1,transcoeff2,transcoeff1mp,transcoeff2mp)
#endif

END SUBROUTINE wigner_su3so3

SUBROUTINE calculate_wigner_su3so3(irrep1,L1,kappa1max,irrep2,L2,kappa2max,irrep3,L3,kappa3max,rhomax,wigner_phys)
!-----------------------------------------------------------------------------------------------------------------
! Calculates SU(3)-SO(3) reduced Wigner coefficients for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3.
!
! Input arguments: irrep1,L1,irrep2,L2,irrep3,L3
!                  rhomax = multipicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
!                  kappa1max = number of occurences of L1 in (lambda1,mu1)
!                  kappa2max,kappa3max
! Output argument: wigner_phys(kappa1,kappa2,kappa3,rho) = the Wigner coefficient
!
! lambda1=irrep1%lambda, mu1=irrep1%mu, lambda2=irrep2%lambda, mu2=irrep2%mu, lambda3=irrep3%lambda, mu3=irrep3%mu
!-----------------------------------------------------------------------------------------------------------------
!#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!USE mpmodule
!#endif
IMPLICIT NONE
TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
INTEGER,INTENT(IN) :: L1,kappa1max,L2,kappa2max,L3,kappa3max,rhomax
INTEGER :: I1,J1,I2,J2,I3,J3,numb,indicator1,indicator2,sum1,sum2,pqdim
REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys ! The indeces are kappa1,kappa2,kappa3,rho, respectively.
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_can
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: matrix1,matrix2,matrix3
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
TYPE(mp_real),ALLOCATABLE,DIMENSION(:,:) :: matrix1mp,matrix2mp
#endif
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1a,p2a,q2a

!INTERFACE
!  SUBROUTINE wigner_canonical_extremal(irrep1,irrep2,irrep3,I3,rhomax,i2,wigner,p1a,p2a,q2a)
!    IMPLICIT NONE
!    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
!    INTEGER,INTENT(IN) :: I3,rhomax
!    INTEGER,INTENT(OUT) :: i2
!    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a
!    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
!  END SUBROUTINE wigner_canonical_extremal
!  SUBROUTINE wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2,&
!                           I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
!    IMPLICIT NONE
!#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
!    REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
!#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!    CLASS(*),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2
!#endif
!    REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix3
!    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
!    INTEGER,INTENT(IN) :: I1,J1,L1,kappa1max,I2,J2,L2,kappa2max,I3,L3,kappa3max,rhomax,numb
!    INTEGER,DIMENSION(:),INTENT(IN) :: p1a,p2a,q2a
!    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wigner_can
!    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
!  END SUBROUTINE wigner_su3so3
!  SUBROUTINE orthonormalization_matrix(I,J,irrep,L,kappamax,matrix)
!    IMPLICIT NONE
!#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
!    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: matrix
!#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
!    CLASS(*),DIMENSION(:,:),INTENT(OUT) :: matrix
!#endif
!    TYPE(su3irrep),INTENT(IN) :: irrep
!    INTEGER,INTENT(IN) :: I,J,L,kappamax
!  END SUBROUTINE orthonormalization_matrix
!END INTERFACE

I2=MAX(irrep2%lambda,irrep2%mu)
pqdim=(MAX(irrep1%lambda,irrep1%mu)+1)*((I2+1)**2)
ALLOCATE(wigner_can(0:MAX(irrep1%lambda,irrep1%mu),0:I2,0:I2,1:rhomax),&
         matrix3(kappa3max,kappa3max),p1a(pqdim),p2a(pqdim),q2a(pqdim))

IF(irrep1%lambda<irrep1%mu)THEN
  I1=1
  J1=0
ELSE
  I1=0
  J1=1
END IF
IF(irrep2%lambda<irrep2%mu)THEN
  I2=1
  J2=0
ELSE
  I2=0
  J2=1
END IF
IF(irrep3%lambda<irrep3%mu)THEN
  I3=1
  J3=0
ELSE
  I3=0
  J3=1
END IF

CALL wigner_canonical_extremal(irrep1,irrep2,irrep3,I3,rhomax,numb,wigner_can,p1a,p2a,q2a)
CALL orthonormalization_matrix(I3,J3,irrep3,L3,kappa3max,matrix3)

#if (defined(NDSU3LIB_DBL) || defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))

ALLOCATE(matrix1(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
CALL orthonormalization_matrix(I1,J1,irrep1,L1,kappa1max,matrix1)
CALL orthonormalization_matrix(I2,J2,irrep2,L2,kappa2max,matrix2)
CALL wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2,&
       I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
DEALLOCATE(wigner_can,matrix1,matrix2,matrix3,p1a,p2a,q2a)

#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))

sum1=irrep1%lambda+irrep1%mu+L1+2
IF(sum1<62)THEN
  indicator1=1
ELSE IF(sum1<=77)THEN
  IF(irrep1%lambda<=sum1/3.AND.irrep1%mu<=(sum1-1)/2-3*irrep1%lambda/4+1)THEN
    indicator1=2
  ELSE
    indicator1=1
  END IF
ELSE
  IF(irrep1%lambda<=sum1/3)THEN
    IF(irrep1%mu<=(sum1-1)/2-irrep1%lambda/2+MAX(0,(sum1-82)/3))THEN
      indicator1=2
    ELSE
      indicator1=1
    END IF
  ELSE IF(irrep1%mu<=(sum1-1)/2-sum1/6-2*(irrep1%lambda-sum1/3)+MAX(0,2*(sum1-86)/3))THEN
    indicator1=2
  ELSE
    indicator1=1
  END IF
END IF

sum2=irrep2%lambda+irrep2%mu+L2+2
IF(sum2<62)THEN
  indicator2=1
ELSE IF(sum2<=77)THEN
  IF(irrep2%lambda<=sum2/3.AND.irrep2%mu<=(sum2-1)/2-3*irrep2%lambda/4+1)THEN
    indicator2=2
  ELSE
    indicator2=1
  END IF
ELSE
  IF(irrep2%lambda<=sum2/3)THEN
    IF(irrep2%mu<=(sum2-1)/2-irrep2%lambda/2+MAX(0,(sum2-82)/3))THEN
      indicator2=2
    ELSE
      indicator2=1
    END IF
  ELSE IF(irrep2%mu<=(sum2-1)/2-sum2/6-2*(irrep2%lambda-sum2/3)+MAX(0,2*(sum2-86)/3))THEN
    indicator2=2
  ELSE
    indicator2=1
  END IF
END IF

IF(indicator1==1.AND.indicator2==1)THEN
  ALLOCATE(matrix1(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,irrep1,L1,kappa1max,matrix1)
  CALL orthonormalization_matrix(I2,J2,irrep2,L2,kappa2max,matrix2)
  CALL wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2,&
         I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1,matrix2,matrix3,p1a,p2a,q2a)
ELSE IF(indicator1==2.AND.indicator2==1)THEN
  ALLOCATE(matrix1mp(kappa1max,kappa1max),matrix2(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,irrep1,L1,kappa1max,matrix1mp)
  CALL orthonormalization_matrix(I2,J2,irrep2,L2,kappa2max,matrix2)
  CALL wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1mp,I2,J2,irrep2,L2,kappa2max,matrix2,&
         I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1mp,matrix2,matrix3,p1a,p2a,q2a)
ELSE IF(indicator1==1.AND.indicator2==2)THEN
  ALLOCATE(matrix1(kappa1max,kappa1max),matrix2mp(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,irrep1,L1,kappa1max,matrix1)
  CALL orthonormalization_matrix(I2,J2,irrep2,L2,kappa2max,matrix2mp)
  CALL wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1,I2,J2,irrep2,L2,kappa2max,matrix2mp,&
         I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1,matrix2mp,matrix3,p1a,p2a,q2a)
ELSE
  ALLOCATE(matrix1mp(kappa1max,kappa1max),matrix2mp(kappa2max,kappa2max))
  CALL orthonormalization_matrix(I1,J1,irrep1,L1,kappa1max,matrix1mp)
  CALL orthonormalization_matrix(I2,J2,irrep2,L2,kappa2max,matrix2mp)
  CALL wigner_su3so3(I1,J1,irrep1,L1,kappa1max,matrix1mp,I2,J2,irrep2,L2,kappa2max,matrix2mp,&
         I3,irrep3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
  DEALLOCATE(wigner_can,matrix1mp,matrix2mp,matrix3,p1a,p2a,q2a)
END IF

#endif

END SUBROUTINE calculate_wigner_su3so3

SUBROUTINE wigner_su3so3_wrapper(irrep1,L1,irrep2,L2,irrep3,L3,&
                                 kappa1max,kappa2max,kappa3max,rhomax,dimen,wigner_phys_block) BIND(C)
!----------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine calculating reduced SU(3)-SO(3) Wigner coefficients
! <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
!
! Input arguments: irrep1,L1,irrep2,L2,irrep3,L3,kappa1max,kappa2max,kappa3max,rhomax,dimen
! Output argument: wigner_phys_block
!
! lambda1=irrep1%lambda, mu1=irrep1%mu, lambda2=irrep2%lambda, mu2=irrep2%mu, lambda3=irrep3%lambda, mu3=irrep3%mu
! kappa1max is the inner multiplicity of L1 within (lambda1,mu1).
! kappa2max is the inner multiplicity of L2 within (lambda2,mu2).
! kappa3max is the inner multiplicity of L3 within (lambda3,mu3).
! rhomax is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3).
! dimen is the size of the array wigner_phys_block. It must be at least kappa1max*kaapa2max*kaapa3max*rhomax.
! wigner_phys_block(ind) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! ind = kappa1+kappa1max*(kappa2-1)+kappa1max*kappa2max*(kappa3-1)+kappa1max*kappa2max*kappa3max*(rho-1)
!-----------------------------------------------------------------------------------------------------------------
USE iso_c_binding
IMPLICIT NONE
TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
INTEGER(C_INT),INTENT(IN) :: L1,L2,L3,kappa1max,kappa2max,kappa3max,rhomax,dimen
REAL(C_DOUBLE),DIMENSION(dimen),INTENT(OUT) :: wigner_phys_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_phys
INTEGER :: ind,rho,kappa1,kappa2,kappa3
!INTERFACE
!  SUBROUTINE calculate_wigner_su3so3(irrep1,L1,kappa1max,irrep2,L2,kappa2max,irrep3,L3,kappa3max,rhomax,wigner_phys)
!    IMPLICIT NONE
!    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
!    INTEGER,INTENT(IN) :: L1,kappa1max,L2,kappa2max,L3,kappa3max,rhomax
!    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
!  END SUBROUTINE calculate_wigner_su3so3
!END INTERFACE
ALLOCATE(wigner_phys(kappa1max,kappa2max,kappa3max,rhomax))
CALL calculate_wigner_su3so3(irrep1,L1,kappa1max,irrep2,L2,kappa2max,irrep3,L3,kappa3max,rhomax,wigner_phys)
ind=0
DO rho=1,rhomax
  DO kappa3=1,kappa3max
    DO kappa2=1,kappa2max
      DO kappa1=1,kappa1max
        ind=ind+1
        wigner_phys_block(ind)=wigner_phys(kappa1,kappa2,kappa3,rho)
      END DO
    END DO
  END DO
END DO
DEALLOCATE(wigner_phys)
END SUBROUTINE wigner_su3so3_wrapper

END MODULE ndsu3lib_wigner_su3so3
