SUBROUTINE wigner_canonical_extremal(lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax,&
                                     i2,wignerfinal,Lambda12a,Lambda22a,epsilon2a)
!--------------------------------------------------------------------------------------------------------------------------------
! Calculates the extremal reduced SU(3)-SU(2)xU(1) Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho for given lambda1,mu1,lambda2,mu2,lambda3,mu3
! using Eq.(20),(21),(17),(18),(8),(1),(32),(35,2B) in Ref.[1]. The phase convention of Ref.[2] is adopted.
!
! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!             [2] K.T.Hecht, Nucl.Phys. 62 (1965) 1
!
! Input arguments: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
! Output arguments: i2,wignerfinal,Lambda12a,Lambda22a,epsilon2a
!
! lambda1x=lambda1
! mu1x=mu1
! lambda2x=lambda2
! mu2x=mu2
! lambda3x=lambda3
! mu3x=mu3
! I3=1 for E=HW, I3=0 for E=LW
! rhomax=the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! 
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho=wignerfinal(2*Lambda1,k,2*Lambda2,rho)
!   where k=(epsilon2+lambda2+2*mu2)/3 for E=HW (I3=1), k=(-epsilon2+mu2+2*lambda2)/3 for E=LW (I3=0)
!         2*Lambda1=Lambda12a(i)
!         2*Lambda2=Lambda22a(i)
!         epsilon2=epsilon2a(i)
!           where 1<=i<=i2
!-------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,p2,q2,Lambda22,i1,j2,j1,p2tilde,q2tilde,i,j,n,a,b,c,d,i2,&
           epsilon2,steps21,eps2,Lambda12,epsilon1,Sq2,Rp2,eta,p1,q1,rho,rhomax,i4,Lambda22max,ABCD,&
           lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,phiprhomax
INTEGER(KIND=8) :: F,G,H,prod
INTEGER,DIMENSION(1) :: Lambda12a,Lambda22a,epsilon2a ! Dimension is (lambda1+mu1+1)*((lambda2+mu2+1)^2)

!INTEGER :: k1,k2,k3,k4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VYMAZAT
!INTEGER,DIMENSION((lambda1x+mu1x+1)*((lambda2x+mu2x+1)**2),rhomax) :: Lambda12at,Lambda22at,epsilon2at !!!!!!!!!!! VYMAZAT

! BACHA, PODLA PATRIKA POLIA, KTORE SU LEN V PODPROGRAME, BY MALI BYT ALOKOVATELNE A NIE AUTOMATICKE ABY NEBOL SEGMENTATION FAULT

REAL(KIND=8),DIMENSION(0:lambda1x+mu1x,0:lambda2x+mu2x) :: wigner ! Dimensions should be (0:lambda1+mu1,0:lambda2+mu2)
REAL(KIND=8),DIMENSION(0:20,0:20,0:20,1:9) :: wignerfinal ! First index is 2*Lambda1, second is (epsilon2-epsilon2HW)/3, third is 2*Lambda2, fourth is rho
! Preco tu nefungure assumed-shape array???
REAL(KIND=8),DIMENSION(rhomax) :: norm ! dimension is rhomax
REAL(KIND=8),DIMENSION(2:rhomax,rhomax-1) :: scalprod ! dimension is (2:rhomax,rhomax-1)

!rhomax=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda3,mu3)

IF(rhomax==0)THEN
  PRINT*,"Coupling not allowed"
  RETURN
END IF

IF(I3==1)THEN ! E=HW
  lambda1=lambda1x
  lambda2=lambda2x
  lambda3=lambda3x
  mu1=mu1x
  mu2=mu2x
  mu3=mu3x
ELSE ! E=LW.
! Coefficients with E=LW are obtained from those with E=HW and flipped lambdas and mus - see Eq.(35,2B),(32) and the text below Eq.(32)
  lambda1=mu1x
  lambda2=mu2x
  lambda3=mu3x
  mu1=lambda1x
  mu2=lambda2x
  mu3=lambda3x
END IF

wignerfinal=0.D0

Lambda22max=0

DO rho=1,rhomax

wigner=0.D0

eta=0
DO
  IF(outer_multiplicity(lambda1,mu1,lambda2-1,mu2-1,lambda3,mu3)==rho-1)EXIT
  lambda2=lambda2-1
  mu2=mu2-1
  eta=eta+1
END DO
! lambda2 is \bar{lambda2}, mu2 is \bar{mu2}

! Beginning of (20)

n=(lambda1+lambda2-lambda3+2*(mu1+mu2-mu3))/3
a=(lambda2+lambda3-lambda1-n)/2 ! POZOR: podla (20) by v A,B,c,D malo byt lambda2 a  nie lambda2bar
b=(lambda3+lambda1-lambda2+n+2)/2
c=(lambda1+lambda2-lambda3-n)/2
d=(lambda1+lambda2+lambda3+2-n)/2
epsilon2=lambda1+2*mu1-lambda3-2*mu3
i1=0
DO p2=0,lambda2
  q2=(2*lambda2+mu2-epsilon2)/3-p2
  IF((q2>=0).AND.(q2<=mu2))THEN
    Lambda22=mu2+p2-q2
    IF((ABS(lambda1-Lambda22)<=lambda3).AND.(lambda3<=lambda1+Lambda22))THEN
      i1=i1+1
      IF(i1==1)THEN
        j2=lambda2-p2
!        Lambda22min=Lambda22
      END IF
      j1=lambda2-p2
!      Lambda22max=Lambda22
      Lambda22a(i1)=Lambda22
    END IF
  END IF
END DO
!print*,"v (20): i1=",i1
IF(i1>1)THEN
  DO i2=1,i1
    Lambda22=Lambda22a(i2)
    p2=p(lambda2,mu2,epsilon2,Lambda22)
    p2tilde=mu2-q(lambda2,mu2,epsilon2,Lambda22)
    q2tilde=lambda2-p2
    IF(p2tilde==0)THEN
      F=1
    ELSE
      F=0
      DO i=0,p2tilde
        prod=1
        DO j=0,p2tilde-1
          IF(j<i)THEN
            prod=prod*(p2+j+1)*(mu1+lambda2+mu2-n+j+2)
          ELSE
            prod=prod*(a+j+1)*(b-j-1)
          END IF
        END DO
        F=F+binom(p2tilde,i)*prod
      END DO
      IF(2*(p2tilde/2)/=p2tilde)F=-F
    END IF
    G=1
    IF(j2>j1)THEN
      DO j=j1,j2-1
        IF(j<q2tilde)THEN
          G=G*(a+n-j)*(b-n+j)*(c+n-j)*(d+n-j)*(lambda2+mu2-j+1)
        ELSE
          G=G*(mu2-n+j+1)
        END IF
      END DO
      G=G*(lambda2-2*q2tilde+n+1)*binom(n,q2tilde)
    END IF
    H=binom(n+1+lambda2-q2tilde,lambda2-q2tilde)
    wigner(lambda1,Lambda22)=DFLOAT(F)*DSQRT(DFLOAT(G)/DFLOAT(H))
  END DO
ELSE
  wigner(lambda1,Lambda22a(1))=1.D0
END IF

! End of (20)

!print*,"-----------------------------------------------------------------------------------------------"
!print*,"po (20):"
!do k1=1,i1
!print*,Lambda22a(k1),wigner(lambda1,Lambda22a(k1))
!end do

! Beginning of (21)

epsilon1=-lambda1-2*mu1-3
steps21=(epsilon2+lambda2+2*mu2)/3 ! steps21 is the number of iterations in Eq.(21)
!print*,steps21
IF(steps21>0)THEN
  DO i2=1,steps21
    epsilon2=epsilon2-3 ! epsilon2 is epsilon2 in Eq.(21)
    epsilon1=epsilon1+3 ! epsilon1 is epsilon1 in Eq.(21)
    DO Lambda22=lambda2-steps21+i2,lambda2+steps21-i2,2 ! Lambda22 is 2*Lambda2 in Eq.(21)
      IF(Lambda22<0)CYCLE
      IF(Lambda22>lambda2+mu2)EXIT
  
      p2=p(lambda2,mu2,epsilon2,Lambda22)
      IF((p2<0).OR.(p2>lambda2))CYCLE
      q2=q(lambda2,mu2,epsilon2,Lambda22)
      IF((q2<0).OR.(q2>mu2))CYCLE
  
      Sq2=S(lambda2,mu2,q2) ! Bacha, aby tu p2 a q2 nevychadzali polociselne, teda aby som neuvazoval invalidne hodnoty Lambda22!
      Rp2=R(lambda2,mu2,p2)
  
!      IF(lambda1-i2>=0)wigner(lambda1-i2,Lambda22)&
!      =DSQRT(DFLOAT(lambda1-i2+1)/DFLOAT((lambda1+2-i2)*(Lambda22+1)*R(lambda1,mu1,p(lambda1,mu1,epsilon1,lambda1+1-i2))))&
!      *(-DSQRT(DFLOAT(Sq2*(lambda1-i2+Lambda22-lambda3+2)&
!      *(lambda1-i2+Lambda22+lambda3+4))/DFLOAT(Lambda22+2))*wigner(lambda1+1-i2,Lambda22+1)/2.D0&
!      -DSQRT(DFLOAT(Rp2*(Lambda22+lambda3-lambda1+i2)&
!      *(lambda3+lambda1-i2-Lambda22+2))/DFLOAT(Lambda22))*wigner(lambda1+1-i2,Lambda22-1)/2.D0)
  
      IF(lambda1-i2>=0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
        j1=p(lambda1,mu1,epsilon1+3,lambda1-i2)
        IF((j1>=0).AND.(j1<=lambda1))THEN

          j2=q(lambda1,mu1,epsilon1+3,lambda1-i2)
          IF((j2>=0).AND.(j2<=mu1))THEN

            p1=p(lambda1,mu1,epsilon1,lambda1+1-i2)
            IF((p1>0).AND.(p1<=lambda1))THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! p1 NESMIE BYT 0 ABY SA NEDELILO NULOU. MOZNO ZE BY SA TU NEMALO CYCLOVAT ALE NULOVAT ELEMENT POLA wigner
              q1=q(lambda1,mu1,epsilon1,lambda1+1-i2)
              IF((q1>=0).AND.(q1<=mu1))THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
                IF((Lambda22+1<=lambda2+mu2))THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
                  p2tilde=p(lambda2,mu2,epsilon2+3,Lambda22+1)
                  IF((p2tilde>=0).AND.(p2tilde<=lambda2))THEN
                    q2tilde=q(lambda2,mu2,epsilon2+3,Lambda22+1)
                    IF((q2tilde>=0).AND.(q2tilde<=mu2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
!print*,(lambda1-i2+Lambda22-lambda3+2)*(lambda1-i2+Lambda22+lambda3+4) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kontrolny vypis

                      ABCD=(lambda1-i2+Lambda22-lambda3+2)*(lambda1-i2+Lambda22+lambda3+4) ! PREMENNU ABCD MOZNO UPLNE ZRUSIT, TJ. NEPOCITAT JU A VSADE JU NAHRADIT PRISLUSNYM VYRAZOM.

!     	              IF(ABCD>0)THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEODMOCNOVALO ZAPORNE CISLO
                      wigner(lambda1-i2,Lambda22)=-DSQRT(DFLOAT(Sq2*ABCD)&
                        /DFLOAT(Lambda22+2))*wigner(lambda1+1-i2,Lambda22+1)/2.D0
!                     ELSE
!                       wigner(lambda1-i2,Lambda22)=0.D0
!                     END IF
         
                    END IF
                  END IF

!               ELSE
!                 wigner(lambda1-i2,Lambda22)=0.D0 !!!!!!!!!!!!!!!!!!!!!! TOTO MOZNO NIE JE NUTNE A MOZNO JE
                END IF

!print*,(Lambda22+lambda3-lambda1+i2)*(lambda3+lambda1-i2-Lambda22+2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kontrolny vypis

                IF(Lambda22-1>=0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
                  p2tilde=p(lambda2,mu2,epsilon2+3,Lambda22-1)
                  IF((p2tilde>=0).AND.(p2tilde<=lambda2))THEN
                    q2tilde=q(lambda2,mu2,epsilon2+3,Lambda22-1)
                    IF((q2tilde>=0).AND.(q2tilde<=mu2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
                      ABCD=(Lambda22+lambda3-lambda1+i2)*(lambda3+lambda1-i2-Lambda22+2)
!                     IF(ABCD>0)
                      wigner(lambda1-i2,Lambda22)=wigner(lambda1-i2,Lambda22)-DSQRT(DFLOAT(Rp2*& !!!!!!!!!!!!! PODMIENKA NA ABCD PRIDANA ABY SA NEODMOCNOVALO ZAPORNE CISLO
                      ABCD)/DFLOAT(Lambda22))*wigner(lambda1+1-i2,Lambda22-1)/2.D0

                    END IF
                  END IF

                END IF
                wigner(lambda1-i2,Lambda22)=wigner(lambda1-i2,Lambda22)&
                  *DSQRT(DFLOAT(lambda1-i2+1)/DFLOAT((lambda1+2-i2)*(Lambda22+1)*R(lambda1,mu1,p1)))


              END IF
            END IF
          END IF
        END IF

      END IF


      DO Lambda12=lambda1+1-i2,lambda1-1+i2,2
        IF(Lambda12<0)CYCLE
        IF(Lambda12>lambda1+mu1)EXIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
        p1=p(lambda1,mu1,epsilon1,Lambda12)
        IF((p1<0).OR.(p1>lambda1))CYCLE
        q1=q(lambda1,mu1,epsilon1,Lambda12)
        IF((q1<=0).OR.(q1>mu1))CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! q1 NESMIE BYT 0 ABY SA NEDELILO NULOU. MOZNO ZE BY SA TU NEMALO CYCLOVAT ALE NULOVAT ELEMENT POLA wigner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
!        IF(Lambda12+1<=lambda1+mu1)wigner(Lambda12+1,Lambda22)&
!        =DSQRT(DFLOAT(Lambda12+2)/DFLOAT((Lambda12+1)*(Lambda22+1)*S(lambda1,mu1,q(lambda1,mu1,epsilon1,Lambda12))))& ! Bacha, aby tu q1 nevychadzalo polociselne
!        *(-DSQRT(DFLOAT(Sq2*(Lambda22+lambda3-Lambda12+1)*(lambda3+Lambda12-Lambda22+1))/DFLOAT(Lambda22+2))&
!        *wigner(Lambda12,Lambda22+1)/2.D0&
!        +DSQRT(DFLOAT(Rp2*(Lambda12+Lambda22-lambda3+1)*(Lambda12+Lambda22+lambda3+3))/DFLOAT(Lambda22))&
!        *wigner(Lambda12,Lambda22-1)/2.D0)

        IF(Lambda12+1<=lambda1+mu1)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
          j1=p(lambda1,mu1,epsilon1+3,Lambda12+1)
          IF((j1<0).OR.(j1>lambda1))CYCLE
          j1=q(lambda1,mu1,epsilon1+3,Lambda12+1)
          IF((j1<0).OR.(j1>mu1))CYCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
          IF(Lambda22+1<=lambda2+mu2)THEN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
            p2tilde=p(lambda2,mu2,epsilon2+3,Lambda22+1)
            IF((p2tilde>=0).AND.(p2tilde<=lambda2))THEN
              q2tilde=q(lambda2,mu2,epsilon2+3,Lambda22+1)
              IF((q2tilde>=0).AND.(q2tilde<=mu2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
!print*,(Lambda22+lambda3-Lambda12+1)*(lambda3+Lambda12-Lambda22+1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kontrolny vypis

                ABCD=(Lambda22+lambda3-Lambda12+1)*(lambda3+Lambda12-Lambda22+1)

!                IF(ABCD>0)THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEODMOCNOVALO ZAPORNE CISLO
                wigner(Lambda12+1,Lambda22)=-DSQRT(DFLOAT(Sq2*ABCD)&
                  /DFLOAT(Lambda22+2))*wigner(Lambda12,Lambda22+1)/2.D0
!                ELSE
!                  wigner(Lambda12+1,Lambda22)=0.D0
!                END IF

              END IF
            END IF

!         ELSE
!           wigner(Lambda12+1,Lambda22)=0.D0 !!!!!!!!!!!!!!!!!!!!!!! TOTO MOZNO NIE JE NUTNE A MOZNO JE
          END IF

!print*,(Lambda12+Lambda22-lambda3+1)*(Lambda12+Lambda22+lambda3+3) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kontrolny vypis

          IF(Lambda22-1>=0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST
            p2tilde=p(lambda2,mu2,epsilon2+3,Lambda22-1)
            IF((p2tilde>=0).AND.(p2tilde<=lambda2))THEN
              q2tilde=q(lambda2,mu2,epsilon2+3,Lambda22-1)
              IF((q2tilde>=0).AND.(q2tilde<=mu2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONIEC TESTU
                ABCD=(Lambda12+Lambda22-lambda3+1)*(Lambda12+Lambda22+lambda3+3)
!                IF(ABCD>0)
                wigner(Lambda12+1,Lambda22)=wigner(Lambda12+1,Lambda22)+DSQRT(DFLOAT(Rp2*& !!!!!!!!!!!!! PODMIENKA NA ABCD PRIDANA ABY SA NEODMOCNOVALO ZAPORNE CISLO
                  ABCD)/DFLOAT(Lambda22))*wigner(Lambda12,Lambda22-1)/2.D0

              END IF
            END IF

          END IF
          wigner(Lambda12+1,Lambda22)=wigner(Lambda12+1,Lambda22)&
            *DSQRT(DFLOAT(Lambda12+2)/DFLOAT((Lambda12+1)*(Lambda22+1)*S(lambda1,mu1,q1)))
        END IF

!	 IF((Lambda12==lambda1+1-i2).AND.(Lambda12-1>=0))wigner(Lambda12-1,Lambda22)=FACTOR*wigner(Lambda12,Lambda22+1)+FACTOR*wigner(Lambda12,Lambda22-1)

      END DO
    END DO
  END DO
END IF

! End of (21)

!print*,"-------------------------------------------------------------------"
!print*,"po (21):"
!do k1=0,lambda1+mu1
!do k2=0,lambda2+mu2
!if(wigner(k1,k2)/=0.d0)print*,k1,k2,wigner(k1,k2)
!end do
!end do

! Beginning of (17)

epsilon1=epsilon1+3
IF(eta>0)THEN
  DO i1=1,eta
    lambda2=lambda2+1
    mu2=mu2+1
    epsilon1=epsilon1+3

!	! Zaciatok bloku, ktory mozno nie je nutny
!    DO i2=0,lambda1+mu1
!      wigner(i2,lambda2)=0.D0
!    END DO
!    ! Koniec bloku, ktory mozno nie je nutny
!    DO p2=0,lambda1 ! p2 is p1
!      q2=(2*lambda1+mu1-epsilon1)/3-p2 ! q2 is q1
!      IF((q2>=0).AND.(q2<=mu1))THEN
!        Lambda12=mu1+p2-q2
!   	    IF((ABS(Lambda12-lambda2)<=lambda3).AND.(lambda3<=Lambda12+lambda2))THEN
!	      wigner(Lambda12,lambda2)=-DSQRT(DFLOAT(R(lambda1,mu1,p2+1)*(lambda2+lambda3-Lambda12)*(lambda3+Lambda12-lambda2+2))&
!		  /DFLOAT(4*(Lambda12+2)*(Lambda12+1)))*wigner(Lambda12+1,lambda2-1)&
!		  +DSQRT(DFLOAT(S(lambda1,mu1,q2+1)*(Lambda12+lambda2-lambda3)*(Lambda12+lambda2+lambda3+2))&
!		  /DFLOAT(4*(Lambda12+1)*Lambda12))*wigner(Lambda12-1,lambda2-1)		  
!	    END IF
!      END IF
!    END DO
    
    IF(wigner(1,lambda2-1)/=0.D0)THEN
      wigner(0,lambda2)=-DSQRT(DFLOAT(R(lambda1,mu1,p(lambda1,mu1,epsilon1,0)+1)*(lambda2+lambda3)&
        *(lambda3-lambda2+2))/8.D0)*wigner(1,lambda2-1)
    ELSE
      wigner(0,lambda2)=0.D0
    END IF
    DO Lambda12=1,lambda1+mu1-1
      IF(wigner(Lambda12+1,lambda2-1)/=0.D0)THEN
        wigner(Lambda12,lambda2)=-DSQRT(DFLOAT(R(lambda1,mu1,p(lambda1,mu1,epsilon1,Lambda12)+1)*(lambda2+lambda3-Lambda12)&
          *(lambda3+Lambda12-lambda2+2))/DFLOAT(4*(Lambda12+2)*(Lambda12+1)))*wigner(Lambda12+1,lambda2-1)
      ELSE
        wigner(Lambda12,lambda2)=0.D0
      END IF
      IF(wigner(Lambda12-1,lambda2-1)/=0.D0)wigner(Lambda12,lambda2)=wigner(Lambda12,lambda2)&
        +DSQRT(DFLOAT(S(lambda1,mu1,q(lambda1,mu1,epsilon1,Lambda12)+1)*(Lambda12+lambda2-lambda3)&
        *(Lambda12+lambda2+lambda3+2))/DFLOAT(4*(Lambda12+1)*Lambda12))*wigner(Lambda12-1,lambda2-1)
    END DO
    Lambda12=lambda1+mu1
    IF(wigner(Lambda12-1,lambda2-1)/=0.D0)THEN
      wigner(Lambda12,lambda2)=DSQRT(DFLOAT(S(lambda1,mu1,q(lambda1,mu1,epsilon1,Lambda12)+1)*(Lambda12+lambda2-lambda3)&
        *(Lambda12+lambda2+lambda3+2))/DFLOAT(4*(Lambda12+1)*Lambda12))*wigner(Lambda12-1,lambda2-1)
    ELSE
      wigner(Lambda12,lambda2)=0.D0
    END IF
  
  END DO
END IF

! End of (17)

!print*,"----------------------------------------------------------------------------------------------------------------"
!print*,"po (17):"
!do k1=0,lambda1+mu1
!do k2=0,lambda2+mu2
!if(wigner(k1,k2)/=0.D0)print*,k1,k2,wigner(k1,k2)
!end do
!end do

! This block fills the array wignerfinal with coefficients <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)HW||(lambda3,mu3)HW>
epsilon2=-lambda2-2*mu2
i2=0
DO p1=0,lambda1 ! V tomto bloku by asi bolo dobre vyhladavat validne hodnoty Lambda1 len v prvej iteracii cyklu cez rho, v ktorej sa naplni pole Lambda12a, a v ostatnych iteraciach vyuzivat hodnoty v tomto poli.
  q1=(2*lambda1+mu1-epsilon1)/3-p1
  IF((q1>=0).AND.(q1<=mu1))THEN
    Lambda12=mu1+p1-q1
    IF((ABS(Lambda12-lambda2)<=lambda3).AND.(lambda3<=Lambda12+lambda2))THEN
      i2=i2+1
      IF(rho==rhomax)THEN
        Lambda12a(i2)=Lambda12
        Lambda22a(i2)=lambda2
        epsilon2a(i2)=epsilon2
        IF((epsilon2==lambda1+2*mu1-lambda3-2*mu3).AND.(Lambda12==lambda1))Lambda22max=MAX(Lambda22max,Lambda22) ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW. This is needed for phase convention.
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zaciatok bloku na vymazanie
!      Lambda12at(i2,rho)=Lambda12
!      Lambda22at(i2,rho)=lambda2
!      epsilon2at(i2,rho)=epsilon2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! koniec bloku na vymazanie
      wignerfinal(Lambda12,0,lambda2,rho)=wigner(Lambda12,lambda2)
    END IF
  END IF
END DO

! Beginning of (18)

DO i1=1,lambda2+mu2 ! This is loop over epsilon2=epsilon2HW+3,...,2*lambda2+mu2; i1 is (epsilon2-epsilon2HW)/3
  epsilon1=epsilon1-3 ! epsilon1 is epsilon1 in Eq.(18)
  IF(epsilon1<-lambda1-2*mu1)EXIT
  epsilon2=epsilon2+3 ! epsilon2 is epsilon2+3 in Eq.(18)
  DO Lambda22=lambda2-i1,lambda2+i1,2 ! Vid vceli diagram! ! Lambda22 is 2*Lambda2' in Eq.(18)
    IF(Lambda22<0)CYCLE
    IF(Lambda22>lambda2+mu2)EXIT

    p2tilde=p(lambda2,mu2,epsilon2,Lambda22)
    IF((p2tilde<0).OR.(p2tilde>lambda2))CYCLE
    q2tilde=q(lambda2,mu2,epsilon2,Lambda22)
    IF((q2tilde<0).OR.(q2tilde>mu2))CYCLE

    DO p1=0,lambda1
      q1=(2*lambda1+mu1-epsilon1)/3-p1 ! 2*lambda1+mu1 by asi bolo dobre dat do nejakej premennej, aby sa to furt nepocitalo.
      IF((q1>=0).AND.(q1<=mu1))THEN
        Lambda12=mu1+p1-q1
        IF((ABS(Lambda12-Lambda22)<=lambda3).AND.(lambda3<=Lambda12+Lambda22))THEN

          IF(lambda2>=mu2)THEN
            IF((Lambda22==lambda2-i1).OR.(Lambda22==0))THEN

              p2=p(lambda2,mu2,epsilon2-3,Lambda22+1)
              IF(p2==0)CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEDELILO NULOU

              IF(Lambda12-1>=0)THEN

!print*,(Lambda22+lambda3-Lambda12+2)*(lambda3+Lambda12-Lambda22) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS

                wignerfinal(Lambda12,i1,Lambda22,rho)=-DSQRT(DFLOAT((R(lambda1,mu1,p1)*(Lambda22+lambda3-Lambda12+2)&
                  *(lambda3+Lambda12-Lambda22)))/DFLOAT(4*Lambda12))*wignerfinal(Lambda12-1,i1-1,Lambda22+1,rho)
              ELSE
                wignerfinal(Lambda12,i1,Lambda22,rho)=0.D0
              END IF
              IF(Lambda12+1<=lambda1+mu1)THEN
!print*,(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS
                wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                  +DSQRT(DFLOAT(S(lambda1,mu1,q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
                  /DFLOAT(4*(Lambda12+2)))*wignerfinal(Lambda12+1,i1-1,Lambda22+1,rho)
              END IF
              wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*(Lambda22+2)*R(lambda2,mu2,p2)))
            ELSE

              q2=q(lambda2,mu2,epsilon2-3,Lambda22-1)
              IF(q2==0)CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEDELILO NULOU

              IF(Lambda12-1>=0)THEN

!print*,(Lambda12+Lambda22-lambda3)*(Lambda12+Lambda22+lambda3+2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS

                wignerfinal(Lambda12,i1,Lambda22,rho)=-DSQRT(DFLOAT(R(lambda1,mu1,p1)*(Lambda12+Lambda22-lambda3)&
                  *(Lambda12+Lambda22+lambda3+2))/DFLOAT(4*Lambda12))*wignerfinal(Lambda12-1,i1-1,Lambda22-1,rho)
              ELSE
                wignerfinal(Lambda12,i1,Lambda22,rho)=0.D0
              END IF
              IF(Lambda12+1<=lambda1+mu1)THEN
!print*,(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS
                wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                  -DSQRT(DFLOAT(S(lambda1,mu1,q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
                  /DFLOAT(4*(Lambda12+2)))*wignerfinal(Lambda12+1,i1-1,Lambda22-1,rho)
              END IF
              wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*Lambda22*S(lambda2,mu2,q2)))
            END IF
          ELSE
            IF((Lambda22==lambda2+i1).OR.(Lambda22==lambda2+mu2))THEN

              q2=q(lambda2,mu2,epsilon2-3,Lambda22-1)
              IF(q2==0)CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEDELILO NULOU

              IF(Lambda12-1>=0)THEN

!print*,(Lambda12+Lambda22-lambda3)*(Lambda12+Lambda22+lambda3+2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS

                wignerfinal(Lambda12,i1,Lambda22,rho)=-DSQRT(DFLOAT(R(lambda1,mu1,p1)*(Lambda12+Lambda22-lambda3)&
                  *(Lambda12+Lambda22+lambda3+2))/DFLOAT(4*Lambda12))*wignerfinal(Lambda12-1,i1-1,Lambda22-1,rho)
              ELSE
                wignerfinal(Lambda12,i1,Lambda22,rho)=0.D0
              END IF
              IF(Lambda12+1<=lambda1+mu1)THEN
!print*,(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS
                wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                  -DSQRT(DFLOAT(S(lambda1,mu1,q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
                  /DFLOAT(4*(Lambda12+2)))*wignerfinal(Lambda12+1,i1-1,Lambda22-1,rho)
              END IF
              wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*Lambda22*S(lambda2,mu2,q2)))
            ELSE

              p2=p(lambda2,mu2,epsilon2-3,Lambda22+1)
              IF(p2==0)CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIDANE ABY SA NEDELILO NULOU

              IF(Lambda12-1>=0)THEN

!print*,(Lambda22+lambda3-Lambda12+2)*(lambda3+Lambda12-Lambda22) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS

                wignerfinal(Lambda12,i1,Lambda22,rho)=-DSQRT(DFLOAT(R(lambda1,mu1,p1)*(Lambda22+lambda3-Lambda12+2)&
                  *(lambda3+Lambda12-Lambda22))/DFLOAT(4*Lambda12))*wignerfinal(Lambda12-1,i1-1,Lambda22+1,rho)
              ELSE
                wignerfinal(Lambda12,i1,Lambda22,rho)=0.D0
              END IF
              IF(Lambda12+1<=lambda1+mu1)THEN
!print*,(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KONTROLNY VYPIS
                wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                  +DSQRT(DFLOAT(S(lambda1,mu1,q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
                  /DFLOAT(4*(Lambda12+2)))*wignerfinal(Lambda12+1,i1-1,Lambda22+1,rho)
              END IF
              wignerfinal(Lambda12,i1,Lambda22,rho)=wignerfinal(Lambda12,i1,Lambda22,rho)&
                *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*(Lambda22+2)*R(lambda2,mu2,p2)))
            END IF
          END IF

          i2=i2+1
          IF(rho==rhomax)THEN
            Lambda12a(i2)=Lambda12
            Lambda22a(i2)=Lambda22
            epsilon2a(i2)=epsilon2
            IF((epsilon2==lambda1+2*mu1-lambda3-2*mu3).AND.(Lambda12==lambda1))Lambda22max=MAX(Lambda22max,Lambda22) ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW. This is needed for phase convention.
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zaciatok bloku na vymazanie
!          Lambda12at(i2,rho)=Lambda12
!          Lambda22at(i2,rho)=Lambda22
!          epsilon2at(i2,rho)=epsilon2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! koniec bloku na vymazanie

        END IF
      END IF
    END DO
  END DO
END DO

! End of (18)

!print*,"-------------------------------------------------------------------------------------"
!print*,"po (18):"
!do k1=0,20
!do k2=0,20
!do k3=0,20
!do k4=1,9
!if(wignerfinal(k1,k2,k3,k4)/=0.d0)print*,k1,3*k2-lambda2-2*mu2,k3,k4,wignerfinal(k1,k2,k3,k4)
!end do
!end do
!end do
!end do

! Valid values of Lambda1, epsilon2 and Lambda2 are elements of arrays Lambda12a, epsilon2a and Lambda22a with indeces from 1 to i2.

END DO

! kontrolny vypis
!do rho=1,rhomax
!  do i1=1,i2
!    if(Lambda12at(i1,rho)/=Lambda12a(i1))print*,"Lambda12at/=Lambda12a pre rho=",rho,"a i1=",i1
!    if(epsilon2at(i1,rho)/=epsilon2a(i1))print*,"epsilon2at/=epsilon2a pre rho=",rho,"a i1=",i1
!    if(Lambda22at(i1,rho)/=Lambda22a(i1))print*,"Lambda22at/=Lambda22a pre rho=",rho,"a i1=",i1
!  end do
!end do

! Beginning of orthonormalization according to (8) with alpha3=(epsilon3,Lambda3)=HW.

! Gram-Schmidt orthonormalization is performed, where Wigner coefficients for given rho represent a vector,
! whose components are indexed by alpha1 and alpha2, and the scalar product is defined in Eq.(8).

norm(1)=0.D0
DO i1=1,i2
  Lambda12=Lambda12a(i1)
  epsilon2=epsilon2a(i1)
  Lambda22=Lambda22a(i1)
  norm(1)=norm(1)+wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,1)&
  *wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,1)
END DO

IF(rhomax>1)THEN
 DO rho=2,rhomax
   DO i4=1,rho-1
     scalprod(rho,i4)=0.D0
     DO i1=1,i2
       Lambda12=Lambda12a(i1)
       epsilon2=epsilon2a(i1)
       Lambda22=Lambda22a(i1)
       scalprod(rho,i4)=scalprod(rho,i4)+wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)&
         *wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,i4)
     END DO
   END DO
   DO i4=1,rho-1
     DO i1=1,i2
       Lambda12=Lambda12a(i1)
       epsilon2=epsilon2a(i1)
       Lambda22=Lambda22a(i1)
       wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)&
         =wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)&
         -scalprod(rho,i4)*wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,i4)/norm(i4)
     END DO
   END DO
   norm(rho)=0.D0
   DO i1=1,i2
     Lambda12=Lambda12a(i1)
     epsilon2=epsilon2a(i1)
     Lambda22=Lambda22a(i1)
     norm(rho)=norm(rho)+wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)&
       *wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)
   END DO
 END DO
END IF

DO rho=1,rhomax
  DO i1=1,i2
    Lambda12=Lambda12a(i1)
    epsilon2=epsilon2a(i1)
    Lambda22=Lambda22a(i1)
    wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)&
    =wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)/DSQRT(norm(rho))
  END DO
END DO

! End of the orthonormalization

! Beginning of phase convention (K.T.Hecht, Nucl.Phys. 62 (1965) 1)

! The phase is such that <(lambda1,mu1)LW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)LW>_rho > 0, which
! according to formula 2B from Eq.(35) with rho_max instead of eta_max (see text therein) means that
! <(lambda1,mu1)HW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)HW>_rho*(-1)^E>0
! where E=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax-rho+(lambda1+2*Lambda2_max-lambda3)/2

!a=MAXVAL(Lambda22a) ! a is 2*Lambda2_max
i=lambda1+lambda2-lambda3+mu1+mu2-mu3
i4=i+rhomax+(lambda1+Lambda22max-lambda3)/2 ! i3 is E+rho
i=(i+mu1+mu2-mu3)/3 ! i is (epsilon2-epsilon2_HW)/3 for epsilon2=epsilon3_HW-epsilon1_HW
DO rho=1,rhomax
  IF(((2*((i4-rho)/2)==i4-rho).AND.(wignerfinal(lambda1,i,Lambda22max,rho)<0.D0))&
  .OR.((2*((i4-rho)/2)/=i4-rho).AND.(wignerfinal(lambda1,i,Lambda22max,rho)>0.D0)))THEN ! i3-rho is E
    DO i1=1,i2
      Lambda12=Lambda12a(i1)
      epsilon2=epsilon2a(i1)
      Lambda22=Lambda22a(i1)
      wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)=-wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)
    END DO
  END IF
END DO

! End of phase convention

! kontrolny vypis
!DO i1=1,i2
!  Lambda12=Lambda12a(i1)
!  epsilon2=epsilon2a(i1)
!  Lambda22=Lambda22a(i1)
!  PRINT*,Lambda12,epsilon2,Lambda22,(wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho),rho=1,rhomax)
!END DO

IF(I3==0)THEN ! E=LW
  phiprhomax=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax ! phiprhomax is phi+rhomax
  DO i1=1,i2
    Lambda12=Lambda12a(i1)
    epsilon2=epsilon2a(i1)
    Lambda22=Lambda22a(i1)
    DO rho=1,rhomax
      ! See Eq.(35,2B)
      IF(BTEST(phiprhomax-rho+(Lambda12+Lambda22-lambda3)/2,0))&
      wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)=-wignerfinal(Lambda12,(epsilon2+lambda2+2*mu2)/3,Lambda22,rho)
    END DO
  END DO
  epsilon2a(1:i2)=-epsilon2a(1:i2) ! For I3=0 (E=LW), values of epsilon2 must change sign (see Eq.(32)).
END IF

CONTAINS
 FUNCTION p(lambda,mu,epsilon,Lam2) RESULT(rp)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: lambda,mu,epsilon,Lam2
  INTEGER :: rp
  rp=(2*lambda+mu-epsilon)/3-mu+Lam2
  IF(2*(rp/2)==rp)THEN
    rp=rp/2
  ELSE
    rp=-1 ! If the values of epsilon and Lam2 are invalid, p is negative to indicate that.
  END IF
 END FUNCTION p
 
 FUNCTION q(lambda,mu,epsilon,Lam2) RESULT(rq)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: lambda,mu,epsilon,Lam2
  INTEGER :: rq
  rq=(2*lambda+mu-epsilon)/3+mu-Lam2
  IF(2*(rq/2)==rq)THEN
    rq=rq/2
  ELSE
    rq=-1 ! If the values of epsilon and Lam2 are invalid, q is negative to indicate that.
  END IF
 END FUNCTION q

 FUNCTION R(lambda,mu,p) RESULT(rR)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: lambda,mu,p
  INTEGER :: rR
  rR=p*(lambda+1-p)*(mu+1+p)
 END FUNCTION R
 
 FUNCTION S(lambda,mu,q) RESULT(rS)
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: lambda,mu,q
  INTEGER :: rS
  rS=q*(mu+1-q)*(lambda+mu+2-q)
 END FUNCTION S
 
 function factorial (n) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer(KIND=8) :: res
    integer :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function binom (n, k) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer(KIND=8) :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function binom

END SUBROUTINE wigner_canonical_extremal
