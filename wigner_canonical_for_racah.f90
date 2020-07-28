SUBROUTINE wigner_canonical_for_racah(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32max,numbLambda3,rhomax,wignerex,wigner)
!--------------------------------------------------------------------------------------------------------------------------------
! Calculates the reduced SU(3)-SU(2)xU(1) Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3 and 2*Lambda3=Lambda32max,Lambda32max-2,...,Lambda32max-2*(numbLambda3-1)
! using Eq.(19) in [1] and Table 9.1 on page 311 in [2].
!
! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!             [2] Varshalovich, Quantum theory of angular momentum
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32max,numbLambda32,rhomax,wignerex
! Output arguments: wigner
!
! rhomax=the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
! wignerex is an array containing the highest-weight coefficients calculated by subroutine wigner_canonical_extremal with I3=1.
! 
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho=wigner(p1,p2,q2,Lambda3ind,rho)
!   where epsilon2=2*lambda2+mu2-3*(p2+q2)
!         epsilon1=epsilon3-epsilon2
!         Lambda1=(mu1+p1-q1)/2
!         Lambda2=(mu2+p2-q2)/2
!         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
!         Lambda32ind=1,2,...,numbLambda3 correspond to 2*Lambda3=Lambda32max,Lambda32max-2,...,Lambda32max-2*(numbLambda3-1)
!
! Note: There is a typo in Eq.(19). In the very last line there should be p instead of q.
!-------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32max,numbLambda3,rhomax
INTEGER :: Lambda32,eps3,rho,p1,q1,p2,q2,Lam32,epsilon2max,epsilon3min,Lambda22,Lambda12,Lam32prime,p3,q3,&
           q2min,p1min,pq1,noname1,noname2,Lambda3ind,numb,p3start,q3start,epsilon3mid,Rp2,Sq2,s2,lm1,lm2,lm3,mu1p2,mu2p2,mu3p2
REAL(KIND=8) :: N3
REAL(KIND=8),DIMENSION(0:,0:,0:,1:) :: wignerex
REAL(KIND=8),DIMENSION(0:,0:,0:,1:,1:),INTENT(OUT) :: wigner

epsilon3min=-lambda3-2*mu3
IF(epsilon3==epsilon3min)THEN
  DO rho=1,rhomax
    DO q2=0,mu2
      DO p2=0,lambda2
        wigner(0:lambda1,p2,q2,1,rho)=wignerex(0:lambda1,p2,q2,rho)
      END DO
    END DO
  END DO
  RETURN
ELSE

lm1=lambda1+mu1+1
lm2=lambda2+mu2+1
lm3=lambda3+mu3+1
mu1p2=mu1+2
mu2p2=mu2+2
mu3p2=mu3+2

Lambda32=Lambda32max-numbLambda3+1 ! This is the middle value of 2*Lambda3
epsilon2max=2*lambda2+mu2
noname1=(epsilon2max-lambda1-2*mu1-epsilon3min)/3
noname2=(epsilon2max+2*lambda1+mu1-epsilon3min)/3
epsilon3min=epsilon3min+3
epsilon3mid=epsilon3-3*(numbLambda3-1)

Lam32=lambda3
p3=lambda3
q3=mu3

DO eps3=epsilon3min,epsilon3mid,3 ! eps3 is epsilon3 in Eq.(19)
  noname1=noname1-1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
  noname2=noname2-1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3
    
  Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
  IF(Lam32<Lambda32)THEN
    Lam32=Lam32+1
    q3=q3-1 ! p3 and q3 correspond to eps3 and Lam32
!    N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((q3+1)*(mu3-q3)*(lambda3+mu3+1-q3)))
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((q3+1)*(mu3-q3)*(lm3-q3)*(Lam32prime+2)*4))
  ELSE
    Lam32=Lam32-1
    p3=p3-1 ! p3 and q3 correspond to eps3 and Lam32
!    N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((p3+1)*(lambda3-p3)*(mu3+2+p3)))
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((p3+1)*(lambda3-p3)*(mu3p2+p3)*Lam32prime*4))
  END IF ! Lam32 is 2*Lambda3 in Eq.(19) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.

  DO p2=MAX(0,noname1-mu2,(Lam32-mu1+noname2-mu2)/2-lambda1),MIN(lambda2,noname2)
    Rp2=(p2+1)*(lambda2-p2)*(mu2p2+p2)
    q2min=MAX(0,noname1-p2,(mu2+noname2-Lam32-mu1)/2-lambda1)
    pq1=noname2-p2-q2min ! pq1 is p1+q1
    Lambda22=mu2+p2-q2min ! Lambda22 is 2*Lambda2 in Eq.(19)
    DO q2=q2min,MIN(mu2,noname2-p2,(mu2+noname2+Lam32-mu1)/2)
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
!     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
!     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
      Sq2=(q2+1)*(mu2-q2)*(lm2-q2)
      p1min=MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2)
      q1=pq1-p1min
      Lambda12=mu1+p1min-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
!      if(p1min>MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2))print*,"!!!!!!!!!!!!!"
      DO p1=p1min,MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
  
        s2=Lambda12+Lambda22+Lam32prime+1
        IF(Lam32>Lam32prime)THEN

          IF(q1/=mu1)THEN
            wignerex(p1,p2,q2,1:rhomax)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
             *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wignerex(p1,p2,q2,1:rhomax)
          ELSE
            wignerex(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=lambda1)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p2+p1)&
            *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wignerex(p1+1,p2,q2,1:rhomax)

          IF(p2/=lambda2)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)-DSQRT(DFLOAT(Rp2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wignerex(p1,p2+1,q2,1:rhomax)

          IF(q2/=mu2)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wignerex(p1,p2,q2+1,1:rhomax)

        ELSE

          IF(q1/=mu1)THEN
            wignerex(p1,p2,q2,1:rhomax)=-DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wignerex(p1,p2,q2,1:rhomax)
          ELSE
            wignerex(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=lambda1)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p2+p1)&
            *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wignerex(p1+1,p2,q2,1:rhomax)

          IF(p2/=lambda2)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Rp2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wignerex(p1,p2+1,q2,1:rhomax)

          IF(q2/=mu2)wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wignerex(p1,p2,q2+1,1:rhomax)

        END IF

        wignerex(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)*N3

        q1=q1-1
        Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
      END DO
      pq1=pq1-1 ! pq1 is p1+q1
      Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
    END DO
  END DO
END DO

!----------------------------------------------------------------------------------------
  DO rho=1,rhomax
    DO q2=0,mu2
      DO p2=0,lambda2
        wigner(0:lambda1,p2,q2,1,rho)=wignerex(0:lambda1,p2,q2,rho)
      END DO
    END DO
  END DO

! Now p3 and q3 correspond to epsilon3mid and Lambda32 (the middle value of 2*Lambda3)
  p3start=p3
  q3start=q3

!  DO eps3=epsilon3mid+3,epsilon3,3 ! eps3 is epsilon3 in Eq.(19)
  DO numb=2,numbLambda3

    noname1=noname1-1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
    noname2=noname2-1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3



    Lambda3ind=numb-1
    Lam32=Lambda32-Lambda3ind
    Lam32prime=Lam32+1
    p3=p3start-numb+1
!    N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((p3+1)*(lambda3-p3)*(mu3+2+p3)))
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((p3+1)*(lambda3-p3)*(mu3p2+p3)*Lam32prime*4))

!*****************************************************************************

    DO p2=MAX(0,noname1-mu2,(Lam32-mu1+noname2-mu2)/2-lambda1),MIN(lambda2,noname2)
      Rp2=(p2+1)*(lambda2-p2)*(mu2p2+p2)
      q2min=MAX(0,noname1-p2,(mu2+noname2-Lam32-mu1)/2-lambda1)
      pq1=noname2-p2-q2min ! pq1 is p1+q1
      Lambda22=mu2+p2-q2min ! Lambda22 is 2*Lambda2 in Eq.(19)
      DO q2=q2min,MIN(mu2,noname2-p2,(mu2+noname2+Lam32-mu1)/2)
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
        Sq2=(q2+1)*(mu2-q2)*(lm2-q2)
        p1min=MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2)
        q1=pq1-p1min
        Lambda12=mu1+p1min-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
        DO p1=p1min,MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

!          IF(q1/=mu1)THEN
!            wigner(p1,p2,q2,numb,1:rhomax)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lambda1+mu1+1-q1)))&
!              *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12-1,Lam32)*wigner(p1,p2,q2,Lambda3ind,1:rhomax)
!          ELSE
!            wigner(p1,p2,q2,numb,1:rhomax)=0.D0
!          END IF

!          IF(p1/=lambda1)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)&
!            +DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1+2+p1)))&
!            *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12+1,Lam32)*wigner(p1+1,p2,q2,Lambda3ind,1:rhomax)

!          IF(p2/=lambda2)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)&
!            +DSQRT(DFLOAT((p2+1)*(lambda2-p2)*(mu2+2+p2)))&
!            *DRR3(Lambda12,Lambda22+1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2+1,q2,Lambda3ind,1:rhomax)

!          IF(q2/=mu2)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)&
!            +DSQRT(DFLOAT((q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)))&
!            *DRR3(Lambda12,Lambda22-1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2,q2+1,Lambda3ind,1:rhomax)

          s2=Lambda12+Lambda22+Lam32prime+1

          IF(q1/=mu1)THEN
            wigner(p1,p2,q2,numb,1:rhomax)=-DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,Lambda3ind,1:rhomax)
          ELSE
            wigner(p1,p2,q2,numb,1:rhomax)=0.D0
          END IF

          IF(p1/=lambda1)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p2+p1)&
            *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,Lambda3ind,1:rhomax)

          IF(p2/=lambda2)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)+DSQRT(DFLOAT(Rp2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,Lambda3ind,1:rhomax)

          IF(q2/=mu2)wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,Lambda3ind,1:rhomax)


          wigner(p1,p2,q2,numb,1:rhomax)=wigner(p1,p2,q2,numb,1:rhomax)*N3

          q1=q1-1
          Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
        END DO
        pq1=pq1-1 ! pq1 is p1+q1
        Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
      END DO
    END DO

!*****************************************************************************

    q3start=q3start-1
    Lam32prime=Lambda32+numb
    Lam32=Lam32prime+1
    q3=q3start
    DO Lambda3ind=1,numb-1
      Lam32prime=Lam32prime-2
      Lam32=Lam32-2
!      N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((q3+1)*(mu3-q3)*(lambda3+mu3+1-q3)))
      N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((q3+1)*(mu3-q3)*(lm3-q3)*(Lam32prime+2)*4))

!***************************************************************************

      DO p2=MAX(0,noname1-mu2,(Lam32-mu1+noname2-mu2)/2-lambda1),MIN(lambda2,noname2)
        Rp2=(p2+1)*(lambda2-p2)*(mu2p2+p2)
        q2min=MAX(0,noname1-p2,(mu2+noname2-Lam32-mu1)/2-lambda1)
        pq1=noname2-p2-q2min ! pq1 is p1+q1
        Lambda22=mu2+p2-q2min ! Lambda22 is 2*Lambda2 in Eq.(19)
        DO q2=q2min,MIN(mu2,noname2-p2,(mu2+noname2+Lam32-mu1)/2)
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
          Sq2=(q2+1)*(mu2-q2)*(lm2-q2)
          p1min=MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2)
          q1=pq1-p1min
          Lambda12=mu1+p1min-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
          DO p1=p1min,MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

!            IF(q1/=mu1)THEN
!              wigner(p1,p2,q2,Lambda3ind,1:rhomax)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lambda1+mu1+1-q1)))&
!                *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12-1,Lam32)*wigner(p1,p2,q2,Lambda3ind,1:rhomax)
!            ELSE
!              wigner(p1,p2,q2,Lambda3ind,1:rhomax)=0.D0
!            END IF

!            IF(p1/=lambda1)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)&
!              +DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1+2+p1)))&
!              *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12+1,Lam32)*wigner(p1+1,p2,q2,Lambda3ind,1:rhomax)

!            IF(p2/=lambda2)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)&
!              +DSQRT(DFLOAT((p2+1)*(lambda2-p2)*(mu2+2+p2)))&
!              *DRR3(Lambda12,Lambda22+1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2+1,q2,Lambda3ind,1:rhomax)

!            IF(q2/=mu2)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)&
!              +DSQRT(DFLOAT((q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)))&
!              *DRR3(Lambda12,Lambda22-1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2,q2+1,Lambda3ind,1:rhomax)

            s2=Lambda12+Lambda22+Lam32prime+1

            IF(q1/=mu1)THEN
              wigner(p1,p2,q2,Lambda3ind,1:rhomax)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
               *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,Lambda3ind,1:rhomax)
            ELSE
              wigner(p1,p2,q2,Lambda3ind,1:rhomax)=0.D0
            END IF

            IF(p1/=lambda1)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)&
              +DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p2+p1)&
              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,Lambda3ind,1:rhomax)

            IF(p2/=lambda2)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)-DSQRT(DFLOAT(Rp2&
              *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,Lambda3ind,1:rhomax)

            IF(q2/=mu2)wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)+DSQRT(DFLOAT(Sq2&
              *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,Lambda3ind,1:rhomax)


            wigner(p1,p2,q2,Lambda3ind,1:rhomax)=wigner(p1,p2,q2,Lambda3ind,1:rhomax)*N3

            q1=q1-1
            Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
          END DO
          pq1=pq1-1 ! pq1 is p1+q1
          Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
        END DO
      END DO

!***************************************************************************

      q3=q3+1
    END DO

  END DO

!END DO

END IF

END SUBROUTINE wigner_canonical_for_racah
