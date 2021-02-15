!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wigner_canonical.f90 -- SU(3)-SU(2)xU(1) coupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,I3,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
!--------------------------------------------------------------------------------------------------------------------------------
! Calculates the reduced SU(3)-SU(2)xU(1) Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3 using Eq.(19) in [1] and its conjugate and Table 9.1 on page 311 in [2].
!
! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!             [2] Varshalovich, Quantum theory of angular momentum
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,I3,rhomax,wignerex
! Input and output arguments: numb,p1a,p2a,q2a,wigner
!
! Lambda32=2*Lambda3 (epsilon3 and Lambda32 must be valid)
! rhomax=the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
! wignerex is an array containing the extremal coefficients calculated by subroutine wigner_canonical_extremal.
! If I3=1, coefficients are calculated from the highest-weight (HW) coefficients.
! If I3=0, coefficients are calculated from the lowest-weight (LW) coefficients.
! 
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho=wigner(p1,p2,q2,rho)
!   where epsilon2=2*lambda2+mu2-3*(p2+q2)
!         epsilon1=epsilon3-epsilon2
!         Lambda1=(mu1+p1-q1)/2
!         Lambda2=(mu2+p2-q2)/2
!         p1=p1a(i) 
!         p2=p2a(i)
!         q2=q2a(i)
!         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
!         1<=i<=numb
!   unless I3=0 and epsilon3=epsilon3LW=2*lambda3+mu3, in this case:
!         q1=mu1-p1a(i)
!         p2=lambda2-q2a(i)
!         q2=mu2-p2a(i)
!         p1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-q1-p2-q2
!
! Note: There is a typo in Eq.(19). In the very last line there should be p instead of q.
!-------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,I3,rhomax
!INTEGER,INTENT(OUT) :: numb
INTEGER :: numb,eps3,rho,p1,q1,p2,q2,Lam32,epsilon2max,epsilon3ex,Lambda22,Lambda12,Lam32prime,p3,q3,&
           pq1,noname1,noname2,Rp2,Sq2,s2,lm1,lm2,lm3,mu1p,mu2p,mu3p,q2ex,p1ex,lambda1p,lambda2p,lambda3p
REAL(KIND=8) :: N3
!INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a ! Dimension is at least (lambda1+1)*(lambda2+1)*(mu2+1)
INTEGER,DIMENSION(:) :: p1a,p2a,q2a ! Dimension is at least (lambda1+1)*(lambda2+1)*(mu2+1)
REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wignerex
REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner

IF(I3==1)THEN

DO rho=1,rhomax
  DO q2=0,mu2
    DO p2=0,lambda2
      wigner(0:lambda1,p2,q2,rho)=wignerex(0:lambda1,p2,q2,rho)
    END DO
  END DO
END DO
!do noname1=1,numb
!p1=p1a(noname1)
!p2=p2a(noname1)
!q2=q2a(noname1)
!wigner(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)
!end do

epsilon3ex=-lambda3-2*mu3
IF(epsilon3==epsilon3ex)RETURN

lm1=lambda1+mu1+1
lm2=lambda2+mu2+1
lm3=lambda3+mu3+1
mu1p=mu1+2
mu2p=mu2+2
mu3p=mu3+2

epsilon2max=2*lambda2+mu2
noname1=(epsilon2max-lambda1-2*mu1-epsilon3ex)/3
noname2=(epsilon2max+2*lambda1+mu1-epsilon3ex)/3
epsilon3ex=epsilon3ex+3

numb=0

Lam32=lambda3
p3=lambda3
q3=mu3

DO eps3=epsilon3ex,epsilon3,3 ! eps3 is epsilon3 in Eq.(19)
  noname1=noname1-1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
  noname2=noname2-1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3

  Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
  IF(Lam32<Lambda32)THEN
    Lam32=Lam32+1
    q3=q3-1 ! p3 and q3 correspond to eps3 and Lam32
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((q3+1)*(mu3-q3)*(lm3-q3)*(Lam32prime+2)*4))
  ELSE
    Lam32=Lam32-1
    p3=p3-1 ! p3 and q3 correspond to eps3 and Lam32
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT((p3+1)*(lambda3-p3)*(mu3p+p3)*Lam32prime*4))
  END IF ! Lam32 is 2*Lambda3 in Eq.(19) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.

  DO p2=MAX(0,noname1-mu2,(Lam32-mu1+noname2-mu2)/2-lambda1),MIN(lambda2,noname2)
    Rp2=(p2+1)*(lambda2-p2)*(mu2p+p2)
    q2ex=MAX(0,noname1-p2,(mu2+noname2-Lam32-mu1)/2-lambda1)
    pq1=noname2-p2-q2ex ! pq1 is p1+q1
    Lambda22=mu2+p2-q2ex ! Lambda22 is 2*Lambda2 in Eq.(19)
    DO q2=q2ex,MIN(mu2,noname2-p2,(mu2+noname2+Lam32-mu1)/2)
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
!     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
!     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
      Sq2=(q2+1)*(mu2-q2)*(lm2-q2)
      p1ex=MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2)
      q1=pq1-p1ex
      Lambda12=mu1+p1ex-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
!      if(p1min>MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2))print*,"!!!!!!!!!!!!!"
      DO p1=p1ex,MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

        IF(eps3==epsilon3)THEN
          numb=numb+1
          p1a(numb)=p1
          p2a(numb)=p2
          q2a(numb)=q2
        END IF

        s2=Lambda12+Lambda22+Lam32prime+1
        IF(Lam32>Lam32prime)THEN

          IF(q1/=mu1)THEN
            wigner(p1,p2,q2,1:rhomax)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
             *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,1:rhomax)
          ELSE
            wigner(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=lambda1)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p+p1)&
            *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,1:rhomax)

          IF(p2/=lambda2)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)-DSQRT(DFLOAT(Rp2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,1:rhomax)

          IF(q2/=mu2)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,1:rhomax)

        ELSE

          IF(q1/=mu1)THEN
            wigner(p1,p2,q2,1:rhomax)=-DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lm1-q1)&
              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,1:rhomax)
          ELSE
            wigner(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=lambda1)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1p+p1)&
            *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,1:rhomax)

          IF(p2/=lambda2)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Rp2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,1:rhomax)

          IF(q2/=mu2)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,1:rhomax)

        END IF

        wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)*N3

        q1=q1-1
        Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
      END DO
      pq1=pq1-1 ! pq1 is p1+q1
      Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
    END DO
  END DO
END DO

ELSE

epsilon3ex=2*lambda3+mu3
epsilon2max=2*lambda2+mu2
noname1=(epsilon2max-lambda1-2*mu1-epsilon3ex)/3
noname2=(epsilon2max+2*lambda1+mu1-epsilon3ex)/3
DO p1=0,lambda1
  DO p2=0,lambda2
    DO q2=0,mu2
!      q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3max)/3-p1-p2-q2
      q1=noname2-p1-p2-q2
      wigner(p1,p2,q2,1:rhomax)=wignerex(mu1-q1,mu2-q2,lambda2-p2,1:rhomax)
    END DO
  END DO
END DO
!do noname1=1,numb
!p1=p1a(noname1)
!p2=p2a(noname1)
!q2=q2a(noname1)
!wigner(p1,p2,q2,1:rhomax)=wignerex(p1,p2,q2,1:rhomax)
!end do

IF(epsilon3==epsilon3ex)RETURN

lambda1p=lambda1+1
lambda2p=lambda2+1
lambda3p=lambda3+1
mu1p=mu1+1
mu2p=mu2+1
mu3p=mu3+1
lm1=lambda1p+mu1p
lm2=lambda2p+mu2p
lm3=lambda3p+mu3p

epsilon3ex=epsilon3ex-3

numb=0

Lam32=mu3
p3=0
q3=0

DO eps3=epsilon3ex,epsilon3,-3 ! eps3 is epsilon3 in Eq.(19)
  noname1=noname1+1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
  noname2=noname2+1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3
    
  Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
  IF(Lam32<Lambda32)THEN
    Lam32=Lam32+1
    p3=p3+1 ! p3 and q3 correspond to eps3 and Lam32
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(p3*(lambda3p-p3)*(mu3p+p3)*(Lam32prime+2)*4)) ! N3 is sqrt((2*Lambda3+1)/(4*(2*Lambda3'+1)))/N3
  ELSE
    Lam32=Lam32-1
    q3=q3+1 ! p3 and q3 correspond to eps3 and Lam32
    N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(q3*(mu3p-q3)*(lm3-q3)*Lam32prime*4)) ! N3 is sqrt((2*Lambda3+1)/(8*Lambda3'))/N3
  END IF ! Lam32 is 2*Lambda3 in Eq.(19).

  DO p2=MIN(lambda2,noname2),MAX(0,noname1-mu2,(Lam32-mu1+noname2-mu2)/2-lambda1),-1
    Rp2=p2*(lambda2p-p2)*(mu2p+p2) ! Rp2 is R(p2)
    q2ex=MIN(mu2,noname2-p2,(mu2+noname2+Lam32-mu1)/2)
    pq1=noname2-p2-q2ex ! pq1 is p1+q1
    Lambda22=mu2+p2-q2ex ! Lambda22 is 2*Lambda2 in Eq.(19)
    DO q2=q2ex,MAX(0,noname1-p2,(mu2+noname2-Lam32-mu1)/2-lambda1),-1
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
!     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
!     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
      Sq2=q2*(mu2p-q2)*(lm2-q2)
      p1ex=MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
      q1=pq1-p1ex
      Lambda12=mu1+p1ex-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
!      if(p1min>MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2))print*,"!!!!!!!!!!!!!"
      DO p1=p1ex,MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2),-1
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
  
        IF(eps3==epsilon3)THEN
          numb=numb+1
          p1a(numb)=p1
          p2a(numb)=p2
          q2a(numb)=q2
        END IF
  
        s2=Lambda12+Lambda22+Lam32prime+1
        IF(Lam32>Lam32prime)THEN

          IF(q1/=0)THEN
            wigner(p1,p2,q2,1:rhomax)=-DSQRT(DFLOAT(q1*(mu1p-q1)*(lm1-q1)&
            *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1,p2,q2,1:rhomax) ! clen s Lambda1p=Lambda1+1/2
          ELSE
            wigner(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(p1*(lambda1p-p1)*(mu1p+p1)&
            *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1-1,p2,q2,1:rhomax) ! clen s Lambda1p=Lambda1-1/2

          IF(q2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
           *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2,q2-1,1:rhomax) ! clen s Lambda2p=Lambdda2+1/2

          IF(p2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Rp2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2-1,q2,1:rhomax) ! clen s Lambda2p=Lambda2-1/2

        ELSE

          IF(q1/=0)THEN
            wigner(p1,p2,q2,1:rhomax)=DSQRT(DFLOAT(q1*(mu1p-q1)*(lm1-q1)&
            *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1,p2,q2,1:rhomax) ! clen s Lambda1p=Lambda1+1/2
          ELSE
            wigner(p1,p2,q2,1:rhomax)=0.D0
          END IF

          IF(p1/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(p1*(lambda1p-p1)*(mu1p+p1)&
            *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1-1,p2,q2,1:rhomax) ! clen s Lambda1p=Lambda1-1/2

          IF(q2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
            *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2,q2-1,1:rhomax) ! clen s Lambda2p=Lambda2+1/2

          IF(p2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)-DSQRT(DFLOAT(Rp2&
            *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2-1,q2,1:rhomax) ! clen s Lambda2p=Lambda2-1/2

        END IF

        wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)*N3
        
        q1=q1+1
        Lambda12=Lambda12-2 ! Lambda12 is 2*Lambda1 in Eq.(19)
      END DO
      pq1=pq1+1 ! pq1 is p1+q1
      Lambda22=Lambda22+1 ! Lambda22 is 2*Lambda2 in Eq.(19)
    END DO
  END DO
END DO

END IF

END SUBROUTINE wigner_canonical
