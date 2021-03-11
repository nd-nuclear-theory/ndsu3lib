SUBROUTINE wigner_canonical_extremal(lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax,&
                                     i2,wigner,p1a,p2a,q2a)!,Lambda12a,Lambda22a,epsilon2a)
!--------------------------------------------------------------------------------------------------------------------------------
! Calculates the extremal reduced SU(3)-SU(2)xU(1) Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho for given lambda1,mu1,lambda2,mu2,lambda3,mu3
! using Eq.(20),(21),(17),(18),(8),(1),(32),(35,2B) in Ref.[1]. The phase convention of Ref.[2] is adopted.
!
! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!             [2] K.T.Hecht, Nucl.Phys. 62 (1965) 1
!             [3] D.Goldberg, ACM Computing Surveys, Vol.23, No.1 (1991) 5
!
! Input arguments: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
! Output arguments: i2,wigner,p1a,p2a,q2a
!
! lambda1x=lambda1
! mu1x=mu1
! lambda2x=lambda2
! mu2x=mu2
! lambda3x=lambda3
! mu3x=mu3
! I3=1 for E=HW, I3=0 for E=LW
! rhomax=the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
! 
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho=wigner(p1,p2,q2,rho)
!   where epsilon2=2*lambda2+mu2-3*(p2+q2) for I3=1, epsilon2=-2*mu2-lambda2+3*(p2+q2) for I3=0
!         epsilon1=epsilon3E-epsilon2
!         Lambda1=(mu1+p1-q1)/2 for I3=1, Lambda1=(lambda1+p1-q1)/2 for I3=0
!         Lambda2=(mu2+p2-q2)/2 for I3=1, Lambda2=(lambda2+p2-q2)/2 for I3=0
!         p1=p1a(i)
!         p2=p2a(i)
!         q2=q2a(i)
!         q1=(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3-p1-p2-q2 for I3=1
!         q1=(2*(mu1+mu2+lambda3)+lambda1+lambda2+mu3)/3-p1-p2-q2 for I3=0
!         1<=i<=i2
!
! Note: There are typos in Eq.(18),(20). There are 4 expressions for X in Eq.(18). The 3rd one is for X(Lambda1+1/2,Lambda2-1/2),
!       not X(Lambda1-1/2,Lambda2-1/2). In Eq.(20), there should be \bar{lambda2} instead of lambda 2 in a,b,c,d.
!-------------------------------------------------------------------------------------------------------------------------------
USE binomial_coeff
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER,INTENT(IN) :: lambda1x,mu1x,lambda2x,mu2x,lambda3x,mu3x,I3,rhomax
INTEGER,INTENT(OUT) :: i2
INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,p2,q2,Lambda22,i1,j2,j1,p2tilde,q2tilde,i,j,n,a,b,c,d,&
           epsilon2,steps21,eps2,Lambda12,epsilon1,eta,p1,q1,rho,i4,Lambda22max,ABCD,&
           phiprhomax,p2min,p2max,noname1,p1min,noname2,p1max
INTEGER(KIND=8) :: Sq2,Rp2
!INTEGER(KIND=8) :: prod!,F,G,H
REAL(KIND=8) :: F,G,H,scalprod,CC,Y,T!,prod!,norm
!INTEGER,DIMENSION(4096) :: Lambda12a,Lambda22a,epsilon2a ! Dimension is (lambda1+1)*(lambda2+1)*(mu2+1)
INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a ! Dimension is (lambda1+1)*(lambda2+1)*(mu2+1)
REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner ! First index is p1, second is p2, third is q2, fourth is rho

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

!wigner=0
wigner(0:lambda1,0:lambda2,0:mu2,1:rhomax)=0.D0
Lambda22max=0

DO rho=1,rhomax

eta=0
DO
  IF(outer_multiplicity(lambda1,mu1,lambda2-1,mu2-1,lambda3,mu3)==rho-1)EXIT
  lambda2=lambda2-1
  mu2=mu2-1
  eta=eta+1
END DO
! lambda2 is \bar{lambda2}, mu2 is \bar{mu2}
!**************************************************************************************
! Beginning of (20)
!**************************************************************************************
!n=(lambda1+lambda2-lambda3+2*(mu1+mu2-mu3))/3
!a=(lambda2+eta+lambda3-lambda1-n)/2
!b=(lambda3+lambda1-lambda2+eta+n+2)/2
!c=(lambda1+lambda2+eta-lambda3-n)/2
!d=(lambda1+lambda2+eta+lambda3+2-n)/2
epsilon2=lambda1+2*mu1-lambda3-2*mu3
noname1=(2*lambda2+mu2-epsilon2)/3
p2max=MIN(lambda2,noname1,(lambda1+lambda3-mu2+noname1)/2)
!p2min=MAX(0,noname1-mu2,(lambda3-lambda1-mu2+noname1)/2,(lambda1-lambda3-mu2+noname1)/2)
p2min=MAX(0,noname1-mu2,(lambda3-lambda1-mu2+noname1+1)/2,(lambda1-lambda3-mu2+noname1+1)/2)
!if(p2min/=MAX(0,noname1-mu2,(lambda3-lambda1-mu2+noname1)/2,(lambda1-lambda3-mu2+noname1)/2))then
!        print*,"p2min/=p2min"
!        stop
!end if
q2=noname1-p2min
IF(p2min/=p2max)THEN
  n=(lambda1+lambda2-lambda3+2*(mu1+mu2-mu3))/3
  a=(lambda2+lambda3-lambda1-n)/2
  b=(lambda3+lambda1-lambda2+n+2)/2
  c=(lambda1+lambda2-lambda3-n)/2
  d=(lambda1+lambda2+lambda3+2-n)/2
  j1=lambda2-p2max
  j2=lambda2-p2min
  DO p2=p2min,p2max
! The upper and lower bound on p2 are such that:
! 1) 0<=p2<=lambda2
! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
! 3) ABS(lambda1-Lambda22)<=lambda3<=lambda1+Lambda22, where Lambda22=mu2-(2*lambda2+mu2-epsilon2)/3+2*p2
    p2tilde=mu2-q2
    q2tilde=lambda2-p2
    IF(p2tilde==0)THEN
      F=1.D0
    ELSE
      F=DFLOAT((a+1)*(b-1))
      DO j=1,p2tilde-1
        F=F*DFLOAT((a+j+1)*(b-j-1))
      END DO
      i1=p2tilde*(p2tilde+1)/2
      DO i=1,p2tilde
        i1=i1+1
        scalprod=DFLOAT((p2+1)*(mu1+lambda2+mu2-n+2))
        DO j=1,p2tilde-1
          IF(j<i)THEN
            scalprod=scalprod*DFLOAT((p2+j+1)*(mu1+lambda2+mu2-n+j+2))
          ELSE
            scalprod=scalprod*DFLOAT((a+j+1)*(b-j-1))
          END IF
        END DO
        F=F+binom(i1)*scalprod
      END DO
!      F=DFLOAT((a+1)*(b-1))
!      DO j=1,p2tilde-1
!        F=F*DFLOAT((a+j+1)*(b-j-1))
!      END DO
!      CC=0.D0
!      DO i=1,p2tilde
!        scalprod=DFLOAT((p2+1)*(mu1+lambda2+mu2-n+2))
!        DO j=1,p2tilde-1
!          IF(j<i)THEN
!            scalprod=scalprod*DFLOAT((p2+j+1)*(mu1+lambda2+mu2-n+j+2))
!          ELSE
!            scalprod=scalprod*DFLOAT((a+j+1)*(b-j-1))
!          END IF
!        END DO
!        scalprod=binom(p2tilde,i)*scalprod
!        Y=scalprod-CC
!        T=F+Y
!        CC=T-F
!        CC=CC-Y
!        F=T
!      END DO
      IF(BTEST(p2tilde,0))F=-F
    END IF
    IF(j1<q2tilde)THEN
      G=DFLOAT(INT8(a+n-j1)*(b-n+j1)*(c+n-j1)*(d+n-j1)*(lambda2+mu2-j1+1))
    ELSE
      G=DFLOAT(mu2-n+j1+1)
    END IF
    DO j=j1+1,j2-1
      IF(j<q2tilde)THEN
        G=G*DFLOAT(INT8(a+n-j)*(b-n+j)*(c+n-j)*(d+n-j)*(lambda2+mu2-j+1))
      ELSE
        G=G*DFLOAT(mu2-n+j+1)
      END IF
    END DO
    G=G*DFLOAT(lambda2-2*q2tilde+n+1)*binom(n*(n+1)/2+q2tilde)
    i1=n+1+lambda2-q2tilde
    H=binom(i1*(i1+1)/2+lambda2-q2tilde)
    wigner(lambda1,p2,q2,rho)=F*DSQRT(G/H)
    q2=q2-1
  END DO
ELSE
  wigner(lambda1,p2min,q2,rho)=1.D0
END IF
!**************************************************************************************
! End of (20) and  beginning of (21)
!**************************************************************************************
!epsilon1=-lambda1-2*mu1-3
!noname2=(2*lambda1+mu1-epsilon1)/3
noname2=lambda1+mu1+1
steps21=(epsilon2+lambda2+2*mu2)/3 ! steps21 is the number of iterations in Eq.(21)
DO i2=1,steps21
! epsilon2=epsilon2-3 ! epsilon2 is epsilon2 in Eq.(21)
! epsilon1=epsilon1+3 ! epsilon1 is epsilon1 in Eq.(21)
  noname1=noname1+1
  noname2=noname2-1
! p1min=MAX(0,noname2-mu1,(lambda1+1-i2-mu1+noname2)/2)
  p1min=MAX(0,noname2-mu1,(lambda1-i2-mu1+noname2+2)/2)
  p1max=MIN(lambda1,noname2-1,(lambda1-1+i2-mu1+noname2)/2)
! p2min=MAX(0,noname1-mu2,(lambda2-steps21+i2-mu2+noname1)/2)
  p2min=MAX(0,noname1-mu2,(lambda2-steps21+i2-mu2+noname1+1)/2)
  q2=noname1-p2min
  Lambda22=mu2+p2min-q2-2
  DO p2=p2min,MIN(lambda2,noname1,(lambda2+steps21-i2-mu2+noname1)/2)
! The upper and lower bound on p2 are such that:
! 1) 0<=p2<=lambda2
! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
! 3) lambda2-steps21+i2<=Lambda22<=lambda2+steps21-i2, where Lambda22=mu2+p2-q2
    Lambda22=Lambda22+2 ! Lambda22 is 2*Lambda2 in Eq.(21)
    Sq2=INT8(q2)*(mu2+1-q2)*(lambda2+mu2+2-q2) ! Sq2 is S(q2)
    Rp2=INT8(p2)*(lambda2+1-p2)*(mu2+1+p2) ! Rp2 is R(p2)
    
    IF(lambda1>=i2)THEN

!if(p1min/=(lambda1-i2-mu1+noname2+1)/2)then
!  print*,"wigner_canonical_extremal: PROBLEM in (21)!"
!  stop
!end if

      IF(Lambda22+1<=lambda2+mu2)THEN
        IF(q2>0)THEN
          ABCD=(lambda1-i2+Lambda22-lambda3+2)*(lambda1-i2+Lambda22+lambda3+4)
          IF(ABCD>0)THEN
!!            wigner(p1min-1,p2,q2,rho)=-DSQRT(DFLOAT(Sq2*ABCD)/DFLOAT(Lambda22+2))*wigner(p1min,p2,q2-1,rho)/2.D0
wigner(p1min-1,p2,q2,rho)=-DSQRT(DFLOAT((lambda1-i2+1)*Sq2*ABCD)&
/DFLOAT(INT8(lambda1+2-i2)*(Lambda22+1)*p1min*(lambda1+1-p1min)*(mu1+1+p1min)*(Lambda22+2)))*wigner(p1min,p2,q2-1,rho)!/2.D0
!         ELSE
!           wigner(p1min-1,p2,q2,rho)=0.D0
          END IF
        END IF
!     ELSE
!       wigner(p1min-1,p2,q2,rho)=0.D0
      END IF
      IF((Lambda22>=1).AND.(p2>0))THEN
        ABCD=(Lambda22+lambda3-lambda1+i2)*(lambda3+lambda1-i2-Lambda22+2)
!!        IF(ABCD>0)wigner(p1min-1,p2,q2,rho)=wigner(p1min-1,p2,q2,rho)&
!!                  -DSQRT(DFLOAT(Rp2*ABCD)/DFLOAT(Lambda22))*wigner(p1min,p2-1,q2,rho)/2.D0
IF(ABCD>0)wigner(p1min-1,p2,q2,rho)=wigner(p1min-1,p2,q2,rho)&
-DSQRT(DFLOAT((lambda1-i2+1)*Rp2*ABCD)/DFLOAT(INT8(lambda1+2-i2)*(Lambda22+1)*p1min*(lambda1+1-p1min)*(mu1+1+p1min)&
*Lambda22))*wigner(p1min,p2-1,q2,rho)!/2.D0
      END IF
!!      wigner(p1min-1,p2,q2,rho)=wigner(p1min-1,p2,q2,rho)&
!!        *DSQRT(DFLOAT(lambda1-i2+1)/DFLOAT((lambda1+2-i2)*(Lambda22+1)*p1min*(lambda1+1-p1min)*(mu1+1+p1min)))
!wigner(p1min-1,p2,q2,rho)=wigner(p1min-1,p2,q2,rho)/2.D0
    END IF
    q1=noname2-p1min
    Lambda12=mu1+p1min-q1-2
    DO p1=p1min,p1max
! The upper and lower bound on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1 cannot be 0 to avoid division by 0
! 3) lambda1+1-i2<=Lambda12<=lambda1-1+i2, where Lambda12=mu1+p1-q1
      Lambda12=Lambda12+2
!     IF(Lambda12+1<=lambda1+mu1)THEN ! This IF is redundant
        IF(Lambda22+1<=lambda2+mu2)THEN
          IF(q2>0)THEN
            ABCD=(Lambda22+lambda3-Lambda12+1)*(lambda3+Lambda12-Lambda22+1)
            IF(ABCD>0)THEN
!!              wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(Sq2*ABCD)/DFLOAT(Lambda22+2))*wigner(p1,p2,q2-1,rho)/2.D0
wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT((Lambda12+2)*Sq2*ABCD)&
/DFLOAT(INT8(Lambda12+1)*(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda22+2)))*wigner(p1,p2,q2-1,rho)!/2.D0
!           ELSE
!             wigner(p1,p2,q2,rho)=0.D0
            END IF
          END IF
!       ELSE
!         wigner(p1,p2,q2,rho)=0.D0
        END IF
        IF((Lambda22>=1).AND.(p2>0))THEN
          ABCD=(Lambda12+Lambda22-lambda3+1)*(Lambda12+Lambda22+lambda3+3)
!!          IF(ABCD>0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!                    +DSQRT(DFLOAT(Rp2*ABCD)/DFLOAT(Lambda22))*wigner(p1,p2-1,q2,rho)/2.D0
IF(ABCD>0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
+DSQRT(DFLOAT((Lambda12+2)*Rp2*ABCD)/DFLOAT(INT8(Lambda12+1)*(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)&
*Lambda22))*wigner(p1,p2-1,q2,rho)!/2.D0
        END IF
!!        wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!          *DSQRT(DFLOAT(Lambda12+2)/DFLOAT((Lambda12+1)*(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)))
!wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)*0.5D0
!     END IF
      q1=q1-1
    END DO
    q2=q2-1
  END DO
END DO
!**************************************************************************************
! End of (21) and beginning of (17)
!**************************************************************************************
noname2=noname2-1
DO i1=1,eta
  lambda2=lambda2+1
  mu2=mu2+1
  noname2=noname2-1
  IF(noname2==mu1)THEN
    IF(wigner(1,lambda2-1,mu2-1,rho)/=0.D0)THEN ! This IF is needed.
      wigner(0,lambda2,mu2,rho)=-DSQRT(DFLOAT(INT8(lambda1)*(mu1+2)*(lambda2+lambda3)&
        *(lambda3-lambda2+2))/2.D0)*wigner(1,lambda2-1,mu2-1,rho)
!   ELSE
!     wigner(0,lambda2,mu2,rho)=0.D0
    END IF
  END IF
  p1min=MAX(0,noname2-mu1,(2-mu1+noname2)/2)
  q1=noname2-p1min
  Lambda12=mu1+p1min-q1-2
  DO p1=p1min,MIN(lambda1,noname2,(lambda1-1+noname2)/2)
! The upper and lower bound on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) 1<=Lambda12<=lambda1+mu1-1, where Lambda12=mu1+p1-q1
    Lambda12=Lambda12+2
    IF(wigner(p1+1,lambda2-1,mu2-1,rho)/=0.D0)THEN ! This IF is nedeed.

!IF(p1==noname2-mu1.OR.p1==lambda1)then
!        print*,"p1==noname2-mu1.OR.p1==lambda1"
!        stop
!end if
!epsilon1=-lambda3-2*mu3+lambda2+2*mu2

      wigner(p1,lambda2,mu2,rho)=-DSQRT(DFLOAT(INT8(p1+1)*(lambda1-p1)*(mu1+2+p1)*(lambda2+lambda3-Lambda12)&
        *(lambda3+Lambda12-lambda2+2))/DFLOAT((Lambda12+2)*(Lambda12+1)))*wigner(p1+1,lambda2-1,mu2-1,rho)
!   ELSE
!     wigner(p1,lambda2,mu2,rho)=0.D0
    END IF
    IF(wigner(p1,lambda2-1,mu2-1,rho)/=0.D0)then
            
!IF(p1==noname2-mu1.OR.p1==lambda1)then
!        print*,"p1==noname2-mu1.OR.p1==lambda1"
!        stop
!end if
            
      wigner(p1,lambda2,mu2,rho)=wigner(p1,lambda2,mu2,rho)&
      +DSQRT(DFLOAT(INT8(q1+1)*(mu1-q1)*(lambda1+mu1+1-q1)*(Lambda12+lambda2-lambda3)&
      *(Lambda12+lambda2+lambda3+2))/DFLOAT((Lambda12+1)*Lambda12))*wigner(p1,lambda2-1,mu2-1,rho)
    ! This IF is needed.
    end if
    q1=q1-1
  END DO
  IF(noname2==lambda1)THEN ! This IF is necessary.
    IF(wigner(lambda1,lambda2-1,mu2-1,rho)/=0.D0)THEN ! This IF is needed.
       Lambda12=lambda1+mu1
       wigner(lambda1,lambda2,mu2,rho)=DSQRT(DFLOAT(INT8(mu1)*(lambda1+mu1+1)*(Lambda12+lambda2-lambda3)&
        *(Lambda12+lambda2+lambda3+2))/DFLOAT((Lambda12+1)*Lambda12))*wigner(lambda1,lambda2-1,mu2-1,rho)
!   ELSE
!     wigner(p1,lambda2,mu2,rho)=0.D0
    END IF
  END IF
END DO
!epsilon1=epsilon1+3*(1+eta)!mozno netreba
!epsilon1=-lambda1-2*mu1+3*(steps21+eta)
!**************************************************************************************
! End of (17)
!**************************************************************************************
! This block fills the arrays p1a,p2a,q2a with values corresponding to coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)HW||(lambda3,mu3)HW>

! It might be possible to put this block at the very beginning before the loop over rho so tha the condition rho==rhomax is not tested.

IF(rho==rhomax)THEN
  i2=0
  p1min=MAX(0,noname2-mu1,(lambda3-mu1+noname2-lambda2+1)/2,(lambda2-lambda3-mu1+noname2+1)/2)
! q1=noname2-p1min
! Lambda12=mu1+p1min-q1-2
  DO p1=p1min,MIN(lambda1,noname2,(lambda2+lambda3-mu1+noname2)/2)
! The upper and lower bound on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-lambda2)<=lambda3<=Lambda12+lambda2, where Lambda12=mu1+p1-q1
!   Lambda12=Lambda12+2
    i2=i2+1
!   Lambda12a(i2)=Lambda12
!   Lambda22a(i2)=lambda2
!   epsilon2a(i2)=epsilon2
!   IF((epsilon2==lambda1+2*mu1-lambda3-2*mu3).AND.(Lambda12==lambda1))Lambda22max=MAX(Lambda22max,Lambda22)
    p1a(i2)=p1
    p2a(i2)=lambda2
    q2a(i2)=mu2
!   q1=q1-1
  END DO
! IF(epsilon1==-lambda1-2*mu1)Lambda22max=lambda2
  IF(-lambda1-2*mu1+3*(steps21+eta)==-lambda1-2*mu1)Lambda22max=lambda2
  ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW. This is needed for setting the phase.
END IF
!**************************************************************************************
! Beginning of (18)
!**************************************************************************************
!epsilon2=-lambda2-2*mu2
noname1=lambda2+mu2
DO i1=1,lambda2+mu2 ! This is loop over epsilon2=epsilon2HW+3,...,2*lambda2+mu2; i1 is (epsilon2-epsilon2HW)/3
! epsilon1=epsilon1-3 ! epsilon1 is epsilon1 in Eq.(18)
  noname2=noname2+1
! IF(epsilon1<-lambda1-2*mu1)EXIT
  IF(noname2>lambda1+mu1)EXIT
! epsilon2=epsilon2+3 ! epsilon2 is epsilon2+3 in Eq.(18)
  noname1=noname1-1
  p2min=MAX(0,noname1-mu2,(lambda2-i1-mu2+noname1+1)/2)
  q2=noname1-p2min
  Lambda22=mu2+p2min-q2-2 ! Lambda22 is 2*Lambda2' in Eq.(18)
  DO p2=p2min,MIN(lambda2,noname1,(lambda2+i1-mu2+noname1)/2)
    Lambda22=Lambda22+2
    p1min=MAX(0,noname2-mu1,(lambda3-mu1+noname2-Lambda22+1)/2,(Lambda22-lambda3-mu1+noname2+1)/2)
    q1=noname2-p1min
    Lambda12=mu1+p1min-q1-2
    DO p1=p1min,MIN(lambda1,noname2,(Lambda22+lambda3-mu1+noname2)/2)
! The upper and lower bound on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=lambda3<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
      Lambda12=Lambda12+2

!      IF(lambda2>=mu2)THEN
        IF((Lambda22==lambda2-i1).OR.(Lambda22==0))THEN
          IF(Lambda12>=1.and.p1/=0)THEN
!!            wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT((p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda22+lambda3-Lambda12+2)&
!!              *(lambda3+Lambda12-Lambda22)))/DFLOAT(4*Lambda12))*wigner(p1-1,p2+1,q2,rho)
wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8(Lambda22+1)*(p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda22+lambda3-Lambda12+2)&
*(lambda3+Lambda12-Lambda22)))&
/DFLOAT(INT8(Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)*4*Lambda12))*wigner(p1-1,p2+1,q2,rho)

!if(p1==0)then
!        print*,"p1=0"
!        stop
!end if

          ELSE
            wigner(p1,p2,q2,rho)=0.D0
          END IF
!!          IF(Lambda12+1<=lambda1+mu1)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            +DSQRT(DFLOAT(q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
!!            /DFLOAT(4*(Lambda12+2)))*wigner(p1,p2+1,q2,rho)

IF(Lambda12+1<=lambda1+mu1.and.q1/=0)THEN

wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
+DSQRT(DFLOAT(INT8(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
/DFLOAT(INT8(Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)*4*(Lambda12+2)))*wigner(p1,p2+1,q2,rho)

END IF

!!          wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)))
        ELSE
          IF(Lambda12>=1.and.p1/=0)THEN
!!            wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda12+Lambda22-lambda3)&
!!              *(Lambda12+Lambda22+lambda3+2))/DFLOAT(4*Lambda12))*wigner(p1-1,p2,q2+1,rho)
wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8(Lambda22+1)*p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda12+Lambda22-lambda3)&
*(Lambda12+Lambda22+lambda3+2))&
/DFLOAT(INT8(Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)*4*Lambda12))*wigner(p1-1,p2,q2+1,rho)
          ELSE
            wigner(p1,p2,q2,rho)=0.D0
          END IF
!!          IF(Lambda12+1<=lambda1+mu1)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            -DSQRT(DFLOAT(q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
!!             /DFLOAT(4*(Lambda12+2)))*wigner(p1,p2,q2+1,rho)

IF(Lambda12+1<=lambda1+mu1.and.q1/=0)THEN
wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
-DSQRT(DFLOAT(INT8(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
/DFLOAT(INT8(Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)*4*(Lambda12+2)))*wigner(p1,p2,q2+1,rho)
END IF

!!          wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)))
        END IF

!      ELSE ! This is not necassary. It seems to slightly increase the maximal error and slightly decrease the mean error.
           ! Since the effect on the maximal error seems to be more significant, the code after this ELSE should probably
           ! be removed along with IF(lambda2>=mu2)THEN above (and END IF down).

go to 10

        IF((Lambda22==lambda2+i1).OR.(Lambda22==lambda2+mu2))THEN
          IF(Lambda12>=1.and.p1/=0)THEN
!!            wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda12+Lambda22-lambda3)&
!!              *(Lambda12+Lambda22+lambda3+2))/DFLOAT(4*Lambda12))*wigner(p1-1,p2,q2+1,rho)
wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8(Lambda22+1)*p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda12+Lambda22-lambda3)&
*(Lambda12+Lambda22+lambda3+2))&
/DFLOAT(INT8(Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)*4*Lambda12))*wigner(p1-1,p2,q2+1,rho)
          ELSE
            wigner(p1,p2,q2,rho)=0.D0
          END IF
!!          IF(Lambda12+1<=lambda1+mu1)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            -DSQRT(DFLOAT(q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
!!            /DFLOAT(4*(Lambda12+2)))*wigner(p1,p2,q2+1,rho)
IF(Lambda12+1<=lambda1+mu1.and.q1/=0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
-DSQRT(DFLOAT(INT8(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda22+lambda3-Lambda12)*(lambda3+Lambda12-Lambda22+2))&
/DFLOAT(INT8(Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)*4*(Lambda12+2)))*wigner(p1,p2,q2+1,rho)
!!          wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*Lambda22*(q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)))
        ELSE
          IF(Lambda12>=1.and.p1/=0)THEN
!!            wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda22+lambda3-Lambda12+2)&
!!              *(lambda3+Lambda12-Lambda22))/DFLOAT(4*Lambda12))*wigner(p1-1,p2+1,q2,rho)
wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8(Lambda22+1)*p1*(lambda1+1-p1)*(mu1+1+p1)*(Lambda22+lambda3-Lambda12+2)&
*(lambda3+Lambda12-Lambda22))&
/DFLOAT(INT8(Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)*4*Lambda12))*wigner(p1-1,p2+1,q2,rho)
          ELSE
            wigner(p1,p2,q2,rho)=0.D0
          END IF
!!          IF(Lambda12+1<=lambda1+mu1)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            +DSQRT(DFLOAT(q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
!!            /DFLOAT(4*(Lambda12+2)))*wigner(p1,p2+1,q2,rho)
IF(Lambda12+1<=lambda1+mu1.and.q1/=0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
+DSQRT(DFLOAT(INT8(Lambda22+1)*q1*(mu1+1-q1)*(lambda1+mu1+2-q1)*(Lambda12+Lambda22-lambda3+2)*(Lambda12+Lambda22+lambda3+4))&
/DFLOAT(INT8(Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)*4*(Lambda12+2)))*wigner(p1,p2+1,q2,rho)
!!          wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
!!            *DSQRT(DFLOAT(Lambda22+1)/DFLOAT((Lambda12+1)*(Lambda22+2)*(p2+1)*(lambda2-p2)*(mu2+2+p2)))
        END IF

10 continue

!      END IF

      IF(rho==rhomax)THEN
        i2=i2+1
!       Lambda12a(i2)=Lambda12
!       Lambda22a(i2)=Lambda22
!       epsilon2a(i2)=epsilon2
!       IF((epsilon2==lambda1+2*mu1-lambda3-2*mu3).AND.(Lambda12==lambda1))Lambda22max=MAX(Lambda22max,Lambda22)
        IF((p1==lambda1).AND.(q1==mu1))Lambda22max=MAX(Lambda22max,Lambda22)        
        ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW. This is needed for setting the phase.
        p1a(i2)=p1
        p2a(i2)=p2
        q2a(i2)=q2
      END IF

      q1=q1-1
    END DO
    q2=q2-1
  END DO
END DO
!**************************************************************************************
! End of (18)
!**************************************************************************************
! Valid values of p1, p2 and q2 are elements of arrays p1a, p2a and q2a with indeces from 1 to i2.
END DO ! end of the loop over rho

!do i1=1,i2
!p1=p1a(i1)
!p2=p2a(i1)
!q2=q2a(i1)
!print*,p1,p2,q2,(wigner(p1,p2,q2,rho),rho=1,rhomax)
!end do

!**************************************************************************************
! Beginning of orthonormalization according to (8) with alpha3=(epsilon3,Lambda3)=HW
!**************************************************************************************
! Gram-Schmidt orthonormalization is performed, where Wigner coefficients for given rho represent a vector,
! whose components are indexed by alpha1 and alpha2, and the scalar product is defined in Eq.(8).
!!DO rho=1,rhomax
!!  DO i4=1,rho-1
!!    scalprod=0.D0
!!    DO i1=1,i2
!!      p1=p1a(i1)
!!      p2=p2a(i1)
!!      q2=q2a(i1)
!!      scalprod=scalprod+wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,i4)
!!    END DO
!!    DO i1=1,i2
!!      p1=p1a(i1)
!!      p2=p2a(i1)
!!      q2=q2a(i1)
!!      wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)-scalprod*wigner(p1,p2,q2,i4)
!!    END DO
!!  END DO
!!  norm=0.D0
!!  DO i1=1,i2
!!    p1=p1a(i1)
!!    p2=p2a(i1)
!!    q2=q2a(i1)
!!    norm=norm+wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,rho)
!!  END DO
!!  DO i1=1,i2
!!    p1=p1a(i1)
!!    p2=p2a(i1)
!!    q2=q2a(i1)
!!    wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)/DSQRT(norm)
!!  END DO
!!END DO

!--------------------------------------------------------------------------------------
DO rho=1,rhomax
  a=1
!  IF((etamax==rhomax-1).AND.rho==rhomax)a=rhomax
  DO i4=a,rho
    p1=p1a(1)
    p2=p2a(1)
    q2=q2a(1)
    scalprod=wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,i4)
    F=0.D0 ! F->CC
    DO i1=2,i2 ! Kahan summation fomula (see Ref.[3]) is used to calculate the scalar product.
      p1=p1a(i1)
      p2=p2a(i1)
      q2=q2a(i1)
      G=wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,i4)-F ! G->Y
      H=scalprod+G ! H->T
      F=H-scalprod
      F=F-G
      scalprod=H
    END DO
    IF(i4/=rho)THEN
!      IF(DABS(scalprod)<1.D-12)CYCLE
      DO i1=1,i2
        p1=p1a(i1)
        p2=p2a(i1)
        q2=q2a(i1)
        wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)-scalprod*wigner(p1,p2,q2,i4)
      END DO
    ELSE
      scalprod=1.D0/DSQRT(scalprod)
      DO i1=1,i2
        p1=p1a(i1)
        p2=p2a(i1)
        q2=q2a(i1)
        wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)*scalprod
      END DO
    END IF
  END DO
END DO
!--------------------------------------------------------------------------------------

!**************************************************************************************
! End of the orthonormalization and beginning of setting the phase according to Ref.[2]
!**************************************************************************************
! The phase is such that <(lambda1,mu1)LW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)LW>_rho>0,
! which, according to formula 2B from Eq.(35) with rho_max instead of eta_max (see text therein), means that
! <(lambda1,mu1)HW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)HW>_rho*(-1)^e>0,
! where e=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax-rho+(lambda1+2*Lambda2_max-lambda3)/2

phiprhomax=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax ! phiprhomax is phi+rhomax
i4=phiprhomax+(lambda1+Lambda22max-lambda3)/2 ! i4 is e+rho

!if(noname1/=(2*lambda2+mu2-lambda1-2*mu1+lambda3+2*mu3)/3)then
!        print*,"nemame noname1"
!        stop
!end if

!nonamae1=(2*lambda2+mu2-lambda1-2*mu1+lambda3+2*mu3)/3 ! epsilon2=epsilon3_HW-epsilon1_HW
p2=(noname1-mu2+Lambda22max)/2
q2=noname1-p2
DO rho=1,rhomax
  IF(((BTEST(i4-rho+1,0)).AND.(wigner(lambda1,p2,q2,rho)<0.D0))&
     .OR.((BTEST(i4-rho,0)).AND.(wigner(lambda1,p2,q2,rho)>0.D0)))THEN
    DO i1=1,i2
      p1=p1a(i1)
      p2=p2a(i1)
      q2=q2a(i1)
      wigner(p1,p2,q2,rho)=-wigner(p1,p2,q2,rho)
    END DO
  END IF
END DO
!**************************************************************************************
! End of setting the phase
!**************************************************************************************
IF(I3==0)THEN ! E=LW
  DO i1=1,i2
    p1=p1a(i1)
    p2=p2a(i1)
!   q2=q2a(i1)
    DO rho=1,rhomax
      ! See Eq.(35,2B)
      IF(BTEST(phiprhomax-rho+p1+p2+(mu1-(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3+mu2-lambda3)/2,0))THEN
        q2=q2a(i1)
        wigner(p1,p2,q2,rho)=-wigner(p1,p2,q2,rho)
      END IF
    END DO
  END DO
END IF
END SUBROUTINE wigner_canonical_extremal
