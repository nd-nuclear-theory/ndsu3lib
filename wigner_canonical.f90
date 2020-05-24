SUBROUTINE wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,&
                            rhomax,numb,wignerex,wigner,p1a,p2a,q2a)!Lambda12a,Lambda22a,epsilon2a)
!--------------------------------------------------------------------------------------------------------------------------------
! Calculates the reduced SU(3)-SU(2)xU(1) Wigner coefficients
! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3 using Eq.(19) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,rhomax,wignerex
! Output arguments: numb,p1a,p2a,q2a,wigner
!
! Lambda32=2*Lambda3 (epsilon3 and Lambda32 must be valid)
! rhomax=the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
! wignerex is an array containing the highest-weight coefficients calculated by subroutine wigner_canonical_extremal with I3=1.
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
!
! Note: There is a typo in Eq.(19). In the very last line there should be p instead of q.
!-------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: su2racah,DRR3 ! TO DO: Replace DRR3 with su2racah
INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,rhomax,eps3,rho,p1,q1,p2,q2,Lam32,epsilon2max,&
           epsilon3min,Lambda22,Lambda12,numb,Lam32prime,p3,q3,q2min,p1min,pq1,noname1,noname2,n1,n2
REAL(KIND=8) :: N3
!INTEGER,DIMENSION(1) :: Lambda12a,Lambda22a,epsilon2a
INTEGER,DIMENSION(1) :: p1a,p2a,q2a ! Dimension is at least (lambda1+1)*(lambda2+1)*(mu2+1)
REAL(KIND=8),DIMENSION(0:15,0:15,0:15,1:9) :: wignerex,wigner
! Why doesn't assumed-shaped array work?

!p3=p(lambda3,mu3,epsilon3,Lambda32)
!IF((p3<0).OR.(p3>lambda3))THEN
!  PRINT*,"INVALID VALUES OF epsilon3 AND Lambda3 "
!  RETURN
!END IF
!q3=q(lambda3,mu3,epsilon3,Lambda32)
!IF((q3<0).OR.(q3>mu3))THEN
!  PRINT*,"INVALID VALUES OF epsilon3 AND Lambda3 "
!  RETURN
!END IF

!wigner=wignerex ! This should be done more effectively, so that unused elements of the array are not copied. 
                 ! However, with automatic arrays this will be the best way of doing it.
DO p1=0,lambda1
  DO p2=0,lambda2
    DO q2=0,mu2
      DO rho=1,rhomax
        wigner(p1,p2,q2,rho)=wignerex(p1,p2,q2,rho)
      END DO
    END DO
  END DO
END DO

epsilon3min=-lambda3-2*mu3
IF(epsilon3==epsilon3min)RETURN

epsilon2max=2*lambda2+mu2
n1=(epsilon2max-lambda1-2*mu1-epsilon3min)/3
n2=(epsilon2max+2*lambda1+mu1-epsilon3min)/3
epsilon3min=epsilon3min+3

numb=0
DO rho=1,rhomax

  Lam32=lambda3
  p3=lambda3
  q3=mu3
  noname1=n1
  noname2=n2

  DO eps3=epsilon3min,epsilon3,3 ! eps3 is epsilon3 in Eq.(19)
    noname1=noname1-1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
    noname2=noname2-1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3
    
    Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
    IF(Lam32<Lambda32)THEN
      Lam32=Lam32+1
      q3=q3-1 ! p3 and q3 correspond to eps3 and Lam32
      N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((q3+1)*(mu3-q3)*(lambda3+mu3+1-q3)))
    ELSE
      Lam32=Lam32-1
      p3=p3-1 ! p3 and q3 correspond to eps3 and Lam32
      N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT((p3+1)*(lambda3-p3)*(mu3+2+p3)))
    END IF ! Lam32 is 2*Lambda3 in Eq.(19) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.

    DO p2=0,lambda2
      q2min=MAX(0,noname1-p2)
      pq1=noname2-p2-q2min+1
      Lambda22=mu2+p2-q2min+1
      DO q2=q2min,MIN(mu2,noname2-p2)
! The lower and upper bounds on q2 are such that:
! 1) 0<=q2<=mu2
! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
!       epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
!       epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
        pq1=pq1-1 ! pq1 is p1+q1
        Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
        p1min=MAX(0,pq1-mu1,(Lam32-mu1+pq1-Lambda22)/2,(Lambda22-Lam32-mu1+pq1)/2)
        q1=pq1-p1min+1
        Lambda12=mu1+p1min-q1-1
        DO p1=p1min,MIN(lambda1,pq1,(Lambda22+Lam32-mu1+pq1)/2)
! The lower and upper bounds on p1 are such that:
! 1) 0<=p1<=lambda1
! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
          q1=q1-1
          Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
  
          IF((rho==1).AND.(eps3==epsilon3))THEN
            numb=numb+1
!           Lambda12a(numb)=Lambda12
!           epsilon2a(numb)=epsilon2
!           Lambda22a(numb)=Lambda22
            p1a(numb)=p1
            p2a(numb)=p2
            q2a(numb)=q2
          END IF
  
          ! Before SU(2) recoupling coefficients are calculated, it might be good to test triangular inequalities.
  
          IF(q1/=mu1)THEN
            wigner(p1,p2,q2,rho)=DSQRT(DFLOAT((q1+1)*(mu1-q1)*(lambda1+mu1+1-q1)))&
              *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12-1,Lam32)*wigner(p1,p2,q2,rho)
          ELSE
            wigner(p1,p2,q2,rho)=0.D0
          END IF

          IF(p1/=lambda1)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)+DSQRT(DFLOAT((p1+1)*(lambda1-p1)*(mu1+2+p1)))&
            *DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12+1,Lam32)*wigner(p1+1,p2,q2,rho)
  
          IF(p2/=lambda2)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)+DSQRT(DFLOAT((p2+1)*(lambda2-p2)*(mu2+2+p2)))&
            *DRR3(Lambda12,Lambda22+1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2+1,q2,rho)
  
          IF(q2/=mu2)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)+DSQRT(DFLOAT((q2+1)*(mu2-q2)*(lambda2+mu2+1-q2)))&
            *DRR3(Lambda12,Lambda22-1,Lam32,1,Lam32prime,Lambda22)*wigner(p1,p2,q2+1,rho)
  
          wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)*N3
    
        END DO
      END DO
    END DO
  END DO
END DO

END SUBROUTINE wigner_canonical
