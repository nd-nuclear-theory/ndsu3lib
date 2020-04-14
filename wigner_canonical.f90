SUBROUTINE wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,&
                rhomax,numb,wignerex,wigner,Lambda12a,Lambda22a,epsilon2a)
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: su2racah,DRR3
INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda32,rhomax,eps3,steps,rho,p3,q3,p1,q1,p2,q2,Lam32,epsilon2max,&
           epsilon2min,epsilon3min,epsilon2,Lambda22,epsilon1min,epsilon1max,Lambda12,numb,epsilon1,Lam32prime,i2
REAL(KIND=8) :: N3
INTEGER,DIMENSION(1) :: Lambda12a,Lambda22a,epsilon2a
REAL(KIND=8),DIMENSION(0:20,0:20,0:20,1:9) :: wignerex,wigner

p3=p(lambda3,mu3,epsilon3,Lambda32)
IF((p3<0).OR.(p3>lambda3))THEN
  PRINT*,"INVALID VALUES OF epsilon3 AND Lambda3 "
  RETURN
END IF
q3=q(lambda3,mu3,epsilon3,Lambda32)
IF((q3<0).OR.(q3>mu3))THEN
  PRINT*,"INVALID VALUES OF epsilon3 AND Lambda3 "
  RETURN
END IF

wigner=wignerex ! Toto by asi bolo treba urobit efektivnejsie, aby sa neprepisovali nepouzivane elementy pola.

IF(epsilon3==-lambda3-2*mu3)RETURN

epsilon2max=2*lambda2+mu2
epsilon2min=-lambda2-2*mu2
epsilon3min=-lambda3-2*mu3+3
epsilon1min=-lambda1-2*mu1
epsilon1max=2*lambda1+mu1

numb=0
DO rho=1,rhomax

  Lam32=lambda3

  DO eps3=epsilon3min,epsilon3,3 ! eps3 is epsilon3 in Eq.(19)
    
	Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
    IF(Lam32<Lambda32)THEN
      Lam32=Lam32+1
	  N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT(S(lambda3,mu3,q(lambda3,mu3,eps3,Lam32)+1)))
    ELSE
      Lam32=Lam32-1
	  N3=DSQRT(DFLOAT((Lam32prime+1)*(Lam32+1))/DFLOAT(R(lambda3,mu3,p(lambda3,mu3,eps3,Lam32)+1)))
    END IF ! Lam32 is 2*Lambda3 in Eq.(19) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.
	
	DO p2=0,lambda2
	  DO q2=0,mu2
	    epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
		epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
		IF((epsilon1<epsilon1min).OR.(epsilon1>epsilon1max))CYCLE
		Lambda22=mu2+p2-q2 ! Lambda22 is 2*Lambda2 in Eq.(19)
		i2=(epsilon2-epsilon2min)/3
		DO p1=0,lambda1
		  q1=(epsilon1max-epsilon1)/3-p1
		  IF((q1<0).OR.(q1>mu1))CYCLE
		  Lambda12=mu1+p1-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
		  IF((ABS(Lambda12-Lambda22)>Lam32).OR.(Lambda12+Lambda22<Lam32))CYCLE
		  
		  IF((rho==1).AND.(eps3==epsilon3))THEN
		    numb=numb+1
			Lambda12a(numb)=Lambda12
			epsilon2a(numb)=epsilon2
			Lambda22a(numb)=Lambda22
		  END IF
		  
		  ! Tam, kde sa pocitaju SU(2) recoupling koeficienty, by mozno bolo dobre testovat trojuholnikove nerovnosti.
		  
		  IF(p1/=lambda1)THEN
		    wigner(Lambda12,i2,Lambda22,rho)=DSQRT(DFLOAT(R(lambda1,mu1,p1+1)))*DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12+1,Lam32)&
			*wigner(Lambda12+1,i2,Lambda22,rho)
		  ELSE
		    wigner(Lambda12,i2,Lambda22,rho)=0.D0
		  END IF
		  
		  IF(q1/=mu1)wigner(Lambda12,i2,Lambda22,rho)=wigner(Lambda12,i2,Lambda22,rho)&
	+DSQRT(DFLOAT(S(lambda1,mu1,q1+1)))*DRR3(Lambda22,Lam32prime,Lambda12,1,Lambda12-1,Lam32)*wigner(Lambda12-1,i2,Lambda22,rho)
		  
		  IF(p2/=lambda2)wigner(Lambda12,i2,Lambda22,rho)=wigner(Lambda12,i2,Lambda22,rho)&
	+DSQRT(DFLOAT(R(lambda2,mu2,p2+1)))*DRR3(Lambda12,Lambda22+1,Lam32,1,Lam32prime,Lambda22)*wigner(Lambda12,i2-1,Lambda22+1,rho)
		  
		  IF(q2/=mu2)wigner(Lambda12,i2,Lambda22,rho)=wigner(Lambda12,i2,Lambda22,rho)&
	+DSQRT(DFLOAT(S(lambda2,mu2,q2+1)))*DRR3(Lambda12,Lambda22-1,Lam32,1,Lam32prime,Lambda22)*wigner(Lambda12,i2-1,Lambda22-1,rho)
		  
		  wigner(Lambda12,i2,Lambda22,rho)=wigner(Lambda12,i2,Lambda22,rho)*N3
	    
		END DO
	  END DO
	END DO
  
  END DO

END DO

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

END SUBROUTINE wigner_canonical
