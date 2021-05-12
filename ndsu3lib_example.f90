!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ndsu3lib_example.f90 -- simple program demonstrating usage of ndsu3lib
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM ndsu3lib_example
!------------------------------------------------------------------------------------------------------------
! Simple program calculating SU(3) reduced Wigner and recoupling coefficients and printing them on the screen
!------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity,inner_multiplicity
INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3_twice,dimpq,dimw,rhomax,numb,i,p1,q1,p2,q2,&
           epsilon1,Lambda1_twice,epsilon2,Lambda2_twice,rho,L1,L2,L3,kappa1max,kappa2max,kappa3max,kappa1,&
           kappa2,kappa3,lambda,mu,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,info,rhoa,rhob,&
           rhoc,rhod,n,lambda13,mu13,lambda4,mu4,lambda24,mu24,lambda34,mu34,rhomax12,rhomax34,rhomax1234,&
           rhomax13,rhomax24,rhomax1324,rho12,rho34,rho1234,rho13,rho24,rho1324
INTEGER,ALLOCATABLE,DIMENSION(:) :: p1a,p2a,q2a
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: wigner_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wigner_phys
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: rac,Zcoeff
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:,:,:) :: ninelm

INTERFACE ! The following subroutines require INTERFACE because they accept assumed-shape arrays.
  SUBROUTINE calculate_wigner_su3so3(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,&
                                     kappa3max,rhomax,wigner_phys)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,kappa3max,rhomax
    REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
  END SUBROUTINE calculate_wigner_su3so3
  SUBROUTINE calculate_u_coeff(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                               rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                          rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: rac
  END SUBROUTINE calculate_u_coeff
  SUBROUTINE calculate_z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                               rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                          rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff
  END SUBROUTINE calculate_z_coeff
  SUBROUTINE calculate_9_lambda_mu(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                                   lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,ninelm,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                          lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:,:,:,:,:),INTENT(OUT) :: ninelm
  END SUBROUTINE calculate_9_lambda_mu
END INTERFACE

CALL ndsu3lib_init(.TRUE.,400)

!************************************************************************************************************
! SU(3)-SU(2)xU(1) reduced Wigner coefficients <(6,1)epsilon1,Lambda1;(2,1)epsilon2,Lambda2||(5,2)-9,5/2>_rho
!************************************************************************************************************

lambda1=6
mu1=1
lambda2=2
mu2=1
lambda3=5
mu3=2
epsilon3=-9
Lambda3_twice=5 ! Lambda3_twice is 2*Lambda3

rhomax=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda3,mu3)

dimpq=(MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
dimw=rhomax*dimpq

ALLOCATE(wigner_block(dimw),p1a(dimpq),p2a(dimpq),q2a(dimpq))

CALL calculate_wigner_canonical(lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3_twice,dimpq,dimw,&
                                rhomax,numb,wigner_block,p1a,p2a,q2a)

WRITE(*,"(A)")"SU(3)-SU(2)xU(1) reduced Wigner coefficients <(6,1)epsilon1,Lambda1;(2,1)epsilon2,Lambda2||(5,2)-9,5/2>_rho:"
WRITE(*,"(A,I1)")"epsilon1  2*Lambda1  epsilon2  2*Lambda2  coefficients for rho=1,...,rhomax=",rhomax
DO i=1,numb
  p1=p1a(i)
  p2=p2a(i)
  q2=q2a(i)
  q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
  epsilon2=2*lambda2+mu2-3*(p2+q2)
  epsilon1=epsilon3-epsilon2
  Lambda1_twice=mu1+p1-q1 ! Lambda1_twice is 2*Lambda1
  Lambda2_twice=mu2+p2-q2 ! Lambda2_twice is 2*Lambda2
  WRITE(*,"(X,I4,7X,I3,7X,I4,7X,I3,3X,2D25.16)")epsilon1,Lambda1_twice,epsilon2,Lambda2_twice,&
                                                (wigner_block(i+numb*(rho-1)),rho=1,rhomax)
END DO

DEALLOCATE(wigner_block,p1a,p2a,q2a)

!*****************************************************************************************
! SU(3)-SO(3) reduced Wigner coefficients <(6,1)kappa1,2;(2,1)kappa2,3||(5,2)kappa3,3>_rho
!*****************************************************************************************

L1=2
L2=3
L3=3
kappa1max=inner_multiplicity(lambda1,mu1,L1)
kappa2max=inner_multiplicity(lambda2,mu2,L2)
kappa3max=inner_multiplicity(lambda3,mu3,L3)

ALLOCATE(wigner_phys(kappa1max,kappa2max,kappa3max,rhomax))

CALL calculate_wigner_su3so3(lambda1,mu1,L1,kappa1max,lambda2,mu2,L2,kappa2max,lambda3,mu3,L3,&
                             kappa3max,rhomax,wigner_phys)

WRITE(*,*)
WRITE(*,"(A)")"SU(3)-SO(3) reduced Wigner coefficients <(6,1)kappa1,2;(2,1)kappa2,3||(5,2)kappa3,3>_rho:"
WRITE(*,"(A,I1)")"kappa1  kappa2  kappa3  coefficients for rho=1,...,rhomax=",rhomax
DO kappa1=1,kappa1max
  DO kappa2=1,kappa2max
    DO kappa3=1,kappa3max
      WRITE(*,"(2X,I2,6X,I2,6X,I2,X,2D25.16)")kappa1,kappa2,kappa3,(wigner_phys(kappa1,kappa2,kappa3,rho),rho=1,rhomax)
    END DO
  END DO
END DO

DEALLOCATE(wigner_phys)

!***********************************************************************************
! SU(3) recoupling coefficients U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
!***********************************************************************************

lambda1=9
mu1=3
lambda2=1
mu2=1
lambda=6
mu=6
lambda3=2
mu3=2
lambda12=9
mu12=3
lambda23=3
mu23=3

rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
rhomaxc=outer_multiplicity(lambda2,mu2,lambda3,mu3,lambda23,mu23)
rhomaxd=outer_multiplicity(lambda1,mu1,lambda23,mu23,lambda,mu)

ALLOCATE(rac(rhomaxd,rhomaxa*rhomaxb*rhomaxc))

CALL calculate_u_coeff(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                       rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,rhomaxd,info)

WRITE(*,*)
IF(info==0)THEN
  WRITE(*,"(A)")"SU(3) recoupling coefficients U[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]:"
  WRITE(*,"(A)")"rhoa  rhob  rhoc  rhod  coefficient"
  DO rhoa=1,rhomaxa
    DO rhob=1,rhomaxb
      DO rhoc=1,rhomaxc
        DO rhod=1,rhomaxd
          n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
          WRITE(*,"(X,I2,4X,I2,4X,I2,4X,I2,D25.16)")rhoa,rhob,rhoc,rhod,rac(rhod,n)
        END DO
      END DO
    END DO
  END DO
ELSE
  WRITE(*,"(A,I2)")"calculate_u_coeff: MKL subroutine dgesv ran with error: info=",info
END IF

DEALLOCATE(rac)

!***********************************************************************************
! SU(3) recoupling coefficients Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]
!***********************************************************************************

lambda2=9
mu2=3
lambda1=1
mu1=1
lambda=6
mu=6
lambda3=2
mu3=2
lambda12=9
mu12=3
lambda13=3
mu13=3

rhomaxa=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
rhomaxb=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda,mu)
rhomaxc=outer_multiplicity(lambda1,mu1,lambda3,mu3,lambda13,mu13)
rhomaxd=outer_multiplicity(lambda13,mu13,lambda2,mu2,lambda,mu)

ALLOCATE(Zcoeff(rhomaxd,rhomaxa*rhomaxb*rhomaxc))

CALL calculate_z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                       rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,rhomaxd,info)

WRITE(*,*)
IF(info==0)THEN
  WRITE(*,"(A)")"SU(3) recoupling coefficients Z[(9,3)(1,1)(6,6)(2,2);(9,3)rhoa,rhob(3,3)rhoc,rhod]:"
  WRITE(*,"(A)")"rhoa  rhob  rhoc  rhod  coefficient"
  DO rhoa=1,rhomaxa
    DO rhob=1,rhomaxb
      DO rhoc=1,rhomaxc
        DO rhod=1,rhomaxd
          n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
          WRITE(*,"(X,I2,4X,I2,4X,I2,4X,I2,D25.16)")rhoa,rhob,rhoc,rhod,Zcoeff(rhod,n)
        END DO
      END DO
    END DO
  END DO
ELSE
  WRITE(*,"(A,I2)")"calculate_z_coeff: MKL subroutine dgesv ran with error: info=",info
END IF

DEALLOCATE(Zcoeff)

!********************************************************
!                            |(1,1) (0,0)  (1,1)  rho12 |
! 9-(lambda,mu) coefficients |(0,4) (4,2)  (2,4)  rho34 |
!                            |(2,3) (4,2)  (2,4) rho1324|
!                            |rho13 rho24 rho1234       |
!********************************************************

lambda1=1
mu1=1
lambda2=0
mu2=0
lambda12=1
mu12=1
lambda3=0
mu3=4
lambda4=4
mu4=2
lambda34=2
mu34=4
lambda13=2
mu13=3
lambda24=4
mu24=2
lambda=2
mu=4

rhomax12=outer_multiplicity(lambda1,mu1,lambda2,mu2,lambda12,mu12)
rhomax34=outer_multiplicity(lambda3,mu3,lambda4,mu4,lambda34,mu34)
rhomax1234=outer_multiplicity(lambda12,mu12,lambda34,mu34,lambda,mu)
rhomax13=outer_multiplicity(lambda1,mu1,lambda3,mu3,lambda13,mu13)
rhomax24=outer_multiplicity(lambda2,mu2,lambda4,mu4,lambda24,mu24)
rhomax1324=outer_multiplicity(lambda13,mu13,lambda24,mu24,lambda,mu)

ALLOCATE(ninelm(rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324))

CALL calculate_9_lambda_mu(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                           lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,ninelm,info)

WRITE(*,*)
IF(info==0)THEN
  WRITE(*,"(27X,A)")"|(1,1) (0,0)  (1,1)  rho12 |"
  WRITE(*,"(A)")"9-(lambda,mu) coefficients |(0,4) (4,2)  (2,4)  rho34 |:"
  WRITE(*,"(27X,A)")"|(2,3) (4,2)  (2,4) rho1324|"
  WRITE(*,"(27X,A)")"|rho13 rho24 rho1234       |"
  WRITE(*,"(A)")"rho12  rho34  rho1234  rho13  rho24  rho1324  coefficient"
  DO rho12=1,rhomax12
    DO rho34=1,rhomax34
      DO rho1234=1,rhomax1234
        DO rho13=1,rhomax13
          DO rho24=1,rhomax24
            DO rho1324=1,rhomax1324
              WRITE(*,"(X,I2,5X,I2,6X,I2,6X,I2,5X,I2,6X,I2,2X,D25.16)")rho12,rho34,rho1234,rho13,rho24,rho1324,&
                                                                       ninelm(rho12,rho34,rho1234,rho13,rho24,rho1324)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
ELSE
  WRITE(*,"(A,I2)")"calculate_9_lambda_mu: MKL subroutine dgesv ran with error: info=",info
END IF

DEALLOCATE(ninelm)

CALL ndsu3lib_free(.TRUE.)

END PROGRAM ndsu3lib_example
