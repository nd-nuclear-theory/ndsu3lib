!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wigner_physical.f90 -- SU(3)-SO(3) coupling coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wigner_physical(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
 I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,p1a,p2a,q2a,wigner_phys)
!---------------------------------------------------------------------------------------------------------------------------------
! Calculates reduced SU(3)-SO(3) Wigner coefficients <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3 using equations (31),(25),(5) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,I3,lambda3,mu3,L3,kappa3max,matrix3,
!                  rhomax,numb,wigner_can,p1a,p2a,q2a
! Output argument: wigner_phys
!
! (I,J)=(1,0) for lambda<mu => E=HW', (I,J)=(0,1) for lambda>=mu => E=LW'
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
!---------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: DWR3,transformation_coeff ! DWR3 need to be replaced
INTEGER,INTENT(IN) :: I1,J1,lambda1,mu1,L1,kappa1max,I2,J2,lambda2,mu2,L2,kappa2max,I3,lambda3,mu3,L3,kappa3max,rhomax,numb
INTEGER :: kappa1,kappa2,kappa3,rho,K1,K2,K3,K32,M1p,M2p,L12,L22,L32,i,M1pmin,epsilon1,epsilon2,epsilon3,&
           Lambda12,Lambda22,MLambda12,MLambda22,j,K1min,K2min,MLambda12min,MLambda32,Lambda32,K3minm2,MLambda12max,&
           p1,p2,q1,q2,hovadina,Lambda12pM1p,Lambda22pM2p,phase331a,phase332a,phase331b,phase332b
INTEGER,DIMENSION(:),INTENT(IN) :: p1a,p2a,q2a
REAL(KIND=8),DIMENSION(:,:),INTENT(IN) :: matrix1,matrix2,matrix3
REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wigner_can
REAL(KIND=8),DIMENSION(:,:,:,:),INTENT(OUT) :: wigner_phys
!REAL(KIND=8),DIMENSION(kappa1max) :: transcoeff1
!REAL(KIND=8),DIMENSION(kappa2max) :: transcoeff2
REAL(KIND=8) :: cg1,cg2,fact,aux

REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: transcoeff1,transcoeff2
ALLOCATE(transcoeff1(kappa1max),transcoeff2(kappa2max))

wigner_phys(1:kappa1max,1:kappa2max,1:kappa3max,1:rhomax)=0.D0

IF(I3==1)THEN ! G_3E=G_3HW'
  epsilon3=-lambda3-2*mu3
  K3minm2=Kmin(lambda3,mu3,L3)-2
  MLambda32=-lambda3 ! MLambda32 is 2*M_Lambda_3
  Lambda32=lambda3 ! Lambda32 is 2*Lambda_3
  hovadina=(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3
ELSE ! G_3E=G_3LW'
  epsilon3=2*lambda3+mu3
  K3minm2=Kmin(mu3,lambda3,L3)-2
  MLambda32=mu3 ! MLambda32 is 2*M_Lambda_3
  Lambda32=mu3 ! Lambda32 is 2*Lambda_3
  hovadina=(2*(mu1+mu2+lambda3)+lambda1+lambda2+mu3)/3
END IF

L12=2*L1
L22=2*L2
L32=2*L3

IF(I1==1)THEN ! G_1E=G_1HW'
  K1min=Kmin(lambda1,mu1,L1)
  phase331a=2*(lambda1+2*mu1)+6*(L1+K1min)
ELSE ! G_1E=G_1LW'
  K1min=Kmin(mu1,lambda1,L1)
  phase331a=-2*(2*lambda1+mu1)+6*(L1+K1min)
END IF
IF(I2==1)THEN ! G_2E=G_2HW'
  K2min=Kmin(lambda2,mu2,L2)
  phase332a=2*(lambda2+2*mu2)+6*(L2+K2min)
ELSE ! G_2E=G_2LW'
  K2min=Kmin(mu2,lambda2,L2)
  phase332a=-2*(2*lambda2+mu2)+6*(L2+K2min)
END IF

DO i=1,numb ! sum over epsilon1,Lambda_1,epsilon2,Lambda_2
  p1=p1a(i)
  p2=p2a(i)
  q2=q2a(i)
  IF(I3==1)THEN
    q1=hovadina-p1-p2-q2
    epsilon2=2*lambda2+mu2-3*(p2+q2)
    Lambda12=mu1+p1-q1
    Lambda22=mu2+p2-q2
  ELSE
    q1=hovadina-p1-p2-q2
    epsilon2=-2*mu2-lambda2+3*(p2+q2)
    Lambda12=lambda1+p1-q1
    Lambda22=lambda2+p2-q2
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
!    IF(BTEST(Lambda22+M2p,0))EXIT ! This might not be necessary, because the condition is probably never satisfied.

    M2p=K3-M1pmin ! M2p is M'_2
    Lambda12pM1p=Lambda12+M1pmin
    Lambda22pM2p=Lambda22+M2p
    DO M1p=M1pmin,MIN(L1,K3+L2,Lambda12,K3+Lambda22),2 ! M1p is M'_1.
    ! For <G_1|(G_1E)K1,L1,M'1> to be nonzero, M1p has to have the same parity
    ! as Lambda12 and abs(M1p) has to be less than or equal to Lambda12.
    ! This is taken care of in the lower and upper bounds on M1p.

!    M2p=K3-M1p ! M2p is M'_2.
    ! For <G_2|(G_2E)K2,L2,M'2> to be nonzero, M2p has to have the same parity
    ! as Lambda22 and abs(M2p) has to be less than or equal to Lambda22.
    ! This is taken care of in the lower and upper bounds on M1p.

    !Lambda12pM1p=Lambda12+M1p
    !Lambda22pM2p=Lambda22+M2p

      cg1=DWR3(L12,L22,L32,2*M1p,2*M2p,K32) ! cg1 is <L1 M'_1;L2 M'_2|L3 K3>

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
        DO kappa1=1,kappa1max
          transcoeff1(kappa1)=transformation_coeff(I1,J1,lambda1,mu1,epsilon1,Lambda12,MLambda12,K1,L1,M1p)
          K1=K1+2
        END DO
        DO kappa1=kappa1max,1,-1
          transcoeff1(kappa1)=matrix1(kappa1,kappa1)*transcoeff1(kappa1)
          DO j=1,kappa1-1
            transcoeff1(kappa1)=transcoeff1(kappa1)+matrix1(kappa1,j)*transcoeff1(j)
          END DO
        END DO

        ! Calculation of <G_2|(G_2E)K2,L2,M'2> = transcoeff2(kappa2)
        K2=K2min
        DO kappa2=1,kappa2max
          transcoeff2(kappa2)=transformation_coeff(I2,J2,lambda2,mu2,epsilon2,Lambda22,MLambda22,K2,L2,M2p)
          K2=K2+2
        END DO
        DO kappa2=kappa2max,1,-1
          transcoeff2(kappa2)=matrix2(kappa2,kappa2)*transcoeff2(kappa2)
          DO j=1,kappa2-1
            transcoeff2(kappa2)=transcoeff2(kappa2)+matrix2(kappa2,j)*transcoeff2(j)
          END DO
        END DO

        cg2=cg1*DWR3(Lambda12,Lambda22,Lambda32,-MLambda12,-MLambda22,-MLambda32) ! Neviem, preco su v DWR3 tie minusy

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

DEALLOCATE(transcoeff1,transcoeff2)

CONTAINS
  FUNCTION Kmin(lambda,mu,L) RESULT(res) ! This function calculates the lowest K for given lambda,mu,L as given by Eq.(4a) in the reference.
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda,mu,L
    INTEGER :: res
    res=MAX(0,L-mu)
    res=res+MOD(res+lambda,2) ! K and lambda must have the same parity.
    IF(res==0)res=res+MOD(L+mu,2)*2 ! If K=0, L and mu must have the same parity.
  END FUNCTION Kmin
END SUBROUTINE wigner_physical
