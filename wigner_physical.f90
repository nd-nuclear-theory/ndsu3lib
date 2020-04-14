SUBROUTINE wigner_physical(I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,&
 I3,lambda3,mu3,L3,kappa3max,matrix3,rhomax,numb,wigner_can,Lambda12a,Lambda22a,epsilon2a,wigner_phys)
!---------------------------------------------------------------------------------------------------------------------------------
! Calculates reduced SU(3)-SO(3) Wigner coefficients <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
! for given lambda1,mu1,L1,lambda2,mu2,L2,lambda3,mu3,L3 using equations (31),(25),(5) in the reference.
!
! Reference: J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
!
! Input arguments: I1,J1,lambda1,mu1,L1,kappa1max,matrix1,I2,J2,lambda2,mu2,L2,kappa2max,matrix2,I3,lambda3,mu3,L3,kappa3max,matrix3,
!                  rhomax,numb,wigner_can,Lambda12a,Lambda22a,epsilon2a
! Output argument: wigner_phys
!
! (I,J)=(1,0) for lambda<mu => E=HW', (I,J)=(0,1) for lambda>=mu => E=LW'
! kappa1max = the number of occurences of L1 in SU(3) irrep (lambda1,mu1) (similarly kappa2max,kappa3max)
! matrix1(i,j) = element O_ij of the orthonormalization matrix for lambda1,mu1,L1 (similarly matrix2,matrix3)
! rhomax = the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3)
! numb = the number of extremal canonical reduced Wigner coefficients
!        <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)HW> excluding the outer multiplicity
! wigner_can = array containing the extremal canonical reduced Wigner coefficients (see wigner_canonical_extremal.f90 for details)
! Lambda12a = array containing the values of 2*Lambda1
! Lambda22a = array containing the values of 2*Lambda2
! epsilon2a = array containing the values of epsilon2
! wigner_phys(kappa1,kappa2,kappa3,rho) = <(lambda1,mu1)kappa1,L1;(lambda2,mu2)kappa2,L2||(lambda3,mu3)kappa3,L3>_rho
!
! Note: triangular inequality for L1,L2,L3 is not checked.
!---------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8),EXTERNAL :: transformation_coeff,DWR3 ! DWR3 needs to be replaced
INTEGER,DIMENSION(1) :: Lambda12a,Lambda22a,epsilon2a
REAL(KIND=8),DIMENSION(9,9) :: matrix1,matrix2,matrix3
REAL(KIND=8),DIMENSION(0:20,0:20,0:20,1:9) :: wigner_can
REAL(KIND=8),DIMENSION(9,9,9,9) :: wigner_phys
REAL(KIND=8),DIMENSION(9) :: transcoeff2
REAL(KIND=8) :: cg,transcoeff1,fact!,transcoeff2
INTEGER :: I1,J1,lambda1,mu1,L1,kappa1max,I2,J2,lambda2,mu2,L2,kappa2max,I3,lambda3,mu3,L3,kappa3max,rhomax,numb,&
           kappa1,kappa2,kappa3,rho,K1,K2,K3,K32,M1p,M2p,L12,L22,L32,i,M1pmin,epsilon1,epsilon2,epsilon3,&
           Lambda12,Lambda22,MLambda12,MLambda22,epsilon2min,epsilon2ind,j,K1min,K2min,MLambda12min,&
           MLambda32,Lambda32,K3minm2,MLambda12max

wigner_phys(1:kappa1max,1:kappa2max,1:kappa3max,1:rhomax)=0.D0

IF(I3==1)THEN ! G_3E=G_3HW'
  epsilon3=-lambda3-2*mu3
  K3minm2=Kmin(lambda3,mu3,L3)-2
  MLambda32=-lambda3 ! MLambda32 is 2*M_Lambda_3
  Lambda32=lambda3 ! Lambda32 is 2*Lambda_3
  epsilon2min=-lambda2-2*mu2
ELSE ! G_3E=G_3LW'
  epsilon3=2*lambda3+mu3
  K3minm2=Kmin(mu3,lambda3,L3)-2
  MLambda32=mu3 ! MLambda32 is 2*M_Lambda_3
  Lambda32=mu3 ! Lambda32 is 2*Lambda_3
  epsilon2min=-mu2-2*lambda2
END IF

L12=2*L1
L22=2*L2
L32=2*L3

IF(I1==1)THEN ! G_1E=G_1HW'
  K1min=Kmin(lambda1,mu1,L1)
ELSE ! G_1E=G_1LW'
  K1min=Kmin(mu1,lambda1,L1)
END IF
IF(I2==1)THEN ! G_2E=G_2HW'
  K2min=Kmin(lambda2,mu2,L2)
ELSE ! G_2E=G_2LW'
  K2min=Kmin(mu2,lambda2,L2)
END IF

DO i=1,numb ! sum over epsilon1,Lambda_1,epsilon2,Lambda_2
  epsilon2=epsilon2a(i)
  epsilon1=epsilon3-epsilon2
  Lambda12=Lambda12a(i)
  Lambda22=Lambda22a(i)
  IF(I3==1)THEN
    epsilon2ind=(epsilon2-epsilon2min)/3
  ELSE
    epsilon2ind=(-epsilon2-epsilon2min)/3
  END IF
  K3=K3minm2
  K32=2*K3
  MLambda12min=MAX(-Lambda12,MLambda32-Lambda22)
  MLambda12max=MIN(Lambda12,MLambda32+Lambda22)

  DO kappa3=1,kappa3max
    K3=K3+2 ! K3 is K_3
    K32=K32+4 ! K32 is 2*K_3
    M1pmin=MAX(-L1,K3-L2) ! M1pmin is the minimal M'_1
    IF(BTEST(Lambda12+M1pmin,0))M1pmin=M1pmin+1 ! M'_1 has to have the same parity as Lambda12 for <G_1|(G_1E)K1,L1,M'1> to be nonzero.
    M2p=K3-M1pmin+2
    IF(BTEST(Lambda22+M2p,0))EXIT ! This might not be necessary, because the condition is probably never satisfied.

    DO M1p=M1pmin,MIN(L1,K3+L2),2 ! M1p is M'_1. It has to have the same parity as Lambda12 for <G_1|(G_1E)K1,L1,M'1> to be nonzero.
       M2p=M2p-2 ! M2p is M'_2. It has to have the same parity as Lambda22 for <G_2|(G_2E)K2,L2,M'2> to be nonzero.
       cg=DWR3(L12,L22,L32,2*M1p,2*M2p,K32) ! cg is <L1 M'_1;L2 M'_2|L3 K3>
       MLambda22=MLambda32-MLambda12min+2

       DO MLambda12=MLambda12min,MLambda12max,2 ! MLambda12 is 2*M_Lambda_1
          MLambda22=MLambda22-2 ! MLambda22 is 2*M_Lambda_2

          DO kappa2=1,kappa2max
            ! Calculation of <G_2|(G_2E)K2,L2,M'2> = transcoeff2(kappa2)
            K2=K2min
            transcoeff2(kappa2)=matrix2(kappa2,1)*transformation_coeff(I2,J2,lambda2,mu2,epsilon2,Lambda22,MLambda22,K2,L2,M2p)
            DO j=2,kappa2
              K2=K2+2
!              IF(2*((Lambda22+M2p)/2)==Lambda22+M2p)
              transcoeff2(kappa2)=transcoeff2(kappa2)&
                +matrix2(kappa2,j)*transformation_coeff(I2,J2,lambda2,mu2,epsilon2,Lambda22,MLambda22,K2,L2,M2p)
            END DO
          END DO

          DO kappa1=1,kappa1max

            ! Calculation of <G_1|(G_1E)K1,L1,M'1> = transcoeff1
            K1=K1min
            transcoeff1=matrix1(kappa1,1)*transformation_coeff(I1,J1,lambda1,mu1,epsilon1,Lambda12,MLambda12,K1,L1,M1p)
            DO j=2,kappa1
              K1=K1+2
!              IF(2*((Lambda12+M1p)/2)==Lambda12+M1p)
              transcoeff1=transcoeff1&
                +matrix1(kappa1,j)*transformation_coeff(I1,J1,lambda1,mu1,epsilon1,Lambda12,MLambda12,K1,L1,M1p)
            END DO

            DO kappa2=1,kappa2max

!            ! Calculation of <G_2|(G_2E)K2,L2,M'2> = transcoeff2
!            transcoeff2=0.D0
!            K2=K2minm2
!            DO j=1,kappa2
!              K2=K2+2
!              IF(2*((Lambda22+M2p)/2)==Lambda22+M2p)transcoeff2=transcoeff2&
!                +matrix2(kappa2,j)*transformation_coeff(I2,J2,lambda2,mu2,epsilon2,Lambda22,MLambda22,K2,L2,M2p)
!            END DO
             
            fact=cg*transcoeff1*transcoeff2(kappa2)*DWR3(Lambda12,Lambda22,Lambda32,-MLambda12,-MLambda22,-MLambda32)
            ! Neviem, preco su v DWR3 tie minusy

            DO rho=1,rhomax
              wigner_phys(kappa1,kappa2,kappa3,rho)=wigner_phys(kappa1,kappa2,kappa3,rho)&
                                                    +fact*wigner_can(Lambda12,epsilon2ind,Lambda22,rho)
            END DO

          END DO
        END DO
      END DO
    END DO
  END DO
END DO

! Orthonormalization (w.r.t. K3)
DO rho=1,rhomax
  DO kappa1=1,kappa1max
    DO kappa2=1,kappa2max
      DO kappa3=kappa3max,1,-1
        wigner_phys(kappa1,kappa2,kappa3,rho)=matrix3(kappa3,kappa3)*wigner_phys(kappa1,kappa2,kappa3,rho)
        DO j=kappa3-1,1,-1
          wigner_phys(kappa1,kappa2,kappa3,rho)=wigner_phys(kappa1,kappa2,kappa3,rho)&
                                                +matrix3(kappa3,j)*wigner_phys(kappa1,kappa2,j,rho)
        END DO
      END DO
    END DO
  END DO
END DO

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
