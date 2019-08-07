SUBROUTINE racah_wrap(su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,&
                      su3irrep23,KR0A,KR0B,KR0C,KR0D,RHS,INFO)
!--------------------------------------------------------------------------------------
! SUBROUTINE TO CALCULATE SU(3) RACAH COEFFICIENTS
! U((LAM1,MU1)(LAM2,MU2)(LAM,MU)(LAM3,MU3);(LAM12,MU12)RHOA,RHOB(LAM23,MU23)RHOC,RHOD)
!
! INPUT ARGUMENTS: su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,su3irrep23,KR0A,
!                  KR0B,KR0C,KR0D
! OUTPUT ARGUMENTS: RHS,INFO
!
! su3irrep1%lambda=LAM1
! su3irrep1%mu=MU1
! su3irrep2%lambda=LAM2
! su3irrep2%mu=MU2
! su3irrepx%lambda=LAM
! su3irrepx%mu=MU
! su3irrep3%lambda=LAM3
! su3irrep3%mu=MU3
! su3irrep12%lambda=LAM12
! su3irrep12%mu=MU12
! su3irrep23%lambda=LAM23
! su3irrep23%mu=MU23
! KR0A=MULTIPLICITY OF COUPLING (LAM1,MU1)x(LAM2,MU2)->(LAM12,MU12)
! KR0B=MULTIPLICITY OF COUPLING (LAM12,MU12)x(LAM3,MU3)->(LAM,MU)
! KR0C=MULTIPLICITY OF COUPLING (LAM2,MU2)x(LAM3,MU3)->(LAM23,MU23)
! KR0D=MULTIPLICITY OF COUPLING (LAM1,MU1)x(LAM23,MU23)->(LAM,MU)
! RHS(RHOD,N)=U((LAM1,MU1)(LAM2,MU2)(LAM,MU)(LAM3,MU3);(LAM12,MU12)RHOA,RHOB(LAM23,MU23)RHOC,RHOD)
!  WHERE N=RHOA+KR0A*(RHOB-1)+KR0A*KR0B*(RHOC-1)
! INFO=0 IFF MKL SUBROUTINE dgesv CALLED IN SUBROUTINE racah RAN WITHOUT ERRORS
!--------------------------------------------------------------------------------------
USE derived_types
IMPLICIT NONE
INTEGER :: KR0A,KR0B,KR0C,KR0D,INFO,KIMAX2,X1,X2,X3,NCE2A,NCE2B,NCE2D,NABC
REAL(KIND=8), DIMENSION(KR0D,KR0A*KR0B*KR0C) :: RHS
TYPE(su3irrep) :: su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,su3irrep23

KIMAX2=KR0C*IDM(su3irrep3)
X1=su3irrep3%lambda+su3irrep3%mu+1
X2=su3irrep2%lambda+su3irrep2%mu+1
X3=su3irrep23%lambda+su3irrep23%mu+1
NCE2A=X2*(X2+1)*(X2+2)/6
NCE2B=X1*(X1+1)*(X1+2)/6
NCE2D=X3*(X3+1)*(X3+2)/6
NABC=KR0A*KR0B*KR0C

CALL racah(su3irrep1,su3irrep2,su3irrepx,su3irrep3,su3irrep12,&
           su3irrep23,KR0A,KR0B,KR0C,KR0D,RHS,INFO,KIMAX2,&
           X1,X2,X3,NCE2A,NCE2B,NCE2D,NABC)
RETURN
CONTAINS
 FUNCTION IDM(su3irrepx) RESULT(IIDM)
  IMPLICIT NONE
  TYPE(su3irrep),INTENT (IN) :: su3irrepx
  INTEGER :: IIDM
  IIDM=(su3irrepx%lambda+1)*(su3irrepx%mu+1)*(su3irrepx%lambda+su3irrepx%mu+2)/2
 END FUNCTION IDM
END SUBROUTINE racah_wrap
