SUBROUTINE nine_lambda_mu(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                          lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,ninelm,info)
!------------------------------------------------------------------------------------------------------------------------------------
! Calculates 9-(lambda,mu) coefficients
!
! | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
! | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
! |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
! |     rho13           rho24          rho1234            |
!
! for given lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,lambda24,mu24,lambda,mu
! using Eq.(3) in the reference.
!
! Reference: D.J.Millener, J.Math.Phys., Vol.19, No.7 (1978) 1513
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,lambda24,mu24,lambda,mu,
!                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
! Output arguments: ninelm,info
!
! rhomax12 is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12)
! rhomax34 is the multiplicity of coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34)
! rhomax1234 is the multiplicity of coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu)
! rhomax13 is the multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13)
! rhomax24 is the multiplicity of coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24)
! rhomax1324 is the multiplicity of coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu)
! ninelm(rho12,rho34,rho1234,rho13,rho24,rho1324) is the 9-(lambda,mu)  coefficient for given rho12,rho34,rho1234,rho13,rho24,rho1324
! info=0 if MKL subroutine dgesv in subroutines racah and Z_coeff ran without errors.
!------------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,EXTERNAL :: outer_multiplicity
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                      lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
INTEGER,INTENT(OUT) :: info
REAL(KIND=8),DIMENSION(:,:,:,:,:,:),INTENT(OUT) :: ninelm
INTEGER :: lambda0,mu0,rho132,rho04,rho123,rhomax132,rhomax04,rhomax123,rho12,rho34,rho13,rho24,rho1324,nU1,nZ,nU2,rhomax12304,i1,i2
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: U1,Z,U2

INTERFACE
  SUBROUTINE racah(lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,&
                   rhomaxa,rhomaxb,rhomaxc,rhomaxd,rac,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda23,mu23,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: rac
  END SUBROUTINE racah
  SUBROUTINE Z_coeff(lambda2,mu2,lambda1,mu1,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                     rhomaxa,rhomaxb,rhomaxc,rhomaxd,Zcoeff,ldb,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda,mu,lambda3,mu3,lambda12,mu12,lambda13,mu13,rhomaxa,rhomaxb,rhomaxc,rhomaxd,ldb
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:),INTENT(OUT) :: Zcoeff
  END SUBROUTINE Z_coeff
END INTERFACE

ninelm=0.D0

! (lambda13,mu13)x(lambda2,mu2)->(lambda0,mu0)
! (lambda12,mu12)x(lambda3,mu3)->(lambda0,mu0)

DO lambda0=0,MIN(lambda13+lambda2+MIN(mu2,lambda13+mu13),lambda12+lambda3+MIN(mu3,lambda12+mu12))
  DO mu0=0,MIN(mu13+mu2+MIN(lambda13,lambda2),mu12+mu3+MIN(lambda12,lambda3))
    rhomax132=outer_multiplicity(lambda13,mu13,lambda2,mu2,lambda0,mu0)
    IF(rhomax132==0)CYCLE
    rhomax04=outer_multiplicity(lambda0,mu0,lambda4,mu4,lambda,mu)
    IF(rhomax04==0)CYCLE
    rhomax123=outer_multiplicity(lambda12,mu12,lambda3,mu3,lambda0,mu0)
    IF(rhomax123==0)CYCLE
    rhomax12304=rhomax123*rhomax04

    ALLOCATE(U1(rhomax1324,rhomax132*rhomax04*rhomax24),Z(rhomax132,rhomax12*rhomax123*rhomax13),&
             U2(rhomax1234,rhomax123*rhomax04*rhomax34))

    CALL racah(lambda13,mu13,lambda2,mu2,lambda,mu,lambda4,mu4,lambda0,mu0,lambda24,mu24,&
               rhomax132,rhomax04,rhomax24,rhomax1324,U1,rhomax1324,info)
    IF(info/=0)THEN
      DEALLOCATE(U1,Z,U2)
      RETURN
    END IF
    CALL Z_coeff(lambda2,mu2,lambda1,mu1,lambda0,mu0,lambda3,mu3,lambda12,mu12,lambda13,mu13,&
                 rhomax12,rhomax123,rhomax13,rhomax132,Z,rhomax132,info)
    IF(info/=0)THEN
      DEALLOCATE(U1,Z,U2)
      RETURN
    END IF
    CALL racah(lambda12,mu12,lambda3,mu3,lambda,mu,lambda4,mu4,lambda0,mu0,lambda34,mu34,&
               rhomax123,rhomax04,rhomax34,rhomax1234,U2,rhomax1234,info)
    IF(info/=0)THEN
      DEALLOCATE(U1,Z,U2)
      RETURN
    END IF

!-------------------------------------------------------------
!      n=0
!      DO rhoc=1,rhomaxc
!        DO rhob=1,rhomaxb
!          DO rhoa=1,rhomaxa
!            n=n+1 ! n=rhoa+rhomaxa*(rhob-1)+rhomaxa*rhomaxb*(rhoc-1)
!-------------------------------------------------------------

    nU1=0
    DO rho24=1,rhomax24
      i1=-rhomax123
      DO rho04=1,rhomax04
        i1=i1+rhomax123 ! i1=rhomax123*(rho04-1)
        DO rho132=1,rhomax132
          nU1=nU1+1
! nU1=rho132+rhomax132*(rho04-1)+rhomax132*rhomax04*(rho24-1)
          nZ=0
          DO rho13=1,rhomax13
            i2=i1
            DO rho123=1,rhomax123
              i2=i2+1 ! i2=rho123+i1
              DO rho12=1,rhomax12
                nZ=nZ+1
! nZ=rho12+rhomax12*(rho123-1)+rhomax12*rhomax123*(rho13-1)
                nU2=i2-rhomax12304
                DO rho34=1,rhomax34
                  nU2=nU2+rhomax12304
! nU2=rho123+rhomax123*(rho04-1)+rhomax123*rhomax04*(rho34-1)
                  DO rho1324=1,rhomax1324
                    ninelm(rho12,rho34,1:rhomax1234,rho13,rho24,rho1324)=ninelm(rho12,rho34,1:rhomax1234,rho13,rho24,rho1324)&
                      +U1(rho1324,nU1)*Z(rho132,nZ)*U2(1:rhomax1234,nU2)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(U1,Z,U2)

  END DO
END DO

END SUBROUTINE nine_lambda_mu
