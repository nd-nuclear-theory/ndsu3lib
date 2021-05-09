!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nine_lambda_mu_wrapper.F90 -- wrapper for 9-(lambda,mu) coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nine_lambda_mu_wrapper(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,&
           mu13,lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,dimen,ninelm_block,info)
!------------------------------------------------------------------------------------------------------------------------------------
! Wrapper of the subroutine calculating 9-(lambda,mu) coefficients
!
! | (lambda1,mu1)   (lambda2,mu2)  (lambda12,mu12)  rho12 |
! | (lambda3,mu3)   (lambda4,mu4)  (lambda34,mu34)  rho34 |
! |(lambda13,mu13) (lambda24,mu24)   (lambda,mu)   rho1324|
! |     rho13           rho24          rho1234            |
!
! Input arguments: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,lambda24,mu24,lambda,mu,
!                  rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,dimen
! Output arguments: ninelm_block,info
!
! rhomax12 is the multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda12,mu12).
! rhomax34 is the multiplicity of coupling (lambda3,mu3)x(lambda4,mu4)->(lambda34,mu34).
! rhomax1234 is the multiplicity of coupling (lambda12,mu12)x(lambda34,mu34)->(lambda,mu).
! rhomax13 is the multiplicity of coupling (lambda1,mu1)x(lambda3,mu3)->(lambda13,mu13).
! rhomax24 is the multiplicity of coupling (lambda2,mu2)x(lambda4,mu4)->(lambda24,mu24).
! rhomax1324 is the multiplicity of coupling (lambda13,mu13)x(lambda24,mu24)->(lambda,mu).
! dimen is the size of the array ninelm_block. It must be at least rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*rhomax1324.
! ninelm_block(ind) is the 9-(lambda,mu)  coefficient for given rho12,rho34,rho1234,rho13,rho24,rho1324.
! ind = rho12+rhomax12*(rho34-1)+rhomax12*rhomax34*(rho1234-1)+rhomax12*rhomax34*rhomax1234*(rho13-1)
!       +rhomax12*rhomax34*rhomax1234*rhomax13*(rho24-1)+rhomax12*rhomax34*rhomax1234*rhomax13*rhomax24*(rho1324-1)
! info=0 if MKL subroutine dgesv in subroutines racah and Z_coeff ran without errors.
!------------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                      lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,dimen
INTEGER,INTENT(OUT) :: info
REAL(KIND=8),DIMENSION(dimen),INTENT(OUT) :: ninelm_block
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:,:,:) :: ninelm
INTEGER :: ind,rho12,rho34,rho1234,rho13,rho24,rho1324

INTERFACE
  SUBROUTINE calculate_9_lambda_mu(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                                lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,ninelm,info)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                          lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324
    INTEGER,INTENT(OUT) :: info
    REAL(KIND=8),DIMENSION(:,:,:,:,:,:),INTENT(OUT) :: ninelm
  END SUBROUTINE calculate_9_lambda_mu
END INTERFACE

ALLOCATE(ninelm(rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324))
CALL calculate_9_lambda_mu(lambda1,mu1,lambda2,mu2,lambda12,mu12,lambda3,mu3,lambda4,mu4,lambda34,mu34,lambda13,mu13,&
                        lambda24,mu24,lambda,mu,rhomax12,rhomax34,rhomax1234,rhomax13,rhomax24,rhomax1324,ninelm,info)
ind=0
DO rho1324=1,rhomax1324
  DO rho24=1,rhomax24
    DO rho13=1,rhomax13
      DO rho1234=1,rhomax1234
        DO rho34=1,rhomax34
          DO rho12=1,rhomax12
            ind=ind+1
            ninelm_block(ind)=ninelm(rho12,rho34,rho1234,rho13,rho24,rho1324)
          END DO
        END DO
      END DO
    END DO
  END DO
END DO
DEALLOCATE(ninelm)

END SUBROUTINE nine_lambda_mu_wrapper
