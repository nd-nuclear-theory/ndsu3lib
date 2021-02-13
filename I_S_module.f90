MODULE I_S_module
!----------------------------------------------------------------------------------------
! I(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose sigma-n)*(q choose n)
! S(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose n)/(p+q choose sigma+n)
!
! Indexing: I(p,q,sigma)=Ia((p^3+(p+q)^2+q)/2+sigma): p=0,...,upbound
!                                                     q=0,...,p
!                                                     sigma=0,...,p+q
!           S(p,q,sigma)=Sa(upbound*p*(p+upbound+2)+(p+q)*(p+q+1)+sigma): p=0,...,upbound
!                                                                         q=0,...,ubbound
!                                                                         sigma=0,...,p+q
!
! For p<q: I(p,q,sigma)=(-1)^{sigma}*I(q,p,sigma)
!----------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: upbound
REAL(KIND=8),DIMENSION(0:4080500) :: Ia
REAL(KIND=8),DIMENSION(0:8120600) :: Sa
CONTAINS
  SUBROUTINE calculate_I_S
  !------------------------------------------------------------------------------------------
  ! This subroutine calculates values I(p,q,sigma) and S(p,q,sigma) using recursion relations
  !
  ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
  ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1) unless q=upbound
  ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q) iff q=upbound
  !
  ! with starting conditions
  !
  ! I(p,0,sigma)=(p choose sigma)
  ! S(0,q,sigma)=1/(q choose sigma)
  !
  ! and stores them in arrays Ia,Sa. This subroutine should be called at the begining of the
  ! main program calculating SU(3)-SO(3) wigner coefficients (after "USE binomial_coeff",
  ! "USE I_S_module" and "CALL calculate_binomial_coeff").
  !------------------------------------------------------------------------------------------
  USE binomial_coeff
  IMPLICIT NONE
  INTEGER :: p,q,sigma,ind1,ind2,aux,ind3
  upbound=200
  !*****************************
  ! Calculation of I(p,q,sigma)
  !*****************************
  aux=0
  ind2=0
  DO p=0,upbound
    aux=aux+p ! aux1=p*(p+1)/2
    ind1=p*aux ! ind1=(p^3+p^2)/2
    DO sigma=0,p
      Ia(ind1)=binom(ind2) ! I(p,0,sigma)=(p choose sigma)
      ind1=ind1+1 ! ind1=(p^3+p^2)/2+sigma
      ind2=ind2+1 ! ind2=p*(p+1)/2+sigma
    END DO
  END DO
  Ia(3)=Ia(0)
  Ia(4)=0.D0
  Ia(5)=-Ia(0)
  ind1=5
  ind2=1
  DO p=2,upbound
    ind1=ind1+p+1
    DO q=1,p
      DO sigma=0,1
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=Ia(ind2) ! I(p,q,sigma)=I(p-1,q-1,sigma)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
      END DO
      DO sigma=2,p+q-2
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=Ia(ind2)-Ia(ind2-2) ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
      END DO
      ind2=ind2-2
      DO sigma=p+q-1,p+q
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=-Ia(ind2) ! I(p,q,sigma)=-I(p-1,q-1,sigma-2)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma-2
      END DO
    END DO
  END DO
  !*****************************
  ! Calculation of S(p,q,sigma)
  !*****************************
  ind1=0
  DO q=0,upbound
    DO sigma=0,q
      Sa(ind1)=1.D0/binom(ind1)
      ind1=ind1+1 ! ind1=q*(q+1)/2+sigma
    END DO
  END DO
  aux=0
  ind1=upbound*(upbound+3)/2+1
  DO p=1,upbound
    ind2=aux+p
    aux=ind1
    DO q=0,upbound-1
      DO sigma=0,p+q-1
        Sa(ind1)=Sa(ind2)-Sa(ind2+1) ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1)
        ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
        ind2=ind2+1 ! ind2=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q+1)*(p-1+q+1+1))/2+sigma
      END DO
      Sa(ind1)=1.D0 ! S(p,q,p+q)=1
      ind1=ind1+1
      ind2=ind2+1
    END DO
    ! q=upbound
    ind3=ind2-p-upbound
    ind2=ind1-p-upbound
    DO sigma=0,p+upbound-1
      Sa(ind1)=(DFLOAT(upbound-sigma)*Sa(ind2)+DFLOAT(p)*Sa(ind3))/DFLOAT(p+upbound)
      ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q)
      ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
      ind2=ind2+1 ! ind2=(upbound*p*(p+upbound+2)+(p+q-1)*(p+q-1+1))/2+sigma
      ind3=ind3+1 ! ind3=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q)*(p-1+q+1))/2+sigma
    END DO
    Sa(ind1)=1.D0 ! S(p,q,p+q)=1
    ind1=ind1+1
  END DO
  END SUBROUTINE calculate_I_S
END MODULE I_S_module
