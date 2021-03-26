MODULE I_S_module
!--------------------------------------------------------------------------------------------------
! I(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose sigma-n)*(q choose n)
! S(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose n)/(p+q choose sigma+n)
!
! Indexing: I(p,q,sigma)=Ia((p^3+(p+q)^2+q)/2+sigma): p=0,...,upbound_I
!                                                     q=0,...,p
!                                                     sigma=0,...,p+q
!           S(p,q,sigma)=Sa((upbound_S*p*(p+upbound_S+2)+(p+q)*(p+q+1))/2+sigma): p=0,...,upbound_S
!                                                                                 q=0,...,ubbound_S
!                                                                                 sigma=0,...,p+q
!
! For p<q: I(p,q,sigma)=(-1)^{sigma}*I(q,p,sigma)
!--------------------------------------------------------------------------------------------------
USE mpmodule
IMPLICIT NONE
INTEGER :: upbound_I,upbound_S
#if defined(SU3DBL)
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Ia
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Sa
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Ia_quad
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Sa_quad
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Ia_mp
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Sa_mp
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Ia
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Sa
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Ia_quad
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Sa_quad
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Ia_mp
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Sa_mp
#elif (defined(SU3MP) || defined(SU3MP_GNU))
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Ia
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Sa
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Ia_quad
REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: Sa_quad
TYPE(mp_real),ALLOCATABLE,DIMENSION(:) :: Ia_mp
TYPE(mp_real),ALLOCATABLE,DIMENSION(:) :: Sa_mp
#endif
CONTAINS
  SUBROUTINE allocate_I(upboundI)
  !--------------------------------------------------------------------------------------
  ! This subroutine calculates values I(p,q,sigma) using recursion relation
  !
  ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
  !
  ! with starting condition
  !
  ! I(p,0,sigma)=(p choose sigma)
  !
  ! and stores them in array Ia. This subroutine should be called at the begining of the
  ! main program calculating SU(3)-SO(3) Wigner coefficients (after "USE binomial_coeff",
  ! "USE I_S_module" and "CALL allocate_binom").
  !--------------------------------------------------------------------------------------
  USE binomial_coeff
  USE mpmodule
  USE precision_level
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: upboundI
  INTEGER :: p,q,sigma,ind1,ind2,aux
  upbound_I=upboundI
  IF(upbound_I>upbound_binom)THEN
    CALL reallocate_binom(upbound_I-upbound_binom)
  END IF
  ind1=(upbound_I**3+(2*upbound_I)**2+upbound_I)/2+2*upbound_I
  ALLOCATE(Ia(0:ind1),Ia_quad(0:ind1),Ia_mp(0:ind1))
  aux=0
  ind2=0
  DO p=0,upbound_I
    aux=aux+p ! aux1=p*(p+1)/2
    ind1=p*aux ! ind1=(p^3+p^2)/2
    DO sigma=0,p
      Ia(ind1)=binom(ind2) ! I(p,0,sigma)=(p choose sigma)
      Ia_quad(ind1)=binom_quad(ind2)
      Ia_mp(ind1)=binom_mp(ind2)
      ind1=ind1+1 ! ind1=(p^3+p^2)/2+sigma
      ind2=ind2+1 ! ind2=p*(p+1)/2+sigma
    END DO
  END DO
  Ia(3)=Ia(0)
  Ia_quad(3)=Ia_quad(0)
  Ia_mp(3)=Ia_mp(0)
  Ia(4)=0.D0
#if defined(SU3DBL)
  Ia_quad(4)=0.D0
  Ia_mp(4)=0.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
  Ia_quad(4)=0.Q0
  Ia_mp(4)=0.Q0
#elif (defined(SU3MP) || defined(SU3MP_GNU))
  Ia_quad(4)=0.Q0
  Ia_mp(4)=mpreal(0.D0,nwds)
#endif
  Ia(5)=-Ia(0)
  Ia_quad(5)=-Ia_quad(0)
  Ia_mp(5)=-Ia_mp(0)
  ind1=5
  ind2=1
  DO p=2,upbound_I
    ind1=ind1+p+1
    DO q=1,p
      DO sigma=0,1
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=Ia(ind2) ! I(p,q,sigma)=I(p-1,q-1,sigma)
        Ia_quad(ind1)=Ia_quad(ind2)
        Ia_mp(ind1)=Ia_mp(ind2)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
      END DO
      DO sigma=2,p+q-2
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=Ia(ind2)-Ia(ind2-2) ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
        Ia_quad(ind1)=Ia_quad(ind2)-Ia_quad(ind2-2)
        Ia_mp(ind1)=Ia_mp(ind2)-Ia_mp(ind2-2)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
      END DO
      ind2=ind2-2
      DO sigma=p+q-1,p+q
        ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
        Ia(ind1)=-Ia(ind2) ! I(p,q,sigma)=-I(p-1,q-1,sigma-2)
        Ia_quad(ind1)=-Ia_quad(ind2)
        Ia_mp(ind1)=-Ia_mp(ind2)
        ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma-2
      END DO
    END DO
  END DO
  END SUBROUTINE allocate_I
  SUBROUTINE allocate_S(upboundS)
  !-------------------------------------------------------------------------------------
  ! This subroutine calculates values S(p,q,sigma) using recursion relations
  !
  ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1) unless q=upbound_S
  ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q) iff q=upbound_S
  !
  ! with starting condition
  !
  ! S(0,q,sigma)=1/(q choose sigma)
  !
  ! and stores them in array Sa. This subroutine should be called at the begining of the
  ! main program calculating SU(3)-SO(3) Wigner coefficients (after "USE binomial_coeff",
  ! "USE I_S_module" and "CALL allocate_binom").
  !--------------------------------------------------------------------------------------
  USE binomial_coeff
  USE mpmodule
  USE precision_level
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: upboundS
  INTEGER :: p,q,sigma,ind1,ind2,aux,ind3,n,ind4,ind5,nmax,aux2,aux3
  upbound_S=upboundS
  IF(2*upbound_S>upbound_binom)THEN
    CALL reallocate_binom(2*upbound_S-upbound_binom)
  END IF
  n=(upbound_S*upbound_S*(2*upbound_S+2)+2*upbound_S*(2*upbound_S+1))/2+2*upbound_S
  ALLOCATE(Sa(0:n),Sa_quad(0:n),Sa_mp(0:n))
  ind1=0
  DO q=0,upbound_S
    DO sigma=0,q
      Sa(ind1)=1.D0/binom(ind1)
#if defined(SU3DBL)
      Sa_quad(ind1)=1.D0/binom_quad(ind1)
      Sa_mp(ind1)=1.D0/binom_mp(ind1)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
      Sa_quad(ind1)=1.Q0/binom_quad(ind1)
      Sa_mp(ind1)=1.Q0/binom_mp(ind1)
#elif (defined(SU3MP) || defined(SU3MP_GNU))
      Sa_quad(ind1)=1.Q0/binom_quad(ind1)
      Sa_mp(ind1)=1.D0/binom_mp(ind1)
#endif
      ind1=ind1+1 ! ind1=q*(q+1)/2+sigma
    END DO
  END DO
  aux=0
  ind1=upbound_S*(upbound_S+3)/2+1
  aux2=1
  DO p=1,upbound_S
    aux2=aux2+p ! aux2=p*(p+1)/2+1
    aux3=(p+upbound_S)*(p+upbound_S+1)/2+1
    ind2=aux+p
    aux=ind1
    DO q=0,upbound_S-1
      DO sigma=0,p+q-1
        Sa(ind1)=Sa(ind2)-Sa(ind2+1) ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1)
        Sa_quad(ind1)=Sa_quad(ind2)-Sa_quad(ind2+1)
        Sa_mp(ind1)=Sa_mp(ind2)-Sa_mp(ind2+1)
        ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
        ind2=ind2+1 ! ind2=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q+1)*(p-1+q+1+1))/2+sigma
      END DO
      Sa(ind1)=1.D0 ! S(p,q,p+q)=1
#if defined(SU3DBL)
      Sa_quad(ind1)=1.D0
      Sa_mp(ind1)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
      Sa_quad(ind1)=1.Q0
      Sa_mp(ind1)=1.Q0
#elif (defined(SU3MP) || defined(SU3MP_GNU))
      Sa_quad(ind1)=1.Q0
      Sa_mp(ind1)=mpreal(1.D0,nwds)
#endif
      ind1=ind1+1
      ind2=ind2+1
    END DO
    ! q=upbound
    ind3=ind2-p-upbound_S
    ind2=ind1-p-upbound_S
    DO sigma=0,upbound_S-1
      Sa(ind1)=(DFLOAT(upbound_S-sigma)*Sa(ind2)+DFLOAT(p)*Sa(ind3))/DFLOAT(p+upbound_S)
#if defined(SU3DBL)
      Sa_quad(ind1)=(DFLOAT(upbound_S-sigma)*Sa_quad(ind2)+DFLOAT(p)*Sa_quad(ind3))/DFLOAT(p+upbound_S)
      Sa_mp(ind1)=(DFLOAT(upbound_S-sigma)*Sa_mp(ind2)+DFLOAT(p)*Sa_mp(ind3))/DFLOAT(p+upbound_S)
#elif defined(SU3QUAD)
      Sa_quad(ind1)=(QFLOAT(upbound_S-sigma)*Sa_quad(ind2)+QFLOAT(p)*Sa_quad(ind3))/QFLOAT(p+upbound_S)
      Sa_mp(ind1)=(QFLOAT(upbound_S-sigma)*Sa_mp(ind2)+QFLOAT(p)*Sa_mp(ind3))/QFLOAT(p+upbound_S)
#elif defined(SU3QUAD_GNU)
      Sa_quad(ind1)=(REAL(upbound_S-sigma,16)*Sa_quad(ind2)+REAL(p,16)*Sa_quad(ind3))/REAL(p+upbound_S,16)
      Sa_mp(ind1)=(REAL(upbound_S-sigma,16)*Sa_mp(ind2)+REAL(p,16)*Sa_mp(ind3))/REAL(p+upbound_S,16)
#elif defined(SU3MP)
      Sa_quad(ind1)=(QFLOAT(upbound_S-sigma)*Sa_quad(ind2)+QFLOAT(p)*Sa_quad(ind3))/QFLOAT(p+upbound_S)
      Sa_mp(ind1)=(DFLOAT(upbound_S-sigma)*Sa_mp(ind2)+DFLOAT(p)*Sa_mp(ind3))/DFLOAT(p+upbound_S)
#elif defined(SU3MP_GNU)
      Sa_quad(ind1)=(REAL(upbound_S-sigma,16)*Sa_quad(ind2)+REAL(p,16)*Sa_quad(ind3))/REAL(p+upbound_S,16)
      Sa_mp(ind1)=(DFLOAT(upbound_S-sigma)*Sa_mp(ind2)+DFLOAT(p)*Sa_mp(ind3))/DFLOAT(p+upbound_S)
#endif
      ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q)
      ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
      ind2=ind2+1 ! ind2=(upbound*p*(p+upbound+2)+(p+q-1)*(p+q-1+1))/2+sigma
      ind3=ind3+1 ! ind3=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q)*(p-1+q+1))/2+sigma
    END DO
    nmax=p-1
    DO sigma=upbound_S,p+upbound_S-1
      nmax=nmax-1
      ind4=aux2
      ind5=aux3+sigma
      Sa(ind1)=1.D0/binom(ind5-1)
#if defined(SU3DBL)
      Sa_quad(ind1)=1.D0/binom_quad(ind5-1)
      Sa_mp(ind1)=1.D0/binom_mp(ind5-1)
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
      Sa_quad(ind1)=1.Q0/binom_quad(ind5-1)
      Sa_mp(ind1)=1.Q0/binom_mp(ind5-1)
#elif (defined(SU3MP) || defined(SU3MP_GNU))
      Sa_quad(ind1)=1.Q0/binom_quad(ind5-1)
      Sa_mp(ind1)=1.D0/binom_mp(ind5-1)
#endif
      DO n=1,nmax,2
        Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
        Sa(ind1)=Sa(ind1)+binom(ind4+1)/binom(ind5+1)
        Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
        Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
        Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
        Sa_mp(ind1)=Sa_mp(ind1)+binom_mp(ind4+1)/binom_mp(ind5+1)
        ind4=ind4+2
        ind5=ind5+2
      END DO
      IF(nmax>=-1)THEN
        IF(BTEST(nmax,0))THEN ! nmax is odd
          Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
          Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
          Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
        ELSE ! nmax is even
          Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
          Sa(ind1)=Sa(ind1)+binom(ind4+1)/binom(ind5+1)
          Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
          Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
          Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
          Sa_mp(ind1)=Sa_mp(ind1)+binom_mp(ind4+1)/binom_mp(ind5+1)
        END IF
      END IF
      ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
    END DO
    Sa(ind1)=1.D0 ! S(p,q,p+q)=1
#if defined(SU3DBL)
    Sa_quad(ind1)=1.D0
    Sa_mp(ind1)=1.D0
#elif (defined(SU3QUAD) || defined(SU3QUAD_GNU))
    Sa_quad(ind1)=1.Q0
    Sa_mp(ind1)=1.Q0
#elif (defined(SU3MP) || defined(SU3MP_GNU))
    Sa_quad(ind1)=1.Q0
    Sa_mp(ind1)=mpreal(1.D0,nwds)
#endif
    ind1=ind1+1
  END DO
  END SUBROUTINE allocate_S
  SUBROUTINE deallocate_I_S
  DEALLOCATE(Ia,Sa,Ia_quad,Sa_quad,Ia_mp,Sa_mp)
  END SUBROUTINE deallocate_I_S
  SUBROUTINE reallocate_I(incr)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: incr
  DEALLOCATE(Ia,Ia_quad,Ia_mp)
  CALL allocate_I(upbound_I+incr)
  END SUBROUTINE reallocate_I
  SUBROUTINE reallocate_S(incr)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: incr
  DEALLOCATE(Sa,Sa_quad,Sa_mp)
  CALL allocate_S(upbound_S+incr)
  END SUBROUTINE reallocate_S
END MODULE I_S_module
