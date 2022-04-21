!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ndsu3lib_wigner_canonical.F90 -- module for SU(3)-SU(2)xU(1) reduced Wigner coefficients
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE ndsu3lib_wigner_canonical
  USE,INTRINSIC :: iso_c_binding
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  USE mpmodule
#endif
  IMPLICIT NONE
  TYPE,BIND(C) :: su3irrep
     INTEGER(C_INT) :: lambda,mu
  END TYPE su3irrep
#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  INTEGER :: ndig,nwds
  PARAMETER(ndig=37,nwds=INT(ndig/mpdpw+2)) ! 37 is the minimal ndig
#endif
  !----------------------------------------------------------------------------------------------------------------------------
  ! Binomial coefficients: (n choose k)=binom(n*(n+1)/2+k)=binom_quad(n*(n+1)/2+k)=binom_mp(n*(n+1)/2+k): n=0,...,upbound_binom
  !                                                                                                       k=0,...,n
  ! I(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose sigma-n)*(q choose n)=Ia((p^3+(p+q)^2+q)/2+sigma): p=0,...,upbound_I
  !                                                                                            q=0,...,p
  !                                                                                            sigma=0,...,p+q
  ! S(p,q,sigma)=\sum_{n}(-1)^{n}*(p choose n)/(p+q choose sigma+n)
  !             =Sa((upbound_S*p*(p+upbound_S+2)+(p+q)*(p+q+1))/2+sigma): p=0,...,upbound_S
  !                                                                       q=0,...,ubbound_S
  !                                                                       sigma=0,...,p+q
  ! For p<q: I(p,q,sigma)=(-1)^{sigma}*I(q,p,sigma)
  !-----------------------------------------------------------------------------------------------------------------------------
  INTEGER :: upbound_binom,upbound_I,upbound_S
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: binom,Ia,Sa
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
  REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: binom_quad,Ia_quad,Sa_quad
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
  REAL(KIND=16),ALLOCATABLE,DIMENSION(:) :: binom_quad,Ia_quad,Sa_quad
  TYPE(mp_real),ALLOCATABLE,DIMENSION(:) :: binom_mp,Ia_mp,Sa_mp
#endif

CONTAINS

  SUBROUTINE allocate_binom(upbound)
    !-----------------------------------------------------------------------------------------------------------------
    ! Allocates arrays binom, binom_quad and binom_mp, calculates binomial coefficients and stores them in the arrays.
    !-----------------------------------------------------------------------------------------------------------------
    !#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    !USE mpmodule
    !#endif
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: upbound
    INTEGER :: n,k,ind
    upbound_binom=upbound
    ind=upbound_binom*(upbound_binom+1)/2+upbound_binom
#if defined(NDSU3LIB_DBL)
    ALLOCATE(binom(0:ind))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    ALLOCATE(binom(0:ind),binom_quad(0:ind))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    ALLOCATE(binom(0:ind),binom_quad(0:ind),binom_mp(0:ind))
#endif
    binom(0)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    binom_quad(0)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    binom_quad(0)=1.Q0
    binom_mp(0)=mpreal(1.D0,nwds)
#endif
    ind=1
    DO n=1,upbound_binom
       binom(ind)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
       binom_quad(ind)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
       binom_quad(ind)=1.Q0
       binom_mp(ind)=mpreal(1.D0,nwds)
#endif
       ind=ind+1
       DO k=1,n-1
          binom(ind)=binom(ind-n-1)+binom(ind-n)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
          binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
          binom_quad(ind)=binom_quad(ind-n-1)+binom_quad(ind-n)
          binom_mp(ind)=binom_mp(ind-n-1)+binom_mp(ind-n)
#endif
          ind=ind+1
       END DO
       binom(ind)=1.D0
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
       binom_quad(ind)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
       binom_quad(ind)=1.Q0
       binom_mp(ind)=mpreal(1.D0,nwds)
#endif
       ind=ind+1
    END DO
  END SUBROUTINE allocate_binom

  SUBROUTINE deallocate_binom
    !---------------------------------------------------
    ! Deallocates arrays binom, binom_quad and binom_mp.
    !---------------------------------------------------
    IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
    DEALLOCATE(binom)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    DEALLOCATE(binom,binom_quad)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    DEALLOCATE(binom,binom_quad,binom_mp)
#endif
  END SUBROUTINE deallocate_binom

  SUBROUTINE reallocate_binom(incr)
    !------------------------------------------------------------------------
    ! Deallocates arrays binom, binom_quad and binom_mp, allocates them again
    ! with upbound_binom increased by incr and calculates the entries.
    !------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: incr
    CALL deallocate_binom
    CALL allocate_binom(upbound_binom+incr)
  END SUBROUTINE reallocate_binom

  SUBROUTINE allocate_I(upboundI)
    !--------------------------------------------------------
    ! Calculates values I(p,q,sigma) using recursion relation
    !
    ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
    !
    ! with starting condition
    !
    ! I(p,0,sigma)=(p choose sigma)
    !
    ! and stores them in arrays Ia, Ia_quad and Ia_mp.
    !---------------------------------------------------------
    !#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    !USE mpmodule
    !#endif
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: upboundI
    INTEGER :: p,q,sigma,ind1,ind2,aux
    upbound_I=upboundI
    IF(upbound_I>upbound_binom)THEN
       CALL reallocate_binom(upbound_I-upbound_binom)
    END IF
    ind1=(upbound_I**3+(2*upbound_I)**2+upbound_I)/2+2*upbound_I
#if defined(NDSU3LIB_DBL)
    ALLOCATE(Ia(0:ind1))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    ALLOCATE(Ia(0:ind1),Ia_quad(0:ind1))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    ALLOCATE(Ia(0:ind1),Ia_quad(0:ind1),Ia_mp(0:ind1))
#endif
    aux=0
    ind2=0
    DO p=0,upbound_I
       aux=aux+p ! aux1=p*(p+1)/2
       ind1=p*aux ! ind1=(p^3+p^2)/2
       DO sigma=0,p
          Ia(ind1)=binom(ind2) ! I(p,0,sigma)=(p choose sigma)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
          Ia_quad(ind1)=binom_quad(ind2)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
          Ia_quad(ind1)=binom_quad(ind2)
          Ia_mp(ind1)=binom_mp(ind2)
#endif
          ind1=ind1+1 ! ind1=(p^3+p^2)/2+sigma
          ind2=ind2+1 ! ind2=p*(p+1)/2+sigma
       END DO
    END DO
    Ia(3)=Ia(0)
    Ia(4)=0.D0
    Ia(5)=-Ia(0)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    Ia_quad(3)=Ia_quad(0)
    Ia_quad(4)=0.Q0
    Ia_quad(5)=-Ia_quad(0)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    Ia_quad(3)=Ia_quad(0)
    Ia_quad(4)=0.Q0
    Ia_quad(5)=-Ia_quad(0)
    Ia_mp(3)=Ia_mp(0)
    Ia_mp(4)=mpreal(0.D0,nwds)
    Ia_mp(5)=-Ia_mp(0)
#endif
    ind1=5
    ind2=1
    DO p=2,upbound_I
       ind1=ind1+p+1
       DO q=1,p
          DO sigma=0,1
             ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
             Ia(ind1)=Ia(ind2) ! I(p,q,sigma)=I(p-1,q-1,sigma)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
             Ia_quad(ind1)=Ia_quad(ind2)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
             Ia_quad(ind1)=Ia_quad(ind2)
             Ia_mp(ind1)=Ia_mp(ind2)
#endif
             ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
          END DO
          DO sigma=2,p+q-2
             ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
             Ia(ind1)=Ia(ind2)-Ia(ind2-2) ! I(p,q,sigma)=I(p-1,q-1,sigma)-I(p-1,q-1,sigma-2)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
             Ia_quad(ind1)=Ia_quad(ind2)-Ia_quad(ind2-2)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
             Ia_quad(ind1)=Ia_quad(ind2)-Ia_quad(ind2-2)
             Ia_mp(ind1)=Ia_mp(ind2)-Ia_mp(ind2-2)
#endif
             ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma
          END DO
          ind2=ind2-2
          DO sigma=p+q-1,p+q
             ind1=ind1+1 ! ind1=(p^3+(p+q)^2+q)/2+sigma
             Ia(ind1)=-Ia(ind2) ! I(p,q,sigma)=-I(p-1,q-1,sigma-2)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
             Ia_quad(ind1)=-Ia_quad(ind2)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
             Ia_quad(ind1)=-Ia_quad(ind2)
             Ia_mp(ind1)=-Ia_mp(ind2)
#endif
             ind2=ind2+1 ! ind2=((p-1)^3+(p-1+q-1)^2+q-1)/2+sigma-2
          END DO
       END DO
    END DO
  END SUBROUTINE allocate_I

  SUBROUTINE allocate_S(upboundS)
    !--------------------------------------------------------------------------------
    ! Calculates values S(p,q,sigma) using recursion relations
    !
    ! S(p,q,sigma)=S(p-1,q+1,sigma)-S(p-1,q+1,sigma+1) unless q=upbound_S
    ! S(p,q,sigma)=((q-sigma)*S(p,q-1,sigma)+p*S(p-1,q,sigma))/(p+q) iff q=upbound_S
    !
    ! with starting condition
    !
    ! S(0,q,sigma)=1/(q choose sigma)
    !
    ! and stores them in array Sa, Sa_quad and Sa_mp.
    !--------------------------------------------------------------------------------
    !#if (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    !USE mpmodule
    !#endif  
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: upboundS
    INTEGER :: p,q,sigma,ind1,ind2,aux,ind3,n,ind4,ind5,nmax,aux2,aux3
    upbound_S=upboundS
    IF(2*upbound_S>upbound_binom)THEN
       CALL reallocate_binom(2*upbound_S-upbound_binom)
    END IF
    n=(upbound_S*upbound_S*(2*upbound_S+2)+2*upbound_S*(2*upbound_S+1))/2+2*upbound_S
#if defined(NDSU3LIB_DBL)
    ALLOCATE(Sa(0:n))
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    ALLOCATE(Sa(0:n),Sa_quad(0:n))
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    ALLOCATE(Sa(0:n),Sa_quad(0:n),Sa_mp(0:n))
#endif
    ind1=0
    DO q=0,upbound_S
       DO sigma=0,q
          Sa(ind1)=1.D0/binom(ind1)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
          Sa_quad(ind1)=1.Q0/binom_quad(ind1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
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
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
             Sa_quad(ind1)=Sa_quad(ind2)-Sa_quad(ind2+1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
             Sa_quad(ind1)=Sa_quad(ind2)-Sa_quad(ind2+1)
             Sa_mp(ind1)=Sa_mp(ind2)-Sa_mp(ind2+1)
#endif
             ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
             ind2=ind2+1 ! ind2=(upbound*(p-1)*(p-1+upbound+2)+(p-1+q+1)*(p-1+q+1+1))/2+sigma
          END DO
          Sa(ind1)=1.D0 ! S(p,q,p+q)=1
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
          Sa_quad(ind1)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
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
#if defined(NDSU3LIB_QUAD)
          Sa_quad(ind1)=(QFLOAT(upbound_S-sigma)*Sa_quad(ind2)+QFLOAT(p)*Sa_quad(ind3))/QFLOAT(p+upbound_S)
#elif defined(NDSU3LIB_QUAD_GNU)
          Sa_quad(ind1)=(REAL(upbound_S-sigma,16)*Sa_quad(ind2)+REAL(p,16)*Sa_quad(ind3))/REAL(p+upbound_S,16)
#elif defined(NDSU3LIB_MP)
          Sa_quad(ind1)=(QFLOAT(upbound_S-sigma)*Sa_quad(ind2)+QFLOAT(p)*Sa_quad(ind3))/QFLOAT(p+upbound_S)
          Sa_mp(ind1)=(DFLOAT(upbound_S-sigma)*Sa_mp(ind2)+DFLOAT(p)*Sa_mp(ind3))/DFLOAT(p+upbound_S)
#elif defined(NDSU3LIB_MP_GNU)
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
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
          Sa_quad(ind1)=1.Q0/binom_quad(ind5-1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
          Sa_quad(ind1)=1.Q0/binom_quad(ind5-1)
          Sa_mp(ind1)=1.D0/binom_mp(ind5-1)
#endif
          DO n=1,nmax,2
             Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
             Sa(ind1)=Sa(ind1)+binom(ind4+1)/binom(ind5+1)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
             Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
             Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
             Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
             Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
             Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
             Sa_mp(ind1)=Sa_mp(ind1)+binom_mp(ind4+1)/binom_mp(ind5+1)
#endif
             ind4=ind4+2
             ind5=ind5+2
          END DO
          IF(nmax>=-1)THEN
             IF(BTEST(nmax,0))THEN ! nmax is odd
                Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
                Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
                Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
                Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
#endif
             ELSE ! nmax is even
                Sa(ind1)=Sa(ind1)-binom(ind4)/binom(ind5)
                Sa(ind1)=Sa(ind1)+binom(ind4+1)/binom(ind5+1)
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
                Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
                Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
                Sa_quad(ind1)=Sa_quad(ind1)-binom_quad(ind4)/binom_quad(ind5)
                Sa_quad(ind1)=Sa_quad(ind1)+binom_quad(ind4+1)/binom_quad(ind5+1)
                Sa_mp(ind1)=Sa_mp(ind1)-binom_mp(ind4)/binom_mp(ind5)
                Sa_mp(ind1)=Sa_mp(ind1)+binom_mp(ind4+1)/binom_mp(ind5+1)
#endif
             END IF
          END IF
          ind1=ind1+1 ! ind1=(upbound*p*(p+upbound+2)+(p+q)*(p+q+1))/2+sigma
       END DO
       Sa(ind1)=1.D0 ! S(p,q,p+q)=1
#if (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
       Sa_quad(ind1)=1.Q0
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
       Sa_quad(ind1)=1.Q0
       Sa_mp(ind1)=mpreal(1.D0,nwds)
#endif
       ind1=ind1+1
    END DO
  END SUBROUTINE allocate_S

  SUBROUTINE deallocate_I
    !-----------------------------------------
    ! Deallocates arrays Ia, Ia_quad and Ia_mp
    !-----------------------------------------
    IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
    DEALLOCATE(Ia)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    DEALLOCATE(Ia,Ia_quad)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    DEALLOCATE(Ia,Ia_quad,Ia_mp)
#endif
  END SUBROUTINE deallocate_I

  SUBROUTINE deallocate_S
    !----------------------------------------
    ! Deallocates array Sa, Sa_quad and Sa_mp
    !----------------------------------------
    IMPLICIT NONE
#if defined(NDSU3LIB_DBL)
    DEALLOCATE(Sa)
#elif (defined(NDSU3LIB_QUAD) || defined(NDSU3LIB_QUAD_GNU))
    DEALLOCATE(Sa,Sa_quad)
#elif (defined(NDSU3LIB_MP) || defined(NDSU3LIB_MP_GNU))
    DEALLOCATE(Sa,Sa_quad,Sa_mp)
#endif
  END SUBROUTINE deallocate_S

  SUBROUTINE reallocate_I(incr)
    !---------------------------------------------------------------------------------------------------------------
    ! Deallocates arrays Ia, Ia_quad and Ia_mp, reallocates them with size increased by incr and calculates entries.
    !---------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: incr
    CALL deallocate_I
    CALL allocate_I(upbound_I+incr)
  END SUBROUTINE reallocate_I

  SUBROUTINE reallocate_S(incr)
    !---------------------------------------------------------------------------------------------------------------
    ! Deallocates arrays Sa, Sa_quad and Sa_mp, reallocates them with size increased by incr and calculates entries.
    !---------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: incr
    CALL deallocate_S
    CALL allocate_S(upbound_S+incr)
  END SUBROUTINE reallocate_S

  SUBROUTINE ndsu3lib_init(wso3,j2max) BIND(C)
    !----------------------------------------------------------------------------
    ! ndsu3lib initialization subroutine
    ! This subroutine must be called by the main program before calling ndsu3lib
    ! subroutines for SU(3) Wigner or recoupling coefficients.
    !
    ! Input arguments: wso3,j2max
    !
    ! wso3 must be .TRUE. if SU(3)-SO(3) Wigner coefficients are going to be
    ! calculated.
    ! If WIGXJPF is not going to be utilized, j2max is not used. Otherwise j2max
    ! must be greater than or equal to two times the maximal angular momentum
    ! expected in ordinary Clebsch-Gordan or SU(2) recoupling coefficients.
    ! j2max should be at least the maximal expected value of lambda+mu if
    ! SU(3)-SO(3) Wigner coefficients are not going to be calculated. If
    ! SU(3)-SO(3) Wigner coefficients are going to be calculated, j2max should be
    ! at least two times the maximal expected value of lambda+mu. If this j2max
    ! is insufficient, WIGXJPF will terminate the program and display an error
    ! message.
    !----------------------------------------------------------------------------
    !USE iso_c_binding
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
    USE fwigxjpf
#endif
    IMPLICIT NONE
    LOGICAL,INTENT(IN) :: wso3
    INTEGER(C_INT),INTENT(IN) :: j2max
    CALL allocate_binom(100)
    IF(wso3)THEN
       CALL allocate_I(50)
       CALL allocate_S(50)
    END IF
#if defined(NDSU3LIB_RACAH_WIGXJPF)
    CALL fwig_table_init(j2max,6)
    CALL fwig_temp_init(j2max)
#elif defined(NDSU3LIB_WSO3_WIGXJPF)
    CALL fwig_table_init(j2max,3)
    CALL fwig_temp_init(j2max)
#endif
  END SUBROUTINE ndsu3lib_init

  SUBROUTINE ndsu3lib_free(wso3) BIND(C)
    !-------------------------------------------------------------------------------
    ! This subroutine can be called by the main program once SU(3) Wigner or
    ! recoupling coefficients are not going to be calculated anymore to free memory.
    !
    ! Input argument: wso3
    !
    ! wso3 should be .TRUE. if ndsu3lib_init was called with the first argument
    ! being .TRUE.
    !-------------------------------------------------------------------------------
    !USE iso_c_binding
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
    USE fwigxjpf
#endif
    IMPLICIT NONE
    LOGICAL,INTENT(IN) :: wso3
    IF(wso3)THEN
       CALL deallocate_I
       CALL deallocate_S
    END IF
    CALL deallocate_binom
#if (defined(NDSU3LIB_RACAH_WIGXJPF) || defined(NDSU3LIB_WSO3_WIGXJPF))
    CALL fwig_temp_free();
    CALL fwig_table_free();
#endif
  END SUBROUTINE ndsu3lib_free

  FUNCTION outer_multiplicity(irrep1,irrep2,irrep3) RESULT(rhomax) BIND(C)
    !--------------------------------------------------------------------------------------
    ! Calculates multiplicity of SU(3) coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3),
    ! where lambda1=irrep1%lambda, mu1=irrep1%mu,
    !       lambda2=irrep2%lambda, mu2=irrep2%mu,
    !       lambda3=irrep3%lambda, mu3=irrep3%mu
    ! Reference: M.F.O'Reilly, J.Math.Phys. 23 (1982) 2022: Section 5, Proposition 7(a)
    !--------------------------------------------------------------------------------------
    !USE iso_c_binding
    IMPLICIT NONE
    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
    INTEGER :: rhomax,C3,C,D
    C3=irrep1%lambda-irrep1%mu+irrep2%lambda-irrep2%mu-irrep3%lambda+irrep3%mu ! C3 equals 3 times the C in the reference
    C=C3/3 ! C equals the C in the reference if 3 divides C3
    IF((3*C==C3).AND.(C>=-MIN(irrep2%mu,irrep1%mu)).AND.(C<=MIN(irrep1%lambda,irrep2%lambda)))THEN
       D=C+irrep1%mu+irrep2%mu-irrep3%mu ! D equals the D in the reference
       IF((D>=0).AND.(D<=MIN(irrep1%mu+irrep2%mu,irrep2%lambda+irrep2%mu,irrep1%lambda+irrep1%mu))&
            .AND.(D+C>=0).AND.(D+C<=MIN(irrep1%lambda+irrep1%mu,irrep2%lambda+irrep1%lambda,irrep2%lambda+irrep2%mu)))THEN
          rhomax=1+MIN(irrep2%mu,irrep1%lambda+irrep1%mu,D,irrep1%lambda-C)&
               -MAX(0,D-irrep1%mu,D-irrep2%lambda,-C,D-C-irrep1%mu,D+C-irrep2%lambda)
       ELSE
          rhomax=0
       END IF
    ELSE
       rhomax=0
    END IF
    RETURN
  END FUNCTION outer_multiplicity

  SUBROUTINE wigner_canonical_extremal(irrep1,irrep2,irrep3,I3,rhomax,i2,wigner,p1a,p2a,q2a)
    !------------------------------------------------------------------------------------------------------------
    ! Calculates the extremal reduced SU(3)-SU(2)xU(1) Wigner coefficients
    ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho
    ! for given lambda1,mu1,lambda2,mu2,lambda3,mu3 using Eq.(20),(21),(17),(18),(8),(1),(32),(35,2B) in Ref.[1].
    ! The phase convention of Ref.[2] is adopted.
    !
    ! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
    !             [2] K.T.Hecht, Nucl.Phys. 62 (1965) 1
    !             [3] D.Goldberg, ACM Computing Surveys, Vol.23, No.1 (1991) 5
    !
    ! Input arguments: irrep1x,irrep2x,irrep3x,I3,rhomax
    ! Output arguments: i2,wigner,p1a,p2a,q2a
    !
    ! lambda1=irrep1%lambda, mu1=irrep1%mu
    ! lambda2=irrep2%lambda, mu2=irrep2%mu
    ! lambda3=irrep3%lambda, mu3=irrep3%mu
    ! I3=1 for E=HW, I3=0 for E=LW
    ! rhomax=multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
    ! 
    ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)E>_rho=wigner(p1,p2,q2,rho)
    !   where epsilon2=2*lambda2+mu2-3*(p2+q2) for I3=1, epsilon2=-2*mu2-lambda2+3*(p2+q2) for I3=0
    !         epsilon1=epsilon3E-epsilon2
    !         Lambda1=(mu1+p1-q1)/2 for I3=1, Lambda1=(lambda1+p1-q1)/2 for I3=0
    !         Lambda2=(mu2+p2-q2)/2 for I3=1, Lambda2=(lambda2+p2-q2)/2 for I3=0
    !         p1=p1a(i)
    !         p2=p2a(i)
    !         q2=q2a(i)
    !         q1=(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3-p1-p2-q2 for I3=1
    !         q1=(2*(mu1+mu2+lambda3)+lambda1+lambda2+mu3)/3-p1-p2-q2 for I3=0
    !         1<=i<=i2
    !
    ! Note: There are typos in Eq.(18),(20) in [1]. There are 4 expressions for X in Eq.(18).
    !       The 3rd one is for X(Lambda1+1/2,Lambda2-1/2), not X(Lambda1-1/2,Lambda2-1/2).
    !       In Eq.(20), there should be \bar{lambda2} instead of lambda 2 in a,b,c,d.
    !------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
    INTEGER,INTENT(IN) :: I3,rhomax
    INTEGER,INTENT(OUT) :: i2
    TYPE(su3irrep) :: irrep1a,irrep2a,irrep3a
    INTEGER :: lambda1,mu1,lambda2,mu2,lambda3,mu3,p2,q2,Lambda22,i1,j2,j1,p2tilde,q2tilde,i,j,n,a,b,c,d,&
         epsilon2,steps21,Lambda12,eta,p1,q1,rho,i4,Lambda22max,ABCD,phiprhomax,p2min,p2max,&
         noname1,p1min,noname2,p1max
    INTEGER(KIND=8) :: Sq2,Rp2
    REAL(KIND=8) :: F,G,H,scalprod
    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a ! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
    ! Sizes of wigner are at least (0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax).

    IF(I3==1)THEN ! E=HW
       lambda1=irrep1%lambda
       lambda2=irrep2%lambda
       lambda3=irrep3%lambda
       mu1=irrep1%mu
       mu2=irrep2%mu
       mu3=irrep3%mu
       irrep1a=irrep1
       irrep3a=irrep3
    ELSE ! E=LW.
       ! Coefficients with E=LW are obtained from those with E=HW and flipped lambdas and mus.
       ! See Eq.(35,2B),(32) and the text below Eq.(32) in [1].
       lambda1=irrep1%mu
       lambda2=irrep2%mu
       lambda3=irrep3%mu
       mu1=irrep1%lambda
       mu2=irrep2%lambda
       mu3=irrep3%lambda
       irrep1a%lambda=lambda1
       irrep1a%mu=mu1
       irrep3a%lambda=lambda3
       irrep3a%mu=mu3
    END IF


    wigner(0:lambda1,0:lambda2,0:mu2,1:rhomax)=0.D0
    Lambda22max=0

    DO
       IF(MAX(mu2,(lambda1+lambda2-lambda3+2*(mu1+mu2-mu3))/3+1+lambda2)>upbound_binom)THEN
          CALL reallocate_binom(50)
       ELSE
          EXIT
       END IF
    END DO

    DO rho=1,rhomax

       eta=0
       DO
          irrep2a%lambda=lambda2-1
          irrep2a%mu=mu2-1
          IF(outer_multiplicity(irrep1a,irrep2a,irrep3a)==rho-1)EXIT
          lambda2=lambda2-1
          mu2=mu2-1
          eta=eta+1
       END DO
       ! lambda2 is \bar{lambda2}, mu2 is \bar{mu2}
       !*********************************************************************************************************************
       ! Calculation of <(lambda1,mu1)HW;(\bar{lambda2},\bar{mu2})epsilon2,Lambda2||(lambda3,mu3)HW>_rho using Eq.(20) in [1]
       !*********************************************************************************************************************
       epsilon2=lambda1+2*mu1-lambda3-2*mu3
       noname1=(2*lambda2+mu2-epsilon2)/3
       p2max=MIN(lambda2,noname1,(lambda1+lambda3-mu2+noname1)/2)
       p2min=MAX(0,noname1-mu2,(lambda3-lambda1-mu2+noname1+1)/2,(lambda1-lambda3-mu2+noname1+1)/2)
       q2=noname1-p2min
       IF(p2min/=p2max)THEN
          n=(lambda1+lambda2-lambda3+2*(mu1+mu2-mu3))/3
          a=(lambda2+lambda3-lambda1-n)/2
          b=(lambda3+lambda1-lambda2+n+2)/2
          c=(lambda1+lambda2-lambda3-n)/2
          d=(lambda1+lambda2+lambda3+2-n)/2
          j1=lambda2-p2max
          j2=lambda2-p2min
          DO p2=p2min,p2max
             ! The upper and lower bound on p2 are such that:
             ! 1) 0<=p2<=lambda2
             ! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
             ! 3) ABS(lambda1-Lambda22)<=lambda3<=lambda1+Lambda22, where Lambda22=mu2-(2*lambda2+mu2-epsilon2)/3+2*p2
             p2tilde=mu2-q2
             q2tilde=lambda2-p2
             IF(p2tilde==0)THEN
                F=1.D0
             ELSE
                F=DFLOAT((a+1)*(b-1))
                DO j=1,p2tilde-1
                   F=F*DFLOAT((a+j+1)*(b-j-1))
                END DO
                i1=p2tilde*(p2tilde+1)/2
                DO i=1,p2tilde
                   i1=i1+1
                   scalprod=DFLOAT((p2+1)*(mu1+lambda2+mu2-n+2))
                   DO j=1,p2tilde-1
                      IF(j<i)THEN
                         scalprod=scalprod*DFLOAT((p2+j+1)*(mu1+lambda2+mu2-n+j+2))
                      ELSE
                         scalprod=scalprod*DFLOAT((a+j+1)*(b-j-1))
                      END IF
                   END DO
                   F=F+binom(i1)*scalprod
                END DO
                IF(BTEST(p2tilde,0))F=-F
             END IF
             IF(j1<q2tilde)THEN
                G=DFLOAT(INT8((a+n-j1)*(b-n+j1))*(c+n-j1)*(d+n-j1)*(lambda2+mu2-j1+1))
             ELSE
                G=DFLOAT(mu2-n+j1+1)
             END IF
             DO j=j1+1,j2-1
                IF(j<q2tilde)THEN
                   G=G*DFLOAT(INT8((a+n-j)*(b-n+j))*(c+n-j)*(d+n-j)*(lambda2+mu2-j+1))
                ELSE
                   G=G*DFLOAT(mu2-n+j+1)
                END IF
             END DO
             G=G*DFLOAT(lambda2-2*q2tilde+n+1)*binom(n*(n+1)/2+q2tilde)
             i1=n+1+lambda2-q2tilde
             H=binom(i1*(i1+1)/2+lambda2-q2tilde)
             wigner(lambda1,p2,q2,rho)=F*DSQRT(G/H)
             q2=q2-1
          END DO
       ELSE
          wigner(lambda1,p2min,q2,rho)=1.D0
       END IF
       !******************************************************************************************************
       ! Calculation of <(lambda1,mu1)epsilon1,Lambda1;(\bar{lambda2},\bar{mu2})HW||(lambda3,mu3)HW>_rho from
       ! <(lambda1,mu1)HW;(\bar{lambda2},\bar{mu2})epsilon2,Lambda2||(lambda3,mu3)HW>_rho using Eq.(21) in [1]
       !******************************************************************************************************
       noname2=lambda1+mu1+1
       steps21=(epsilon2+lambda2+2*mu2)/3 ! steps21 is the number of iterations in Eq.(21)
       DO i2=1,steps21
          noname1=noname1+1
          noname2=noname2-1
          p1min=MAX(0,noname2-mu1,(lambda1-i2-mu1+noname2+2)/2)
          p1max=MIN(lambda1,noname2-1,(lambda1-1+i2-mu1+noname2)/2)
          p2min=MAX(0,noname1-mu2,(lambda2-steps21+i2-mu2+noname1+1)/2)
          q2=noname1-p2min
          Lambda22=mu2+p2min-q2-2
          DO p2=p2min,MIN(lambda2,noname1,(lambda2+steps21-i2-mu2+noname1)/2)
             ! The upper and lower bound on p2 are such that:
             ! 1) 0<=p2<=lambda2
             ! 2) 0<=q2<=mu2, where q2=(2*lambda2+mu2-epsilon2)/3-p2
             ! 3) lambda2-steps21+i2<=Lambda22<=lambda2+steps21-i2, where Lambda22=mu2+p2-q2
             Lambda22=Lambda22+2 ! Lambda22 is 2*Lambda2 in Eq.(21)
             Sq2=q2*(mu2+1-q2)*INT8(lambda2+mu2+2-q2) ! Sq2 is S(q2)
             Rp2=p2*(lambda2+1-p2)*INT8(mu2+1+p2) ! Rp2 is R(p2)
             IF(lambda1>=i2)THEN
                IF(Lambda22+1<=lambda2+mu2.AND.q2>0)THEN
                   ABCD=(lambda1-i2+Lambda22-lambda3+2)*(lambda1-i2+Lambda22+lambda3+4)
                   IF(ABCD>0)wigner(p1min-1,p2,q2,rho)=-DSQRT(DFLOAT((lambda1-i2+1)*Sq2*ABCD)&
                        /DFLOAT(INT8((lambda1+2-i2)*(Lambda22+1)*p1min)*(lambda1+1-p1min)&
                        *(mu1+1+p1min)*(Lambda22+2)))*wigner(p1min,p2,q2-1,rho)
                END IF
                IF(Lambda22>=1.AND.p2>0)THEN
                   ABCD=(Lambda22+lambda3-lambda1+i2)*(lambda3+lambda1-i2-Lambda22+2)
                   IF(ABCD>0)wigner(p1min-1,p2,q2,rho)=wigner(p1min-1,p2,q2,rho)&
                        -DSQRT(DFLOAT((lambda1-i2+1)*Rp2*ABCD)/DFLOAT(INT8((lambda1+2-i2)*(Lambda22+1)*p1min)&
                        *(lambda1+1-p1min)*(mu1+1+p1min)*Lambda22))*wigner(p1min,p2-1,q2,rho)
                END IF
             END IF
             q1=noname2-p1min
             Lambda12=mu1+p1min-q1-2
             DO p1=p1min,p1max
                ! The upper and lower bound on p1 are such that:
                ! 1) 0<=p1<=lambda1
                ! 2) 0<q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1 cannot be 0 to avoid division by 0
                ! 3) lambda1+1-i2<=Lambda12<=lambda1-1+i2, where Lambda12=mu1+p1-q1
                Lambda12=Lambda12+2
                IF(Lambda22+1<=lambda2+mu2.AND.q2>0)THEN
                   ABCD=(Lambda22+lambda3-Lambda12+1)*(lambda3+Lambda12-Lambda22+1)
                   IF(ABCD>0)wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT((Lambda12+2)*Sq2*ABCD)&
                        /DFLOAT(INT8((Lambda12+1)*(Lambda22+1)*q1)*(mu1+1-q1)*(lambda1+mu1+2-q1)&
                        *(Lambda22+2)))*wigner(p1,p2,q2-1,rho)
                END IF
                IF(Lambda22>=1.AND.p2>0)THEN
                   ABCD=(Lambda12+Lambda22-lambda3+1)*(Lambda12+Lambda22+lambda3+3)
                   IF(ABCD>0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
                        +DSQRT(DFLOAT((Lambda12+2)*Rp2*ABCD)/DFLOAT(INT8((Lambda12+1)*(Lambda22+1)&
                        *q1)*(mu1+1-q1)*(lambda1+mu1+2-q1)*Lambda22))*wigner(p1,p2-1,q2,rho)
                END IF
                q1=q1-1
             END DO
             q2=q2-1
          END DO
       END DO
       !******************************************************************************************************
       ! Calculation of <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)HW||(lambda3,mu3)HW>_rho from
       ! <(lambda1,mu1)epsilon1,Lambda1;(\bar{lambda2},\bar{mu2})HW||(lambda3,mu3)HW>_rho using Eq.(17) in [1]
       !******************************************************************************************************
       noname2=noname2-1
       DO i1=1,eta
          lambda2=lambda2+1
          mu2=mu2+1
          noname2=noname2-1
          IF(noname2==mu1.AND.wigner(1,lambda2-1,mu2-1,rho)/=0.D0)wigner(0,lambda2,mu2,rho)&
               =-DSQRT(DFLOAT(INT8(lambda1*(mu1+2))&
               *(lambda2+lambda3)*(lambda3-lambda2+2))/2.D0)*wigner(1,lambda2-1,mu2-1,rho)
          p1min=MAX(0,noname2-mu1,(2-mu1+noname2)/2)
          q1=noname2-p1min
          Lambda12=mu1+p1min-q1-2
          DO p1=p1min,MIN(lambda1,noname2,(lambda1-1+noname2)/2)
             ! The upper and lower bound on p1 are such that:
             ! 1) 0<=p1<=lambda1
             ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
             ! 3) 1<=Lambda12<=lambda1+mu1-1, where Lambda12=mu1+p1-q1
             Lambda12=Lambda12+2
             IF(wigner(p1+1,lambda2-1,mu2-1,rho)/=0.D0)wigner(p1,lambda2,mu2,rho)=-DSQRT(DFLOAT(INT8((p1+1)&
                  *(lambda1-p1))*(mu1+2+p1)*(lambda2+lambda3-Lambda12)*(lambda3+Lambda12-lambda2+2))&
                  /DFLOAT((Lambda12+2)*(Lambda12+1)))*wigner(p1+1,lambda2-1,mu2-1,rho)
             IF(wigner(p1,lambda2-1,mu2-1,rho)/=0.D0)wigner(p1,lambda2,mu2,rho)=wigner(p1,lambda2,mu2,rho)&
                  +DSQRT(DFLOAT(INT8((q1+1)*(mu1-q1))*(lambda1+mu1+1-q1)*(Lambda12+lambda2-lambda3)&
                  *(Lambda12+lambda2+lambda3+2))/DFLOAT((Lambda12+1)*Lambda12))*wigner(p1,lambda2-1,mu2-1,rho)
             q1=q1-1
          END DO
          IF(noname2==lambda1.AND.wigner(lambda1,lambda2-1,mu2-1,rho)/=0.D0)THEN
             Lambda12=lambda1+mu1
             wigner(lambda1,lambda2,mu2,rho)=DSQRT(DFLOAT(INT8(mu1*(lambda1+mu1+1))*(Lambda12+lambda2-lambda3)&
                  *(Lambda12+lambda2+lambda3+2))/DFLOAT((Lambda12+1)*Lambda12))*wigner(lambda1,lambda2-1,mu2-1,rho)
          END IF
       END DO

       ! Now the arrays p1a,p2a,q2a are filled with values corresponding to coefficients
       ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)HW||(lambda3,mu3)HW>.
       IF(rho==rhomax)THEN
          i2=0
          p1min=MAX(0,noname2-mu1,(lambda3-mu1+noname2-lambda2+1)/2,(lambda2-lambda3-mu1+noname2+1)/2)
          DO p1=p1min,MIN(lambda1,noname2,(lambda2+lambda3-mu1+noname2)/2)
             ! The upper and lower bound on p1 are such that:
             ! 1) 0<=p1<=lambda1
             ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
             ! 3) ABS(Lambda12-lambda2)<=lambda3<=Lambda12+lambda2, where Lambda12=mu1+p1-q1
             i2=i2+1
             p1a(i2)=p1
             p2a(i2)=lambda2
             q2a(i2)=mu2
          END DO
          IF(-lambda1-2*mu1+3*(steps21+eta)==-lambda1-2*mu1)Lambda22max=lambda2
          ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW.
          ! This is needed for setting the phase.
       END IF
       !**************************************************************************************************
       ! Calculation of <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)HW>_rho
       ! from <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)HW||(lambda3,mu3)HW>_rho using Eq.(18) in [1]
       !**************************************************************************************************
       noname1=lambda2+mu2
       DO i1=1,lambda2+mu2 ! This is loop over epsilon2=epsilon2HW+3,...,2*lambda2+mu2; i1 is (epsilon2-epsilon2HW)/3
          noname2=noname2+1
          IF(noname2>lambda1+mu1)EXIT
          noname1=noname1-1
          p2min=MAX(0,noname1-mu2,(lambda2-i1-mu2+noname1+1)/2)
          q2=noname1-p2min
          Lambda22=mu2+p2min-q2-2 ! Lambda22 is 2*Lambda2' in Eq.(18)
          DO p2=p2min,MIN(lambda2,noname1,(lambda2+i1-mu2+noname1)/2)
             Lambda22=Lambda22+2
             p1min=MAX(0,noname2-mu1,(lambda3-mu1+noname2-Lambda22+1)/2,(Lambda22-lambda3-mu1+noname2+1)/2)
             q1=noname2-p1min
             Lambda12=mu1+p1min-q1-2
             DO p1=p1min,MIN(lambda1,noname2,(Lambda22+lambda3-mu1+noname2)/2)
                ! The upper and lower bound on p1 are such that:
                ! 1) 0<=p1<=lambda1
                ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                ! 3) ABS(Lambda12-Lambda22)<=lambda3<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1
                Lambda12=Lambda12+2
                IF(Lambda22==lambda2-i1.OR.Lambda22==0)THEN
                   IF(Lambda12>=1.AND.p1/=0)THEN
                      wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8((Lambda22+1)*p1*(lambda1+1-p1))*(mu1+1+p1)&
                           *(Lambda22+lambda3-Lambda12+2)*(lambda3+Lambda12-Lambda22))/DFLOAT(INT8((Lambda12+1)&
                           *(Lambda22+2)*(p2+1))*(lambda2-p2)*(mu2+2+p2)*4*Lambda12))*wigner(p1-1,p2+1,q2,rho)
                   ELSE
                      wigner(p1,p2,q2,rho)=0.D0
                   END IF
                   IF(Lambda12+1<=lambda1+mu1.AND.q1/=0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
                        +DSQRT(DFLOAT(INT8((Lambda22+1)*q1*(mu1+1-q1))*(lambda1+mu1+2-q1)*(Lambda12+Lambda22-lambda3+2)&
                        *(Lambda12+Lambda22+lambda3+4))/DFLOAT(INT8((Lambda12+1)*(Lambda22+2)*(p2+1))*(lambda2-p2)&
                        *(mu2+2+p2)*4*(Lambda12+2)))*wigner(p1,p2+1,q2,rho)
                ELSE
                   IF(Lambda12>=1.AND.p1/=0)THEN
                      wigner(p1,p2,q2,rho)=-DSQRT(DFLOAT(INT8((Lambda22+1)*p1*(lambda1+1-p1))*(mu1+1+p1)&
                           *(Lambda12+Lambda22-lambda3)*(Lambda12+Lambda22+lambda3+2))/DFLOAT(INT8((Lambda12+1)&
                           *Lambda22*(q2+1))*(mu2-q2)*(lambda2+mu2+1-q2)*4*Lambda12))*wigner(p1-1,p2,q2+1,rho)
                   ELSE
                      wigner(p1,p2,q2,rho)=0.D0
                   END IF
                   IF(Lambda12+1<=lambda1+mu1.and.q1/=0)wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)&
                        -DSQRT(DFLOAT(INT8((Lambda22+1)*q1*(mu1+1-q1))*(lambda1+mu1+2-q1)*(Lambda22+lambda3-Lambda12)&
                        *(lambda3+Lambda12-Lambda22+2))/DFLOAT(INT8((Lambda12+1)*Lambda22*(q2+1))*(mu2-q2)&
                        *(lambda2+mu2+1-q2)*4*(Lambda12+2)))*wigner(p1,p2,q2+1,rho)
                END IF

                IF(rho==rhomax)THEN
                   i2=i2+1
                   IF(p1==lambda1.AND.q1==mu1)Lambda22max=MAX(Lambda22max,Lambda22)
                   ! Lambda22max is the greatest value of 2*Lambda2 for which (epsilon1,Lambda1)=HW.
                   ! This is needed for setting the phase.
                   p1a(i2)=p1
                   p2a(i2)=p2
                   q2a(i2)=q2
                END IF

                q1=q1-1
             END DO
             q2=q2-1
          END DO
       END DO

       ! Valid values of p1, p2 and q2 are elements of arrays p1a, p2a and q2a with indeces from 1 to i2.
    END DO ! End of the loop over rho
    !**********************************************************************
    ! Orthonormalization according to (8) with alpha3=(epsilon3,Lambda3)=HW
    !**********************************************************************
    ! Gram-Schmidt orthonormalization is performed, where Wigner coefficients for given rho represent a vector,
    ! whose components are indexed by alpha1 and alpha2, and the scalar product is defined in Eq.(8) in [1].
    DO rho=1,rhomax
       DO i4=1,rho
          p1=p1a(1)
          p2=p2a(1)
          q2=q2a(1)
          scalprod=wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,i4)
          F=0.D0
          DO i1=2,i2 ! Kahan summation fomula (see Ref.[3]) is used to calculate the scalar product.
             p1=p1a(i1)
             p2=p2a(i1)
             q2=q2a(i1)
             G=wigner(p1,p2,q2,rho)*wigner(p1,p2,q2,i4)-F
             H=scalprod+G
             F=H-scalprod
             F=F-G
             scalprod=H
          END DO
          IF(i4/=rho)THEN
             DO i1=1,i2
                p1=p1a(i1)
                p2=p2a(i1)
                q2=q2a(i1)
                wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)-scalprod*wigner(p1,p2,q2,i4)
             END DO
          ELSE
             scalprod=1.D0/DSQRT(scalprod)
             DO i1=1,i2
                p1=p1a(i1)
                p2=p2a(i1)
                q2=q2a(i1)
                wigner(p1,p2,q2,rho)=wigner(p1,p2,q2,rho)*scalprod
             END DO
          END IF
       END DO
    END DO
    !***************************************
    ! Setting the phase according to Ref.[2]
    !***************************************
    ! The phase is such that <(lambda1,mu1)LW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)LW>_rho>0,
    ! which, according to formula 2B from Eq.(35) in [1] with rho_max instead of eta_max (see text therein),
    ! means that <(lamda1,mu1)HW;(lambda2,mu2)epsilon2,Lambda2_max||(lambda3,mu3)HW>_rho*(-1)^e>0,
    ! where e=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax-rho+(lambda1+2*Lambda2_max-lambda3)/2.
    phiprhomax=lambda1+lambda2-lambda3+mu1+mu2-mu3+rhomax ! phiprhomax is phi+rhomax
    i4=phiprhomax+(lambda1+Lambda22max-lambda3)/2 ! i4 is e+rho
    p2=(noname1-mu2+Lambda22max)/2
    q2=noname1-p2
    DO rho=1,rhomax
       IF((BTEST(i4-rho+1,0).AND.wigner(lambda1,p2,q2,rho)<0.D0)&
            .OR.(BTEST(i4-rho,0).AND.wigner(lambda1,p2,q2,rho)>0.D0))THEN
          DO i1=1,i2
             p1=p1a(i1)
             p2=p2a(i1)
             q2=q2a(i1)
             wigner(p1,p2,q2,rho)=-wigner(p1,p2,q2,rho)
          END DO
       END IF
    END DO

    IF(I3==0)THEN ! E=LW
       DO i1=1,i2
          p1=p1a(i1)
          p2=p2a(i1)
          DO rho=1,rhomax
             ! See Eq.(35,2B)
             IF(BTEST(phiprhomax-rho+p1+p2+(mu1-(2*(lambda1+lambda2+mu3)+mu1+mu2+lambda3)/3+mu2-lambda3)/2,0))THEN
                q2=q2a(i1)
                wigner(p1,p2,q2,rho)=-wigner(p1,p2,q2,rho)
             END IF
          END DO
       END DO
    END IF
  END SUBROUTINE wigner_canonical_extremal

  SUBROUTINE wigner_canonical(irrep1,irrep2,irrep3,epsilon3,Lambda32,I3,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
    !----------------------------------------------------------------------------------------------------------------------
    ! Calculates the reduced SU(3)-SU(2)xU(1) Wigner coefficients
    ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
    ! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3 using Eq.(19) in [1] and its
    ! conjugate and Table 9.1 on page 311 in [2].
    !
    ! References: [1] J.P.Draayer, Y.Akiyama, J.Math.Phys., Vol.14, No.12 (1973) 1904
    !             [2] Varshalovich, Quantum theory of angular momentum
    !
    ! Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,I3,rhomax,wignerex
    ! Input and output arguments: numb,p1a,p2a,q2a,wigner
    !
    ! lambda1=irrep1%lambda, mu1=irrep1%mu
    ! lambda2=irrep2%lambda, mu2=irrep2%mu
    ! lambda3=irrep3%lambda, mu3=irrep3%mu
    ! Lambda32=2*Lambda3 (epsilon3 and Lambda32 must be valid)
    ! rhomax=multiplicity of coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3) which must be greater than 0
    ! wignerex is an array containing the extremal coefficients calculated by subroutine wigner_canonical_extremal.
    ! If I3=1, coefficients are calculated from the highest-weight (HW) coefficients.
    ! If I3=0, coefficients are calculated from the lowest-weight (LW) coefficients.
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
    !   unless I3=0 and epsilon3=epsilon3LW=2*lambda3+mu3, in this case:
    !         q1=mu1-p1a(i)
    !         p2=lambda2-q2a(i)
    !         q2=mu2-p2a(i)
    !         p1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-q1-p2-q2
    !
    ! Note: There is a typo in Eq.(19). In the very last line there should be p instead of q.
    !----------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
    INTEGER,INTENT(IN) :: epsilon3,Lambda32,I3,rhomax
    INTEGER :: numb,eps3,rho,p1,q1,p2,q2,Lam32,epsilon2max,epsilon3ex,Lambda22,Lambda12,Lam32prime,p3,q3,&
         pq1,noname1,noname2,s2,lm1,lm2,lm3,mu1p,mu2p,mu3p,q2ex,p1ex,lambda1p,lambda2p,lambda3p
    INTEGER(KIND=8) :: Rp2,Sq2
    REAL(KIND=8) :: N3
    INTEGER,DIMENSION(:) :: p1a,p2a,q2a ! Size is at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1)
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wignerex
    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
    ! Sizes of wignerex and wigner are at least (0:MAX(lambda1,mu1),0:MAX(lambda2,mu2),0:MAX(lambda2,mu2),1:rhomax).

    IF(I3==1)THEN

       DO rho=1,rhomax
          DO q2=0,irrep2%mu
             DO p2=0,irrep2%lambda
                wigner(0:irrep1%lambda,p2,q2,rho)=wignerex(0:irrep1%lambda,p2,q2,rho)
             END DO
          END DO
       END DO

       epsilon3ex=-irrep3%lambda-2*irrep3%mu
       IF(epsilon3==epsilon3ex)RETURN

       lm1=irrep1%lambda+irrep1%mu+1
       lm2=irrep2%lambda+irrep2%mu+1
       lm3=irrep3%lambda+irrep3%mu+1
       mu1p=irrep1%mu+2
       mu2p=irrep2%mu+2
       mu3p=irrep3%mu+2

       epsilon2max=2*irrep2%lambda+irrep2%mu
       noname1=(epsilon2max-irrep1%lambda-2*irrep1%mu-epsilon3ex)/3
       noname2=(epsilon2max+2*irrep1%lambda+irrep1%mu-epsilon3ex)/3
       epsilon3ex=epsilon3ex+3

       numb=0

       Lam32=irrep3%lambda
       p3=irrep3%lambda
       q3=irrep3%mu

       DO eps3=epsilon3ex,epsilon3,3 ! eps3 is epsilon3 in Eq.(19)
          noname1=noname1-1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
          noname2=noname2-1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3

          Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
          IF(Lam32<Lambda32)THEN
             Lam32=Lam32+1
             q3=q3-1 ! p3 and q3 correspond to eps3 and Lam32
             N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(INT8((q3+1)*(irrep3%mu-q3))*(lm3-q3)*(Lam32prime+2)*4))
          ELSE
             Lam32=Lam32-1
             p3=p3-1 ! p3 and q3 correspond to eps3 and Lam32
             N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(INT8((p3+1)*(irrep3%lambda-p3))*(mu3p+p3)*Lam32prime*4))
          END IF ! Lam32 is 2*Lambda3 in Eq.(19) and N3 is sqrt((2*Lambda3'+1)*(2*Lambda3+1))/N3.

          DO p2=MAX(0,noname1-irrep2%mu,(Lam32-irrep1%mu+noname2-irrep2%mu)/2-irrep1%lambda),MIN(irrep2%lambda,noname2)
             Rp2=INT8(p2+1)*(irrep2%lambda-p2)*(mu2p+p2)
             q2ex=MAX(0,noname1-p2,(irrep2%mu+noname2-Lam32-irrep1%mu)/2-irrep1%lambda)
             pq1=noname2-p2-q2ex ! pq1 is p1+q1
             Lambda22=irrep2%mu+p2-q2ex ! Lambda22 is 2*Lambda2 in Eq.(19)
             DO q2=q2ex,MIN(irrep2%mu,noname2-p2,(irrep2%mu+noname2+Lam32-irrep1%mu)/2)
                ! The lower and upper bounds on q2 are such that:
                ! 1) 0<=q2<=mu2
                ! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
                !     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
                !     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
                Sq2=INT8(q2+1)*(irrep2%mu-q2)*(lm2-q2)
                p1ex=MAX(0,pq1-irrep1%mu,(Lam32-irrep1%mu+pq1-Lambda22)/2,(Lambda22-Lam32-irrep1%mu+pq1)/2)
                q1=pq1-p1ex
                Lambda12=irrep1%mu+p1ex-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
                DO p1=p1ex,MIN(irrep1%lambda,pq1,(Lambda22+Lam32-irrep1%mu+pq1)/2)
                   ! The lower and upper bounds on p1 are such that:
                   ! 1) 0<=p1<=lambda1
                   ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                   ! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

                   IF(eps3==epsilon3)THEN
                      numb=numb+1
                      p1a(numb)=p1
                      p2a(numb)=p2
                      q2a(numb)=q2
                   END IF

                   s2=Lambda12+Lambda22+Lam32prime+1
                   IF(Lam32>Lam32prime)THEN

                      IF(q1/=irrep1%mu)THEN
                         wigner(p1,p2,q2,1:rhomax)=DSQRT(DFLOAT(INT8((q1+1)*(irrep1%mu-q1))*(lm1-q1)&
                              *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,1:rhomax)
                      ELSE
                         wigner(p1,p2,q2,1:rhomax)=0.D0
                      END IF

                      IF(p1/=irrep1%lambda)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)&
                           +DSQRT(DFLOAT(INT8((p1+1)*(irrep1%lambda-p1))*(mu1p+p1)&
                           *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,1:rhomax)

                      IF(p2/=irrep2%lambda)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)-DSQRT(DFLOAT(Rp2&
                           *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,1:rhomax)

                      IF(q2/=irrep2%mu)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
                           *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,1:rhomax)

                   ELSE

                      IF(q1/=irrep1%mu)THEN
                         wigner(p1,p2,q2,1:rhomax)=-DSQRT(DFLOAT(INT8((q1+1)*(irrep1%mu-q1))*(lm1-q1)&
                              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1,p2,q2,1:rhomax)
                      ELSE
                         wigner(p1,p2,q2,1:rhomax)=0.D0
                      END IF

                      IF(p1/=irrep1%lambda)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)&
                           +DSQRT(DFLOAT(INT8((p1+1)*(irrep1%lambda-p1))*(mu1p+p1)&
                           *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1+1,p2,q2,1:rhomax)

                      IF(p2/=irrep2%lambda)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Rp2&
                           *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2+1,q2,1:rhomax)

                      IF(q2/=irrep2%mu)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
                           *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2,q2+1,1:rhomax)

                   END IF

                   wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)*N3

                   q1=q1-1
                   Lambda12=Lambda12+2 ! Lambda12 is 2*Lambda1 in Eq.(19)
                END DO
                pq1=pq1-1 ! pq1 is p1+q1
                Lambda22=Lambda22-1 ! Lambda22 is 2*Lambda2 in Eq.(19)
             END DO
          END DO
       END DO

    ELSE

       epsilon3ex=2*irrep3%lambda+irrep3%mu
       epsilon2max=2*irrep2%lambda+irrep2%mu
       noname1=(epsilon2max-irrep1%lambda-2*irrep1%mu-epsilon3ex)/3
       noname2=(epsilon2max+2*irrep1%lambda+irrep1%mu-epsilon3ex)/3
       DO p1=0,irrep1%lambda
          DO p2=0,irrep2%lambda
             DO q2=0,irrep2%mu
                !       q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3max)/3-p1-p2-q2
                q1=noname2-p1-p2-q2
                wigner(p1,p2,q2,1:rhomax)=wignerex(irrep1%mu-q1,irrep2%mu-q2,irrep2%lambda-p2,1:rhomax)
             END DO
          END DO
       END DO

       IF(epsilon3==epsilon3ex)RETURN

       lambda1p=irrep1%lambda+1
       lambda2p=irrep2%lambda+1
       lambda3p=irrep3%lambda+1
       mu1p=irrep1%mu+1
       mu2p=irrep2%mu+1
       mu3p=irrep3%mu+1
       lm1=lambda1p+mu1p
       lm2=lambda2p+mu2p
       lm3=lambda3p+mu3p

       epsilon3ex=epsilon3ex-3

       numb=0

       Lam32=irrep3%mu
       p3=0
       q3=0

       DO eps3=epsilon3ex,epsilon3,-3 ! eps3 is epsilon3 in Eq.(19)
          noname1=noname1+1 ! noname1 is (epsilon1min-eps3+epsilon2max)/3
          noname2=noname2+1 ! noname2 is (epsilon1max-eps3+epsilon2max)/3

          Lam32prime=Lam32 ! Lam32prime is 2*Lambda3' in Eq.(19)
          IF(Lam32<Lambda32)THEN
             Lam32=Lam32+1
             p3=p3+1 ! p3 and q3 correspond to eps3 and Lam32
             N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(INT8(p3*(lambda3p-p3))*(mu3p+p3)*(Lam32prime+2)*4))
             ! N3 is sqrt((2*Lambda3+1)/(4*(2*Lambda3'+1)))/N3
          ELSE
             Lam32=Lam32-1
             q3=q3+1 ! p3 and q3 correspond to eps3 and Lam32
             N3=DSQRT(DFLOAT(Lam32+1)/DFLOAT(INT8(q3*(mu3p-q3))*(lm3-q3)*Lam32prime*4))
             ! N3 is sqrt((2*Lambda3+1)/(8*Lambda3'))/N3
          END IF ! Lam32 is 2*Lambda3 in Eq.(19).

          DO p2=MIN(irrep2%lambda,noname2),MAX(0,noname1-irrep2%mu,(Lam32-irrep1%mu+noname2-irrep2%mu)/2-irrep1%lambda),-1
             Rp2=INT8(p2)*(lambda2p-p2)*(mu2p+p2) ! Rp2 is R(p2)
             q2ex=MIN(irrep2%mu,noname2-p2,(irrep2%mu+noname2+Lam32-irrep1%mu)/2)
             pq1=noname2-p2-q2ex ! pq1 is p1+q1
             Lambda22=irrep2%mu+p2-q2ex ! Lambda22 is 2*Lambda2 in Eq.(19)
             DO q2=q2ex,MAX(0,noname1-p2,(irrep2%mu+noname2-Lam32-irrep1%mu)/2-irrep1%lambda),-1
                ! The lower and upper bounds on q2 are such that:
                ! 1) 0<=q2<=mu2
                ! 2) -lambda1-2*mu1<=epsilon1<=2*lambda1+mu1, where epsilon1=eps3-epsilon2max+3*p2+3*q2
                !     epsilon2=epsilon2max-3*(p2+q2) ! epsilon2 is epsilon2 in Eq.(19)
                !     epsilon1=eps3-epsilon2 ! epsilon1 is epsilon1 in Eq.(19)
                Sq2=INT8(q2)*(mu2p-q2)*(lm2-q2)
                p1ex=MIN(irrep1%lambda,pq1,(Lambda22+Lam32-irrep1%mu+pq1)/2)
                q1=pq1-p1ex
                Lambda12=irrep1%mu+p1ex-q1 ! Lambda12 is 2*Lambda1 in Eq.(19)
                DO p1=p1ex,MAX(0,pq1-irrep1%mu,(Lam32-irrep1%mu+pq1-Lambda22)/2,(Lambda22-Lam32-irrep1%mu+pq1)/2),-1
                   ! The lower and upper bounds on p1 are such that:
                   ! 1) 0<=p1<=lambda1
                   ! 2) 0<=q1<=mu1, where q1=(2*lambda1+mu1-epsilon1)/3-p1
                   ! 3) ABS(Lambda12-Lambda22)<=Lam32<=Lambda12+Lambda22, where Lambda12=mu1+p1-q1

                   IF(eps3==epsilon3)THEN
                      numb=numb+1
                      p1a(numb)=p1
                      p2a(numb)=p2
                      q2a(numb)=q2
                   END IF

                   s2=Lambda12+Lambda22+Lam32prime+1
                   IF(Lam32>Lam32prime)THEN

                      ! term with Lambda1p=Lambda1+1/2
                      IF(q1/=0)THEN
                         wigner(p1,p2,q2,1:rhomax)=-DSQRT(DFLOAT(INT8(q1*(mu1p-q1))*(lm1-q1)&
                              *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1,p2,q2,1:rhomax)
                      ELSE
                         wigner(p1,p2,q2,1:rhomax)=0.D0
                      END IF

                      ! term with Lambda1p=Lambda1-1/2
                      IF(p1/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(INT8(p1*(lambda1p-p1))*(mu1p+p1)&
                           *(s2+2)*(s2-2*Lambda22))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1-1,p2,q2,1:rhomax)

                      ! term with Lambda2p=Lambdda2+1/2
                      IF(q2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
                           *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2,q2-1,1:rhomax)

                      ! term with Lambda2p=Lambda2-1/2
                      IF(p2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Rp2&
                           *(s2+2)*(s2-2*Lambda12))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2-1,q2,1:rhomax)

                   ELSE

                      ! term with Lambda1p=Lambda1+1/2
                      IF(q1/=0)THEN
                         wigner(p1,p2,q2,1:rhomax)=DSQRT(DFLOAT(INT8(q1*(mu1p-q1))*(lm1-q1)&
                              *(s2+2)*(s2-2*Lambda22))/DFLOAT((Lambda12+1)*(Lambda12+2)))*wigner(p1,p2,q2,1:rhomax)
                      ELSE
                         wigner(p1,p2,q2,1:rhomax)=0.D0
                      END IF

                      ! term with Lambda1p=Lambda1-1/2
                      IF(p1/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(INT(p1)*(lambda1p-p1)*(mu1p+p1)&
                           *(s2-2*Lambda12)*(s2-2*Lam32prime))/DFLOAT(Lambda12*(Lambda12+1)))*wigner(p1-1,p2,q2,1:rhomax)

                      ! term with Lambda2p=Lambda2+1/2
                      IF(q2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)+DSQRT(DFLOAT(Sq2&
                           *(s2+2)*(s2-2*Lambda12))/DFLOAT((Lambda22+1)*(Lambda22+2)))*wigner(p1,p2,q2-1,1:rhomax)

                      ! term with Lambda2p=Lambda2-1/2
                      IF(p2/=0)wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)-DSQRT(DFLOAT(Rp2&
                           *(s2-2*Lam32prime)*(s2-2*Lambda22))/DFLOAT(Lambda22*(Lambda22+1)))*wigner(p1,p2-1,q2,1:rhomax)

                   END IF

                   wigner(p1,p2,q2,1:rhomax)=wigner(p1,p2,q2,1:rhomax)*N3

                   q1=q1+1
                   Lambda12=Lambda12-2 ! Lambda12 is 2*Lambda1 in Eq.(19)
                END DO
                pq1=pq1+1 ! pq1 is p1+q1
                Lambda22=Lambda22+1 ! Lambda22 is 2*Lambda2 in Eq.(19)
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE wigner_canonical

  SUBROUTINE calculate_wigner_canonical(irrep1,irrep2,irrep3,epsilon3,Lambda32,&
       dimpq,dimw,rhomax,numb,wigner_block,p1a,p2a,q2a) BIND(C)
    !-------------------------------------------------------------------------------------------------------------------
    ! Calculates all SU(3)-SU(2)xU(1) reduced Wigner coefficients
    ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho
    ! for given lambda1,mu1,lambda2,mu2,lambda3,mu3,epsilon3,Lambda3=Lambda32/2. The quantum numbers epsilon3,Lambda32
    ! must be valid, i.e. there must be integers p3,q3 satisfying:
    !   0<=p3<=lambda3
    !   0<=q3<=mu3
    !   epsilon3=2*lambda3+mu3-3*(p3+q3)
    !   Lambda3=(mu3+p3-q3)/2
    !
    ! Input arguments: irrep1,irrep2,irrep3,epsilon3,Lambda32,dimpq,dimw,rhomax
    ! Output arguments: numb,wigner_block,p1a,p2a,q2a
    !
    ! lambda1=irrep1%lambda, mu1=irrep1%mu
    ! lambda2=irrep2%lambda, mu2=irrep2%mu
    ! lambda3=irrep3%lambda, mu3=irrep3%mu
    ! rmomax is the multiplicity of the coupling (lambda1,mu1)x(lambda2,mu2)->(lambda3,mu3).
    ! dimpq is the size of the arrays p1a,p2a,q2a, which should be at least (MAX(lambda1,mu1)+1)*(lambda2+1)*(mu2+1).
    ! dimw is the size of the array wigner_block, which should be at least rhomax*dimpq.
    !
    ! <(lambda1,mu1)epsilon1,Lambda1;(lambda2,mu2)epsilon2,Lambda2||(lambda3,mu3)epsilon3,Lambda3>_rho=wigner_block(ind)
    !   where epsilon2=2*lambda2+mu2-3*(p2+q2)
    !         epsilon1=epsilon3-epsilon2
    !         Lambda1=(mu1+p1-q1)/2
    !         Lambda2=(mu2+p2-q2)/2
    !         p1=p1a(i) 
    !         p2=p2a(i)
    !         q2=q2a(i)
    !         q1=(2*(lambda1+lambda2)+mu1+mu2-epsilon3)/3-p1-p2-q2
    !         ind=i+numb*(rho-1)
    !         1<=i<=numb
    !-------------------------------------------------------------------------------------------------------------------
    !USE iso_c_binding
    IMPLICIT NONE
    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
    INTEGER(C_INT),INTENT(IN) :: epsilon3,Lambda32,dimpq,dimw,rhomax
    INTEGER(C_INT),INTENT(OUT) :: numb
    INTEGER(C_INT),DIMENSION(dimpq),INTENT(OUT) :: p1a,p2a,q2a
    REAL(C_DOUBLE),DIMENSION(dimw),INTENT(OUT) :: wigner_block
    INTEGER :: i,rho,ind
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:,:) :: wignerex,wigner

    !INTERFACE
    !  SUBROUTINE wigner_canonical_extremal(irrep1,irrep2,irrep3,I3,rhomax,i2,wigner,p1a,p2a,q2a)
    !    IMPLICIT NONE
    !    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3
    !    INTEGER,INTENT(IN) :: I3,rhomax
    !    INTEGER,INTENT(OUT) :: i2
    !    INTEGER,DIMENSION(:),INTENT(OUT) :: p1a,p2a,q2a
    !    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
    !  END SUBROUTINE wigner_canonical_extremal
    !  SUBROUTINE wigner_canonical(irrep1,irrep2,irrep3,epsilon3,Lambda32,I3,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
    !    IMPLICIT NONE
    !    TYPE(su3irrep),INTENT(IN) :: irrep1,irrep2,irrep3 
    !    INTEGER,INTENT(IN) :: epsilon3,Lambda32,I3,rhomax
    !    INTEGER :: numb
    !    INTEGER,DIMENSION(:) :: p1a,p2a,q2a
    !    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(IN) :: wignerex
    !    REAL(KIND=8),DIMENSION(0:,0:,0:,1:),INTENT(OUT) :: wigner
    !  END SUBROUTINE wigner_canonical
    !END INTERFACE

    i=MAX(irrep2%lambda,irrep2%mu)
    ALLOCATE(wignerex(0:MAX(irrep1%lambda,irrep1%mu),0:i,0:i,1:rhomax),wigner(0:MAX(irrep1%lambda,irrep1%mu),0:i,0:i,1:rhomax))

    IF(2*epsilon3<=irrep3%lambda-irrep3%mu)THEN
       CALL wigner_canonical_extremal(irrep1,irrep2,irrep3,1,rhomax,numb,wignerex,p1a,p2a,q2a)
       CALL wigner_canonical(irrep1,irrep2,irrep3,epsilon3,Lambda32,1,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
    ELSE
       CALL wigner_canonical_extremal(irrep1,irrep2,irrep3,0,rhomax,numb,wignerex,p1a,p2a,q2a)
       CALL wigner_canonical(irrep1,irrep2,irrep3,epsilon3,Lambda32,0,rhomax,numb,wignerex,wigner,p1a,p2a,q2a)
       IF(epsilon3==2*irrep3%lambda+irrep3%mu)THEN
          rho=(2*(irrep1%lambda-irrep1%mu-irrep2%mu)-irrep2%lambda-epsilon3)/3
          DO i=1,numb
             ind=p2a(i)
             p1a(i)=rho+p1a(i)+q2a(i)+ind
             p2a(i)=irrep2%lambda-q2a(i)
             q2a(i)=irrep2%mu-ind
          END DO
       END IF
    END IF

    ind=0
    DO rho=1,rhomax
       DO i=1,numb
          ind=ind+1 ! ind=i+numb*(rho-1)
          wigner_block(ind)=wigner(p1a(i),p2a(i),q2a(i),rho)
       END DO
    END DO

    DEALLOCATE(wignerex,wigner)

  END SUBROUTINE calculate_wigner_canonical

END MODULE ndsu3lib_wigner_canonical
