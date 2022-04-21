! This code is a mini version of ndsu3lib and a program that calls it.
! It contains a module representing ndsu3lib and an OpenMP parallel program representing
! a program that calls ndsu3lib. The module contains a global array and 2 subroutines.
! The 1st subroutine allocates the array and fills it with numbers. The 2nd subroutine
! does some calculations using the array, but first it checks whether the size of the
! array is sufficient and resizes the array if needed. It uses the rwlock.

MODULE modul
  ! This module represents ndsu3lib. It contains an array, a subroutine allocating
  ! the array and filling it with non-negative integers, and a subroutine using the array.
  USE rwlock
  IMPLICIT NONE
  INTEGER,ALLOCATABLE,DIMENSION(:) :: array
  INTEGER :: upbound_array ! the size of `array`
  TYPE(ReadWriteLock_t) :: lock
CONTAINS

  SUBROUTINE allocate_array(upbound)
    ! allocates `array` and fills it with non-negative integers
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: upbound
    INTEGER :: i
    upbound_array=upbound
    ALLOCATE(array(upbound_array))
    DO i=1,upbound_array
       array(i)=i
    END DO
  END SUBROUTINE allocate_array

  SUBROUTINE routine(numb,suma)
    ! calculates the sum of the first `numb` elements of `array`
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: numb
    INTEGER,INTENT(OUT) :: suma
    INTEGER :: i

    ! acquire a reader lock
    CALL lock%reader_lock()

    ! if we need to resize `array`, try to get a writer lock
    DO WHILE(upbound_array<numb)
       IF(.NOT.lock%writer_lock(.TRUE.))THEN
          ! if we didn't get the writer lock, some other thread is resizing
          ! the array -- wait for it to finish, and re-check if resizing is
          ! still required
          CALL lock%reader_unlock()
          CALL lock%reader_lock()
          ! flush, to get update from other thread which did update
          !$omp flush acquire
          CYCLE
       END IF
       ! flush, to get any updates from other writers
       !$omp flush acquire

       ! resize the array
       DEALLOCATE(array)
       CALL allocate_array(upbound_array+10)
       ! flush, to publish update to other threads
       !$omp flush release

       ! release the writer lock
       CALL lock%writer_unlock(.TRUE.)
    END DO

    suma=0
    DO i=1,numb
       suma=suma+array(i)
    END DO

    ! release the reader lock
    CALL lock%reader_unlock()

  END SUBROUTINE routine

END MODULE modul

PROGRAM rwlock_example
  ! This program represents a program that calls ndsu3lib.
  USE modul
  USE omp_lib
  IMPLICIT NONE
  INTEGER :: numb,suma

  CALL allocate_array(10)

  ! initialize the lock
  CALL lock%init()

  !$omp parallel do schedule(static) default(none) &
  !$omp     private(numb, suma) &
  !$omp     shared(lock, upbound_array, array)
  DO numb=1,50
     CALL routine(numb,suma)
  END DO
  !$omp end parallel do

  DEALLOCATE(array)

END PROGRAM rwlock_example
