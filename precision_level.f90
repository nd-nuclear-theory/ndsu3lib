!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! precision_level.f90 -- level of arbitrary precision
!
! Jakub Herko
! University of Notre Dame
!
! SPDX-License-Identifier: MIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE precision_level
USE mpmodule
IMPLICIT NONE
INTEGER :: ndig,nwds
PARAMETER(ndig=37,nwds=INT(ndig/mpdpw+2)) ! 37 is the minimal ndig
END MODULE precision_level
