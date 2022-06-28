MODULE Constants

  implicit none

  integer :: Nl(3)
  double precision :: t

!#include "convert.F90" 
  
! integer, parameter :: nb=3 ! number of bands      

  ! Lattice and frequency grid parameters. (all must be 2^n!)
! integer, parameter :: m=8192
  
! integer, parameter :: lcx =4 
! integer, parameter :: lcy = 4 
! integer, parameter :: lcz = 1
! integer, parameter :: nc = lcx*lcy*lcz  

! integer, parameter :: nl = llx*lly*llz

  ! MPI related constants.
! integer, parameter :: mp = m/np  ! must be an integer
! integer, parameter :: ncp = nc/np  ! must be an integer

  ! Convenience constants
! integer, parameter :: mp1 = mp - 1
! integer, parameter :: m1 = m - 1

! integer, parameter :: llx1 = llx - 1
! integer, parameter :: lly1 = lly - 1
! integer, parameter :: llz1 = llz - 1

! integer, parameter :: lcx1 = lcx - 1
! integer, parameter :: lcy1 = lcy - 1
! integer, parameter :: lcz1 = lcz - 1

! integer, parameter :: nc1 = nc - 1
! integer, parameter :: ncp1 = ncp - 1

  ! Mathematical and physical constants
  double precision, parameter :: kb = 8.61734315d-05  ! eV/K
  double precision, parameter :: mub = 5.788381755d-5 ! eV/T
  double precision, parameter :: gs = 2.002319

END MODULE CONSTANTS
