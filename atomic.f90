!>
!! \brief This module contains atomic constants and parameter definitions
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! Author: Garrelt Mellema
!!
!! Date: 2003-06-01
!!
!! Version: cgs units
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<
module atomic

  use precision, only: dp

  implicit none

  private

  !> adiabatic index
  real(kind=dp),public,parameter :: gamma = 5.0_dp/3.0_dp
  !> adiabatic index - 1
  real(kind=dp),public,parameter :: gamma1 = gamma - 1.0_dp

end module atomic
