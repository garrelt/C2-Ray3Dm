!>
!! \brief This module contains radiation-related physical constants
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: cgs units
!!

module cgsphotoconstants

  ! A collection of physical constants and conversion factors for 
  ! photo-ionization calculations
  ! Units: cgs
  
  use precision, only: dp
  use cgsconstants

  !	!> Helium ionization potentials (eV)
  !	real(kind=dp), dimension(0:1),parameter :: ethe=(/24.587,54.416/)
  !> HI cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HI_at_ion_freq=6.346e-18
  !> HeI cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HeI_at_ion_freq=7.430e-18
  !> HeII cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HeII_at_ion_freq=1.589e-18
  !> HI ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HI=ev2fr*eth0
  !> HeI ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HeI=ev2fr*ethe(0)
  !> HeII ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HeII=ev2fr*ethe(1)

  ! HI cross-section at HeI ionization threshold
  real(kind=dp) :: sigma_H_heth = 1.238e-18_dp 

  ! HI cross-section at HeII Lya
  real(kind=dp) :: sigma_H_heLya = 9.907e-22_dp

  ! HeI cross-section at HeII Lya
  real(kind=dp) :: sigma_He_heLya = 1.301e-20_dp

  ! HeI cross-section at HeII ionization threshold
  real(kind=dp) :: sigma_He_he2 = 1.690780687052975e-18_dp

  ! HI cross-section at HeII ionization threshold
  real(kind=dp) :: sigma_H_he2 =1.230695924714239e-19_dp

end module cgsphotoconstants




