module cosmology_parameters

  ! This module defines essential cosmological parameters

  ! Version: WMAP3

  use precision, only: dp
  use mathconstants, only: pi
  use cgsconstants, only: G_grav
  use astroconstants, only: Mpc

  implicit none

  ! WMAP3

  real(kind=dp),parameter :: h=0.73 ! Hubble constant (in 100 km/s/Mpc)
  !real(kind=dp),parameter :: Omega0=0.238 ! Total matter density (in critical density)
  real(kind=dp),parameter :: Omega0=0.24 ! Total matter density (in critical density)
  !real(kind=dp),parameter :: Omega_B=0.0418 ! Baryon density (in critical density)
  real(kind=dp),parameter :: Omega_B=0.0223/(h*h) ! Baryon density (in critical density)
  real(kind=dp),parameter :: cmbtemp=2.726 ! CMB temperature

  ! Derived parameters
  real(kind=dp),parameter :: H0=h*100.0*1e5/Mpc ! Hubble constant (cgs)
  real(kind=dp),parameter :: rho_crit_0=3.0*H0*H0/(8.0*pi*G_grav) ! critical density (cgs)

end module cosmology_parameters
