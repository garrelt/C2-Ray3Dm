module cosmology_parameters

  ! This module defines essential cosmological parameters

  ! Version: WMAP5

  use precision, only: dp
  use mathconstants, only: pi
  use cgsconstants, only: G_grav
  use astroconstants, only: Mpc

  implicit none

  ! WMAP5 LG; Note: Omega_lambda=1-Omega_0 and for this model sigma8=0.817, n_s=0.96.

  character(len=10),parameter :: cosmo_id="WMAP5"

  real(kind=dp),parameter :: h=0.7 ! Hubble constant (in 100 km/s/Mpc)

  real(kind=dp),parameter :: Omega0=0.279 ! Total matter density (in critical density)
  real(kind=dp),parameter :: Omega_B=0.046 ! Baryon density (in critical density)
  real(kind=dp),parameter :: cmbtemp=2.726 ! CMB temperature

  real(kind=dp),parameter :: sigma8=0.817 !< sigma8

  real(kind=dp),parameter :: n_s=0.96 !< slope of density power spectrum

  ! Derived parameters
  real(kind=dp),parameter :: H0=h*100.0*1e5/Mpc ! Hubble constant (cgs)
  real(kind=dp),parameter :: rho_crit_0=3.0*H0*H0/(8.0*pi*G_grav) ! critical density (cgs)

end module cosmology_parameters
