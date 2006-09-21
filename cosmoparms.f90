module cosmology_parameters

  use precision, only: dp
  use mathconstants, only: pi
  use cgsconstants, only: G_grav
  use astroconstants, only: Mpc
  implicit none

  logical,parameter :: cosmological=.true.

  ! WMAP3
  real(kind=dp),parameter :: h=0.73 ! Hubble constant (in 100 km/s/Mpc)
  real(kind=dp),parameter :: Omega0=0.24 ! Matter density (in critical density)
  real(kind=dp),parameter :: Omega_B=0.0223/(h*h) ! Baryon density (in critical density)

  real(kind=dp),parameter :: H0=h*100.0*1e5/Mpc ! Hubble constant (cgs)
  real(kind=dp),parameter :: rho_crit_0=3.0*H0*H0/(8.0*pi*G_grav) ! critical density (cgs)
  real,parameter :: cmbtemp=2.726

end module cosmology_parameters
