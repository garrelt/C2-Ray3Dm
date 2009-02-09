!>
!! \brief This module contains definitions of cosmological parameters
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 05-Oct-2008
!!
!! Version: WMAP3+

module cosmology_parameters

  ! This module defines essential cosmological parameters

  ! Version: WMAP3+

  use precision, only: dp
  use mathconstants, only: pi
  use cgsconstants, only: G_grav
  use astroconstants, only: Mpc

  implicit none

  ! WMAP3+; Note: Omega_lambda=1-Omega_0 and for this model sigma8=0.8, n_s=0.96.
  character(len=10),parameter :: cosmo_id="WMAP3+"

  real(kind=dp),parameter :: h=0.7 !< Hubble constant (in 100 km/s/Mpc)
  !real(kind=dp),parameter :: Omega0=0.238 !< Total matter density (in critical density)
  real(kind=dp),parameter :: Omega0=0.27 !< Total matter density (in critical density)
  !real(kind=dp),parameter :: Omega_B=0.0418 !< Baryon density (in critical density)
  real(kind=dp),parameter :: Omega_B=0.044 !< Baryon density (in critical density)
  real(kind=dp),parameter :: cmbtemp=2.726 !< CMB temperature

  ! Derived parameters
  real(kind=dp),parameter :: H0=h*100.0*1e5/Mpc !< Hubble constant (cgs)
  real(kind=dp),parameter :: rho_crit_0=3.0*H0*H0/(8.0*pi*G_grav) !< critical density (cgs)

end module cosmology_parameters
