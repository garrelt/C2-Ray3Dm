!>
!! \brief This module contains SED parameters specific for C2-Ray
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! Note: this file contains parameters for different versions
!!       of C2-Ray (dimensions, single/multiple sources), so 
!!       not all parameters will be used.

module sed_parameters

  ! This module collects parameters needed by C2-Ray

  use precision, only: dp
  use cgsconstants, only: ev2fr
  use cgsphotoconstants, only: ion_freq_HI, ion_freq_HeII
  use astroconstants, only: YEAR

  implicit none

  !> Stellar SED type
  !! Black body: 1
  !! Power law: 2
  integer,parameter :: stellar_SED_type=1
  
  !> Parameters for SED (BB)
  !> Effective temperature (K); if set to zero, the code will ask
  !! for SED parameters
  real(kind=dp),parameter :: bb_Teff=5.0e4
  !> Number of ionizing photons / second (reference)
  real(kind=dp),parameter :: bb_S_star=1e48_dp
  !> nominal minimum and maximum frequency for BB source
  real(kind=dp),parameter :: bb_MinFreq=ion_freq_HI
  real(kind=dp),parameter :: bb_MaxFreq=ion_freq_HeII * 10.00_dp

  !> Parameters for SED (PL)
  !> nominal power law index (for photon number)
  real(kind=dp),parameter :: pl_index=3.0_dp
  !> Number of ionizing photons / second (reference)
  real(kind=dp),parameter :: pl_S_star=1e48_dp
  !> nominal minimum and maximum frequency for power law source
  real(kind=dp),parameter :: pl_MinFreq=ion_freq_HI
  real(kind=dp),parameter :: pl_MaxFreq=ion_freq_HeII! * 100.00_dp

  !> nominal Eddington efficiency
  real(kind=dp),parameter :: EddLeff_nominal=1.0_dp
  !> nominal black hole mass for Eddington luminosity (M0)
  real(kind=dp),parameter :: mass_nominal=1.0e6_dp
  !> Eddington luminosity per mass_nominal solar mass (erg/s)
  real(kind=dp),parameter :: EddLum=1.38e38*mass_nominal

  !> Use an X-ray SED as well?
  !! This is not yet implemented so should be false.
  logical,parameter :: use_xray_SED=.false.
  !> Source properties: X-ray photons per baryon. Mesinger et al. (2012) use
  !! 0.02 as their nominal value. Note that this depends on your integration
  !! limits. Mesinger et al. use 300 eV as lowest energy.
  real,parameter :: xray_phot_per_atom = 0.02

end module sed_parameters
