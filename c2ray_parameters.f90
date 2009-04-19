!>
!! \brief This module contains parameters specific for C2-Ray
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! Note: this file contains parameters for different versions
!!       of C2-Ray (dimensions, single/multiple sources), so 
!!       not all parameters will be used.

module c2ray_parameters

  ! This module collects parameters needed by C2-Ray

  use precision, only: dp
  use astroconstants, only: YEAR

  !> Which fraction of the cells can be left unconverged in order
  !! to improve performance (used in rad_evolve3d)
  real(kind=dp),parameter :: convergence_fraction=1.0e-5

  !> Set to true to let C2-Ray not change the temperature
  logical,parameter :: isothermal=.true.

  !> A really small number
  real(kind=dp),parameter :: epsilon=1e-14_dp

  !> Convergence criterion for per source calculation (evolve0d)
  real(kind=dp),parameter :: convergence1=1.0e-3

  !> Convergence criterion for global calculation (evolve0d)
  real(kind=dp),parameter :: convergence2=1.0e-3

  !> Convergence criterion for neutral fraction (evolve4_periodic)
  real(kind=dp),parameter :: convergence_frac=1.0e-8

  !> Size increase of subboxes around sources (evolve4_periodic)
  !! If 10, we do 10^3, 20^3, 30^3 cubes around a source until
  !! no photons escape or we reach the edge of the (possibly periodic)
  !! grid
  integer,parameter :: subboxsize=10

  !> Parameters for nominal SED
  real(kind=dp),parameter :: teff_nominal=50000.0
  real(kind=dp),parameter :: s_star_nominal=1e48_dp
  !real(kind=dp),parameter :: s_star_nominal=1e50_dp
  
  !> Subgrid clumping\n
  !! 1: constant clumping (with clumping_factor)\n
  !! 2: 3.5Mpc PM, WMAP1 clumping\n
  !! 3: 3.5Mpc PM, WMAP3 clumping\n
  !! 4: 1 Mpc P3M\n
  !! 5: position dependent clumping
  integer,parameter :: type_of_clumping=5
  !> Clumping factor if constant
  real,parameter :: clumping_factor=1.0  

  !> Cosmological cooling
  logical,parameter :: cosmological=.true.

  !> Thermal: minimum temperature
  real(kind=dp),parameter :: minitemp=1.0 ! minimum temperature
  !> Thermal: fraction of the cooling time step below which no iteration is done
  real(kind=dp),parameter :: relative_denergy=0.1

  !> Source properties: Photon per atom for high mass sources
  real,parameter :: phot_per_atom1=10.0
  !> Source properties: Photon per atom for low mass sources
  real,parameter :: phot_per_atom2=150.0
  !> Source properties: Life time of sources (if set at compile time)
  real,parameter :: lifetime=20e6*YEAR
  !> Source properties: Smallest number of particles that makes a reliable halo
  real,parameter :: MinParticleContent=20. ! in solar masses
  !> Source properties: Upper limit for low mass sources (not used)
  real,parameter :: LowMassLimit=1e9 ! in solar masses
  !> Source properties: Lower limit neutral fraction for suppression criterion
  real,parameter :: StillNeutral=0.1 ! lower limit of neutral criterium

end module c2ray_parameters
