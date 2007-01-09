module c2ray_parameters

  ! This module collects parameters needed by C2-Ray

  use precision
  use astroconstants

  ! Which fraction of the cells can be left unconverged in order
  ! to improve performance (used in rad_evolve3d)
  real(kind=dp),parameter :: convergence_fraction=1.5e-5

  ! Set to true to let C2-Ray not change the temperature
  logical,parameter :: isothermal=.true.

  ! A really small number
  real(kind=dp),parameter :: epsilon=1e-40_dp

  ! Convergence criterion for per source calculation (evolve0d)
  real(kind=dp),parameter :: convergence1=1.0e-3

  ! Convergence criterion for global calculation (evolve0d)
  real(kind=dp),parameter :: convergence2=1.0e-2

  ! Parameters for nominal SED
  real(kind=dp),parameter :: teff_nominal=50000.0
  real(kind=dp),parameter :: s_star_nominal=1e50_dp
  
  ! Subgrid clumping
  ! 1: constant clumping (with clumping_factor)
  ! 2: 3.5Mpc PM, WMAP1 clumping
  ! 3: 3.5Mpc PM, WMAP3 clumping
  ! 4: 1 Mpc P3M
  integer,parameter :: type_of_clumping=1
  ! Clumping factor if constant
  real,parameter :: clumping_factor=1.0  

  ! Cosmological cooling
  logical,parameter :: cosmological=.true.

  ! Thermal
  ! minimum temperature
  real(kind=dp),parameter :: minitemp=1.0 ! minimum temperature
  ! fraction of the cooling time step below which no iteration is done
  real(kind=dp),parameter :: relative_denergy=0.1

  ! Source properties
  ! Photon per atom for high mass sources
  real,parameter :: phot_per_atom1=250.0
  ! Photon per atom for low mass sources
  real,parameter :: phot_per_atom2=250.0
  ! Life time of sources
  real,parameter :: lifetime=20e6*YEAR
  ! Upper limit for low mass sources
  real,parameter :: LowMassLimit=1e9 ! in solar masses
  ! Upper limit for suppression criterion
  real,parameter :: StillNeutral=0.9 ! lower limit of neutral criterium

end module c2ray_parameters
