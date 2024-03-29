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
  use cgsconstants, only: ev2fr
  use astroconstants, only: YEAR
  !use sizes, only: mesh

  implicit none

  !> Which fraction of the cells can be left unconverged in order
  !! to improve performance (used in rad_evolve3d)
  real(kind=dp),parameter :: convergence_fraction=1.0e-4

  !> Set to true to let C2-Ray not change the temperature
  logical,parameter :: isothermal=.true.

  !> A really small number
  real(kind=dp),parameter :: epsilon=1e-14_dp

  !> Convergence criterion for per source calculation (evolve0d)
  real(kind=dp),parameter :: minimum_fractional_change=1.0e-3

  !> Convergence criterion for global calculation (evolve0d)
  real(kind=dp),parameter :: convergence2=1.0e-3

  !> Convergence criterion for neutral fraction (evolve4_periodic)
  real(kind=dp),parameter :: minimum_fraction_of_atoms=1.0e-8
  
  !> Logical that determines the use of grey opacities
  logical,parameter :: grey = .false. 
  !> Modification factor for the HI ionization cross section. Should
  !! only be used in combination with grey opacities!! Will in that
  !! case ensure that the average ionization cross sections will be
  !! similar between the grey and non-grey cases.
  real(kind=dp),parameter :: freq_factor=1.0_dp

  !> Size increase of subboxes around sources (evolve4_periodic)
  !! If 10, we do 10^3, 20^3, 30^3 cubes around a source until
  !! no photons escape or we reach the edge of the (possibly periodic)
  !! grid
  integer,parameter :: subboxsize=5

  !> The maximum number	of cells on EITHER side	of the	source for which
  !! ray tracing is done. This	is a very crude	mean free path parameter
  !! which sets	up a photon barrier at exactly this distance. However,
  !! this barrier is cubic in shape (it's a subbox) so we prefer to set
  !! a spherical barrier, which must be done below (type_LLS=3).
  integer,parameter :: max_subbox=1000

  !> Add photon losses back into volume or not
  logical,parameter :: add_photon_losses=.false.

  !> Fraction of remaining photons below we stop ray-tracing
  real(kind=dp) :: loss_fraction=1e-2_dp
  
  !> Subgrid clumping
  !! 1: constant clumping (with clumping_factor)
  !! 2: Globally Averaged Clumping Model GCM, (provide parameters file)
  !! 3: Deterministic Clumping Model DCM (Mao et al. 2019), (provide parameters file)
  !! 4: Stochastic Clumping Model SCM (Bianco et al. 2020), (provide parameters file)
  !! 5: Pre-computed grid of clumping, (provide files)
  integer,parameter :: type_of_clumping=1
  !> Clumping factor if constant
  real,parameter :: clumping_factor=1.0

  !> Include LLS?
  logical,parameter :: use_LLS=.true.
  !> Type of LLS approach
  !! 0: no LLS
  !! 1: homogeneous optical depth due to LLS
  !! 2: position dependent optical depth due to LLS (requires LLS input
  !!     files)
  !! 3: hard barrier at a distance R_max_cMpc (see below)
  integer,parameter :: type_of_LLS=1

  !> Model selection in case of type_of_LLS=1 and 2
  !! 1: standard Worseck et al. (2014) mean free path fit
  !! 2: Worseck et al. (2014) lower mean free path
  !! 3: Worseck et al. (2014) higher mean free path
  !! 4: constant proper Mpc
  !! 5: constant comoving Mpc
  integer,parameter :: LLS_model=5

  !> Value of maximum comoving distance for photons from source in case of
  !! type_of_LLS=3
  real(kind=dp),parameter :: R_max_cMpc=10.0

  !> Should we stop when photon conservation violation is detected?
  logical,parameter :: stop_on_photon_violation = .false.

  !> Cosmology (in C2Ray.F90 and mat_ini) and Cosmological cooling (in cosmology)
  logical,parameter :: cosmological=.true. 

  !> Thermal: minimum temperature
  real(kind=dp),parameter :: minitemp=1.0 ! minimum temperature
  !> Thermal: fraction of the cooling time step below which no iteration is done
  real(kind=dp),parameter :: relative_denergy=0.1
  !> Thermal: initial temperature (if used)
  real(kind=dp),parameter :: initial_temperature=1e4

  !> Source properties: Number of different sources
  integer,parameter :: Number_Sourcetypes=2
  !> Source properties: Photon per atom for different source types (high to low mass)
  real,dimension(Number_Sourcetypes),parameter :: phot_per_atom= (/ 10.0, 150.0 /)
  !real,dimension(Number_Sourcetypes),parameter :: phot_per_atom= (/ 10.0, 150.0 , 0.0 /)

  !> Source properties: zeta parameter for collapsed fraction growth source
  !! model (high to low mass)
  real,dimension(Number_Sourcetypes),parameter :: zeta= (/ 50.0, 0.0 /)

  !> Source properties: X-ray photons per baryon. Mesinger et al. (2012) use
  !! 0.02 as their nominal value. Note that this depends on your integration
  !! limits. Mesinger et al. use 300 eV as lowest energy.
  real,parameter :: xray_phot_per_atom = 0.02
  !> Source properties: Life time of sources (if set at compile time)
  real,parameter :: lifetime=10e6*YEAR
  !> Source properties: Smallest number of particles that makes a reliable halo
  real,parameter :: MinParticleContent=20. ! in solar masses
  !> Source properties: Upper limit for low mass sources (not used)
  real,parameter :: LowMassLimit=1e9 ! in solar masses
  !> Source properties: Lower limit neutral fraction for suppression criterion
  real,parameter :: StillNeutral=0.1 ! lower limit of neutral criterium
  !real,parameter :: StillNeutral=0.9 ! lower limit of neutral criterium
  !real,parameter :: StillNeutral=1.1 ! NEVER suppress
  !real,parameter :: StillNeutral=-0.1 ! ALWAYS suppress

end module c2ray_parameters
