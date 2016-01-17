!>
!! \brief This module contains data and routines for handling the data from
!! the Nbody simulations which provide the basis for C2Ray simulations
!! of Reionization.
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 22-May-2008 (previous version was not dated)
!!
!! \b Version: test simulations
!!
!! The test module has a number of parameters hard coded. These are
!!
!! Starting redshift: 9
!!
!! Number of redshift slices: 5 (separated by 10 million years).
!!
!! Box size: 100/h Mpc

module nbody

  ! This file contains routine having to do with the CubeP3M N-body
  ! simulation.
  ! The routine in here is
  ! - nbody_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.
  ! It calculates the mass (in particles) of the simulation box
  ! It sets an identifier (id_str) for the resolution (used when
  ! reading in source list files and density fields.

  ! Authors: Ilian Iliev, Garrelt Mellema

  ! Date: 22-May-2008 (previous version was not dated)

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, file_input
  use astroconstants, only: Mpc, M_SOLAR, YEAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, h, H0

  implicit none

  character(len=10),parameter :: nbody_type="test" !< ID of Nbody type

  real(kind=dp),parameter :: boxsize=100.0  !< Box size in Mpc/h comoving

  ! redshift sequence information
  integer, public :: NumZred               !< number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array !< array of redshifts 
  integer,dimension(:),allocatable,public :: snap !< array of snapshot numbers (for compatibility)
  character(len=8),public :: id_str       !< resolution dependent string

  character(len=480),public :: dir_dens !< Path to directory with density files
  character(len=480),public :: dir_clump !< Path to directory with clump files
  character(len=480),public :: dir_LLS !< Path to directory with LLS files
  character(len=480),public :: dir_src !< Path to directory with source files

  !> Path to directory containing directory with density files:
  character(len=*),parameter,private :: dir_dens_path = "" 
  !> Name of directory with density files
  character(len=180),parameter,private :: dir_dens_name= ""

  !> Path to directory containing directory with clumping files:
  character(len=*),parameter,private :: dir_clump_path = "" 
  !> Name of directory with files used for clumping
  character(len=*),parameter,private :: dir_clump_name= ""

  !> Path to directory containing directory with LLS files:
  character(len=*),parameter,private :: dir_LLS_path = "" 
  !> Name of directory with files used for LLS
  character(len=*),parameter,private :: dir_LLS_name= ""

  !> Path to directory containing directory with source files:
  character(len=*),parameter,private :: dir_src_path = "./" 
  !> Name of directory with source files
  character(len=*),parameter,private :: dir_src_name= ""

  !> Format of density file (unformatted or binary)
#ifdef IFORT
  ! ifort standard for "binary"
  character(len=*),parameter :: densityformat="binary"
  character(len=*),parameter :: densityaccess="sequential"
#else
  ! Fortran2003 standard for "binary"
  character(len=*),parameter :: densityformat="unformatted"
  character(len=*),parameter :: densityaccess="stream"
#endif
  !> Format of clumping file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: clumpingformat="binary"
  character(len=*),parameter :: clumpingaccess="sequential"
#else
  character(len=15),parameter :: clumpingformat="unformatted"
  character(len=*),parameter :: clumpingaccess="stream"
#endif
  !> Format of LLS file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: LLSformat="binary"
  character(len=*),parameter :: LLSaccess="sequential"
#else
  character(len=15),parameter :: LLSformat="unformatted"
  character(len=*),parameter :: LLSaccess="stream"
#endif
  !> density file with header?
  logical,parameter :: densityheader=.true.
  !> clumping file with header?
  logical,parameter :: clumpingheader=.true.
  !> LLS file with header?
  logical,parameter :: LLSheader=.true.

  !> unit of density in density file
  !! can be "grid", "particle", "M0Mpc3"
  character(len=*),parameter :: density_unit="none"

  ! Parameters of simulations boxes
  ! properties of the box:
  ! M_box      - mass in box
  ! M_particle - mass per particle
  ! M_grid - mean mass per pmfast cell
  real(kind=dp),parameter,public :: M_box=rho_crit_0*Omega0*(boxsize*Mpc/h)**3 !< mass in box
  real(kind=dp),parameter,public :: M_grid=0.0 !< mean mass per grid cell: not used
  real(kind=dp),parameter,public :: M_particle=0.0 !< mass per particle: not used

  !> Conversion factor for comoving gas (number) density (cm^-3)
  real(kind=dp),parameter,public :: density_convert_grid=1.0
  !> Conversion factor for comoving gas (number) density (cm^-3)
  real(kind=dp),parameter,public :: density_convert_particle=1.0
  !> Conversion factor for lenght scales: not used
  real(kind=dp),parameter,public :: lscale=1.0
  !> Conversion factor for time scale: not used
  real(kind=dp),parameter,public :: tscale= 1.0

#ifdef MPI
  integer,private :: mympierror !< MPI error flag variable
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini ()
    
    real(kind=dp) :: t0,timestep
    integer :: nz ! loop counter
    character(len=256) :: value
    integer :: len, status

    ! Set directories: not currently used
    dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
    dir_clump=trim(adjustl(dir_clump_path))//trim(adjustl(dir_clump_name))
    dir_LLS=trim(adjustl(dir_LLS_path))//trim(adjustl(dir_LLS_name))
    dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))

    ! Construct redshift sequence
    if (rank == 0) then

       ! Set the number of redshift slices
       NumZred=5
       allocate(zred_array(NumZred))

       ! Time step
       timestep=1e7*YEAR

       ! Starting redshift
       zred_array(1)=9.

       ! Cosmological time corresponding to (initial) redshift
       ! NOTE: Good only for high-z!!!
       t0 = 2.*(1.+zred_array(1))**(-1.5)/(3.*H0*sqrt(Omega0))
       do nz=2,NumZred
          zred_array(nz)=-1+(1.+zred_array(1))* &
               (t0/(t0+real(nz-1)*timestep))**(2./3.)
       enddo
    endif
       dir_clump=trim(adjustl(dir_clump_path))//trim(adjustl(dir_clump_name))
       dir_LLS=trim(adjustl(dir_LLS_path))//trim(adjustl(dir_LLS_name))
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif

    ! Set id_str for compatibility reasons
    id_str="test"

  end subroutine nbody_ini

end module nbody
