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

  character(len=180),public :: dir_dens !< Path to directory with density files
  character(len=180),public :: dir_src="./" !< Path to directory with source files

#ifdef MPI
  integer,private :: mympierror !< MPI error flag variable
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini (ierror)
    
    integer,intent(out) :: ierror

    real(kind=dp) :: t0,timestep
    integer :: nz ! loop counter
    character(len=256) :: value
    integer :: len, status

    ! Set error flag to zero
    ierror=0

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

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif

  end subroutine nbody_ini

end module nbody
