!>
!! \brief This module contains data and routines for handling the data from
!! the Nbody simulations which provide the basis for C2Ray simulations
!! of Reionization.
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 22-Jan-2016 (22-May-2008 (previous version was not dated)
!!
!! \b Version: PMFAST simulations

module nbody

  ! This file contains routine having to do with the PMFAST N-body
  ! simulation.

  ! The routine in here is
  ! - nbody_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.
  ! It calculates the mass (in particles) of the simulation box
  ! It sets an identifier (id_str) for the resolution (used when
  ! reading in source list files and density fields.

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, file_input
  use cgsconstants, only: m_p
  use astroconstants, only: Mpc, M_SOLAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, Omega_B, h, H0
  use abundances, only: mu

  implicit none

  character(len=10),parameter :: nbody_type="pmfast" !< ID of Nbody type

  ! Boxsize: change here
  real(kind=dp),parameter :: boxsize=35.0  !< Box size in Mpc/h comoving
  integer,parameter :: n_box=3248  !< cells/side (in N-body,fine grid)

  !real(kind=dp),parameter :: boxsize=100.0  !< Box size in Mpc/h comoving
  !integer,parameter,private :: n_box=3248    !< cells/side (in N-body)

  !> Path to directory containing directory with density files:
  character(len=180),parameter,private :: dir_dens_path = "../" 
  !> Name of directory with density files
  character(len=180),parameter,private :: dir_dens_name= "coarser_densities/"

  !> Path to directory containing directory with source files:
  character(len=180),parameter,private :: dir_src_path = "./" 
  !> Name of directory with source files
  character(len=180),parameter,private :: dir_src_name= "sources/"

   !> Path to directory containing directory with clumping files:
  character(len=*),parameter,private :: dir_clump_path = "../" 
  !> Name of directory with files used for clumping
  !character(len=180),parameter,private :: dir_clump_name= "coarser_densities/"
  character(len=*),parameter,private :: dir_clump_name= "coarser_densities/halos_included/"

  !> Path to directory containing directory with LLS files:
  character(len=*),parameter,private :: dir_LLS_path = "../" 
  !> Name of directory with files used for LLS
  character(len=*),parameter,private :: dir_LLS_name= "halos/"

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
  logical,parameter :: densityheader=.false.
  !> clumping file with header?
  logical,parameter :: clumpingheader=.false.
  !> LLS file with header?
  logical,parameter :: LLSheader=.true.

  !> unit of density in density file
  !! can be "grid", "particle", "M0Mpc3"
  character(len=20),parameter :: density_unit="grid"

  ! Parameters of simulations boxes
  ! properties of the box:
  ! M_grid - mean mass per pmfast cell
  real(kind=dp),parameter,public :: M_box=rho_crit_0*Omega0*(boxsize*Mpc/h)**3 !< mass in box
  real(kind=dp),parameter,public :: M_grid=M_box/(real(n_box)**3) !< mean mass per grid cell
  real(kind=dp),parameter,public :: M_particle=8.0*M_grid !< mass per particle

  !> Conversion factor for comoving gas (number) density (cm^-3)
  real(kind=dp),parameter,public :: density_convert_grid=rho_crit_0*Omega_B/(mu*m_p)*(real(mesh(1))/real(n_box))**3
  !> Conversion factor for comoving gas (number) density (cm^-3)
  real(kind=dp),parameter,public :: density_convert_particle=8.0*density_convert_grid
  !> Conversion factor for (comoving) cubep3m lenght scales
  real(kind=dp),parameter,public :: lscale=boxsize*Mpc/h/n_box
  !> Conversion factor for cubep3m time scale (divide by (1+z)^2 to get proper
  !! converison factor for time)
  real(kind=dp),parameter,public :: tscale= 2.d0/(3.d0*sqrt(Omega0)*H0)

  ! redshift sequence information
  integer, public :: NumZred               !< number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array !< array of redshifts
  integer,dimension(:),allocatable,public :: snap !< array of snapshot numbers (for compatibility)
  character(len=8),public :: id_str="unknown" !< resolution dependent string

  character(len=180),public :: dir_dens !< Path to directory with density files
  character(len=180),public :: dir_src !< Path to directory with source files
  character(len=480),public :: dir_clump !< Path to directory with clump files
  character(len=480),public :: dir_LLS !< Path to directory with LLS files

  integer,parameter,public :: tot_nfiles=7 !< number of files at 812^3

#ifdef MPI
  integer,private :: mympierror !< MPI error flag variable
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini (ierror)
    
    integer,intent(out) :: ierror

    ! Set error flag to zero
    ierror=0

    ! Set the base directory names
    call set_directory_names ()

    ! In some cases a special file system is used, and its name is
    ! found from an environment variable. This needs to be added
    ! to the directory names. This behaviour is triggered by a preprocessor
    ! flag, so -DDEISA
#ifdef DEISA
    call set_directory_prefix(ierror)
#endif

    ! Read in the list of redshifts
    call set_list_of_redshifts(ierror)

    ! Determine resolution string (id_str)
    call set_resolution_string (ierror)
    
  end subroutine nbody_ini

  !---------------------------------------------------------------------------

  subroutine set_directory_names
    
    dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
    dir_clump=trim(adjustl(dir_clump_path))//trim(adjustl(dir_clump_name))
    dir_LLS=trim(adjustl(dir_LLS_path))//trim(adjustl(dir_LLS_name))
    dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))
    
  end subroutine set_directory_names

  !---------------------------------------------------------------------------

  subroutine set_directory_prefix(ierror)

    integer,intent(inout) :: ierror
    
    character(len=20) :: dataroot="DEISA_DATA"
    character(len=256) :: value
    integer :: len, status

    ! In some cases a special file system is used, and its name is
    ! found from an environment variable. This needs to be added
    ! to the directory names.
    call get_environment_variable (dataroot, value, len, status, .true.)
    if (status == 0) then
       ! The directory with density files is located in the dataroot
       ! plus the dir_dens_path parameter
       if (len > 0) then
          dir_dens=value(1:len)//trim(adjustl(dir_dens))
          dir_clump=value(1:len)//trim(adjustl(dir_clump))
          dir_LLS=value(1:len)//trim(adjustl(dir_LLS))
          dir_src=value(1:len)//trim(adjustl(dir_src))
       endif
    elseif (status == -1) then
       ! Warning
       write(logf,*) "Data file system name is truncated"
       ierror=5
    endif

  end subroutine set_directory_prefix

  !---------------------------------------------------------------------------

  subroutine set_list_of_redshifts(ierror)

    integer,intent(inout) :: ierror
    character(len=180) :: redshift_file ! name of file with list of redshifts
    integer :: nz ! loop counter
    integer :: open_status

    ! Ask for redshift file
    if (rank == 0) then
       if (.not.file_input) write(*,"(A,$)") "File with redshifts: "
       read(stdinput,*) redshift_file
       
       ! Open and read redshift file
       open(unit=60,file=redshift_file,form="formatted",status="old", &
            iostat=open_status)
       if (open_status == 0) then
          read(unit=60,fmt=*) NumZred
          allocate(zred_array(NumZred))
          do nz=1,NumZred
             read(unit=60,fmt=*) zred_array(nz)
          enddo
          close(60)
       else
          ierror=ierror+open_status
       endif
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif

  end subroutine set_list_of_redshifts

  !---------------------------------------------------------------------------

  subroutine set_resolution_string (ierror)

    integer,intent(inout) :: ierror

    ! Set identifying string (resolution-dependent)
    ! Construct the file name
    select case (mesh(1))
       case(203) 
          id_str="coarsest"
       case(406) 
          id_str="coarser"
       case(812) 
          id_str="coarse"
       end select
    if (rank == 0) write(unit=logf,fmt=*) "Type of resolution: ",id_str

    if (id_str == "unknown" ) ierror=ierror+1

  end subroutine set_resolution_string

end module nbody
