module pmfast

  ! This file contains routine having to do with the PMFAST N-body
  ! simulation.
  ! The routine in here is
  ! - pmfast_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.
  ! It calculates the mass (in particles) of the simulation box
  ! It sets an identifier (id_str) for the resolution (used when
  ! reading in source list files and density fields.

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput
  use astroconstants, only: Mpc, M_SOLAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, h

  implicit none

  integer,parameter,private :: n_box=3248    ! cells/side (in N-body)
  real(kind=dp),parameter,private :: boxsize=100.0  ! Box size in Mpc/h comoving

  ! properties of the box:
  ! M_box      - mass in box
  ! M_particle - mass per particle
  ! M_grid - mean mass per pmfast cell
  real(kind=dp),private :: M_box,M_particle
  real(kind=dp),public :: M_grid
  
  ! redshift sequence information
  integer, public :: NumZred               ! number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array ! array of redshifts
  character(len=8),public :: id_str       ! resolution dependent string

  character(len=180),public :: dir_dens
  character(len=180),parameter,private :: dir_dens_path = &
       "/disk/sn-12/garrelt/Simulations/Reionization/100Mpc_WMAP3/203/"
  integer,parameter,public :: tot_nfiles=7 ! number of files at 812^3

#ifdef MPI
  integer,private :: ierror
#endif

contains

  subroutine pmfast_ini ()
    
    character(len=180) :: redshift_file ! name of file with list of redshifts
    real(kind=dp) :: rho_crit_0_MpcM_sun,rho_matter
    integer :: nz ! loop counter
    character(len=20) :: dataroot='DEISA_DATA'
    character(len=256) :: value
    integer :: len, status

    ! In some cases a special file system is used, and its name is
    ! found from an environment variable.
#ifdef DEISA
    call get_environment_variable (dataroot, value, len, status, .true.)
    if (status == 0) then
       ! The directory with density files is located in the dataroot
       ! plus the dir_dens_path parameter
       if (len > 0) then
          dir_dens=value(1:len)//dir_dens_path
       else
          dir_dens=dir_dens_path
       endif
    elseif (status == 1) then
       ! Assume that the whole path is set in the parameter
       dir_dens=dir_dens_path
    elseif (status == -1) then
       ! Warning
       write(30,*) 'Data file system name is truncated'
    endif
#else
    dir_dens=dir_dens_path
#endif
       
    ! Ask for redshift file
    if (rank == 0) then
       write(*,'(A,$)') 'File with redshifts: '
       read(stdinput,*) redshift_file
       
       ! Open and read redshift file
       open(unit=60,file=redshift_file,form='formatted',status='old')
       read(60,*) NumZred
       allocate(zred_array(NumZred))
       do nz=1,NumZred
          read(60,*) zred_array(nz)
       enddo
       close(20)
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
#endif

    ! Parameters of simulations boxes

    rho_matter = rho_crit_0*Omega0             ! mean matter density (g/cm^3)
    M_box      = rho_matter*(boxsize*Mpc/h)**3 ! mass in box (g, not M0) 
    M_particle = 8.0*M_box/(real(n_box)**3)    ! mass per particle (g, not M0)
    M_grid = M_particle/8.                     ! mass in grid cell (g)

    ! Set identifying string (resolution-dependent)
    ! Construct the file name
    if (mesh(1) == 203) id_str="coarsest"
    if (mesh(1) == 406) id_str="coarser"
    if (mesh(1) == 812) id_str="coarse"
    write(30,*) 'Type of resolution: ',id_str

    return
  end subroutine pmfast_ini

end module pmfast
