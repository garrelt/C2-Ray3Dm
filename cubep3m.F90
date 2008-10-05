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
  use astroconstants, only: Mpc, M_SOLAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, h

  implicit none

  character(len=10),parameter :: nbody_type="cubep3m"

  real(kind=dp),parameter :: boxsize=64.0  ! Box size in Mpc/h comoving
  integer,parameter,private :: n_box=3456   ! cells/side (in N-body,fine grid)
  !real(kind=dp),parameter :: boxsize=114.0  ! Box size in Mpc/h comoving
  !integer,parameter,private :: n_box=6144   ! cells/side (in N-body,fine grid)

  character(len=180),parameter,private :: dir_dens_path = "../" 
  character(len=180),parameter,private :: dir_dens_name= "coarser_densities/"
  character(len=180),parameter,private :: dir_src_path = "./" 
  character(len=180),parameter,private :: dir_src_name= "sources/"

  ! properties of the box:
  ! M_box      - mass in box
  ! M_particle - mass per particle
  ! M_grid - mean mass per pmfast cell
  real(kind=dp),private :: M_box,M_particle
  real(kind=dp),public :: M_grid
  
  ! redshift sequence information
  integer, public :: NumZred               ! number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array ! array of redshifts
  integer,dimension(:),allocatable,public :: snap ! array of snapshot numbers
                                                  ! (for compatibility)
  character(len=8),public :: id_str       ! resolution dependent string

  character(len=180),public :: dir_dens, dir_src

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini ()
    
    character(len=180) :: redshift_file ! name of file with list of redshifts
    real(kind=dp) :: rho_crit_0_MpcM_sun,rho_matter
    integer :: nz ! loop counter
    character(len=20) :: dataroot="DEISA_DATA"
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
          dir_dens=value(1:len)//trim(adjustl(dir_dens_path)) &
               //trim(adjustl(dir_dens_name))
          dir_src=value(1:len)//trim(adjustl(dir_src_path)) &
               //trim(adjustl(dir_src_name))
       else
          dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
          dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))
       endif
    elseif (status == 1) then
       ! Assume that the whole path is set in the parameter
       dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
       dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))
    elseif (status == -1) then
       ! Warning
       write(logf,*) "Data file system name is truncated"
    endif
#else
    dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
    dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))
#endif
       
    ! Ask for redshift file
    if (rank == 0) then
       if (.not.file_input) write(*,"(A,$)") "File with redshifts: "
       read(stdinput,*) redshift_file
       
       ! Open and read redshift file
       open(unit=60,file=redshift_file,form="formatted",status="old")
       read(unit=60,fmt=*) NumZred
       allocate(zred_array(NumZred))
       do nz=1,NumZred
          read(unit=60,fmt=*) zred_array(nz)
       enddo
       close(20)
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif

    ! Parameters of simulations boxes

    rho_matter = rho_crit_0*Omega0             ! mean matter density (g/cm^3)
    M_box      = rho_matter*(boxsize*Mpc/h)**3 ! mass in box (g, not M0) 
    M_particle = 8.0*M_box/(real(n_box)**3)    ! mass per particle (g, not M0)
    M_grid = M_particle/8.                     ! mass in grid cell (g)

    ! Set identifying string (resolution-dependent)
    ! Construct the file name
    select case (mesh(1))
       case(216,256)	
          id_str="coarsest"
       case(432,384) 
          id_str="coarser"
       case(864,512) 
          id_str="coarse"
       end select
    if (rank == 0) write(unit=logf,fmt=*) "Type of resolution: ",id_str

  end subroutine nbody_ini

end module nbody
