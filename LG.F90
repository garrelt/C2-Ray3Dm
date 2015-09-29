!>
!! \brief This module contains data and routines for handling the data from
!! the Nbody simulations which provide the basis for C2Ray simulations
!! of Reionization.
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 25-Jun-2013 (09-Dec-2009 (22-May-2008, previous versions were not dated)
!!
!! \b Version: LG simulations

module nbody

  ! This file contains routine having to do with the LG constrained IC's 
  ! Gadget N-body simulation.

  ! The routine in here is
  ! - nbody_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.

  ! Authors: Ilian Iliev, Garrelt Mellema

  ! Date: 22-May-2008 (previous version was not dated)

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, file_input
  use astroconstants, only: Mpc, M_SOLAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, h

  implicit none

  character(len=10),parameter :: nbody_type="LG"

  real(kind=dp),parameter :: boxsize=64.0  ! Box size in Mpc/h comoving

  !> Path to directory containing directory with density files:
  character(len=*),parameter,private :: dir_dens_path = "../" 
  !> Name of directory with density files
  character(len=*),parameter,private :: dir_dens_name= "coarser_densities/"

  !> Path to directory containing directory with clumping files:
  character(len=*),parameter,private :: dir_clump_path = "../" 
  !> Name of directory with files used for clumping
  !character(len=180),parameter,private :: dir_clump_name= "coarser_densities/"
  character(len=*),parameter,private :: dir_clump_name= "coarser_densities/halos_included/"

  !> Path to directory containing directory with LLS files:
  character(len=*),parameter,private :: dir_LLS_path = "../" 
  !> Name of directory with files used for LLS
  character(len=*),parameter,private :: dir_LLS_name= "halos/"

  !> Path to directory containing directory with source files:
  character(len=*),parameter,private :: dir_src_path = "./" 
  !> Name of directory with source files
  character(len=*),parameter,private :: dir_src_name= "sources/"

#ifdef IFORT
  ! ifort standard for "unformatted"
  character(len=*),parameter :: densityformat="unformatted"
  character(len=*),parameter :: densityaccess="sequential"
#else
  ! Fortran2003 standard for "binary"
  character(len=*),parameter :: densityformat="unformatted"
  character(len=*),parameter :: densityaccess="sequential"
#endif
  !> Format of clumping file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: clumpingformat="unformatted"
  character(len=*),parameter :: clumpingaccess="sequential"
#else
  character(len=15),parameter :: clumpingformat="unformatted"
  character(len=*),parameter :: clumpingaccess="sequential"
#endif
  !> Format of LLS file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: LLSformat="unformatted"
  character(len=*),parameter :: LLSaccess="sequential"
#else
  character(len=15),parameter :: LLSformat="unformatted"
  character(len=*),parameter :: LLSaccess="sequential"
#endif
  !> density file with header?
  logical,parameter :: densityheader=.true.
  !> clumping file with header?
  logical,parameter :: clumpingheader=.true.
  !> LLS file with header?
  logical,parameter :: LLSheader=.true.
  !> unit of density in density file
  !! can be "grid", "particle", "M0Mpc3"
  character(len=20),parameter :: density_unit="M0Mpc3"

  ! redshift sequence information
  integer, public :: NumZred               !< number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array !< array of redshifts 
  integer,dimension(:),allocatable,public :: snap !< array of snapshot numbers (for compatibility)
  character(len=8),public :: id_str       !< resolution dependent string

  character(len=180),public :: dir_dens !< Path to directory with density files
  character(len=480),public :: dir_clump !< Path to directory with clump files
  character(len=480),public :: dir_LLS !< Path to directory with LLS files
  character(len=180),public :: dir_src !< Path to directory with source files

#ifdef MPI
  integer,private :: mympierror !< MPI error flag variable
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini (ierror)
    
    integer,intent(out) :: ierror
    character(len=180) :: redshift_file ! name of file with list of redshifts
    integer :: nz ! loop counter
    character(len=20) :: dataroot="DEISA_DATA"
    character(len=256) :: value
    integer :: len, status
    real :: a ! scale factor?

    ! Set error flag to zero
    ierror=0

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
       allocate(snap(NumZred))
       do nz=1,NumZred
          read(unit=60,fmt=*) snap(nz), zred_array(nz)
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

    ! Set identifying string (resolution-dependent)
    ! Construct the file name
    select case (mesh(1))
       case(128) 
          id_str="coarsest"
       case(256) 
          id_str="coarser"
       case(512) 
          id_str="original"
       end select
    if (rank == 0) write(unit=logf,fmt=*) "Type of resolution: ",id_str

  end subroutine nbody_ini

end module nbody
