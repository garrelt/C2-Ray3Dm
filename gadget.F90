module nbody

  ! This file contains routine having to do with the GADGET N-body
  ! simulation.
  ! The routine in here is
  ! - nbody_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.

  ! Authors: Ilian Iliev, Garrelt Mellema

  ! Date: 22-May-2008 (previous version was not dated)

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, file_input
  use my_mpi

  implicit none

  character(len=10),parameter :: nbody_type="gadget"
  
  ! Boxsize: change here
  real(kind=dp),parameter :: boxsize=12.5  ! Box size in Mpc/h comoving

  ! Paths to directories with density and source files: change here
  character(len=180),parameter,private :: dir_dens_path = "../" 
  character(len=180),parameter,private :: dir_dens_name= "coarser_densities/"
  character(len=180),parameter,private :: dir_src_path = "../" 
  character(len=180),parameter,private :: dir_src_name= "sources/"

  ! redshift sequence information
  integer, public :: NumZred               ! number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array ! array of redshifts
  integer,dimension(:),allocatable,public :: snap ! array of snapshot numbers
                                                  ! (for compatibility)

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini ()
    
    character(len=180) :: redshift_file ! name of file with list of redshifts
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
       NumZred=1
       allocate(zred_array(NumZred))
       if (.not.file_input) write(*,'(A,$)') 'Initial redshift: '
       read(stdinput,*) zred_array(1)
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
