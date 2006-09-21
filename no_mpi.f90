module my_mpi

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01
  ! This module is also accepted by the F compiler (Dec 9, 2003)
 
  ! This is a dummy module for systems where there is no MPI
  ! for F90.
  !
  !----------------------------------------------------------------------------
  

  !implicit none

  !private

  integer,parameter,public :: NPDIM=3 ! dimension of problem

  integer,public :: rank              ! rank of the processor
  integer,public :: npr               ! number of processors
  integer,public :: MPI_COMM_NEW      ! the (new) communicator
  integer,dimension(NPDIM),public :: dims ! number of processors in 
                                             !  each dimension
  integer,dimension(NPDIM),public :: grid_struct ! coordinates of 
                                               !the processors in the grid
  
  integer,public ::  nbrleft,nbrright ! left and right neighbours
  integer,public ::  nbrdown,nbrup    ! up and down neighbours 
  integer,public ::  nbrabove,nbrbelow        ! above and below neighbours 

  integer,parameter,public :: MPI_PROC_NULL=-1

contains

  !----------------------------------------------------------------------------

  subroutine mpi_setup ( )

    character(len=10) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: hostnm
    character(len=100) :: hostname

    call mpi_basic ()

    ! Open processor dependent log file
    write(unit=number,fmt="(I4)") rank
    filename=trim(adjustl("log."//trim(adjustl(number))))
    open(unit=30,file=filename,status="unknown",action="write")

    write(unit=30,fmt=*) "Log file for rank ",rank
    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror=hostnm(hostname)
    write(unit=30,fmt=*) "The Processor is ",hostname
    call flush(30)

    call mpi_topology ()

  end subroutine mpi_setup

  !----------------------------------------------------------------------------

  subroutine mpi_basic ( )

    rank=0 ! Find processor rank

    ! Find total number of processors (npr)
    npr=1

  end subroutine mpi_basic

  !----------------------------------------------------------------------------

  subroutine mpi_topology ( )


    ! Make a new topology
    dims=1

    ! makes MPI_COMM_NEW    
    MPI_COMM_NEW=0

    ! makes grid_struct               
    grid_struct=0
      
    ! Find the neighbours.
    ! My neighbors are now +/- 1 with my rank. Handle the case of the 
    ! boundaries by using MPI_PROC_NULL.
    call fnd2dnbrs ( )

  end subroutine mpi_topology

  !----------------------------------------------------------------------------

  subroutine mpi_end ( )

    ! Close log file
    close(30)

  end subroutine mpi_end

  !----------------------------------------------------------------------------

  subroutine fnd2dnbrs ( )
    
    ! This routine determines the neighbours in a 2-d decomposition of
    ! the grid. This assumes that MPI_Cart_create has already been called 

    ! Single processor version

    nbrleft=MPI_PROC_NULL
    nbrright=MPI_PROC_NULL
    nbrup=MPI_PROC_NULL
    nbrdown=MPI_PROC_NULL
    nbrabove=MPI_PROC_NULL
    nbrbelow=MPI_PROC_NULL
    
  end subroutine fnd2dnbrs

end module my_mpi
