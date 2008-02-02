module my_mpi

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01
  ! This module is also accepted by the F compiler (Dec 9, 2003)
 
  ! This is a dummy module for systems where there is no MPI
  ! for F90.
  !
  !----------------------------------------------------------------------------

  use file_admin, only: log

#ifdef XLF
  USE XLFUTILITY, only: hostnm => hostnm_ , flush => flush_
#endif

#ifdef IFORT
  USE IFPORT, only: hostnm, flush
  !$ USE OMP_LIB, only: omp_get_num_threads, omp_get_thread_num
#endif

  implicit none

  integer,parameter,public :: NPDIM=3 ! dimension of problem

  integer,public :: rank            ! rank of the processor
  integer,public :: npr             ! number of processors
  integer,public :: nthreads        ! number of threads (per processor)
  integer,public :: MPI_COMM_NEW    ! the (new) communicator

  integer,dimension(NPDIM),public :: dims ! number of processors in 
                                             !  each dimension
  integer,dimension(NPDIM),public :: grid_struct ! coordinates of 
                                               !the processors in the grid
  
  integer,public ::  nbrleft,nbrright  ! left and right neighbours
  integer,public ::  nbrdown,nbrup     ! up and down neighbours 
  integer,public ::  nbrabove,nbrbelow ! above and below neighbours 

  integer,parameter,public :: MPI_PROC_NULL=-1

contains

  !----------------------------------------------------------------------------

  subroutine mpi_setup ( )

    character(len=10) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: tn
    character(len=100) :: hostname

    call mpi_basic ()

    ! Open processor dependent log file
    write(unit=number,fmt="(I4)") rank
    filename=trim(adjustl("log."//trim(adjustl(number))))
    open(unit=log,file=filename,status="unknown",action="write")

    write(unit=log,fmt=*) "Log file for rank ",rank

    nthreads=1
    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror=hostnm(hostname)
    if (ierror == 0) then
       write(log,*) "Running on processor named ",hostname
    else 
       write(log,*) "Error establishing identity of processor."
    endif

    ! Report number of OpenMP threads
    !$omp parallel default(shared)
    !$nthreads=omp_get_num_threads()
    !$omp end parallel
    !$write(log,*) ' Number of OpenMP threads is ',nthreads

    ! Let OpenMP threads report
    !$omp parallel default(shared)
    !$tn=omp_get_thread_num()+1
    !$write(log,*) 'Thread number ',tn,' reporting'
    !$omp end parallel

    call flush(log)

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
    close(log)

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
