module my_mpi

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01

  ! This module contains the routines for using the MPI library
  ! mpi_setup :   set up MPI
  ! mpi_basic :   basic MPI initialization
  ! mpi_topology: domain decomposition
  ! mpi_end:      close down the MPI interface
  ! fnd3dnbrs:    find neighbours in 3D domain decomposition
  
  ! Each processor has its own log file, which is called log.0, log.1,
  ! etc. This is opened in mpi_setup, and is connected to unit 30.
  ! The file is close in mpi_end

  ! This is the system module:
  !!include '/beosoft/mpich/include/mpif.h'        ! necessary for MPI
  !use mpi
  use file_admin, only:log

#ifdef XLF
  USE XLFUTILITY, only: hostnm => hostnm_ , flush => flush_
#endif

#ifdef IFORT
  USE IFPORT, only: hostnm, flush
  !$ USE OMP_LIB, only: omp_get_num_threads, omp_get_thread_num
#endif

  implicit none

  include 'mpif.h'

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

contains

  !----------------------------------------------------------------------------

  subroutine mpi_setup ( )

    character(len=10) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: tn
    character(len=100) :: hostname

    call mpi_basic ()

#ifdef MPILOG
    ! Open processor dependent log file
    write(unit=number,fmt="(I4)") rank
    filename=trim(adjustl("log."//trim(adjustl(number))))
    open(unit=log,file=filename,status="unknown",action="write")

    write(unit=log,fmt=*) "Log file for rank ",rank," of a total of ",npr

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
#else
    if (rank == 0) then
       filename=trim(adjustl("C2Ray.log"//trim(adjustl(number))))
       open(unit=log,file=filename,status="unknown",action="write")
       write(unit=log,fmt="(A)") "Log file for C2-Ray run"
    endif
#endif

    call mpi_topology ()

  end subroutine mpi_setup

  !----------------------------------------------------------------------------

  subroutine mpi_basic

    integer :: ierror          ! control variable for MPI

    call MPI_INIT (ierror)  ! Initialize MPI

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror) ! Find processor rank

    ! Find total number of processors (npr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npr,ierror)

  end subroutine mpi_basic

  !----------------------------------------------------------------------------

  subroutine mpi_topology

    logical,dimension(NPDIM) :: periods ! for periodic grid
    logical                  :: reorder ! reorder the MPI_COMM_WORLD
    integer          :: ierror=0

    ! Make a new topology
    dims(:)=0

    call MPI_Dims_create(npr,NPDIM,dims,ierror)

    periods(:)=.FALSE.      ! non-periodic boundaries

    reorder=.TRUE.
    ! makes MPI_COMM_NEW    
    call MPI_Cart_create(MPI_COMM_WORLD,NPDIM,dims,periods,reorder, &
         MPI_COMM_NEW,ierror)
    ! makes grid_struct               
    call MPI_Cart_get(MPI_COMM_NEW,NPDIM,dims, & ! makes grid_struct
         periods,grid_struct,ierror)
      
    ! Find the neighbours.
    ! My neighbors are now +/- 1 with my rank. Handle the case of the 
    ! boundaries by using MPI_PROC_NULL.
    call fnd3dnbrs ()

  end subroutine mpi_topology

  !----------------------------------------------------------------------------

  subroutine mpi_end

    integer :: ierror=0

    ! Close log file
    close(log)

    ! Close MPI
    call MPI_FINALIZE(ierror)

  end subroutine mpi_end

  !----------------------------------------------------------------------------

  subroutine fnd3dnbrs
    
    ! This routine determines the neighbours in a 3-d decomposition of
    ! the grid. This assumes that MPI_Cart_create has already been called 

    integer :: ierror=0

    call MPI_Cart_shift( MPI_COMM_NEW, 0,  1, nbrleft,  nbrright, ierror )
    call MPI_Cart_shift( MPI_COMM_NEW, 1,  1, nbrdown,  nbrup,    ierror )
    call MPI_Cart_shift( MPI_COMM_NEW, 2,  1, nbrbelow, nbrabove, ierror )

  end subroutine fnd3dnbrs

end module my_mpi
