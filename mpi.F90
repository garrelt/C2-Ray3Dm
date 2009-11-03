!>
!! \Brief This module contains data and routines for MPI parallelization
!!
!! Module for C2Ray / Capreole (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-06-01
!!
!! \b Version: True MPI (no dummy). Also reports on OpenMP parallelization.
!! Log files for nodes 1 and higher are called 'log.n', so node 0 it is
!! 'C2Ray.log'.
!! This module is also accepted by the F compiler (Dec 9, 2003)\n

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

  ! rank 0 has a log file called C2Ray.log associated with it. If
  ! the log unit is equal to 6, no file is opened and all log output
  ! is sent to standard output.
  ! If the code is compiled with the option -DMPILOG the other processors get
  ! their own log files, called log.1, log.2, etc.
  ! All these files are opened in mpi_setup and closed in mpi_end.

  ! This is the system module:
  !!include '/beosoft/mpich/include/mpif.h'        ! necessary for MPI
  !use mpi
  use file_admin, only: logf, results_dir

#ifdef XLF
  USE XLFUTILITY, only: hostnm => hostnm_ , flush => flush_
#endif

#ifdef IFORT
  USE IFPORT, only: hostnm
#ifdef _OPENMP
  USE OMP_LIB, only: omp_get_num_threads, omp_get_thread_num
#endif
#endif

  implicit none

  include 'mpif.h'

  integer,parameter,public :: NPDIM=3 ! dimension of problem

  integer,public :: rank            ! rank of the processor
  integer,public :: npr             ! number of processors
  integer,public :: nthreads=1      ! number of threads (per processor)
  integer,public :: MPI_COMM_NEW    ! the (new) communicator
  integer,public,dimension(MPI_STATUS_SIZE) :: mympi_status ! status array

  logical,parameter :: reorder=.false. !< reorder the mpi structure (for hydro)
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

    character(len=512) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: tn
    character(len=100) :: hostname

    call mpi_basic ()

    if (rank == 0) then
       if (logf /= 6) then
          filename=trim(adjustl(trim(adjustl(results_dir))//"C2Ray.log"))
          open(unit=logf,file=filename,status="unknown",action="write",&
               position="append")
       endif
       write(unit=logf,fmt="(A)") "Log file for C2-Ray run"
       write(unit=logf,fmt=*) " Number of MPI ranks used: ",npr
    endif

    ! Find number of OpenMP threads (needed to establish OpenMP character
    ! of run (see evolve)
    !$omp parallel default(shared)
#ifdef _OPENMP
    nthreads=omp_get_num_threads()
#endif
    !$omp end parallel

    ! Report OpenMP usage
#ifdef _OPENMP
    if (rank == 0) write(logf,*) " Running in OpenMP mode"
#endif

#ifdef MPILOG
    ! Open processor dependent log file
    write(unit=number,fmt="(I4)") rank
    filename=trim(adjustl("log."//trim(adjustl(number))))
    if (rank /= 0) open(unit=logf,file=filename,status="unknown",action="write")

    write(unit=logf,fmt=*) "Log file for rank ",rank," of a total of ",npr

    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror=hostnm(hostname)
    if (ierror == 0) then
       write(logf,*) "Running on processor named ",trim(adjustl(hostname))
    else 
       write(logf,*) "Error establishing identity of processor."
    endif

    ! Report number of OpenMP threads
#ifdef _OPENMP
    write(logf,*) ' Number of OpenMP threads is ',nthreads
#endif

    ! Let OpenMP threads report
    !$omp parallel default(private)
#ifdef _OPENMP
    tn=omp_get_thread_num()+1
    write(logf,*) 'Thread number ',tn,' reporting'
#endif
    !$omp end parallel

    write(logf,*) "almost end of mpi setup"
    flush(logf)
#endif

    if (reorder) then
       call mpi_topology ()
    else
       MPI_COMM_NEW=MPI_COMM_WORLD
    endif

#ifdef MPILOG
    write(logf,*) "end of mpi setup"
    flush(logf)
#endif

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
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Dims_create"
    write(logf,*) "MPI_Dims_create"
    flush(logf)
#endif


    periods(:)=.FALSE.      ! non-periodic boundaries

    reorder=.TRUE.
    ! makes MPI_COMM_NEW
    ! Warning: openmpi + gfortran seems to have problems here
    ! GM 091102
    call MPI_Cart_create(MPI_COMM_WORLD,NPDIM,dims,periods,reorder, &
         MPI_COMM_NEW,ierror)
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Cart_create"
    write(logf,*) "MPI_Cart_create"
    flush(logf)
#endif
    ! makes grid_struct               
    call MPI_Cart_get(MPI_COMM_NEW,NPDIM,dims, & ! makes grid_struct
         periods,grid_struct,ierror)
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Cart_get"
    write(logf,*) "MPI_Cart_get"
    flush(logf)
#endif
      
    ! Find the neighbours.
    ! My neighbors are now +/- 1 with my rank. Handle the case of the 
    ! boundaries by using MPI_PROC_NULL.
    call fnd3dnbrs ()

  end subroutine mpi_topology

  !----------------------------------------------------------------------------

  subroutine mpi_end

    integer :: ierror=0
    logical :: openlog

    ! Find out if log file is open
    inquire(unit=logf,opened=openlog)

    ! Close log file
    if (openlog) close(logf)

    ! Close MPI
    call MPI_FINALIZE(ierror)

  end subroutine mpi_end

  !----------------------------------------------------------------------------

  subroutine fnd3dnbrs
    
    ! This routine determines the neighbours in a 3-d decomposition of
    ! the grid. This assumes that MPI_Cart_create has already been called 

    integer :: ierror=0

    call MPI_Cart_shift( MPI_COMM_NEW, 0,  1, nbrleft,  nbrright, ierror )
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Cart_shift",0
    write(logf,*) "MPI_Cart_shift",0
    flush(logf)
#endif
    call MPI_Cart_shift( MPI_COMM_NEW, 1,  1, nbrdown,  nbrup,    ierror )
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Cart_shift",1
    write(logf,*) "MPI_Cart_shift",1
    flush(logf)
#endif
    call MPI_Cart_shift( MPI_COMM_NEW, 2,  1, nbrbelow, nbrabove, ierror )
#ifdef MPILOG
    if (ierror /= 0 ) write(logf,*) "error in MPI_Cart_shift",2
    write(logf,*) "MPI_Cart_shift",2
    flush(logf)
#endif

  end subroutine fnd3dnbrs

end module my_mpi
