!>
!! \brief This module contains routines distributing the work over a number
!! of (MPI) slaves.
!! The module supplies one routine to the outside (do_grid), 
!! but internally has two
!! methods to distribute the work: one which distributes
!! the work statically over the available MPI processes (do_grid_static)
!! and one which works with a master-slave set up (do_grid_master_slave). 
!! In this latter set up rank 0 is the master and distributes the work over 
!! the slaves dynamically.
!! The decision of which one to use depends on the number of MPI processes.
!! Above min_numproc_master_slave (module parameter) the master-slave set up
!! is used.
!!
!! The actual work done is contained in the (external) routine do_source.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-09-05
!!
!! \b Version: 

module master_slave_processing

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,timefile,iterdump, results_dir, dump_dir
  use clocks, only: timestamp_wallclock
  use sourceprops, only: NumSrc, srcpos
  use evolve_source, only: do_source
  
  implicit none

  private

  public :: do_grid

  !> Minimum number of MPI processes for using the master-slave setup 
  integer, parameter ::  min_numproc_master_slave=10

  !> Report about the sending of every 1000th source if check_progress
  !! is 1000
  integer, parameter :: check_progress = 1000000

contains

  ! ===========================================================================

  !> Does the ray-tracing over the sources by distributing
  !! the sources evenly over the available MPI processes-
  subroutine do_grid (dt,niter)

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    ! Ray trace the whole grid for all sources.
    ! We can do this in two ways, depending on
    ! the number of processors. For many processors
    ! the master-slave setup should be more efficient.
    if (npr > min_numproc_master_slave) then
       call do_grid_master_slave (dt,niter)
    else
       call do_grid_static (dt,niter)
    endif
    
  end subroutine do_grid

  ! ===========================================================================

  !> Does the ray-tracing over the sources by distributing
  !! the sources evenly over the available MPI processes-
  subroutine do_grid_static (dt,niter)

    ! Does the ray-tracing over the sources by distributing
    ! the sources evenly over the available MPI processes-
    
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: ns1

    ! Source Loop - distributed for the MPI nodes
    do ns1=1+rank,NumSrc,npr
#ifdef MPILOG
       ! Report
       write(logf,*) 'Processor ',rank,' received: ',ns1
       write(logf,*) ' that is source ',ns1 !SrcSeries(ns1)
       write(logf,*) ' at:',srcpos(:,ns1)
       flush(logf)
#endif
       call do_source(dt,ns1,niter)
    enddo

  end subroutine do_grid_static

  ! ===========================================================================
  
  !> Ray tracing the entire grid for all the sources using the
  !! master-slave model for distributing the sources over the
  !! MPI processes.
  subroutine do_grid_master_slave (dt,niter)

    ! Ray tracing the entire grid for all the sources using the
    ! master-slave model for distributing the sources over the
    ! MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    if (rank == 0) then
       call do_grid_master ()
    else
       call do_grid_slave (dt,niter)
    endif

  end subroutine do_grid_master_slave

  ! ===========================================================================

  !> The master task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_master ()

    ! The master task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.

    integer :: ns1
    integer :: sources_done,whomax,who,answer
    ! counter for master-slave process
    integer,dimension(:),allocatable :: counts
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    ! Source Loop - Master Slave with rank=0 as Master
    sources_done = 0
          
    ns1 = 0
    
    ! Allocate counter for master-slave process
    if (.not.(allocated(counts))) allocate(counts(0:npr-1))

    ! send tasks to slaves 
    
    whomax = min(NumSrc,npr-1)
    do who=1,whomax
       if (ns1 <= NumSrc) then
          ns1=ns1+1
          call MPI_Send (ns1, 1, MPI_INTEGER, who, 1,  &
               MPI_COMM_NEW, mympierror)
       endif
    enddo
    
    do while (sources_done < NumSrc)
       
       ! wait for an answer from a slave. 
       
       call MPI_Recv (answer,     & ! address of receive buffer
            1,		   & ! number of items to receive
            MPI_INTEGER,	   & ! type of data
            MPI_ANY_SOURCE,  & ! can receive from any other
            1,		   & ! tag
            MPI_COMM_NEW,	   & ! communicator
            mympi_status,	   & ! status
            mympierror)
       
       who = mympi_status(MPI_SOURCE) ! find out who sent us the answer
       sources_done=sources_done+1 ! and the number of sources done
       
       ! put the slave on work again,
       ! but only if not all tasks have been sent.
       ! we use the value of num to detect this */
       if (ns1 < NumSrc) then
          ns1=ns1+1

          ! Report on the sending on nth source
          if (mod(ns1,check_progress) == 0) then
             write(logf,"(A,I12,A,I12,A,I6,A,I12)") &
                  "Sending source ",ns1," of ",NumSrc," to processor ",who, &
                  " sources done ", sources_done
             flush(logf)
          endif
          
          call MPI_Send (ns1, 1, MPI_INTEGER, &
               who,		&	
               1,		&	
               MPI_COMM_NEW, &
               mympierror)
       endif
    enddo
    
    ! Now master sends a message to the slaves to signify that they 
    ! should end the calculations. We use a special tag for that:
    
    do who = 1,npr-1
       call MPI_Send (0, 1, MPI_INTEGER, &
            who,			  &
            2,			  & ! tag 
            MPI_COMM_NEW,	          &
            mympierror)
       
       ! the slave will send to master the number of calculations
       ! that have been performed. 
       ! We put this number in the counts array.
       
       call MPI_Recv (counts(who), & ! address of receive buffer
            1,                & ! number of items to receive
            MPI_INTEGER,      & ! type of data 
            who,              & ! receive from process who 
            7,                & ! tag 
            MPI_COMM_NEW,     & ! communicator 
            mympi_status,     & ! status
            mympierror)
    enddo
    
    ! Reporting to the log file
    write(logf,*) 'Mean number of sources per processor: ', &
         real(NumSrc)/real(npr-1)
    write(logf,*) 'Counted mean number of sources per processor: ', &
         real(sum(counts(1:npr-1)))/real(npr-1)
    write(logf,*) 'Minimum and maximum number of sources ', &
         'per processor: ', &
         minval(counts(1:npr-1)),maxval(counts(1:npr-1))
    flush(logf)

#endif

  end subroutine do_grid_master

  ! ===========================================================================

  !> The slave task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_slave(dt,niter)

    ! The slave task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: local_count
    integer :: ns1
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    local_count=0
    call MPI_Recv (ns1,  & ! address of receive buffer
         1,    & ! number of items to receive
         MPI_INTEGER,  & ! type of data
         0,		  & ! can receive from master only
         MPI_ANY_TAG,  & ! can expect two values, so
         ! we use the wildcard MPI_ANY_TAG 
         ! here
         MPI_COMM_NEW, & ! communicator
         mympi_status, & ! status
         mympierror)
    
    ! if tag equals 2, then skip the calculations
    
    if (mympi_status(MPI_TAG) /= 2) then
       do 
#ifdef MPILOG
          ! Report
          write(logf,*) 'Processor ',rank,' received: ',ns1
          write(logf,*) ' that is source ',ns1 !SrcSeries(ns1)
          write(logf,*) ' at:',srcpos(:,ns1)
          flush(logf)
#endif
          ! Do the source at hand
          call do_source(dt,ns1,niter)
          
          ! Update local counter
          local_count=local_count+1
          
          ! Send 'answer'
          call MPI_Send (local_count, 1,  & ! sending one int 
               MPI_INTEGER, 0, & ! to master
               1,              & ! tag
               MPI_COMM_NEW,   & ! communicator
               mympierror)
          
          ! Receive new source number
          call MPI_Recv (ns1,     & ! address of receive buffer
               1,            & ! number of items to receive
               MPI_INTEGER,  & ! type of data
               0,            & ! can receive from master only
               MPI_ANY_TAG,  & !  can expect two values, so
               !  we use the wildcard MPI_ANY_TAG 
               !  here
               MPI_COMM_NEW, & ! communicator
               mympi_status, & ! status
               mympierror)
          
          ! leave this loop if tag equals 2
          if (mympi_status(MPI_TAG) == 2) then
#ifdef MPILOG
             write(logf,*) 'Stop received'
             flush(logf)
#endif
             exit 
          endif
       enddo
    endif
    
    ! this is the point that is reached when a task is received with
    ! tag = 2
    
    ! send the number of calculations to master and return
    
#ifdef MPILOG
    ! Report
    write(logf,*) 'Processor ',rank,' did ',local_count,' sources'
    flush(logf)
#endif
    call MPI_Send (local_count,  &
         1,           & 
         MPI_INTEGER, & ! sending one int
         0,           & ! to master
         7,           & ! tag
         MPI_COMM_NEW,& ! communicator
         mympierror)
#endif

  end subroutine do_grid_slave

end module master_slave_processing
