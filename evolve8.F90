!>
!! \brief This module contains routines for calculating the ionization and temperature evolution of the entire grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 3D, no OpenMP

module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (3D).
     
  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! - evolve3D : step through grid
  ! - evolve2D : step through z-plane
  ! - evolve0D : take care of one grid point
  ! - evolve0D_global: update entire grid

  ! Needs:
  ! doric : ionization calculation for one point + photo-ionization rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,iterdump
  use sizes, only: Ndim, mesh
  use grid, only: x,y,z,vol,dr
  use material, only: ndens, xh, temper, get_temperature_point
  use sourceprops, only: NumSrc, srcpos, NormFlux !SrcSeries
  use radiation, only: NumFreqBnd
  use photonstatistics, only: state_before, calculate_photon_statistics, &
       photon_loss, LLS_loss, report_photonstatistics, state_after, total_rates, &
       total_ionizations, update_grandtotal_photonstatistics
  use c2ray_parameters, only: convergence_fraction, subboxsize, max_subbox, S_star_nominal

  implicit none

  private

  public :: evolve3D, phih_grid, evolve_ini

  !> Periodic boundary conditions, has to be true for this version
  logical,parameter :: periodic_bc = .true.

  !> Minimum number of MPI processes for using the master-slave setup 
  integer, parameter ::  min_numproc_master_slave=10

  ! Grid variables

  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: phih_grid
  !> Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_av
  !real(kind=dp),dimension(:,:,:),allocatable :: xh_av
  !> Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_intermed
  !real(kind=dp),dimension(:,:,:),allocatable :: xh_intermed
  !> H0 Column density (outgoing)
  real(kind=dp),dimension(:,:,:),allocatable :: coldensh_out
  !> Buffer for MPI communication
  real(kind=dp),dimension(:,:,:),allocatable :: buffer
  !> Photon loss from the grid
  real(kind=dp) :: photon_loss_all(1:NumFreqBnd)
  !> Photon loss from one source
  real(kind=dp),dimension(:,:),allocatable :: photon_loss_src_thread
  real(kind=dp) :: photon_loss_src(1:NumFreqBnd)

  ! mesh positions of end points for RT
  integer,dimension(Ndim) :: lastpos_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: lastpos_r !< mesh position of right end point for RT
  integer,dimension(Ndim) :: last_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: last_r !< mesh position of right end point for RT

  integer :: sum_nbox !< sum of all nboxes (on one processor)
  integer :: sum_nbox_all !< sum of all nboxes (on all processors)

  integer :: tn !< thread number

contains

  ! =======================================================================

  !> Allocate the arrays needed for evolve
  subroutine evolve_ini ()
    
    allocate(phih_grid(mesh(1),mesh(2),mesh(3)))
    allocate(xh_av(mesh(1),mesh(2),mesh(3),0:1))
    !allocate(xh_av(mesh(1),mesh(2),mesh(3)))
    allocate(xh_intermed(mesh(1),mesh(2),mesh(3),0:1))
    !allocate(xh_intermed(mesh(1),mesh(2),mesh(3)))
    allocate(coldensh_out(mesh(1),mesh(2),mesh(3)))
    allocate(buffer(mesh(1),mesh(2),mesh(3)))
    allocate(photon_loss_src_thread(1:NumFreqBnd,nthreads))

  end subroutine evolve_ini

  ! ===========================================================================

  !> Evolve the entire grid over a time step dt
  subroutine evolve3D (dt,restart)

    ! Calculates the evolution of the hydrogen ionization state
     
    ! Author: Garrelt Mellema
     
    ! Date: 28-Feb-2008 (21-Aug-2006 (f77/OMP: 13-Jun-2005))

    ! Version: Multiple sources / Using average fractions to converge
    ! loop over sources
    
    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz
    ! 21-Aug-2006 (GM) : MPI parallelization over the sources (static model).
    ! 28-Feb-2008 (GM) : Added master-slave model for distributing
    !                    over the processes. The program decides which
    !                    model to use.

    ! The time step
    real(kind=dp),intent(in) :: dt !< time step
    integer,intent(in) :: restart !< restart flag

    ! Loop variables
    integer :: niter  ! iteration counter

    ! Wall clock counting
    ! 8 bytes to beat the maxcount
    integer(kind=8) :: wallclock1
    integer(kind=8) :: wallclock2
    integer(kind=8) :: countspersec

    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag

    ! Minimum number of cells which are allowed to be non-converged
    integer :: conv_criterion 

#ifdef MPI
    integer :: mympierror
#endif

    ! End of declarations

    ! Initialize wall clock counter (for dumps)
    call system_clock(wallclock1)

     ! Initial state (for photon statistics)
    call state_before (xh)

    ! initialize average and intermediate results to initial values
    if (restart == 0) then
       xh_av(:,:,:,:)=xh(:,:,:,:)
       xh_intermed(:,:,:,:)=xh(:,:,:,:)
       !xh_av(:,:,:)=xh(:,:,:,1)
       !xh_intermed(:,:,:)=xh(:,:,:,1)
       niter=0 ! iteration starts at zero
       conv_flag=mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
    else
       ! Reload xh_av,xh_intermed,photon_loss,niter
       call start_from_dump(niter)
       call global_pass (conv_flag,dt)
    endif

    ! Set the conv_criterion, if there are few sources we should make
    ! sure that things are converged around these sources.
    conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)), &
         (NumSrc-1)/3)

    ! Iterate to reach convergence for multiple sources
    do
       ! Update xh if converged and exit
       if (conv_flag <= conv_criterion) then
          xh(:,:,:,:)=xh_intermed(:,:,:,:)
          !xh(:,:,:)=xh_intermed(:,:,:)
          exit
       else
          if (niter > 100) then
             ! Complain about slow convergence
             if (rank == 0) write(logf,*) 'Multiple sources not converging'
             exit
          endif
       endif
       
       ! Iteration loop counter
       niter=niter+1
       
       call pass_all_sources (niter,dt)
       
       ! Report subbox statistics
       if (rank == 0) &
            write(logf,*) "Average number of subboxes: ",sum_nbox_all/NumSrc

       if (rank == 0) then
          call system_clock(wallclock2,countspersec)
          ! Write iteration dump if more than 15 minutes have passed.
          ! system_clock starts counting at 0 when it reaches
          ! a max value. To catch this, test also for negative
          ! values of wallclock2-wallclock1
          write(logf,*) "Time and limit are: ", &
               wallclock2-wallclock1, 15.0*60.0*countspersec
          if (wallclock2-wallclock1 > 15*60*countspersec .or. &
               wallclock2-wallclock1 < 0 ) then
             call write_iteration_dump(niter)
             wallclock1=wallclock2
          endif
       endif

       call global_pass (conv_flag,dt)
       
    enddo

    ! Calculate photon statistics
    call calculate_photon_statistics (dt,xh,xh_av)
    call report_photonstatistics (dt)
    call update_grandtotal_photonstatistics (dt)

  end subroutine evolve3D

  ! ===========================================================================

  subroutine write_iteration_dump (niter)

    integer,intent(in) :: niter  ! iteration counter

    integer :: ndump=0
    
    character(len=20) :: iterfile

    ndump=ndump+1
    if (mod(ndump,2) == 0) then
       iterfile="iterdump2.bin"
    else
       iterfile="iterdump1.bin"
    endif

    open(unit=iterdump,file=iterfile,form="unformatted", &
         status="unknown")

    write(iterdump) niter
    write(iterdump) photon_loss_all
    write(iterdump) phih_grid
    write(iterdump) xh_av
    write(iterdump) xh_intermed

    close(iterdump)

  end subroutine write_iteration_dump

  ! ===========================================================================

  subroutine start_from_dump(niter)

    integer,intent(out) :: niter  ! iteration counter

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then
       open(unit=iterdump,file="iterdump.bin",form="unformatted", &
            status="old")
       
       read(iterdump) niter
       read(iterdump) photon_loss_all
       read(iterdump) phih_grid
       read(iterdump) xh_av
       read(iterdump) xh_intermed
       
       close(iterdump)
       write(logf,*) "Read iteration ",niter," from dump file"
       write(logf,*) 'photon loss counter: ',photon_loss_all
       write(logf,*) "Intermediate result for mean ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
            !sum(xh_intermed(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(niter,1, &
         MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(photon_loss_all,NumFreqBnd, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(phih_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xh_av,mesh(1)*mesh(2)*mesh(3)*2, &
    !call MPI_BCAST(xh_av,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xh_intermed,mesh(1)*mesh(2)*mesh(3)*2, &
    !call MPI_BCAST(xh_intermed,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine start_from_dump

  ! ===========================================================================

  subroutine pass_all_sources(niter,dt)
    
    ! For random permutation of sources:
    use  m_ctrper, only: ctrper

    integer,intent(in) :: niter  ! iteration counter
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D

    real(kind=dp) :: LLS_loss_all

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) write(logf,*) 'Doing all sources '
    ! reset global rates to zero for this iteration
    phih_grid(:,:,:)=0.0
    
    ! reset photon loss counters
    photon_loss(:) = 0.0
    LLS_loss = 0.0 ! make this a NumFreqBnd vector if needed later (GM/101129)

    ! Reset sum of subboxes counter
    sum_nbox=0

    ! Make a randomized list of sources :: call in serial
    ! disabled / GM110512
    !if ( rank == 0 ) call ctrper (SrcSeries(1:NumSrc),1.0)
    
#ifdef MPI
    ! Distribute the source list to the other nodes
    ! disabled / GM110512
    !call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Ray trace the whole grid for all sources.
    ! We can do this in two ways, depending on
    ! the number of processors. For many processors
    ! the master-slave setup should be more efficient.
    if (npr > min_numproc_master_slave) then
       call do_grid_master_slave (dt,niter)
    else
       call do_grid_static (dt,niter)
    endif
    
#ifdef MPI
    ! accumulate (sum) the MPI distributed photon losses
    call MPI_ALLREDUCE(photon_loss, photon_loss_all, NumFreqBnd, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)

    ! accumulate (sum) the MPI distributed photon losses
    call MPI_ALLREDUCE(LLS_loss, LLS_loss_all, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Put LLS_loss_all back in the LLS variable
    LLS_loss = LLS_loss_all

    ! accumulate (sum) the MPI distributed phih_grid
    call MPI_ALLREDUCE(phih_grid, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    
    ! Overwrite the processor local values with the accumulated value
    phih_grid(:,:,:)=buffer(:,:,:)

    ! accumulate (sum) the MPI distributed sum of number of boxes
    call MPI_ALLREDUCE(sum_nbox, sum_nbox_all, 1, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_NEW, mympierror)

    ! Only on the first iteration does evolve2D (evolve0D) change the
    ! ionization fractions
    if (niter == 1) then
       ! accumulate (max) MPI distributed xh_av
       call MPI_ALLREDUCE(xh_av(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
       !call MPI_ALLREDUCE(xh_av(:,:,:), buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)
       
       ! Overwrite the processor local values with the accumulated value
       xh_av(:,:,:,1) = buffer(:,:,:)
       xh_av(:,:,:,0) = max(0.0_dp,min(1.0_dp,1.0-xh_av(:,:,:,1)))
       !xh_av(:,:,:) = buffer(:,:,:)
       
       ! accumulate (max) MPI distributed xh_intermed
       call MPI_ALLREDUCE(xh_intermed(:,:,:,1), buffer, &
       !call MPI_ALLREDUCE(xh_intermed(:,:,:), buffer, &
            mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_NEW, mympierror)
       
       ! Overwrite the processor local values with the accumulated value
       xh_intermed(:,:,:,1)=buffer(:,:,:)
       xh_intermed(:,:,:,0)=max(0.0_dp,min(1.0_dp,1.0-xh_intermed(:,:,:,1)))
       !xh_intermed(:,:,:)=buffer(:,:,:)
    endif
#else
    photon_loss_all(:)=photon_loss(:)
    sum_nbox_all=sum_nbox
#endif
    
  end subroutine pass_all_sources

  ! ===========================================================================

  subroutine global_pass (conv_flag,dt)

    ! Flag variable (passed back from evolve0D_global)
    integer,intent(out) :: conv_flag
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D

    integer :: i,j,k  ! mesh position

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: pos

    ! Report photon losses over grid boundary 
    if (rank == 0) write(logf,*) 'photon loss counter: ',photon_loss_all(:)
    
    ! Turn total photon loss into a mean per cell (used in evolve0d_global)
    ! GM/110225: Possibly this should be 
    ! photon_loss_all(:)/(real(mesh(1))*real(mesh(2))
    ! Since this is the correct answer for an optically thin volume
    ! Compromise: photon_loss_all(:)/(real(mesh(1))**3)*sum(xh_av(:,:,:,1))
    ! then if the box is fully ionized the optically thin 1/N^2 is used
    ! and if the box is more neutral something close to 1/N^3 is used.
    ! Not sure
    photon_loss(:)=photon_loss_all(:)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
 
    ! Report minimum value of xh_av(0) to check for zeros
    if (rank == 0) then
       write(logf,*) "min xh_av(0): ",minval(xh_av(:,:,:,0))
       !write(logf,*) "min xh_av(0): ",minval(1.0_dp-xh_av(:,:,:))
    endif
    
    ! Apply total photo-ionization rates from all sources (phih_grid)
    conv_flag=0 ! will be used to check for convergence
    
    ! Loop through the entire mesh
    if (rank == 0) write(logf,*) 'Doing global '
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             pos=(/ i,j,k /)
             call evolve0D_global(dt,pos,conv_flag)
          enddo
       enddo
    enddo
    
    ! Report on convergence and intermediate result
    if (rank == 0) then
       write(logf,*) "Number of non-converged points: ",conv_flag
       write(logf,*) "Intermediate result for mean H ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
            !sum(xh_intermed(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
    endif
    
    ! Report on photon conservation
    call calculate_photon_statistics (dt,xh_intermed,xh_av)
    call report_photonstatistics (dt)
    
  end subroutine global_pass

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
          
#ifdef MPILOG
          ! Report ionization fractions
          write(logf,*) sum(xh_intermed(:,:,:,1))/ &
          !write(logf,*) sum(xh_intermed(:,:,:))/ &
               real(mesh(1)*mesh(2)*mesh(3))
          write(logf,*) sum(xh_av(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          !write(logf,*) sum(xh_av(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
          write(logf,*) local_count
#endif
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
       call do_source(dt,ns1,niter)
    enddo

  end subroutine do_grid_static

  ! ===========================================================================
  
  !> Does the ray-tracing over the entire 3D grid for one source.
  !! The number of this source in the current list is ns1.
  subroutine do_source(dt,ns1,niter)

    ! Does the ray-tracing over the entire 3D grid for one source.
    ! The number of this source in the current list is ns1.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer, intent(in) :: ns1 !< number of the source being done
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: naxis,nplane,nquadrant
    integer :: ns
    integer :: k
    integer :: nbox
    integer ::nnt

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: rtpos
      
    ! Pick up source number from the source list
    !ns=SrcSeries(ns1)
    ns=ns1
    
    ! reset column densities for new source point
    ! coldensh_out is unique for each source point
    coldensh_out(:,:,:)=0.0
    
    ! Find the mesh position for the end points of the loop
    ! We trace until we reach max_subbox (set in c2ray_parameters)
    ! or the end of the grid. In the periodic case the end of the
    ! grid is always mesh/2 away from the source. If the grid is
    ! even-sized we trave mesh/2 cells to the left and mesh/2-1
    ! cell to the right. If it is odd, it is mesh/2 in either direction.
    ! The mod(mesh,2) takes care of handling this.
    if (periodic_bc) then
       lastpos_r(:)=srcpos(:,ns)+min(max_subbox,mesh(:)/2-1+mod(mesh(:),2))
       lastpos_l(:)=srcpos(:,ns)-min(max_subbox,mesh(:)/2)
    else
       lastpos_r(:)=min(srcpos(:,ns)+max_subbox,mesh(:))
       lastpos_l(:)=max(srcpos(:,ns)-max_subbox,1)
    endif

    ! Loop through grid in the order required by 
    ! short characteristics
    
    ! Transfer is done in a set of cubes of increasing size.
    ! If the HII region is small we do not waste time calculating
    ! column densities of parts of the grid where no radiation
    ! penetrates. To test whether the current subbox is large
    ! enough we use the photon_loss_src. If this is non-zero,
    ! photons are leaving this subbox and we need to do another
    ! one. We also stop once we have done the whole grid.
    nbox=0 ! subbox counter
    photon_loss_src(:)=NormFlux(ns)*S_star_nominal !-1.0 ! to pass the first while test
    last_r(:)=srcpos(:,ns) ! to pass the first while test
    last_l(:)=srcpos(:,ns) ! to pass the first while test

    ! Loop through boxes of increasing size
    ! NOTE: make this limit on the photon_loss a fraction of
    ! a source flux loss_fraction*NormFlux(ns)*S_star_nominal)
    do while (all(photon_loss_src(:) > 1e-6*NormFlux(ns)*S_star_nominal) &
    !do while (all(photon_loss_src(:) /= 0.0) &
         .and. last_r(3) < lastpos_r(3) &
         .and. last_l(3) > lastpos_l(3))
       nbox=nbox+1 ! increase subbox counter
       photon_loss_src(:) = 0.0 ! reset photon_loss_src to zero
       photon_loss_src_thread(:,:) = 0.0 ! reset photon_loss_src to zero
       last_r(:)=min(srcpos(:,ns)+subboxsize*nbox,lastpos_r(:))
       last_l(:)=max(srcpos(:,ns)-subboxsize*nbox,lastpos_l(:))

       ! OpenMP: if we have multiple OpenMP threads (nthreads > 1) we 
       ! parallelize over the threads by doing independent parts of
       ! the mesh.
       if (nthreads > 1) then ! OpenMP parallelization

          ! First do source point (on first pass)
          if (nbox == 1) then
             rtpos(:)=srcpos(:,ns)
             call evolve0D(dt,rtpos,ns,niter)
          endif

          ! do independent areas of the mesh in parallel using OpenMP
          !$omp parallel default(shared) private(tn)
          !!!reduction(+:photon_loss_src)

          ! Find out your thread number
#ifdef MY_OPENMP
          tn=omp_get_thread_num()+1
#else
          tn=1
#endif
          
          ! Then do the the axes
          !$omp do schedule(dynamic,1)
          do naxis=1,6
             call evolve1D_axis(dt,ns,niter,naxis)
          enddo
          !$omp end do
          
          ! Then the source planes
          !$omp do schedule (dynamic,1)
          do nplane=1,12
             call evolve2D_plane(dt,ns,niter,nplane)
          end do
          !$omp end do
          
          ! Then the quadrants
          !$omp do schedule (dynamic,1)
          do nquadrant=1,8
             call evolve3D_quadrant(dt,ns,niter,nquadrant)
          end do
          !$omp end do
          
          !$omp end parallel
          ! Collect photon losses for each thread
          do nnt=1,nthreads
             photon_loss_src(:)=photon_loss_src(:) + &
                  photon_loss_src_thread(:,nnt)
          enddo

       else ! No OpenMP parallelization

          ! 1. transfer in the upper part of the grid 
          !    (srcpos(3)-plane and above)
          do k=srcpos(3,ns),last_r(3)
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do
          
          ! 2. transfer in the lower part of the grid (below srcpos(3))
          do k=srcpos(3,ns)-1,last_l(3),-1
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do
          ! No OpenMP threads so we use position 1
          photon_loss_src(:)=photon_loss_src_thread(:,1)

       endif
    enddo
    
    ! Record the final photon loss, this is the photon loss that leaves
    ! the grid.
    photon_loss(:)=photon_loss(:) + photon_loss_src(:)

    ! Sum the total number of subboxes used for reporting later
    sum_nbox=sum_nbox+nbox

  end subroutine do_source

  ! ===========================================================================

  !> Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
  !! directions.
  subroutine evolve2D(dt,rtpos,ns,niter)

    ! Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
    ! directions.
    
    real(kind=dp),intent(in) :: dt      !! passed on to evolve0D
    integer,dimension(Ndim),intent(inout) :: rtpos !< mesh position, pos(3) is
                                                 !! intent(in)
    integer,intent(in) :: ns           !< current source
    integer,intent(in) :: niter        !< passed on to evolve0D

    integer :: i,j ! mesh positions

    ! sweep in `positive' j direction
    do j=srcpos(2,ns),last_r(2)
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
       end do
    end do
    
    ! sweep in `negative' j direction
    do j=srcpos(2,ns)-1,last_l(2),-1
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
       end do
    end do

  end subroutine evolve2D

  ! ===========================================================================

  ! Ray tracing for the axes going through the source point
  ! should be called after having done the source point
  subroutine evolve1D_axis(dt,ns,niter,naxis)

    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: naxis        ! axis to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (naxis)
    case(1)
       ! sweep in +i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)+1,last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `negative' k
       end do
    end select
    
  end subroutine evolve1D_axis

  ! ===========================================================================

  !> Ray tracing for planes containing the source point
  !! should be called after evolve1D_axis
  subroutine evolve2D_plane(dt,ns,niter,nplane)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nplane        ! plane to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nplane)
    case(1)
       ! sweep in +i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
       
    end select
    
  end subroutine evolve2D_plane

  ! ===========================================================================

  !> Ray tracing for the 8 octants 
  !! should be called after evolve2D_plane
  subroutine evolve3D_quadrant(dt,ns,niter,nquadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    
    real(kind=dp),intent(in) :: dt     ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nquadrant    ! which quadrant to do    

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nquadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    end select

  end subroutine evolve3D_quadrant

  !=======================================================================

  !> Calculates the photo-ionization rate for one cell due to one source
  !! and adds this contribution to the collective rate.
  subroutine evolve0D(dt,rtpos,ns,niter)
    
    ! Calculates the photo-ionization rate for one cell due to one source
    ! and adds this contribution to the collective rate.
    
    ! Author: Garrelt Mellema
    
    ! Date: 01-Feb-2008 (21-Aug-2006, 20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: multiple sources, fixed temperature
    
    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated.
    ! For the first pass (niter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence speed
    ! in the case of isolated sources.

    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use radiation, only: photoion, photrates
    use material, only: clumping_point !,coldensh_LLS
    use c2ray_parameters, only: epsilon,convergence1,convergence2, &
         type_of_clumping, convergence_frac,use_LLS,type_of_LLS
    use mathconstants, only: pi
    use material, only: coldensh_LLS, LLS_point
    use photonstatistics, only: total_LLS_loss

    ! column density for stopping chemisty
    real(kind=dp),parameter :: max_coldensh=2e19_dp 
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: rtpos ! cell position (for RT)
    integer,intent(in)      :: ns ! source number 
    integer,intent(in)      :: niter ! global iteration number

    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    integer,dimension(Ndim) :: srcpos1
    real(kind=dp) :: dist2,path,vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp) :: ndens_p
    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0
    real(kind=dp) :: convergence
    
    type(photrates) :: phi

    !write(*,*) rtpos
    ! set convergence tolerance
    convergence=convergence1

    ! Map pos to mesh pos, assuming a periodic mesh
    do idim=1,Ndim
       pos(idim)=modulo(rtpos(idim)-1,mesh(idim))+1
    enddo

    ! If coldensh_out is zero, we have not yet done this point
    ! yet, so do it. Otherwise do nothing.
    if (coldensh_out(pos(1),pos(2),pos(3)) == 0.0) then
       ! Initialize local ionization states to the global ones
       yh0(1)=xh(pos(1),pos(2),pos(3),1)
       !yh0(1)=xh(pos(1),pos(2),pos(3))
       yh0(0)=xh(pos(1),pos(2),pos(3),0)
       !yh0(0)=1.0_dp-yh0(1)
       yh_av(1)=xh_av(pos(1),pos(2),pos(3),1)
       !yh_av(1)=xh_av(pos(1),pos(2),pos(3))
       yh_av(0)=xh_av(pos(1),pos(2),pos(3),0)
       !yh_av(0)=1.0_dp-yh_av(1)
       
       ! Initialize local density and temperature
       ndens_p=ndens(pos(1),pos(2),pos(3))
       call get_temperature_point(pos(1),pos(2),pos(3))

       ! Find the column density at the entrance point of the cell (short
       ! characteristics)
       
       if ( all( rtpos(:) == srcpos(:,ns) ) ) then
          ! Do not call cinterp for the source point.
          ! Set coldensh and path by hand
          coldensh_in=0.0
          path=0.5*dr(1)
          
          ! Find the distance to the source (average?)
          !dist=0.5*dr(1) NOT NEEDED         ! this makes vol=dx*dy*dz
          !vol_ph=4.0/3.0*pi*dist**3
          vol_ph=dr(1)*dr(2)*dr(3)
          
       else
          
          ! For all other points call cinterp to find the column density
          !do idim=1,Ndim
          !   srcpos1(idim)=srcpos(idim,ns)
          !enddo
          call cinterp(rtpos,srcpos(:,ns),coldensh_in,path)
          path=path*dr(1)
          
          ! Find the distance to the source
          xs=dr(1)*real(rtpos(1)-srcpos(1,ns))
          ys=dr(2)*real(rtpos(2)-srcpos(2,ns))
          zs=dr(3)*real(rtpos(3)-srcpos(3,ns))
          dist2=xs*xs+ys*ys+zs*zs
          
          ! Find the volume of the shell this cell is part of 
          ! (dilution factor).
          vol_ph=4.0*pi*dist2*path

          ! Add LLS opacity
          ! GM/110224: previously we added this to coldensh_out,
          ! but this gives funny results for the photo-ionization
          ! rate which is based on the difference between the
          ! in and out column density. To mimick a LLS fog, it
          ! may be better to add it here.
          ! Initialize local LLS (if type of LLS is appropriate)
          if (use_LLS) then
             if (type_of_LLS == 2) call LLS_point (pos(1),pos(2),pos(3))
             coldensh_in = coldensh_in + coldensh_LLS * path/dr(1)
          endif
       endif
       
       ! Only ray trace and exit. Do not touch the ionization
       ! fractions. They are updated using phih_grid in evolve0d_global

       ! Only do chemistry if this is the first pass over the sources,
       ! and if column density is below the maximum.
       ! On the first global iteration pass it may be beneficial to assume 
       ! isolated sources, but on later passes the effects of multiple sources 
       ! has to be taken into account. 
       ! Therefore no changes to xh, xh_av, etc. should happen on later passes!
       if (niter == 1 .and. coldensh_in < max_coldensh) then
          call do_chemistry (dt, ndens_p, yh, yh_av, yh0, 0.0_dp, &
               coldensh_in, path, vol_ph, pos, ns, local=.true.)

          ! Copy ion fractions tp global arrays.
          ! This will speed up convergence if
          ! the sources are isolated and only ionizing up.
          ! In other cases it does not make a difference.
          xh_intermed(pos(1),pos(2),pos(3),1)=max(yh(1), &
               xh_intermed(pos(1),pos(2),pos(3),1))
          !xh_intermed(pos(1),pos(2),pos(3))=max(yh(1), &
          !     xh_intermed(pos(1),pos(2),pos(3)))
          !xh_intermed(pos(1),pos(2),pos(3),0)=1.0- &
          !     xh_intermed(pos(1),pos(2),pos(3),1)
          xh_av(pos(1),pos(2),pos(3),1)=max(yh_av(1), &
               xh_av(pos(1),pos(2),pos(3),1))
          !xh_av(pos(1),pos(2),pos(3))=max(yh_av(1), &
          !     xh_av(pos(1),pos(2),pos(3)))
          !xh_av(pos(1),pos(2),pos(3),0)=1.0- &
          !     xh_av(pos(1),pos(2),pos(3),1)
          
       endif
          
       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       ! and add the LLS column density to this.
       ! GM/110224: No! This messes up phi since phi is based
       !  upon the difference between the in and out column density.
       !  Instead add the LLS to coldensh_in, see above
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,yh_av(0),ndens_p) !+ coldensh_LLS * path/dr(1)
       !###################################################################
       
       ! Calculate (photon-conserving) photo-ionization rate
       if (coldensh_in < max_coldensh) then
          call photoion(phi,coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
               vol_ph,ns)
          phi%h=phi%h/(yh_av(0)*ndens_p)
          ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
          ! the photon losses.
          ! GM/110302: Use phi%h_in (not out) here since this is where we draw the
          ! photons from, above.
          call total_LLS_loss(phi%h_in*vol/vol_ph, coldensh_LLS * path/dr(1))
       else
          phi%h=0.0
          phi%h_out=0.0
       endif
       
       ! Add photo-ionization rate to the global array 
       ! (applied in evolve0D_global)
       phih_grid(pos(1),pos(2),pos(3))= &
            phih_grid(pos(1),pos(2),pos(3))+phi%h
       
       ! Photon statistics: register number of photons leaving the grid
       if ( (any(rtpos(:) == last_l(:))) .or. &
            (any(rtpos(:) == last_r(:))) ) then
          photon_loss_src_thread(1,tn)=photon_loss_src_thread(1,tn) + &
               phi%h_out*vol/vol_ph
          !photon_loss_src(1,tn)=photon_loss_src(1,tn) + phi%h_out*vol/vol_ph
       endif

    endif ! end of coldens test

  end subroutine evolve0D

  ! =======================================================================

  !> Calculates the evolution of the hydrogen ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global(dt,pos,conv_flag)

    ! Calculates the evolution of the hydrogen ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.
    
    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping, convergence_frac

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter

    integer :: nx,nit ! loop counters
    real(kind=dp) :: de ! electron density
    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0 ! ionization fractions
    real(kind=dp) :: yh_av0
    real(kind=dp) :: avg_temper ! temperature
    real(kind=dp) :: ndens_p ! local number density

    real(kind=dp) :: phih ! local H photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp) :: phih_total ! local total photo-ionization rate (including
                                ! photon loss term)

    real(kind=dp) :: convergence
    
    ! Set convergence tolerance
    convergence=convergence2

    ! Initialize local ionization states to global ones
    yh0(1)=xh(pos(1),pos(2),pos(3),1)
    !yh0(1)=xh(pos(1),pos(2),pos(3))
    yh0(0)=xh(pos(1),pos(2),pos(3),0)
    !yh0(0)=1.0_dp-yh0(1)
    yh_av(1)=xh_av(pos(1),pos(2),pos(3),1) ! use calculated xh_av
    !yh_av(1)=xh_av(pos(1),pos(2),pos(3))
    yh_av(0)=xh_av(pos(1),pos(2),pos(3),0)
    !yh_av(0)=1.0_dp-yh_av(1)

    ! Initialize local scalars for density and temperature
    ndens_p=ndens(pos(1),pos(2),pos(3))
    call get_temperature_point(pos(1),pos(2),pos(3))
    avg_temper=temper

    ! Use the collected photo-ionization rates
    phih=phih_grid(pos(1),pos(2),pos(3))

    call do_chemistry (dt, ndens_p, yh, yh_av, yh0, phih, &
         0.0_dp, 0.0_dp, 0.0_dp, pos, 0 , local=.false.)

    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence
    yh_av0=xh_av(pos(1),pos(2),pos(3),0) ! use previously calculated xh_av
    !yh_av0=1.0_dp-xh_av(pos(1),pos(2),pos(3)) ! use previously calculated xh_av
    if (abs((yh_av(0)-yh_av0)) > convergence2 .and. &
         (abs((yh_av(0)-yh_av0)/yh_av(0)) > convergence2 .and. &
         (yh_av(0) > convergence_frac))) then
       conv_flag=conv_flag+1
    endif

    ! Copy ion fractions to the global arrays.
    xh_intermed(pos(1),pos(2),pos(3),0)=yh(0)
    xh_intermed(pos(1),pos(2),pos(3),1)=yh(1)
    xh_av(pos(1),pos(2),pos(3),0)=yh_av(0)
    xh_av(pos(1),pos(2),pos(3),1)=yh_av(1)
    !xh_intermed(pos(1),pos(2),pos(3))=yh(1)
    !xh_av(pos(1),pos(2),pos(3))=yh_av(1)

  end subroutine evolve0D_global

  ! ===========================================================================

  subroutine do_chemistry (dt, ndens_p, yh, yh_av, yh0, &
       phih, coldensh_in, path, vol_ph, pos, ns, local)

    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping, & 
         convergence_frac, add_photon_losses
    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use material, only: clumping_point
    use radiation, only: photoion, photrates

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: ndens_p
    real(kind=dp),dimension(0:1),intent(out) :: yh
    real(kind=dp),dimension(0:1),intent(inout) :: yh_av
    real(kind=dp),dimension(0:1),intent(in) :: yh0
    real(kind=dp),intent(in) :: phih !< local photo-ionization rate
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh
    integer,intent(in)      :: ns !< source number 
    logical,intent(in) :: local !< true if doing a non-global calculation.

    real(kind=dp) :: avg_temper
    real(kind=dp) :: phih_cell
    real(kind=dp) :: yh_av0
    real(kind=dp) :: de
    real(kind=dp) :: coldensh_cell

    integer :: nit

    type(photrates) :: phi

    ! Initialize yh to initial value
    yh(:)=yh0(:)

    ! Initialize local temperature
    avg_temper=temper

    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5) call clumping_point (pos(1),pos(2),pos(3))
    
    nit=0
    do 
       nit=nit+1

       ! Save the values of yh_av found in the previous iteration
       yh_av0=yh_av(0)

       ! Copy ionic abundances back to initial values (doric assumes
       ! that it contains this)
       yh(:)=yh0(:)
              
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,yh_av)

       ! Find total photo-ionization rate 
       if (local) then
          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(path,yh_av(0),ndens_p)
          
          ! Calculate (photon-conserving) photo-ionization rate
          call photoion(phi,coldensh_in,coldensh_in+coldensh_cell, &
               vol_ph,ns)
          phih_cell=phi%h/(yh_av(0)*ndens_p)
       else
          ! (direct plus photon losses)
          ! DO THIS HERE, yh_av is changing
          ! (if the cell is ionized, add a fraction of the lost photons)
          !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
          phih_cell=phih 
          !if (add_photon_losses) phih_cell=phih_cell + & 
          !     photon_loss(1)/(vol*yh_av(0)*ndens_p)
          ! GM/110225: New approach to lost photons, taking into
          ! account optical depth of cells. See notes.
          if (add_photon_losses) then
             NormFlux(0)=sum(photon_loss(:))/S_star_nominal
             ! Calculate (time averaged) column density of cell
             coldensh_cell=coldens(dr(1),yh_av(0),ndens_p)
             call photoion(phi,0.0d0,coldensh_cell,vol,0)
             phih_cell=phih_cell + phi%h/(yh_av(0)*ndens_p)
          endif
       endif

       ! Calculate the new and mean ionization states
       call doric(dt,avg_temper,de,ndens_p,yh,yh_av,phih_cell)

       ! Test for convergence on time-averaged neutral fraction
       ! For low values of this number assume convergence
       if ((abs((yh_av(0)-yh_av0)/yh_av(0)) < convergence2 &
             .or. (yh_av(0) < convergence_frac))) exit
                  
       ! Warn about non-convergence and terminate iteration
       if (nit > 5000) then
          if (rank == 0) then
             write(logf,*) 'Convergence failing (global)'
             write(logf,*) 'xh: ',yh_av(0),yh_av0
          endif
          exit
       endif
    enddo

  end subroutine do_chemistry

  ! ===========================================================================

  !> Finds the column density at pos as seen from the source point srcpos
  !! through interpolation. The interpolation
  !! depends on the orientation of the ray. The ray crosses either
  !! a z-plane, a y-plane or an x-plane.
  subroutine cinterp (rtpos,srcpos,cdensi,path)
    
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point srcpos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    integer,dimension(Ndim),intent(in) :: rtpos !< cell position (mesh)
    integer,dimension(Ndim),intent(in) :: srcpos !< source position (mesh)
    real(kind=dp),intent(out) :: cdensi !< column density to cell
    real(kind=dp),intent(out) :: path !< path length over cell

    real(kind=dp),parameter :: sqrt3=sqrt(3.0)
    real(kind=dp),parameter :: sqrt2=sqrt(2.0)

    integer :: i,j,k,i0,j0,k0

    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4
    real(kind=dp) :: di,dj,dk


    !DEC$ ATTRIBUTES FORCEINLINE :: weightf
    ! map to local variables (should be pointers ;)
    i=rtpos(1)
    j=rtpos(2)
    k=rtpos(3)
    i0=srcpos(1)
    j0=srcpos(2)
    k0=srcpos(3)
    
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel=i-i0
    jdel=j-j0
    kdel=k-k0
    idela=abs(idel)
    jdela=abs(jdel)
    kdela=abs(kdel)
    
    ! Find coordinates of points closer to source
    sgni=sign(1,idel)
!      if (idel == 0) sgni=0
    sgnj=sign(1,jdel)
!      if (jdel == 0) sgnj=0
    sgnk=sign(1,kdel)
!      if (kdel == 0) sgnk=0
    im=i-sgni
    jm=j-sgnj
    km=k-sgnk
    di=real(idel)
    dj=real(jdel)
    dk=real(kdel)

    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       alam=(real(km-k0)+sgnk*0.5)/dk
              
       xc=alam*di+real(i0) ! x of crossing point on z-plane 
       yc=alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx=2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1=(1.-dx)*(1.-dy)    ! interpolation weights of
       s2=(1.-dy)*dx         ! corner points to c-point
       s3=(1.-dx)*dy
       s4=dx*dy
       
       ! Map to rtpos to mesh pos, assuming a periodic mesh
       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jp=modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)    !# column densities at the
       c2=coldensh_out(ip,jmp,kmp)     !# four corners
       c3=coldensh_out(imp,jp,kmp)
       c4=coldensh_out(ip,jp,kmp)
       
       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       ! column density at the crossing point
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4) 

       ! Take care of diagonals
       ! if (kdela == idela.or.kdela == jdela) then
       ! if (kdela == idela.and.kdela == jdela) then
       ! cdensi=sqrt3*cdensi
       !else
       !cdensi=sqrt2*cdensi
       !endif
       !endif

       if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
          if (idela == 1.and.jdela == 1) then
             cdensi=sqrt3*cdensi
          else
             cdensi=sqrt2*cdensi
          endif
       endif
       ! if (kdela == 1) then
       ! if ((w3 == 1.0).or.(w2 == 1.0)) cdensi=sqrt(2.0)*cdensi
       ! if (w1 == 1.0) cdensi=sqrt(3.0)*cdensi
       ! write(logf,*) idela,jdela,kdela
       !endif

       ! Path length from c through d to other side cell.
       !dxp=di/dk
       !dyp=dj/dk
       path=sqrt((di*di+dj*dj)/(dk*dk)+1.0) ! pathlength from c to d point  


       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
          
       alam=(real(jm-j0)+sgnj*0.5)/dj
       zc=alam*dk+real(k0)
       xc=alam*di+real(i0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dx=2.0*abs(xc-(real(im)+0.5*sgni))
       s1=(1.-dx)*(1.-dz)
       s2=(1.-dz)*dx
       s3=(1.-dx)*dz
       s4=dx*dz
       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)
       c2=coldensh_out(ip,jmp,kmp)
       c3=coldensh_out(imp,jmp,kp)
       c4=coldensh_out(ip,jmp,kp)

       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       ! Take care of diagonals
       if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
          if (idela == 1.and.kdela == 1) then
             !write(logf,*) 'error',i,j,k
             cdensi=sqrt3*cdensi
          else
             !write(logf,*) 'diagonal',i,j,k
             cdensi=sqrt2*cdensi
          endif
       endif

       !dxp=di/dj
       !dzp=dk/dj
       !path=sqrt(dxp*dxp+1.0+dzp*dzp)
       path=sqrt((di*di+dk*dk)/(dj*dj)+1.0)
       

       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then

       alam=(real(im-i0)+sgni*0.5)/di
       zc=alam*dk+real(k0)
       yc=alam*dj+real(j0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1=(1.-dz)*(1.-dy)
       s2=(1.-dz)*dy
       s3=(1.-dy)*dz
       s4=dy*dz

       imp=modulo(im-1,mesh(1))+1
       jp=modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)
       c2=coldensh_out(imp,jp,kmp)
       c3=coldensh_out(imp,jmp,kp)
       c4=coldensh_out(imp,jp,kp)
       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       if ( idela == 1 .and. ( jdela == 1 .or. kdela == 1 ) ) then
          if ( jdela == 1 .and. kdela == 1 ) then
             cdensi=sqrt3*cdensi
          else
             cdensi=sqrt2*cdensi
          endif
       endif
       
       !dyp=dj/di
       !dzp=dk/di
       !path=sqrt(1.0+dyp*dyp+dzp*dzp)
       path=sqrt(1.0+(dj*dj+dk*dk)/(di*di))
       
    end if
    
  end subroutine cinterp

  ! =========================================================================

  !> Weight function for interpolation in cinterp
  real(kind=dp) function weightf (cd)

    use cgsphotoconstants, only: sigh

    real(kind=dp),intent(in) :: cd

    real(kind=dp),parameter :: minweight=1.0_dp/0.6_dp

    !weightf=1.0
    ! weightf=1.0/max(1.0d0,cd**0.54)
    ! weightf=exp(-min(700.0,cd*0.15*6.3d-18))
    weightf=1.0/max(0.6_dp,cd*sigh)

    ! weightf=1.0/log(max(e_ln,cd))

  end function weightf

end module evolve
