!>
!! \brief Main program for C2Ray-3Dm (3D, multiple sources)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 22-May-2008
!<
Program C2Ray

  ! Authors: Garrelt Mellema, Ilian Iliev

  ! Date: 22-May-2008 (06-Mar-2008 (30-Jan-2008 (8-Dec-2007 (23-Sep-2006))

  ! Goal:
  ! Cosmological reionization simulation using precomputed density fields
  ! and source lists.
  
  ! Does not include hydrodynamics
  ! Assumes constant time step

  ! Needs following `modules'
  ! c2ray_parameters : all the tunable parameters for the code
  ! my_mpi : sets up the MPI parallel environment
  ! output_module : output routines
  ! grid : sets up the grid
  ! radiation : radiation tools
  ! nbody : interface to N-body output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp
  use file_admin, only: stdinput, logf, file_input, flag_for_file_input
  use c2ray_parameters, only: cosmological, type_of_clumping
  use astroconstants, only: YEAR
  use my_mpi !, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini
  use nbody, only: nbody_type, nbody_ini, NumZred, zred_array, snap
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol, &
       time2zred, zred2time, zred
  use material, only: mat_ini, xfrac_ini, dens_ini, set_clumping
  use times, only: time_ini, set_timesteps
  use sourceprops, only: source_properties, NumSrc
  use evolve, only: evolve3D

#ifdef XLF
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

#ifdef PGI
  include 'lib3f.h' ! for iargc, getargc
#endif

  ! Start and end time for CPU report
  real :: tstart !< Start time for CPU report
  real :: tend !< End time for CPU report

  ! Wall clock time variables
  integer :: cntr1 !< Start time wall clock
  integer :: cntr2 !< End time wall clock
  integer :: countspersec !< counts per second (for wall clock time)

  integer :: restart !< restart flag
  integer :: nz !< loop counter for loop over redshift list
  integer :: nz0 !< index of starting redshift from redshift list
  integer :: ierror !< error flag

  ! end_time - end time of the simulation (s)
  ! dt - time step (s)
  ! sim_time - actual time (s)
  ! end_time - end time of calculation (s)
  ! output_time - time interval between outputs (s)
  ! next_output_time - time of next output (s)
  real(kind=dp) :: end_time !< end time of the simulation (s)
  real(kind=dp) :: sim_time !< actual time (s)
  real(kind=dp) :: output_time !< time interval between outputs (s)
  real(kind=dp) :: next_output_time !< time of next output (s)
  real(kind=dp) :: dt !< calculated time step
  real(kind=dp) :: actual_dt !< actual time step (s)
  real(kind=dp) :: zred_interm !< intermediate redshift (for restart)

#ifdef MPI
  integer :: mympierror
#endif

  ! Input file
  character(len=512) :: inputfile !< name of input file

  ! Initialize cpu timer
  call cpu_time(tstart)

  ! Initialize wall cock times
  call system_clock(cntr1)

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (iargc() > 0) then
     call getarg(1,inputfile)
     if (rank == 0) then
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     endif
  endif

  ! Initialize output
  call setup_output ()
  if (rank == 0) call flush(logf)

  ! Initialize grid
  call grid_ini ()

  ! Initialize photo-ionization calculation
  call rad_ini ( )

  ! Initialize the material properties
  call mat_ini (restart, nz0, ierror)

  ! Find the redshifts we are dealing with
  call nbody_ini ()

  ! Initialize time step parameters
  call time_ini ()
  if (rank == 0) call flush(logf)

  ! Set time to zero
  sim_time=0.0

  ! Initialize cosmology
  call cosmology_init(zred_array(nz0),sim_time)

  ! If a restart, read ionization fractions from file
  if (restart == 1) call xfrac_ini(zred_array(nz0))
  if (restart == 2) then
     if (rank == 0) then
        if (.not.file_input) &
             write(*,"(A,$)") "At which redshift to restart x_frac?:"
        read(stdinput,*) zred_interm
     endif
#ifdef MPI
     call MPI_BCAST(zred_interm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
          mympierror)
#endif
     call xfrac_ini(zred_interm)
  end if
  
  ! Loop over redshifts
  do nz=nz0,NumZred-1

     zred=zred_array(nz)
     if (rank == 0) write(logf,*) "Doing redshift: ",zred," to ", &
          zred_array(nz+1)
     
     ! Initialize time parameters
     call set_timesteps(zred,zred_array(nz+1), &
          end_time,dt,output_time)
     if (rank == 0) write(logf,*) &
          "This is time ",sim_time/YEAR," to ",end_time/YEAR
     if (rank == 0 .and. nbody_type == "LG") &
          write(logf,*) "This is snapshot ", snap(nz)
         
     ! Initialize source position
#ifdef MPILOG     
     write(logf,*) 'Calling source_properties'
#endif 
     call source_properties(zred,nz,end_time-sim_time,restart)

     ! Initialize density field
     call dens_ini(zred,nz)
     if (type_of_clumping == 5) call set_clumping(zred)

     ! Set time if restart at intermediate time
     ! Set next output time
     ! Note: this assumes that there are ALWAYS two time steps
     ! between redshifts. Should be made more general.
     if (restart == 2) then
        sim_time=zred2time(zred_interm)
        next_output_time=end_time
     else
        next_output_time=sim_time+output_time
     endif

     ! Reset restart
     restart=0

     ! Loop until end time is reached
     do
        ! Make sure you produce output at the correct time
        actual_dt=min(next_output_time-sim_time,dt)
        
        ! Report time and time step
        if (rank == 0) write(logf,"(A,2(1pe10.3,x),A)") "Time, dt:", &
             sim_time/YEAR,actual_dt/YEAR," (years)"

        ! For cosmological simulations evolve proper quantities
        if (cosmological) then
           call redshift_evol(sim_time+0.5*actual_dt)
           call cosmo_evol()
        endif
        ! Do not call in case of position dependent clumping,
        ! the clumping grid should have been initialized above
        if (type_of_clumping /= 5) call set_clumping(zred)

        ! Take one time step
        if (NumSrc > 0) call evolve3D(actual_dt)

        ! Update time
        sim_time=sim_time+actual_dt
            
        ! Write output
        if (abs(sim_time-next_output_time) <= 1e-6*sim_time) then
           call output(time2zred(sim_time),sim_time,dt)
           next_output_time=next_output_time+output_time
        endif
        if (abs(sim_time-end_time) <= 1e-6*end_time) exit
        if (rank == 0) call flush(logf)
     enddo

     ! stop
     ! Scale to the current redshift
     if (cosmological) then
        call redshift_evol(sim_time)
        call cosmo_evol()
     endif

  enddo

  ! Write final output

  call output(zred,sim_time,dt)

  call close_down ()
  
  ! Find out CPU time
  call cpu_time(tend)
  call system_clock(cntr2,countspersec)

  if (rank == 0) then
     write(logf,*) "CPU time: ",tend-tstart," s"
     write(logf,*) "Wall clock time: ",(cntr2-cntr1)/countspersec," s"
  endif

  ! End the run
  call mpi_end ()

end Program C2Ray
