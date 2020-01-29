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
  ! precision: definition of single and double precision
  ! clocks: data and routines related to clocks (cpu and wall clocks)
  ! file_admin: parameters related to file I/O
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
  use clocks, only: setup_clocks, update_clocks, report_clocks, &
       timestamp_wallclock
  use report_memory_module, only: report_memory
  use file_admin, only: stdinput, logf, timefile, file_input, &
       flag_for_file_input
  use c2ray_parameters, only: cosmological, type_of_clumping, &
       use_LLS, type_of_LLS,stop_on_photon_violation
  use astroconstants, only: YEAR
  use my_mpi !, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini, dr
  use radiation_tables, only: rad_ini
  use nbody, only: nbody_type, nbody_ini, NumZred, zred_array, snap
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol, &
       time2zred, zred2time, zred
  use material, only: material_ini
  use ionfractions_module, only: xfrac_restart_init
  use density_module, only: density_init
  use temperature_module, only: temperature_restart_init
  use clumping_module, only: set_clumping, load_clumping_model
  use lls_module, only: set_LLS
  use times, only: time_ini, set_timesteps
  use sourceprops, only: source_properties_ini, source_properties, NumSrc
  use evolve_data, only: evolve_ini
  use evolve, only: evolve3D
  use c2ray_parameters, only: isothermal ! NB: in some versions set in material

#ifdef XLF
  ! Place for xlf specific statements
#endif

  implicit none

#ifdef PGI
  include 'lib3f.h' ! for iargc, getargc
#endif

  integer :: startup_error=0 ! the use of this flag should be increased.
  ! The idea is to check several things before starting the simulation proper.
  ! If the startup_error is not zero, some problem was detected.
  integer :: restart=0 !< restart flag
  integer :: iter_restart=0 !< restart from iteration flag
  integer :: nz !< loop counter for loop over redshift list
  integer :: nz0 !< index of starting redshift from redshift list
  integer :: ierror !< error flag
  integer :: photcons_flag=0 !< photon conservation flag, non-zero if photon conservation is violated. This stops the simulation

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
  real(kind=dp) :: interm_zred !< calculated intermediate redshift (for restart)
#ifdef MPI
  integer :: mympierror
#endif

  ! Input file
  character(len=512) :: inputfile !< name of input file
  character(len=1) :: answer !< y or n answer

  ! Initialize clocks (cpu and wall)
  call setup_clocks

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (rank == 0) then
     write(logf,*) "input or input?"
     flush(logf)
     if (COMMAND_ARGUMENT_COUNT () > 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     else
        write(logf,*) "reading input from command line"
     endif
     flush(logf)
  endif

  ! Initialize output
  call setup_output ()

  ! Initialize grid
  call grid_ini ()

  if (rank == 0) &
       write(timefile,"(A,F8.1)") "Time after grid_ini: ",timestamp_wallclock ()

#ifdef MPILOG
  write(logf,*) "Before rad_ini"
#endif
  ! Initialize photo-ionization calculation
  call rad_ini ( )

  if (rank == 0) &
       write(timefile,"(A,F8.1)") "Time after rad_ini: ",timestamp_wallclock ()

#ifdef MPILOG
  write(logf,*) "Before material_ini"
#endif
  ! Initialize the material properties
  call material_ini (restart, nz0, ierror)
  startup_error = startup_error + ierror

  if (rank == 0) &
       write(timefile,"(A,F8.1)") "Time after material_ini: ",timestamp_wallclock ()

#ifdef MPILOG
  write(logf,*) "Before nbody_ini"
#endif
  ! Find the redshifts we are dealing with
  if (startup_error == 0) call nbody_ini (ierror)
  startup_error=startup_error+ierror

#ifdef MPILOG
  write(logf,*) "Before source_properties_ini"
#endif
  ! Initialize the source model
  if (startup_error == 0) call source_properties_ini ()

#ifdef MPILOG
  write(logf,*) "Before time_ini"
#endif
  ! Initialize time step parameters
  if (startup_error == 0) call time_ini ()
  if (rank == 0) flush(logf)

#ifdef MPILOG
  write(logf,*) "Before evolve_ini"
#endif
  ! Initialize evolve arrays
  if (startup_error == 0) call evolve_ini ()

  ! Set time to zero
  sim_time=0.0

  ! Initialize cosmology
  ! (always call bacause it initializes some of the things even if 
  ! not doing cosmo. 
  if (startup_error == 0) call cosmology_init(zred_array(nz0),sim_time)

  if (rank == 0) &
       write(timefile,"(A,F8.1)") "Time after cosmology_init: ", &
       timestamp_wallclock ()

  ! If a restart, inquire whether to restart from iteration
  if (startup_error == 0 .and. restart /= 0) then
     if (rank == 0) then
        if (.not.file_input) &
             write(*,"(A,$)") "Restart from iteration dump (y/n) or (0/1/2)? : "
        read(stdinput,*) answer
        write(logf,*) "restart from iteration answer: ",answer
        select case (answer)
        case("y", "Y") 
           ! Set flag, this is passed to evolve3d
           iter_restart=3
           write(logf,*) "Restarting from generic iteration dump."
        case("1") 
           ! Set flag, this is passed to evolve3d
           iter_restart=1
           write(logf,*) "Restarting from iteration dump 1."
        case("2") 
           ! Set flag, this is passed to evolve3d
           iter_restart=2
           write(logf,*) "Restarting from iteration dump 2."
        end select
     endif
#ifdef MPI
     call MPI_BCAST(iter_restart,1,MPI_INTEGER,0,MPI_COMM_NEW, &
          mympierror)
#endif
  endif

  ! MPI_Barrier
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif  

  ! If a restart, read ionization fractions from file
  if (restart == 1) then
     call xfrac_restart_init(zred_array(nz0))
     if (.not.isothermal) call temperature_restart_init(zred_array(nz0))
  endif
  if (restart == 2) then
     if (rank == 0) then
        if (.not.file_input) &
             write(*,"(A,$)") "At which redshift to restart x_frac?: "
        read(stdinput,*) zred_interm
     endif
#ifdef MPI
     call MPI_BCAST(zred_interm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
          mympierror)
#endif
     ! Check the value for consistency
     if (zred_interm > zred_array(nz0) .or. &
          zred_interm < zred_array(nz0+1) ) startup_error=5
     call xfrac_restart_init(zred_interm)
     if (.not.isothermal) call temperature_restart_init(zred_interm)
  end if
  
  if (startup_error == 0) then
     ! Load required clumping model parameters file (passing resolution)
     ! Models are customized on the simulation resolution, 
     ! you should provide parameters for the current simulation resolution
     call load_clumping_model(dr(1))
     
     ! Loop over redshifts
     do nz=nz0,NumZred-1
        
        if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
             "Time before starting redshift evolution step ",nz," :", &
             timestamp_wallclock ()
        
        zred=zred_array(nz)
        if (rank == 0) write(logf,*) "Doing redshift: ",zred," to ", &
             zred_array(nz+1),'of total', NumZred
        
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
        
        if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
             "Time after setting sources for step ",nz," :", &
             timestamp_wallclock ()
        
        ! Initialize density field
        call density_init(zred,nz)
        ! Set clumping and LLS in the case of position dependent values
        ! (read in the pre-computed grid values)
        if (type_of_clumping == 5) call set_clumping(zred)
        if (use_LLS .and. type_of_LLS == 2) call set_LLS(zred)
        
        if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
             "Time after setting material properties for step ",nz," :", &
             timestamp_wallclock ()
        
        ! Set time if restart at intermediate time
        ! Set next output time
        ! Note: this assumes that there are ALWAYS two time steps
        ! between redshifts. Should be made more general.
        if (restart >= 2) then
           interm_zred=time2zred(zred2time(zred)+dt)
           if (abs(interm_zred-zred_interm) < 0.001) then
              sim_time=zred2time(zred)+dt
           else
              sim_time=zred2time(zred_interm)
           endif
           next_output_time=end_time
        else
           next_output_time=sim_time+output_time
        endif
        
#ifdef MPILOG     
        write(logf,*) 'First output'
#endif 
        ! If start of simulation output (and not a restart).
        if (NumSrc > 0 .and. sim_time == 0.0 .and. restart == 0) &
             call output(time2zred(sim_time),sim_time,dt, &
             photcons_flag)
        
        ! Reset restart flag now that everything has been dealt with
        restart=0

#ifdef MPILOG     
        write(logf,*) 'Start of loop'
#endif 
        if (rank == 0) flush(logf)
        ! Loop until end time is reached
        do
           ! Report memory usage
           if (rank == 0) call report_memory(logf)
        
           ! Make sure you produce output at the correct time
           actual_dt=min(next_output_time-sim_time,dt)
           
           ! Report time and time step
           if (rank == 0) write(logf,"(A,2(es10.3,x),A)") "Time, dt:", &
                sim_time/YEAR,actual_dt/YEAR," (years)"
           
           ! For cosmological simulations evolve proper quantities
           ! We choose to scale the proper quantities to the halfway
           ! point between the start and the end of the integration
           ! time step (acutal_dt)
           if (cosmological) then
              call redshift_evol(sim_time+0.5*actual_dt)
              call cosmo_evol()
           endif
           
           ! Do not call in case of the pre-computed clumping data,
           ! the clumping grid should have been initialized above
           ! Same for LLS
           if (type_of_clumping /= 5) call set_clumping(zred)
           if (use_LLS .and. type_of_LLS /= 2) call set_LLS(zred)
           
           ! Take one time step
           if (NumSrc > 0) call evolve3D(sim_time,actual_dt,iter_restart)
           
           ! Reset flag for restart from iteration 
           ! (evolve3D is the last routine affected by this)
           iter_restart=0
           
           ! Update time
           sim_time=sim_time+actual_dt
           
           ! Write output
           if (abs(sim_time-next_output_time) <= 1e-6*sim_time) then
              !if (NumSrc > 0) call output(time2zred(sim_time),sim_time,actual_dt, &
              !photcons_flag)
              call output(time2zred(sim_time),sim_time,actual_dt, &
                   photcons_flag)
              next_output_time=next_output_time+output_time
              if (photcons_flag /= 0 .and. stop_on_photon_violation) then
                 if (rank == 0) write(logf,*) &
                      "Exiting because of photon conservation violation"
                 ! GM (110131): Forgot to check here whether we care about photon
                 ! conservations violations; if the code jumps out of the evolution
                 ! here funny things happen to the time step!
                 exit ! photon conservation violated
              endif
           endif
           ! end time for this redshift interval reached
           if (abs(sim_time-end_time) <= 1e-6*end_time) exit
           if (rank == 0) flush(logf)
        enddo
        
        if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
             "Time after finishing step ",nz," :", &
             timestamp_wallclock ()
        
        ! Get out: photon conservation violated
        if (stop_on_photon_violation .and. photcons_flag /= 0 .and. rank == 0) &
             write(logf,*) "Exiting because of photon conservation violation"
        if (stop_on_photon_violation .and. photcons_flag /= 0) exit ! photon conservation violated
        
        ! Scale to the current redshift
        if (cosmological) then
           call redshift_evol(sim_time)
           call cosmo_evol()
        endif
        
        ! Update clock counters (cpu + wall, to avoid overflowing the counter)
        call update_clocks ()
        
     enddo

     ! Write final output
     ! GM/110414: Bug, this statement previously overwrote the previous
     ! output (from the redshift list, so with multiple outputs per
     ! slice it could be the previous previous one). Fixed by writing
     ! the output for the last redshift: zred_array(NumZRed)
     if (photcons_flag == 0) call output(zred_array(NumZRed), &
          sim_time,actual_dt,photcons_flag)

  endif ! end of startup_error test

  ! End output streams
  call close_down ()
  
  if (rank == 0) write(timefile,"(A,F8.1)") &
       "Time at end of simulation: ", timestamp_wallclock ()

  ! Report clocks (cpu and wall)
  call report_clocks ()

  ! End the run
  call mpi_end ()

end Program C2Ray
