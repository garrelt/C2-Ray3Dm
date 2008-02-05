Program Ifront

  ! Author: Garrelt Mellema

  ! Date: 23-Sep-2006

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
  ! pmfast : interface to PMFAST output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp
  use c2ray_parameters, only: cosmological
  use astroconstants, only: YEAR
  use my_mpi, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini
  use gadget, only: gadget_ini, NumZred, zred_array
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol, &
       time2zred, zred2time, zred
  use material, only: mat_ini, xfrac_ini, dens_ini
  use times, only: time_ini, set_timesteps
  use sourceprops, only: source_properties
  use evolve, only: evolve3D
  use subgrid_clumping, only: set_clumping
  use file_admin, only:stdinput,log

#ifdef XLF
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! Start and end time for CPU report
  real :: tstart,tend
  integer :: cntr1,cntr2,countspersec

  integer :: restart,nz,flag

  ! end_time - end time of the simulation (s)
  ! dt - time step (s)
  ! time - actual time (s)
  ! end_time - end time of calculation (s)
  ! output_time - time interval between outputs (s)
  ! next_output_time - time of next output (s)
  real(kind=dp) :: end_time,time,output_time,next_output_time
  real(kind=dp) :: dt,actual_dt
  real(kind=dp) :: zred_interm
  integer :: ntime
  ! Input file
  character(len=512) :: inputfile

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
        write(*,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
     endif
  endif

  ! Initialize output
  call setup_output()
  call flush(log)
  !Initialize grid
  call grid_ini()

  ! Initialize photo-ionization calculation
  call rad_ini( )

  ! Initialize the material properties
  call mat_ini (restart)

  ! Find the redshifts we are dealing with
  call gadget_ini ()

  ! Initialize time step parameters
  call time_ini ()

  call flush(log)
  ! Set time to zero
  time=0.0

  ! Initialize cosmology
  call cosmology_init(zred_array(1),time)

  ! For cosmological simulations evolve proper quantities
  if (cosmological) then
     call redshift_evol(time)
     call cosmo_evol()
  endif
  call set_clumping(real(1.0/(1.0+zred_array(1)),4))

  ! Loop over redshifts
  do nz=1,1
     
     zred=zred_array(nz)
     end_time=1e7*YEAR
     write(log,*) "Doing redshift: ",zred," to ",time2zred(time+end_time)
     
     ! Initialize time parameters
     write(log,*) "This is time ",time/YEAR," to ",end_time/YEAR
         
     ! Initialize source position
     call source_properties(zred,end_time-time)

     ! print*,"zred before dens_ini=",zred
     ! Initialize density field
     call dens_ini(zred)

     ! Loop until end time is reached
     !dt=end_time/real(100000)
     do ntime=1,1
        ! Make sure you produce output at the correct time
        !dt=10.0**(ntime+3)*YEAR
        !dt=end_time/10.0
        dt=1e7*YEAR
        if (ntime == 4) dt=0.5*dt 
        if (ntime == 5) dt=0.1*dt 
        actual_dt=dt-time
        !actual_dt=dt
        ! Report time and time step
        write(log,"(A,2(1pe10.3,x),A)") "Time, dt:", &
             time/YEAR,actual_dt/YEAR," (years)"

        ! Take one time step
        call evolve3D(actual_dt)

        ! Update time
        time=time+actual_dt
            
        ! Write output
        call output(time2zred(time),time,dt)
        if (abs(time-end_time).lt.1e-6*end_time) exit
        call flush(log)
     enddo
  enddo
  !
  call close_down ()
  
  ! Find out CPU time
  call cpu_time(tend)
  call system_clock(cntr2,countspersec)

  write(log,*) "CPU time: ",tend-tstart," s"
  write(log,*) "Wall clock time: ",(cntr2-cntr1)/countspersec," s"

  ! End the run
  call mpi_end ()

end Program Ifront
