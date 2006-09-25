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
  use my_mpi, only: mpi_setup, mpi_end
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini
  use pmfast, only: pmfast_ini, NumZred, zred_array
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol, &
       time2zred, zred
  use material, only: mat_ini, xfrac_ini, dens_ini
  use times, only: time_ini, set_timesteps
  use sourceprops, only: source_properties
  use evolve, only: evolve3D
  use subgrid_clumping, only: set_clumping

  implicit none

  ! Start and end time for CPU report
  real :: tstart,tend

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

  ! Initialize cpu timer
  call cpu_time(tstart)

  ! Set up MPI structure
  call mpi_setup()

  ! Initialize output
  call setup_output()

  !Initialize grid
  call grid_ini()

  ! Initialize photo-ionization calculation
  call rad_ini( )

  ! Initialize the material properties
  call mat_ini (restart)

  ! Find the redshifts we are dealing with
  call pmfast_ini ()

  ! Initialize time step parameters
  call time_ini ()

  ! Set time to zero
  time=0.0

  ! Initialize cosmology
  call cosmology_init(zred_array(1),time)

  !if restarts read ionization fractions from file
  !if (restart.ne.0) call xfrac_ini(zred_array(1))
  if (restart.eq.1) call xfrac_ini(zred_array(1))
  if (restart.eq.2)then
     write(*,'(A,$)') 'At which redshift to restart x_frac?:'
     read(*,*) zred_interm         
     call xfrac_ini(zred_interm)
  end if
  
  flag=1
  
  ! Loop over redshifts
  do nz=1,NumZred-1
     
     zred=zred_array(nz)
     write(30,*) 'Doing redshift: ',zred,' to ',zred_array(nz+1)
     
     ! Initialize time parameters
     call set_timesteps(zred,zred_array(nz+1), &
          end_time,dt,output_time)
     write(30,*) 'This is time ',time/YEAR,' to ',end_time/YEAR
         
     next_output_time=time+output_time

     ! Initialize source position
     call source_properties(zred,end_time-time)

     ! print*,'zred before dens_ini=',zred
     ! Initialize density field
     call dens_ini(zred)
     
     ! Loop until end time is reached
     do
        ! Make sure you produce output at the correct time
        actual_dt=min(next_output_time-time,dt)
        
        ! Report time and time step
        write(30,'(A,2(1pe10.3,x),A)') 'Time, dt:', &
             time/YEAR,actual_dt/YEAR,' (years)'

        ! For cosmological simulations evolve proper quantities
        if (cosmological) then
           call redshift_evol(time+0.5*actual_dt)
           call cosmo_evol()
        endif
        call set_clumping(real(1.0/(1.0+zred),4))
        ! print*,'check', restart,flag,time
        if(flag.eq.1 .and. restart .eq. 2)then !i.e. do not call the first time around
                                                !if restart is at middle point in time
           flag=2
        else
           ! Take one time step
           call evolve3D(actual_dt)
        end if

        ! Update time
        time=time+actual_dt
            
        ! Write output
        if (abs(time-next_output_time).le.1e-6*time) then
           call output(time2zred(time),time,dt)
           next_output_time=next_output_time+output_time
        endif

        if (abs(time-end_time).lt.1e-6*end_time) exit
     enddo

     ! stop
     ! Scale to the current redshift
     if (cosmological) then
        call redshift_evol(time)
        call cosmo_evol()
     endif
     
  enddo

  ! Write final output

  call output(zred,time,dt)

  call close_down ()
  
  ! Find out CPU time
  call cpu_time(tend)

  write(30,*) 'CPU time: ',tend-tstart,' s'

  ! End the run
  call mpi_end ()

end Program Ifront
