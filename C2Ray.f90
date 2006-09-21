Program Ifront

  ! Author: Garrelt Mellema

  ! Date: 20-Aug-2006

  ! Goal:
  ! Calculates the evolution in time of an ionization front,
  ! trying to make sure photons are conserved.
  
  ! Does not include hydrodynamics
  ! Assumes constant time step

  ! Needs following `modules'
  ! output : output routines
  ! grid_ini : initialize grid
  ! rad_ini : initialize radiative parameters
  ! mat_ini : initialize material properties
  ! time_ini : initialize time properties
  ! evolve : evolve grid in time

  implicit none

  use my_mpi
  use output
  use grid
  use radiation

  integer :: restart,nz,flag

  ! end_time - end time of the simulation (s)
  ! dt - time step (s)
  ! time - actual time (s)
  ! end_time - end time of calculation (s)
  ! output_time - time interval between outputs (s)
  ! next_output_time - time of next output (s)
  real(kind=8) :: end_time,time,output_time,next_output_time
  real(kind=8) :: dt,actual_dt
  real(kind=8) :: zred_interm
  
  external time2zred
  real(kind=8) :: time2zred

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
  call cosmo_ini(zred_array(1),time)

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
     write(*,*) 'Doing redshift: ',zred,' to ',zred_array(nz+1)
     
     ! Initialize time parameters
     call set_timesteps(zred,zred_array(nz+1), &
          end_time,dt,output_time)
     write(*,*) 'This is time ',time/YEAR,' to ',end_time/YEAR
         
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
        write(*,'(A,2(1pe10.3,x),A)') 'Time, dt:', &
             time/YEAR,actual_dt/YEAR,' (years)'

        ! For cosmological simulations evolve proper quantities
        if (cosmological) then
           call redshift_evol(time+0.5*actual_dt)
           call cosmo_evol()
        endif
            
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
  
  ! End the run
  call mpi_end()

end Program Ifront
