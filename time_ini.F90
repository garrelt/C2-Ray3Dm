!>
!! \brief This module handles the time variables
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 20-Aug-2006
!<
module times

  ! This module handles the time variables

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, file_input
  use astroconstants, only: YEAR
  use cosmology, only: t0, zred_t0, zred2time

  implicit none

  integer :: number_timesteps !< Number of time steps between redshift slices
  integer :: number_outputs !< Number of outputs between redshift slices
  real(kind=dp) :: dt_nbody
  
contains

  ! =======================================================================

  !>
  !!  Initializes the module variables
  !<
  subroutine time_ini( )
    
    ! Initializes number of time steps per frame (integration and output)

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 10-Mar-2004)

#ifdef MPI
    integer :: ierror
#endif

    if (rank == 0) then
       ! Ask for number of time steps
       if (.not.file_input) &
            write(*,'(A,$)') 'Enter number of time steps between slices: '
       read(stdinput,*) number_timesteps

       ! Ask for interval between outputs
       if (.not.file_input) &
            write(*,'(A,$)') 'Enter number of outputs between slices: '
       read(stdinput,*) number_outputs
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(number_timesteps,1,MPI_INTEGER,0,&
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(number_outputs,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
    
  end subroutine time_ini

  !=======================================================================

  !>
  !!  Sets the time steps for the interval zred0-zred_end
  !!  The time steps are constant, given by dividing the redshift
  !! range into equal time intervals (number_timesteps)
  !<
  subroutine set_timesteps (zred0,zred_end,end_time,dt,output_dt)

    ! Sets time steps for calculation and outputting

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 19-May-2005)

    real(kind=dp),intent(in) :: zred0 !< starting redshift
    real(kind=dp),intent(in) :: zred_end !< ending redshift
    real(kind=dp),intent(out) :: end_time !< end time
    real(kind=dp),intent(out) :: dt !< time step
    real(kind=dp),intent(out) :: output_dt !< output time step

    real(kind=dp) :: current_time

    ! Convert to time (in seconds)
    current_time=zred2time(zred0) !t0*( ( (1.0+zred_t0)/(1.0+zred0) )**1.5-1.0)
    end_time=zred2time(zred_end) !t0*( ( (1.0+zred_t0)/(1.0+zred_end) )**1.5-1.0)

    ! Set value of time interval between two nbody outputs (density fields
    ! and source lists
    dt_nbody=end_time-current_time

    ! Set value of time step
    dt=dt_nbody/real(number_timesteps)

    ! Convert to time 
    output_dt=dt_nbody/real(number_outputs)

  end subroutine set_timesteps

end module times
