!>
!! \brief This module contains variables and routines for dealing with 
!!   cpu and wall clocks
!!
!! This module uses calls to the Fortran 2003 internal routines cpu_time
!! and system_clock. To get the time passed between two events, one has
!! to call these twice, once at the start and once at the end. 
!! To avoid overflowing the variable (which can happen for codes that run
!! for many hours or days), this module provides "update" routines that
!! update hour and minute counters. These should be called at regular
!! intervals. The module also comes with its own reporting routines.
!!
!! General module
!!
!! Dependencies:
!!   - file_admin module (for the unit of the log file where to write the report)
!!   - my_mpi (for the value of rank, the MPI id of the current processor, reporting is only done for the rank 0 processor)
!!
!! Author: Garrelt Mellema
!!
!! Date: 2010-03-08
!!
!! Version: 1.0
!!
!<

module clocks

  use file_admin, only: logf
  use my_mpi, only: rank

  implicit none

  ! Start and end time for CPU report
  real :: cputime1 !< Start time for call to CPU routine
  real :: cputime2 !< End time for call to CPU routine

  integer :: cpu_hours=0 !< accumulates the CPU hours
  integer :: cpu_minutes=0 !< accumulates the CPU minutes
  real :: cpu_seconds=0.0 !< accumulates the CPU seconds
  
  ! Wall clock time variables
  integer :: cntr1 !< Start time for call to wall clock routine
  integer :: cntr2 !< End time for call to wall clock routine

  integer :: countspersec !< counts per second (for wall clock time)

  integer :: clock_hours=0 !< accumulates the wall clock hours
  integer :: clock_minutes=0 !< accumulates the wall clock minutes
  real :: clock_seconds=0.0 !< accumulates the wall clock seconds
  
contains

  !=======================================================================

  !> Sets up all the clocks (initialization routine)
  subroutine setup_clocks
    
    call setup_cpuclock()
    call setup_wallclock()
    
  end subroutine setup_clocks
  
  !=======================================================================

  !> Sets up cpu clock (initialization routine)
  subroutine setup_cpuclock
    
    ! Initialize cpu timer
    call cpu_time(cputime1)
    
  end subroutine setup_cpuclock
  
  !=======================================================================

  !> Sets up wall clock (initialization routine)
  subroutine setup_wallclock
    
    ! Initialize wall cock timer
    call system_clock(cntr1)
    
  end subroutine setup_wallclock
  
  !=======================================================================

  !> Updates all the clocks
  subroutine update_clocks
    
    call update_cpuclock
    call update_wallclock
    
  end subroutine update_clocks
  
  !=======================================================================

  !> Updates CPU clock
  subroutine update_cpuclock
    
    ! Find out intermediate CPU time (to avoid overflowing the counter)
    call cpu_time(cputime2)
    cpu_seconds=cpu_seconds+real(cputime2-cputime1)
    cputime1=cputime2
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )
    
  end subroutine update_cpuclock
  
  !=======================================================================

  !> Updates wall clock
  subroutine update_wallclock
    
    call system_clock(cntr2,countspersec)
    clock_seconds=clock_seconds+real(cntr2-cntr1)/real(countspersec)
    cntr1=cntr2
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )
    
  end subroutine update_wallclock
  
  !=======================================================================

  !> Reports all the clocks
  subroutine report_clocks
    
    call report_cpuclock
    call report_wallclock
    
  end subroutine report_clocks
  
  !=======================================================================

  !> Reports CPU clock to log file (unit logf)
  subroutine report_cpuclock
    
    call update_cpuclock ()
    if (rank == 0) then
       write(logf,*) "CPU time: ",cpu_hours,' hours',cpu_minutes,' minutes', &
            cpu_seconds,' seconds.'
    endif

  end subroutine report_cpuclock
  
  !=======================================================================

  !> Reports wall clock to log file (unit logf)
  subroutine report_wallclock
    
    call update_wallclock ()
    if (rank == 0) then
       write(logf,*) "Wall clock time: ",clock_hours,' hours', &
            clock_minutes,' minutes',clock_seconds,' seconds.'
    endif

  end subroutine report_wallclock
  
end module clocks
