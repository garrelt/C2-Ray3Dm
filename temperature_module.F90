!>
!! \brief This module contains data and routines for handling the temperature 
!! properties on the grid (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 15-Apr-2014
!!
!! \b Version: 

module temperature_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use c2ray_parameters, only: isothermal

implicit none

  type temperature_states
     real(kind=si) :: current
     real(kind=si) :: average
     real(kind=si) :: intermed
  end type temperature_states

  type temperature_states_dbl
     real(kind=dp) :: current
     real(kind=dp) :: average
     real(kind=dp) :: intermed
  end type temperature_states_dbl

  real(kind=dp) :: temper
  real(kind=dp) :: temper_val
  type(temperature_states),dimension(:,:,:),allocatable :: temperature_grid

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  subroutine temperature_array_init (temper_val_input)
    
    real(kind=dp),intent(in) :: temper_val_input

    ! Initialize global variable temper_val
    temper_val=temper_val_input

    ! Allocate temperature array and initialize if the run is not
    ! isothermal
    if (.not.isothermal) then
       allocate(temperature_grid(mesh(1),mesh(2),mesh(3)))
       temperature_grid(:,:,:)%current=temper_val_input
       temperature_grid(:,:,:)%average=temper_val_input
       temperature_grid(:,:,:)%intermed=temper_val_input
    else
       temper=temper_val_input
    endif

    ! Report on temperature situation
    if (rank == 0) then
       if (isothermal) then
          write(logf,"(A)") "Thermal conditions: isothermal"
       else
          write(logf,"(A)") &
               "Thermal conditions: applying heating and cooling"
       endif
    endif

  end subroutine temperature_array_init

  ! ===========================================================================

  subroutine temperature_restart_init (zred_now)

    ! Initializes temperature on the grid (at redshift zred_now).
    ! They are read from an "Temper3D" file which should have been created
    ! by an earlier run. This file is read when in restart mode.

    ! Author: Garrelt Mellema

    ! Date: 02-June-2011

    real(kind=dp),intent(in) :: zred_now
    
    character(len=512) :: temper_file
    character(len=6) :: zred_str
    integer :: m1,m2,m3

    if (isothermal) then
       if (rank == 0) write(logf,"(A)") &
            "Incorrect call to temper_ini in isothermal case"
    else
       if (rank == 0) then
          write(zred_str,"(f6.3)") zred_now
          temper_file= trim(adjustl(results_dir))// &
               "Temper3D_"//trim(adjustl(zred_str))//".bin"
          
          write(unit=logf,fmt="(2A)") "Reading temperature from ", &
               trim(temper_file)
          ! Open temperature file
          open(unit=20,file=temper_file,form="unformatted",status="old")
          
          ! Read in data
          read(20) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "WARNING: file with temperatures unusable, as"
             write(logf,*) "mesh found in file: ",m1,m2,m3
          else
             read(20) temperature_grid%current
          endif

          ! Fill the other parts of the temperature grid array
          ! See evolve for their use
          temperature_grid(:,:,:)%average=temperature_grid(:,:,:)%current
          temperature_grid(:,:,:)%intermed=temperature_grid(:,:,:)%current
          
          ! close file
          close(20)
       endif
       
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temperature_grid,mesh(1)*mesh(2)*mesh(3)*3,MPI_REAL,0,&
            MPI_COMM_NEW,mympierror)
#endif
    endif

  end subroutine temperature_restart_init

  ! ===========================================================================

  subroutine get_temperature_point (i,j,k,temperature_point)

    ! Puts value of temperature (from grid or initial condition value)
    ! in the module variable temper

    integer,intent(in) :: i,j,k
    type(temperature_states_dbl),intent(out) :: temperature_point

    if (isothermal) then
       temperature_point%current = dble(temper_val)
       temperature_point%average = dble(temper_val)
       temperature_point%intermed =dble(temper_val)
    else
       temperature_point%current = dble(temperature_grid(i,j,k)%current)
       temperature_point%average = dble(temperature_grid(i,j,k)%average)   
       temperature_point%intermed = dble(temperature_grid(i,j,k)%intermed)
    endif

  end subroutine get_temperature_point

  ! ===========================================================================

  subroutine set_temperature_point (i,j,k,temperature_point)
    
    ! Puts value of module variable temper back in temperature grid
    ! (if running not isothermal)

    integer,intent(in) :: i,j,k
    type(temperature_states_dbl),intent(in) :: temperature_point
    
    if (.not.isothermal) then
       temperature_grid(i,j,k)%intermed=real(temperature_point%intermed)
       temperature_grid(i,j,k)%average=real(temperature_point%average)
    endif
    
  end subroutine set_temperature_point

  ! ===========================================================================

  subroutine set_final_temperature_point ()
    
    ! Puts value of module variable temper back in temperature grid
    ! (if running not isothermal)

    !integer,intent(in) :: i,j,k
    !real(kind=dp),intent(in) :: temper
    
    if (.not.isothermal) temperature_grid(:,:,:)%current=temperature_grid(:,:,:)%intermed
    
  end subroutine set_final_temperature_point  

  ! ===========================================================================

end module temperature_module
