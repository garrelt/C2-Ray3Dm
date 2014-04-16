!>
!! \brief This module contains data and routines for handling the material properties on the grid (3D)
!!
!! These properties are; density, temperature, clumping, ionization fractions
!! 
!! \b Author: Garrelt Mellema
!!
!! \b Date: 27-Jun-2013 (20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))
!!
!! \b Version: cubep3m simulations.

module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use file_admin, only: stdinput, logf, file_input
  use my_mpi
  use nbody, only: nbody_type
  use density_module, only: density_array_init
  use ionfractions_module, only: xfrac_array_init
  use temperature_module, only: temperature_array_init
  use clumping_module, only: clumping_fit_file
  use lls_module, only: lls_init

  implicit none

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ============================================================================

  subroutine mat_ini (restart, nz0, ierror)

    ! Initializes material properties on grid

    ! Authors: Garrelt Mellema, Ilian Iliev

    ! Date: 30-Jan-2008 (20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f)))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file)
    ! - Initially completely neutral
    ! - Source points set elsewhere (multiple sources)

    ! History:
    ! - integrated with new cosmological approach developed for
    !    test 4.
    ! - adapted for multiple cosmological sources
    ! - f90 version with MPI

    integer,intent(out) :: restart ! will be /= 0 if a restart is intended
    integer,intent(out) :: nz0    ! nz0 is the number of the starting slice
    integer,intent(out) :: ierror ! will be /=0 if an error occurred

    integer :: i,j,k,n ! loop counters
    character(len=1) :: answer
    real(kind=dp) :: temper_val

    ! Set error flag to zero
    ierror=0
    ! Check consistency with nbody module
    if (nbody_type /= "cubep3m") then
       write(logf,*) "Error: wrong material module was compiled."
       write(logf,*) "       Expected cubep3m, but got", &
            trim(adjustl(nbody_type))
       ierror=1
    else
       ! restart
       restart=0 ! no restart by default

       if (rank == 0) then
          ! Ask for temperature, restart. Read in values
          if (.not.file_input) &
               write(*,"(A,$)") "Enter initial temperature (K): "
          read(stdinput,*) temper_val
          if (.not.file_input) write(*,"(A,$)") "Restart (y/n)? : "
          read(stdinput,*) answer
          if (answer == "y" .or. answer == "Y") then
             restart=1
             ! Apparently confusing if this question is only asked
             ! conditionally
             !if (.not.file_input) &
             !     write(*,"(A,$)") "Restart at midpoint (y/n)? : "
             !read(stdinput,*) answer
             !if (answer == "y" .or. answer == "Y") restart=2
          endif
          if (.not.file_input) &
               write(*,"(A,$)") "Restart at midpoint (y/n)? : "
          read(stdinput,*) answer
          if (answer == "y" .or. answer == "Y") restart=2
          if (.not.file_input) write(*,"(A,$)") "Number of starting slice: "
          read(stdinput,*) nz0
          if (.not.file_input) &
               write(*,"(A,$)") "Clumping fit table file: "
          read(stdinput,*) clumping_fit_file
          !print*,"Number of starting slice: ",nz0
       endif
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(restart,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(nz0,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
       
       ! Initialize density data
       call density_array_init (1.0_dp)
       
       ! Initialize temperature data
       call temperature_array_init(temper_val)

       ! Initialize ionization fractions
       call xfrac_array_init()

       ! Initialize LLS parametets
       call LLS_init ()

    endif

  end subroutine mat_ini

end module material
