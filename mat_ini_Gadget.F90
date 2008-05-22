module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf
  use my_mpi
  use grid, only: vol
  use cgsconstants, only: m_p
  use cosmology_parameters, only: Omega_B,Omega0,rho_crit_0
  use nbody, only: nbody_type
  use abundances, only: mu, abu_he
  use c2ray_parameters, only: type_of_clumping,clumping_factor

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh(1),mesh(2),mesh(3))
  real(kind=dp) :: temper
  real(kind=dp) :: xh(mesh(1),mesh(2),mesh(3),0:1)
  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  public :: set_clumping, clumping_point
  logical isothermal

#ifdef MPI
  integer,private :: mympierror
#endif

contains
  ! ============================================================================
  subroutine mat_ini (restart, ierror)

    ! Initializes material properties on grid

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))

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
    integer,intent(out) :: ierror

    integer :: i,j,k,n ! loop counters
    real(kind=dp) :: temper_val
    character(len=1) :: answer

    ierror=0
    if (nbody_type /= "gadget") then
       write(logf,*) "Error: wrong material module was compiled."
       ierror=1
    else
       ! restart
       restart=0 ! no restart by default
       
       if (rank == 0) then
          ! Ask for temperature, restart. Read in values
          write(*,"(A,$)") "Enter initial temperature (K): "
          read(stdinput,*) temper_val
          write(*,"(A,$)") "Restart (y/n)? : "
          read(stdinput,*) answer
          if (answer.eq."y".or.answer.eq."Y") restart=1
          write(*,"(A,$)") "Restart at midpoint (y/n)? : "
          read(stdinput,*) answer
          if (answer.eq."y".or.answer.eq."Y") restart=2
       endif
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(restart,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
       
       ! Scalar version for constant temperature
       isothermal=.true.
       temper=temper_val
       
       ! Assign dummy density to the grid
       ! This should be overwritten later (in dens_ini)
       ndens(:,:,:)=1.0
       
       ! Assign ionization fractions (completely neutral)
       ! In case of a restart this will be overwritten in xfrac_ini
       xh(:,:,:,0)=1.0
       xh(:,:,:,1)=0.0
    endif

  end subroutine mat_ini

  ! ===========================================================================
  subroutine dens_ini (zred_now,nz)

    ! Initializes density on the grid (at redshift zred_now)

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    integer,intent(in) :: nz ! number in the list of redshifts (for
                             ! compatibility reasons)
    
    integer :: i,j,k,n,nfile ! loop counters
    character(len=512):: dens_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: convert ! conversion factor
    real(kind=dp) :: summed_density
    real(kind=dp) :: avg_dens

    ! density in file is in 4B reals, read in via this array
    real(kind=dp),dimension(:,:,:),allocatable :: ndens_real

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       ! Allocate array needed to read in data
       allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
       write(logf,*) "Reading Gadget input"
       dens_file=trim(adjustl(zred_str))// &
               "rho_gadget.dat"
          
       ! Open density file: note that it is in `binary" form
       open(unit=20,file=dens_file,form="binary",status="old")
       ! Read in data and store it in ndens
       read(20) ndens_real
       ndens(:,:,:)=ndens_real(:,:,:)
       
       ! close file
       close(20)
       
       ! Deallocate array needed for reading in the data.
       deallocate(ndens_real)
    endif
    !write(logf,*) "Distributing the density"
#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    !write(logf,*) "Density distributed"
       
    ! Report on data: min, max, total
    write(logf,*) "minimum: ",minval(ndens)
    write(logf,*) "maximum: ",maxval(ndens)
    write(logf,*) "summed density: ",sum(ndens)

    ! The original values are in terms of mass density
    ! Below is the conversion factor to number density.
    ! This makes the density ndens, the TOTAL number density.
    convert=1.0/(mu*m_p)

    ! Assign density to the grid
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             ndens(i,j,k)=ndens(i,j,k)*convert
             ! take care of empty cells (0.1 particles)
             if (ndens(i,j,k)<= 0.0) ndens(i,j,k)=0.1*convert
          enddo
       enddo
    enddo
    
    ! Find summed_and average density
    avg_dens=sum(ndens)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
    
    ! Report average density
    write(logf,"(A,1pe10.3,A)") "Average density = ",avg_dens," cm^-3"
    write(logf,"(A,1pe10.3,A)") "Theoretical value = ", &
         rho_crit_0*Omega_B/(mu*m_p)*(1.0+zred_now)**3, &
         " cm^-3" 
    write(logf,"(A,1pe10.3,A)") "(at z=0 : ", &
         rho_crit_0/(mu*m_p)*Omega_B, &
         " cm^-3)"

    ! Take only the H part of the particle density
    ndens=ndens*(1.0-abu_he)

  end subroutine dens_ini

  ! ===========================================================================
  subroutine xfrac_ini (zred_now)

    ! Initializes ionization fractions on the grid (at redshift zred_now).
    ! They are read from an "Ifront3" file which should have been created
    ! by an earlier run. This file is read in restart mode.

    ! Author: Garrelt Mellema

    ! Date: 19-May-2005

    real(kind=dp),intent(in) :: zred_now
    
    character(len=180) :: xfrac_file
    character(len=6) :: zred_str
    integer :: m1,m2,m3
    ! Array needed to read in 4B reals
    real(kind=si),dimension(:,:,:),allocatable :: xh1_real

    if (rank == 0) then
       allocate(xh1_real(mesh(1),mesh(2),mesh(3)))
       write(zred_str,"(f6.3)") zred_now
       xfrac_file= "./Ifront3_"//trim(adjustl(zred_str))//".bin"
       
       ! Open density file
       open(unit=20,file=xfrac_file,form="unformatted",status="old")
       
       ! Read in data
       read(20) m1,m2,m3
       if (m1.ne.mesh(1).or.m2.ne.mesh(2).or.m3.ne.mesh(3)) then
          write(*,*) "Warning: file with ionization fractions unusable"
       else
          read(20) xh1_real
          xh(:,:,:,1)=real(xh1_real(:,:,:),8)
          xh(:,:,:,0)=1.0-xh(:,:,:,1)
       endif

       ! close file
       close(20)
       deallocate(xh1_real)
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(xh,mesh(1)*mesh(2)*mesh(3)*2,MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine xfrac_ini

  ! ===========================================================================

  subroutine set_clumping(z)

    ! This routine initializes the clumping factor.
    ! It used position dependent clumping derived from PMFAST simulations
    ! since no Gadget clumping data is available

    real(kind=dp),intent(in) :: z

    select case (type_of_clumping)
    case(1)
       clumping = clumping_factor
    case(2) 
       clumping = 27.466*exp(-0.114*z+0.001328*z*z)
    case(3)
       clumping = 26.2917*exp(-0.1822*z+0.003505*z*z)
    case(4)
       clumping = 17.57*exp(-0.101*z+0.0011*z*z)
    case(5)
       ! This one should not be called!!
       call clumping_init (z)
    end select

    if (rank == 0) write(logf,*) "Setting (mean) global clumping factor to ", &
         clumping,"(type ", type_of_clumping,")"
    
  end subroutine set_clumping

  ! ===========================================================================

  subroutine clumping_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_clumping /= 5) then
       write(logf,*) "Error: using position dependent clumping, "
       write(logf,*) "       but array is not initialized."
    else
       clumping = clumping_grid(i,j,k)
    endif

  end subroutine clumping_point

  ! ===========================================================================

  subroutine clumping_init (zred_now)

    ! Initializes position dependent clumping (at redshift zred_now)

    real(kind=dp),intent(in) :: zred_now
    
    if (rank == 0) then
       write(logf,*) "Position dependent clumping is not available for"
       write(logf,*) " this Gadget material module."
    endif

  end subroutine clumping_init

end module material
