module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use grid, only: vol
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, MPC
  use cosmology_parameters, only: Omega_B,Omega0,rho_crit_0,h
  use nbody, only: nbody_type, id_str, dir_dens, boxsize
  use nbody, only: densityformat, densityheader, density_unit
  use abundances, only: mu, abu_he
  use c2ray_parameters, only: type_of_clumping,clumping_factor

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh(1),mesh(2),mesh(3))
  real(kind=dp) :: temper
  real(kind=dp) :: xh(mesh(1),mesh(2),mesh(3),0:1)
  logical isothermal
  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  public :: set_clumping, clumping_point

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ============================================================================

  subroutine mat_ini (restart, nz0, ierror)

    ! Initializes material properties on grid

    ! Authors: Ilian Iliev, Garrelt Mellema

    ! Date: 28-Dec-2007 (20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f)))

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
    real(kind=dp) :: temper_val
    character(len=1) :: answer

    ierror=0
    ! Check consistency with nbody module
    if (nbody_type /= "LG") then
       write(logf,*) "Error: wrong material module was compiled."
       write(logf,*) "       Expected LG, but got", &
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
          !print*,"Number of starting slice: ",nz0
       endif
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(restart,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(nz0,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
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

    ! Authors: Ilian Iliev, Garrelt Mellema

    ! Date: 26-Dec-2007 (20-Aug-2006 (19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004)))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    integer,intent(in) :: nz ! number in redshift list
    
    integer :: i,j,k,n,nfile ! loop counters
    character(len=512):: dens_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    character(len=3) :: id
    real(kind=dp) :: convert ! conversion factor
    real(kind=dp) :: summed_density
    real(kind=dp) :: avg_dens
    integer :: m1,m2,m3

    ! density in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: ndens_real

    if (rank == 0) then
       write(logf,*) "Reading LG input"
       write(logf,*) "Reading ",nz," input"

       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       write(id,'(i3.3)') nz
       if (id_str /= "dmdens_cic") then
          dens_file=trim(adjustl(dir_dens))// &
               trim(adjustl(id))//"rho_"//trim(adjustl(id_str))//".dat"
       else
          dens_file=trim(adjustl(dir_dens))// &
               trim(adjustl(id))//trim(adjustl(id_str))//".dat"
       endif
       
       ! Open density file: note that it is in "unformatted" form
       open(unit=20,file=dens_file,form=densityformat,status="old")
       write(logf,*) " from file ",trim(adjustl(dens_file))

       ! Read in data
       if (densityheader) then
          write(logf,*) " reading density header"
          read(20) m1,m2,m3
          write(logf,*) " density header: ",m1,m2,m3
          write(logf,*) " should match: ",mesh


          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with densities unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif

       ! Read in data and store it in ndens
       ! Allocate array needed to read in data
       allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
       do k=1,mesh(3)
          read(20) ((ndens_real(i,j,k),i=1,mesh(1)),j=1,mesh(2))
       end do
       write(logf,*) "Data read",m1,m2,m3
       ndens(:,:,:)=ndens_real(:,:,:)

       ! close file
       close(20)

       ! Deallocate array needed for reading in the data.
       deallocate(ndens_real)
 
    endif
#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
       
    if (rank == 0) write(logf,*) "M_box from data [M_sun]= ", &
         sum(ndens)*(boxsize/mesh(1))**3/h,mesh,h   

    ! The original values are in terms of mass density
    ! Below is the conversion factor to number density.
    ! This makes the density ndens, the TOTAL number density.

    ! The original values in terms of the mean density
    ! Below is the conversion factor.
    ! For LG the standard is "M0Mpc3", see LG.F90
    ! vol is redshift dependent, that is why we need to recalculate this
    ! M_particle and M_grid should be in g
    select case(density_unit)
    case ("M0Mpc3")
       convert=M_solar/Mpc**3*h**2*Omega_B/Omega0/(mu*m_p)*(1.0+zred_now)**3 
    end select

    ! Assign density to the grid
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             ndens(i,j,k)=ndens(i,j,k)*convert
             ! take care of empty cells, if any (0.001 M_solar/Mpc^3)
             if (ndens(i,j,k)<= 0.0) ndens(i,j,k)=0.001*convert
          enddo
       enddo
    enddo
    
    ! Find summed and average density
    avg_dens=sum(ndens)/(real(mesh(1),dp)*real(mesh(2),dp)*real(mesh(3),dp))
    
    ! Report density field properties
    if (rank == 0) then
       write(logf,*) "Raw density diagnostics (cm^-3)"
       write(logf,*) "minimum density: ",minval(ndens)
       write(logf,*) "maximum density: ",maxval(ndens)
       write(logf,"(A,1pe10.3,A)") "Average density = ",avg_dens," cm^-3"
       write(logf,"(A,1pe10.3,A)") "Theoretical value = ", &
            rho_crit_0*Omega_B/(mu*m_p)*(1.0+zred_now)**3, &
            " cm^-3" 
       write(logf,"(A,1pe10.3,A)") "(at z=0 : ", &
            rho_crit_0*Omega_B/(mu*m_p), &
            " cm^-3)"
       write(logf,"(A,1pe10.3)")"calculated at z=", zred_now	 
    endif

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
    
    character(len=512) :: xfrac_file
    character(len=6) :: zred_str
    integer :: m1,m2,m3
    ! Array needed to read in 4B reals
    real(kind=dp),dimension(:,:,:),allocatable :: xh1
    !real(kind=si),dimension(:,:,:),allocatable :: xh1_real

    if (rank == 0) then
       allocate(xh1(mesh(1),mesh(2),mesh(3)))
       write(zred_str,"(f6.3)") zred_now
       xfrac_file= trim(adjustl(results_dir))// &
            "xfrac3d_"//trim(adjustl(zred_str))//".bin"
       
       write(unit=logf,fmt="(2A)") "Reading ionization fractions from ", &
            trim(xfrac_file)
       ! Open ionization fractions file
       open(unit=20,file=xfrac_file,form="unformatted",status="old")
       
       ! Read in data
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with ionization fractions unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
       else
          read(20) xh1
          xh(:,:,:,1)=xh1(:,:,:)
          !xh(:,:,:,1)=real(xh1_real(:,:,:),dp)
          xh(:,:,:,0)=1.0-xh(:,:,:,1)
       endif

       ! close file
       close(20)
       deallocate(xh1)
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
    ! since no Gadget/LG clumping data is available

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
       write(logf,*) " this LG material module."
    endif

  end subroutine clumping_init

end module material
