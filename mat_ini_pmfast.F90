!>
!! \brief This module contains data and routines for handling the material properties on the grid (3D)
!!
!! These properties are; density, temperature, clumping, ionization fractions
!! 
!! \b Author: Garrelt Mellema
!!
!! \b Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))
!!
!! \b Version: PMFAST simulations

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
  use astroconstants, only: M_solar, Mpc
  use cosmology_parameters, only: Omega_B, Omega0, rho_crit_0, h
  use nbody, only: nbody_type, M_grid, M_particle, id_str, dir_dens, NumZred, Zred_array
  use nbody, only: densityformat, densityheader, clumpingformat, clumpingheader, density_unit, tot_nfiles
  use nbody, only: density_convert_particle, density_convert_grid
  use abundances, only: mu
  use c2ray_parameters, only: type_of_clumping,clumping_factor

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=si) :: ndens(mesh(1),mesh(2),mesh(3)) !< number density (cm^-3) of a cell
  real(kind=dp) :: temper !< temperature (K) of a cell 
  real(kind=dp) :: xh(mesh(1),mesh(2),mesh(3),0:1) !< ionization fractions for one cell
  logical isothermal !< isothermal?
  real,public :: clumping !< clumping factor of cell
  real,dimension(:,:,:),allocatable :: clumping_grid !< grid of clumping factors
  public :: set_clumping, clumping_point

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ============================================================================

  !> Initializes material properties on grid
  subroutine mat_ini (restart, nz0, ierror)

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
    integer,intent(out) :: nz0    ! nz0 is the number of the starting slice
    integer,intent(out) :: ierror ! will be /=0 if an error occurred

    integer :: i,j,k,n ! loop counters
    real(kind=dp) :: temper_val
    character(len=1) :: answer

    ierror=0
    ! Check consistency with nbody module
    if (nbody_type /= "pmfast") then
       write(logf,*) "Error: wrong material module was compiled."
       ierror=1
       write(logf,*) "       Expected cubep3m, but got", &
            trim(adjustl(nbody_type))
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

  !> Initializes density on the grid (at redshift zred_now)
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
    integer :: m1,m2,m3

    ! density in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: ndens_real

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       if (id_str /= "coarse") then
          dens_file=trim(adjustl(dir_dens))// &
               trim(adjustl(zred_str))// &
               "rho_"//trim(adjustl(id_str))//".dat"
          write(unit=logf,fmt="(4A)") "Reading ",id_str, &
               " density input from ",trim(dens_file)
          
          ! Open density file: note that it is in `binary" form
          open(unit=20,file=dens_file,form=densityformat,status="old")
          if (densityheader) then
             read(20) m1,m2,m3
             if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
                write(logf,*) "Warning: file with densities unusable"
                write(logf,*) "mesh found in file: ",m1,m2,m3
                stop
             endif
          endif
          ! Allocate array needed to read in data
          allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
          ! Read in data and store it in ndens
          read(20) ndens_real
          ndens(:,:,:)=ndens_real(:,:,:)
          
          ! close file
          close(20)

          ! Deallocate array needed for reading in the data.
          deallocate(ndens_real)
       else
          ! For the highest resolution the density is spread out
          ! over tot_nfiles. Otherwise the same.
          do nfile=0,tot_nfiles-1
             write(nfile_str,"(I1)") nfile
             dens_file=trim(adjustl(dir_dens))// &
                  trim(adjustl(zred_str))// &
                  "rho_c"//nfile_str//".dat"
             open(unit=20,file=dens_file,form=densityformat,status="old")
             if (densityheader) then
                read(20) m1,m2,m3
                if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
                   write(logf,*) "Warning: file with densities unusable"
                   write(logf,*) "mesh found in file: ",m1,m2,m3
                   stop
                endif
             endif
             ! Allocate array needed to read in data
             allocate(ndens_real(mesh(1),mesh(2),mesh(3)/tot_nfiles))
             read(20) ndens_real
             ndens(:,:, &
                  1+nfile*(mesh(3)/tot_nfiles): &
                  (1+nfile)*(mesh(3)/tot_nfiles))=ndens_real(:,:,:)
             close(20)
          enddo
          deallocate(ndens_real)
       endif
    endif
#ifdef MPI       
    ! Distribute the density to the other nodes
    !call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
    call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_REAL,0,&
         MPI_COMM_NEW,mympierror)
#endif
       
    ! The original values in terms of the mean density
    ! Below is the conversion factor
    ! vol is redshift dependent, that is why we need to recalculate this
    ! M_particle and M_grid should be in g
    select case(density_unit)
    case ("grid")
       convert=density_convert_grid*(1.0+zred_now)**3 
    case ("particle")
       convert=density_convert_particle*(1.0+zred_now)**3 
    case ("M0Mpc3")
       convert=M_solar/Mpc**3*h**2*Omega_B/Omega0/(mu*m_p)*(1.0+zred_now)**3 
    end select

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
    endif

  end subroutine dens_ini

  ! ===========================================================================

  !> Initializes ionization fractions on the grid (at redshift zred_now).
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
    real(kind=dp),dimension(:,:,:),allocatable :: xh1_real

    if (rank == 0) then
       allocate(xh1_real(mesh(1),mesh(2),mesh(3)))
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
          read(20) xh1_real
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          !xh(:,:,:,1)=real(xh1_real(:,:,:),dp)
          xh(:,:,:,1)=xh1_real(:,:,:)
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

  !> Initialize clumping factor
  subroutine set_clumping(z)

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
       call clumping_init (z)
    end select

    if (rank == 0) write(logf,*) "Setting (mean) global clumping factor to ", &
         clumping,"(type ", type_of_clumping,")"
    
  end subroutine set_clumping

  ! ===========================================================================

  !> set clumping factor for current cell
  subroutine clumping_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_clumping /= 5) then
       write(logf,*) &
         "Error: calling position dependent, but array is not initialized."
    else
       clumping = clumping_grid(i,j,k)
    endif

  end subroutine clumping_point

  ! ===========================================================================

  !> Initializes position dependent clumping (at redshift zred_now)
  subroutine clumping_init (zred_now)

    ! Initializes position dependent clumping (at redshift zred_now)

    ! Author: Ilian Iliev, modified from code by Garrelt Mellema

    ! Date: 03-Mar-2008 (5-Mar-2007( 20-Aug-2006, 19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004))

    ! Version: 
    ! - Three-dimensional. 
    ! - gas clumping field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    
    integer :: i,j,k,n,nfile ! loop counters
    integer :: m1,m2,m3 ! size of mesh in clumping file (header)
    character(len=512):: clump_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: summed_density
    real(kind=dp) :: avg_dens

    ! clumping in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(clumping_grid))) &
         allocate(clumping_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       if (id_str == "coarsest") then
          ! Allocate array needed to read in data
          allocate(clumping_real(mesh(1),mesh(2),mesh(3)))
          write(30,*) "Reading ",id_str," input"
          clump_file=trim(adjustl(dir_dens))// &
               trim(adjustl(zred_str))// &
               "clump_"//trim(adjustl(id_str))//".dat"
          
          ! Open clumping file: note that it is in `binary" form
          open(unit=20,file=clump_file,form=clumpingformat,status="old")
          
          ! Read in data and store it in clumping_grid
          if (clumpingheader) then
             read(20) m1,m2,m3
             ! disabled because of errors in clumping files
             ! GM 090324
             !if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             !   write(logf,*) "Warning: file with clumping factors unusable"
             !   write(logf,*) "mesh found in file: ",m1,m2,m3
             !   stop
             !endif
          endif
          read(20) clumping_real
          clumping_grid(:,:,:)=clumping_real(:,:,:)
          
          ! close file
          close(20)

          ! Deallocate array needed for reading in the data.
          deallocate(clumping_real)

       else if (id_str == "coarser") then
          ! For the highest resolution the density is spread out
          ! over tot_nfiles. Otherwise the same.
          allocate(clumping_real(mesh(1),mesh(2),mesh(3)/tot_nfiles))
          do nfile=0,tot_nfiles-1
             write(nfile_str,"(I1)") nfile
             clump_file=trim(adjustl(dir_dens))// &
                  trim(adjustl(zred_str))// &
                  "clump"//nfile_str//".dat"
             write(30,*) clump_file
             open(unit=20,file=clump_file,form=clumpingformat,status="old")
             if (clumpingheader) then
                read(20) m1,m2,m3
                ! disabled because of errors in clumping files
                ! GM 090324
                !if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
                !   write(logf,*) "Warning: file with clumping factors unusable"
                !   write(logf,*) "mesh found in file: ",m1,m2,m3
                !   stop
                !endif
             endif
             read(20) clumping_real
             clumping_grid(:,:, &
                  1+nfile*(mesh(3)/tot_nfiles): &
                  (1+nfile)*(mesh(3)/tot_nfiles))=clumping_real(:,:,:)
             close(20)
          enddo
          deallocate(clumping_real)
       else!812^3 run
          print*,"Cannot do that, data unavailable"!for 812^3 runs just use the 406^3 values in each 2x2x2? 
       endif
    endif

#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(logf,*) "minimum: ",minval(clumping_grid)
       write(logf,*) "maximum: ",maxval(clumping_grid)
       write(logf,*) "mean clumping: ",clumping
    endif

  end subroutine clumping_init

end module material
