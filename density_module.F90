module density_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use cgsconstants, only: m_p
  use astroconstants, only: M_solar, Mpc
  use c2ray_parameters, only: cosmological
  use cosmology_parameters, only: Omega_B, Omega0, rho_crit_0, h, H0
  use abundances, only: mu
  use nbody, only: nbody_type, M_grid, M_particle, id_str, dir_dens
  use nbody, only: NumZred, Zred_array
  use nbody, only: densityformat, densityaccess, densityheader
  use nbody, only: density_convert_particle, density_convert_grid, density_unit

  implicit none
  
  ! ndens - number density (cm^-3) of a cell
  ! SINGLE PRECISION! Be careful when passing this as argument to
  ! functions and subroutines.
  real(kind=si),dimension(:,:,:),allocatable :: ndens
  real(kind=dp) :: avg_dens !< average density

#ifdef MPI
  integer,private :: mympierror
#endif

contains
  
  ! ============================================================================

  subroutine density_array_init (density_value_input)

    real(kind=dp),intent(in) :: density_value_input

    ! Allocate density array
    allocate(ndens(mesh(1),mesh(2),mesh(3)))

    ! Assign dummy density to the grid
    ! This should be overwritten later (in dens_ini)
    ndens(:,:,:)=1.0
    
  end subroutine density_array_init

  ! ===========================================================================

  subroutine density_init (redshift,nz)

    ! Initializes density on the grid (at redshift zred_now)

    ! Authors: Garrelt Mellema, Ilian Iliev

    ! Date: 30-Jan-2008 (20-Aug-2006 (19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004)))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: redshift
    integer,intent(in) :: nz ! number in the list of redshifts (for
                             ! compatibility reasons)
    
    integer :: i,j,k,n ! loop counters
    real(kind=dp) :: summed_density
    integer :: m1,m2,m3

    ! density in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: ndens_real

    call set_density(redshift,nz)

    call density_diagnostics(redshift)

  end subroutine density_init

  ! ============================================================================

  subroutine set_density(redshift,nz)

    real(kind=dp),intent(in) :: redshift
    integer,intent(in) :: nz
    
    character(len=512):: dens_file

    select case (nbody_type)
       ! test problem: constant average density
       ! if cosmological: corresponding to current redshift
       ! if not: choose the first redshift (non-changing density)
    case("test") 
       if (cosmological) then
          call set_constant_average_density(redshift)
       else
          call set_constant_average_density(Zred_array(1))
       endif

       ! cubep3m and LG: read density field from file
    case("cubep3m","LG") 
       if (rank == 0) then
          
          ! construct filename
          dens_file=construct_densfilename(redshift,nz)
          
          ! read density from file
          call read_density_file(dens_file)
          
       endif
       
#ifdef MPI       
       ! Distribute the density to the other nodes
       !call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
       call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_REAL,0,&
            MPI_COMM_NEW,mympierror)
       call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
       
       ! scale the density to correct units
       call scale_density(redshift)

    end select
    
  end subroutine set_density

  ! ============================================================================

  subroutine set_constant_average_density(redshift)

    real(kind=dp),intent(in) :: redshift

    integer :: i,j,k

    ! Calculate average density of atoms (H+He)
    avg_dens=rho_crit_0*Omega_B/(mu*m_p)*(1.0+redshift)**3

    ! Assign density to the grid (average density at this redshift)
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             ndens(i,j,k)=avg_dens
          enddo
       enddo
    enddo
    
  end subroutine set_constant_average_density

  ! ============================================================================

  function construct_densfilename(redshift,number_of_redshift)
    
    real(kind=dp),intent(in) :: redshift
    integer,intent(in) :: number_of_redshift
    character(len=512) :: construct_densfilename

    character(len=6) :: redshift_str
    character(len=3) :: number_of_redshift_str
    character(len=512) :: dens_file

    select case (nbody_type)
    case("cubep3m")
       write(redshift_str,"(f6.3)") redshift
       dens_file=trim(adjustl(dir_dens))// &
            trim(adjustl(redshift_str))// &
            "n_all.dat"
    case("pmfast")
       ! This case is more complicated in reality. If id_str="coarse"
       ! the density data is spread out over several files which need
       ! to be read sequentially. This is not implemented here as it
       ! is unclear whether this will ever be needed.
       write(redshift_str,"(f6.3)") redshift
       dens_file=trim(adjustl(dir_dens))// &
            trim(adjustl(redshift_str))// &
            "rho_"//trim(adjustl(id_str))//".dat"
    case("LG")
       write(redshift_str,"(f6.3)") redshift
       write(number_of_redshift_str,'(i3.3)') number_of_redshift
       if (id_str /= "dmdens_cic") then
          dens_file=trim(adjustl(dir_dens))// &
               trim(adjustl(number_of_redshift_str))//"rho_"//trim(adjustl(id_str))//".dat"
       else
          dens_file=trim(adjustl(dir_dens))// &
               trim(adjustl(number_of_redshift_str))//trim(adjustl(id_str))//".dat"
       endif
    case("gadget")
       write(redshift_str,"(f6.3)") redshift
       dens_file=trim(adjustl(dir_dens))//trim(adjustl(redshift_str))// &
               "rho_gadget.dat"
    end select

    ! Report to log file
    write(unit=logf,fmt="(4A)") "Reading ",id_str, &
         " density input from ",trim(dens_file)

    ! Copy result to function variable
    construct_densfilename=dens_file
    
  end function construct_densfilename

  ! ============================================================================

  subroutine read_density_file (density_file)

    character(len=512),intent(in):: density_file

    integer :: m1,m2,m3
    integer :: i,j,k

    ! Open density file: note that it is in `binary" form
    open(unit=20,file=density_file,form=densityformat, &
         access=densityaccess,status="old")
       
    
    ! Read in data header
    if (densityheader) then
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with densities unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
          stop
       endif
    endif
    
    ! Read in data and store it in ndens
    ! Allocate array needed to read in data
    ! GM/110308: Disabled now that ndens is single precision
    !allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
    !read(20) ndens_real
    !ndens(:,:,:)=ndens_real(:,:,:)
    select case (nbody_type)
    case("cubep3m","pmfast","gadget")
       read(20) ndens
    case("LG")
       do k=1,mesh(3)
          read(20) ((ndens(i,j,k),i=1,mesh(1)),j=1,mesh(2))
       enddo
    end select
       
    ! close file
    close(20)

  end subroutine read_density_file

  ! ============================================================================

  subroutine scale_density (redshift)
    
    real(kind=dp),intent(in) :: redshift

    real(kind=dp) :: convert ! conversion factor
    integer :: i,j,k

    ! The original values in terms of the mean density
    ! Below is the conversion factor
    ! vol is redshift dependent, that is why we need to recalculate this
    ! M_particle and M_grid should be in g
    select case(density_unit)
    case ("grid")
       convert=density_convert_grid
    case ("particle")
       convert=density_convert_particle
    case ("M0Mpc3")
       convert=M_solar/Mpc**3*h**2*Omega_B/Omega0/(mu*m_p)
    case("mass_density")
       convert=1.0/(mu*m_p)
    end select
    
    ! Check separately for cosmological. We need to set comoving values here
    ! to initialize correctly.
    if (cosmological) then
       convert=convert*(1.0+redshift)**3
    endif
    
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
    
  end subroutine scale_density
  
  ! ============================================================================

  subroutine density_diagnostics (redshift)
    
    real(kind=dp),intent(in) :: redshift

    real(kind=si),dimension(:),allocatable :: partial_sum_ndens
    integer :: k

    ! Find summed and average density
    ! GM/110208: Make this a partial sum (since we are dealing
    ! with bigger meshes now.
    allocate(partial_sum_ndens(mesh(3)))
    do k=1,mesh(3)
       partial_sum_ndens(k)=sum(ndens(:,:,k))
    enddo
    avg_dens=sum(partial_sum_ndens)/ &
         (real(mesh(1),dp)*real(mesh(2),dp)*real(mesh(3),dp))
    deallocate(partial_sum_ndens)
    
    ! Report density field properties
    if (rank == 0) then
       write(logf,*) "Raw density diagnostics (cm^-3)"
       write(logf,*) "minimum density: ",minval(ndens)
       write(logf,*) "maximum density: ",maxval(ndens)
       write(logf,"(A,es10.3,A)") "Average density = ",avg_dens," cm^-3"
       write(logf,"(A,es10.3,A)") "Theoretical value = ", &
            rho_crit_0*Omega_B/(mu*m_p)*(1.0+redshift)**3, &
            " cm^-3" 
       write(logf,"(A,es10.3,A)") "(at z=0 : ", &
            rho_crit_0*Omega_B/(mu*m_p), &
            " cm^-3)"
    endif
    
  end subroutine density_diagnostics
  
end module density_module
