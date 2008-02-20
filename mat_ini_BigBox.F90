module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, log
  use my_mpi
  use grid, only: vol
  use cgsconstants, only: m_p
  use cosmology_parameters, only: Omega_B,Omega0,rho_crit_0
  use pmfast, only: M_grid, id_str, dir_dens, tot_nfiles
  use abundances, only: mu

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh(1),mesh(2),mesh(3))
  real(kind=dp) :: temper
  real(kind=dp) :: xh(mesh(1),mesh(2),mesh(3),0:1)
  logical isothermal

#ifdef MPI
  integer,private :: ierror
#endif

contains
  ! ============================================================================
  subroutine mat_ini (restart)

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

    integer :: i,j,k,n ! loop counters
    real(kind=dp) :: temper_val
    character(len=1) :: answer

    ! restart
    restart=0 ! no restart by default

    if (rank == 0) then
       ! Ask for temperature, restart. Read in values
       write(*,"(A,$)") "Enter initial temperature (K): "
       read(stdinput,*) temper_val
       write(*,"(A,$)") "Restart (y/n)? : "
       read(stdinput,*) answer
       if (answer == "y" .or. answer == "Y") restart=1
       write(*,"(A,$)") "Restart at midpoint (y/n)? : "
       read(stdinput,*) answer
       if (answer == "y".or.answer == "Y") restart=2
    endif
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
       call MPI_BCAST(restart,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
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

    return
  end subroutine mat_ini

  ! ===========================================================================
  subroutine dens_ini (zred_now)

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
    
    integer :: i,j,k,n,nfile ! loop counters
    character(len=512):: dens_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: convert ! conversion factor
    real(kind=dp) :: summed_density
    real(kind=dp) :: avg_dens

    ! density in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: ndens_real

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       if (id_str /= "coarse") then
          ! Allocate array needed to read in data
          allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
          dens_file=trim(adjustl(dir_dens))//"coarser_densities/"// &
               trim(adjustl(zred_str))// &
               "rho_"//trim(adjustl(id_str))//".dat"
          write(unit=log,fmt="(4A)") "Reading ",id_str, &
               " density input from ",trim(dens_file)
          
          ! Open density file: note that it is in `binary" form
          open(unit=20,file=dens_file,form="binary",status="old")
          
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
          allocate(ndens_real(mesh(1),mesh(2),mesh(3)/tot_nfiles))
          do nfile=0,tot_nfiles-1
             write(nfile_str,"(I1)") nfile
             dens_file=trim(adjustl(dir_dens))// &
                  trim(adjustl(zred_str))// &
                  "rho_c"//nfile_str//".dat"
             open(unit=20,file=dens_file,form="binary",status="old")
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
    call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,ierror)
#endif
       
    ! Report on data: min, max, total
    if (rank == 0) then
       write(log,*) "Raw density diagnostics (in simulation units)"
       write(log,*) "minimum density: ",minval(ndens)/8.
       write(log,*) "maximum density: ",maxval(ndens)/8.
       write(log,*) "summed density:  ",sum(ndens)/8.
    endif

    ! The original values in terms of the mean density
    ! Below is the conversion factor
    ! vol is redshift dependent, that is why we need to recalculate this
    ! M_particle should be in g
    !convert=M_particle*M_solar/vol*Omega_B/Omega0/(mu*m_p)
    convert=M_grid*Omega_B/Omega0/(mu*m_p)/vol

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
    if (rank == 0) then
       write(log,"(A,1pe10.3,A)") "Average density = ",avg_dens," cm^-3"
       write(log,"(A,1pe10.3,A)") "Theoretical value = ", &
            rho_crit_0*Omega_B/(mu*m_p)*(1.0+zred_now)**3, &
            " cm^-3" 
       write(log,"(A,1pe10.3,A)") "(at z=0 : ", &
            rho_crit_0/(mu*m_p)*Omega_B, &
            " cm^-3)"
    endif

    return
  end subroutine dens_ini

  ! ===========================================================================
  subroutine xfrac_ini (zred_now)

    ! Initializes ionization fractions on the grid (at redshift zred_now).
    ! They are read from an 'Ifront3' file which should have been created
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
       xfrac_file= "./xfrac3d_"//trim(adjustl(zred_str))//".bin"

       write(unit=log,fmt="(2A)") "Reading ionization fractions from ", &
            trim(xfrac_file)
       ! Open ionization fractions file
       open(unit=20,file=xfrac_file,form="unformatted",status="old")
       
       ! Read in data
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(*,*) "Warning: file with ionization fractions unusable"
       else
          read(20) xh1_real
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          !xh(:,:,:,1)=real(xh1_real(:,:,:),8)
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
         MPI_COMM_NEW,ierror)
#endif
    
    return
  end subroutine xfrac_ini

end module material
