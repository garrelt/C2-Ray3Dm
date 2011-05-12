!>
!! \brief This module contains data and routines for handling the sources
!! of reionization.
!!
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 10-May-2011 (30-Jan-2008)
!!
!! \b Version: Processed output from Gadget simulation of Andreas Pawlik and Joop Schaye. The source lists already contain the UV luminosity.

module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, Number_Sourcetypes 

  implicit none

  !> base name of source list files
  character(len=100),private :: sourcelistfile_base="_gadget_sources.dat"

  integer :: NumSrc=0 !< Number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos !< grid position of sources
  real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources

  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01
  character(len=6),private :: z_str !< string value of redshift

  character(len=512),private :: sourcelistfile,sourcelistfilesuppress

contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  !! Authors: Garrelt Mellema, Ilian Iliev
  !! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))

  subroutine source_properties(zred_now,nz,lifetime2,restart)

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    integer :: ns,ns0

#ifdef MPI
    integer :: mympierror
#endif

    ! Deallocate source arrays
    if (NumSrc0 /= 0) then
       deallocate(srcpos)
       deallocate(rsrcpos)
       deallocate(srcMass)
       deallocate(NormFlux)
       deallocate(srcSeries)
    endif

    ! Rank 0 reads in sources
    if (rank == 0) then

       ! Sources are read from files with redshift in the file name:
       ! construct redshift string
       write(z_str,'(f6.3)') zred_now
 
       ! Construct the file name
       sourcelistfile=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//trim(adjustl(sourcelistfile_base))
      
       call establish_number_of_active_sources (restart)

    endif ! end of rank 0 test
#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif

    ! Allocate arrays for this NumSrc
    if (NumSrc > 0) then
       allocate(srcpos(3,NumSrc))
       allocate(rsrcpos(3,NumSrc))
       allocate(NormFlux(0:NumSrc)) ! 0 will hold lost photons
       allocate(SrcSeries(NumSrc))
       
       ! Fill in the source arrays
       if (rank == 0) then
          call read_in_sources (restart)

          call assign_uv_luminosities (lifetime2,nz)

          write(logf,*) 'Total flux= ',sum(NormFlux(1:NumSrc))*S_star_nominal,' s^-1'
          ! Create array of source numbers for generating random order
          do ns=1,NumSrc
             SrcSeries(ns)=ns
          enddo
          
          ! Make a random order
          call ctrper(SrcSeries(1:NumSrc),1.0)
       endif

#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       ! Distribute the source series to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    endif
    
  end subroutine source_properties

  ! =======================================================================

  subroutine establish_number_of_active_sources (restart)

    integer,intent(in) :: restart !< not used here
    ! the source list can be read in for a restart

    open(unit=50,file=sourcelistfile,status='old')
    ! Number of sources
    read(50,*) NumSrc0
    NumSrc=NumSrc0 ! no suppression in these source lists
    close(50)
    
    write(logf,*) "Number of sources: ",NumSrc

  end subroutine establish_number_of_active_sources

  ! =======================================================================

  subroutine read_in_sources (restart)

    integer,intent(in) :: restart !< not used here

    integer :: ns
    integer :: ns0

    open(unit=50,file=sourcelistfile,status='old')
    ! Number of sources
    read(50,*) NumSrc0
    ! Read in source positions and mass
    ns=0
    do ns0=1,NumSrc0
       read(50,*) srcpos0(1),srcpos0(2),srcpos0(3),NormFlux(ns0)
       ns=ns0
       ! Source positions in file start at zero!
       srcpos(1,ns)=srcpos0(1)
       srcpos(2,ns)=srcpos0(2)
       srcpos(3,ns)=srcpos0(3)
       ! Source is always at cell centre!!
       rsrcpos(1,ns)=x(srcpos(1,ns))
       rsrcpos(2,ns)=y(srcpos(2,ns))
       rsrcpos(3,ns)=z(srcpos(3,ns))
    enddo

  end subroutine read_in_sources

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)

    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    ! Please put your conversion of mass to luminosity (ionizing
    ! photons / second) here.
    ! This version (Pawlik & Schaye output): source list supplies 
    ! luminosity, only normalization is needed
    do ns=1,NumSrc
       NormFlux(ns)=NormFlux(ns)/S_star_nominal
    enddo
    
  end subroutine assign_uv_luminosities

  ! =======================================================================
  
  !> Initialization routine: determine the source model and optionally read 
  !! in source properties
  !! Author: Garrelt Mellema
  
  !! This version: dummy (source model fixed outside of the code)

  subroutine source_properties_ini ()
    
  end subroutine source_properties_ini
  
end module sourceprops
