module sourceprops

  use precision
  use my_mpi, only: rank
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0
  use pmfast, only: id_str, M_grid
  use material, only: xh
  use grid, only: x,y,z
  !use sizes
  !use mathconstants
  !use cgs

  use c2ray_parameters, only: phot_per_atom1, phot_per_atom2, lifetime

  integer :: NumSrc
  real(kind=dp),dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: srcMass
  real(kind=dp),dimension(:),allocatable :: NormFlux
  integer,dimension(:),allocatable :: srcSeries

  integer,private :: NumSrc0
  real(kind=dp),dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01

contains
  
  ! =======================================================================

  subroutine source_properties(zred_now)

    ! Input routine: establish the source properties
    ! Author: Garrelt Mellema
    ! Update: 20-Sep-2006 (3-jan-2005, 15-Apr-2004)

    ! For random permutation of sources
    use  m_ctrper

    real(kind=8),intent(in) :: zred_now ! current redshift

    character(len=180) :: sourcelistfile,sourcelistfilesuppress
    integer :: ns,ns0

    ! Ask for input
    if (rank == 0) then
       ! Sources are read from file

       ! Construct the file name
       if (zred_now.ge.10.0) then
          write(sourcelistfile,'(f6.3,3A)') zred_now, &
               "-",trim(adjustl(id_str)),"_sources.dat"
          write(sourcelistfilesuppress,'(f6.3,3A)') zred_now, &
               "-",trim(adjustl(id_str)),"_sources_used_wfgamma.dat"
       else
          write(sourcelistfile,'(f5.3,3A)') zred_now, &
               "-",trim(adjustl(id_str)),"_sources.dat"
          write(sourcelistfilesuppress,'(f5.3,3A)') zred_now, &
               "-",trim(adjustl(id_str)),"_sources_used_wfgamma.dat"
       endif
      
       open(unit=50,file=sourcelistfile,status='old')
       open(unit=49,file=sourcelistfilesuppress,status='unknown')
      
       ! Number of sources
       read(50,*) NumSrc0

       ! Report
       write(30,*) 'Number of sources, no suppression: ',NumSrc0

       ! Read in source positions and mass to establish number
       ! of non-suppressed sources
       NumSrc = 0
       do ns0=1,NumSrc0
          read(50,*) srcpos0(1),srcpos0(2),srcpos0(3),SrcMass00,SrcMass01
          ! the cell is still neutral, no suppression
          if (SrcMass00 /= 0.0 .or. &
               xh(srcpos0(1),srcpos0(2),srcpos0(3),1).lt.0.1) NumSrc=NumSrc+1
       enddo
    endif
#ifdef MPI
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
    deallocate(srcpos)
    deallocate(rsrcpos)
    deallocate(srcMass)
    deallocate(NormFlux)
    deallocate(srcSeries)
    allocate(srcpos(3,NumSrc))
    allocate(rsrcpos(3,NumSrc))
    allocate(SrcMass(NumSrc))
    allocate(NormFlux(NumSrc))
    allocate(SrcSeries(NumSrc))

    if (rank == 0) then
       ! Read in source positions and mass
       do ns0=1,NumSrc0
          
          read(50,*) srcpos0(1),srcpos0(2),srcpos0(3), &
               SrcMass00,SrcMass01
          
          if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1).lt.0.1) then
             ! the cell is still neutral, no suppression
             ns=ns+1
             ! Source positions in file start at zero!
             srcpos(1,ns)=srcpos0(1)
             srcpos(2,ns)=srcpos0(2)
             srcpos(3,ns)=srcpos0(3)
            
             ! Source is always at cell centre!!
             rsrcpos(1,ns)=x(srcpos(1,ns))
             rsrcpos(2,ns)=y(srcpos(2,ns))
             rsrcpos(3,ns)=z(srcpos(3,ns))
             SrcMass(ns)=SrcMass00*phot_per_atom1  & !massive sources
                  +SrcMass01*phot_per_atom2      !small sources   
          elseif (SrcMass00.gt.0.0d0) then
             !the cell is ionized but source is massive enough to survive
                                !and is assumed Pop. II 
             ns=ns+1
             ! Source positions in file start at zero!
             srcpos(1,ns)=srcpos0(1)
             srcpos(2,ns)=srcpos0(2)
             srcpos(3,ns)=srcpos0(3)
             ! Source is always at cell centre!!
             rsrcpos(1,ns)=x(srcpos(1,ns))
             rsrcpos(2,ns)=y(srcpos(2,ns))
             rsrcpos(3,ns)=z(srcpos(3,ns))
             SrcMass(ns)=SrcMass00*phot_per_atom1
          else
             ! Report
             write(30,*) 'Source dropped: ', &
                  srcpos0(1),srcpos0(2),srcpos0(3)
             write(30,*) 'mass= ',SrcMass01*(M_grid/M_SOLAR), &
                  xh(srcpos0(1),srcpos0(2),srcpos0(3),0)
          endif
       enddo
       ! Report
       write(30,*) 'Number of sources: ',NumSrc

       ! Save new source list, without the suppressed ones
       write(49,*) NumSrc

       do ns0=1,NumSrc
          write(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
               SrcMass(ns0)
       end do

       close(49)
    endif

#ifdef MPI
    ! Distribute the source parameters to the other nodes
    call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(rsrcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(SrcMass,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif

    ! Turn masses into luminosities
    do ns=1,NumSrc
       NormFlux(ns)=SrcMass(ns)*M_grid*  &!note that now photons/atom are included in SrcMass
            Omega_B/(Omega0*m_p*lifetime)!/S_star
    enddo

    if (rank == 0) then
       ! Create array of source numbers for generating random order
       do ns=1,NumSrc
          SrcSeries(ns)=ns
       enddo
       
       ! Make a random order
       call ctrper(SrcSeries(1:NumSrc),1.0)
    endif
#ifdef MPI
    ! Distribute the source series to the other nodes
    call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif

    return
  end subroutine source_properties

end module sourceprops
