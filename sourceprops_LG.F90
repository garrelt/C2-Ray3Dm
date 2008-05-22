module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: id_str, dir_src
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom1, phot_per_atom2, lifetime, &
       S_star_nominal, StillNeutral

  implicit none

  integer :: NumSrc 
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: srcMass
  real(kind=dp),dimension(:),allocatable :: NormFlux
  integer,dimension(:),allocatable :: srcSeries

  integer,private :: NumSrc0=0
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01
  character(len=6) :: z_str 
  integer,private :: NumMassiveSrc,NumSupprbleSrc,NumSupprsdSrc

contains
  
  ! =======================================================================

  subroutine source_properties(zred_now,nz,lifetime2,restart)

    ! Input routine: establish the source properties
    ! Authors: Garrelt Mellema, Ilian Iliev
    ! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    character(len=512) :: sourcelistfile,sourcelistfilesuppress
    integer :: ns,ns0
    character(len=3) :: id !output number string

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(rsrcpos)) deallocate(rsrcpos)
    if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(srcSeries)) deallocate(srcSeries)
    
    ! Ask for input
    if (rank == 0) then
       
       ! Sources are read from file
       write(z_str,'(f6.3)') zred_now
       write(id,'(i3.3)') nz

       ! Construct the file name
       sourcelistfile=trim(adjustl(dir_src))//trim(adjustl(id))//"-"// &
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))//"_sources.dat"
       sourcelistfilesuppress=trim(adjustl(dir_src))//trim(adjustl(id))// &
            "-"//trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            "_sources_used_wfgamma.dat"
       
       if (restart == 0 .or. restart == 1) then
          open(unit=50,file=sourcelistfile,status='old')
          ! Number of sources
          read(50,*) NumSrc0
          
          ! Report
          write(logf,*) "Total number of source locations, no suppression: ", &
               NumSrc0
          
          ! Read in source positions and mass to establish number
          ! of non-suppressed sources
          NumSrc = 0
          NumMassiveSrc = 0
          NumSupprbleSrc = 0
          NumSupprsdSrc = 0
          do ns0=1,NumSrc0
             read(50,*) srcpos0(1),srcpos0(2),srcpos0(3),SrcMass00,SrcMass01
             ! the cell is still neutral, no suppression
             if (SrcMass00 /= 0.0 .or. &
                  xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) &
                  NumSrc=NumSrc+1
             ! Count different types of sources
             if (SrcMass00 /= 0.0) NumMassiveSrc=NumMassiveSrc+1
             if (SrcMass01 /= 0.0) NumSupprbleSrc=NumSupprbleSrc+1
             ! How many suppressed?
             if (SrcMass01 /= 0.0 .and. &
                  xh(srcpos0(1),srcpos0(2),srcpos0(3),1) > StillNeutral) &
                  NumSupprsdSrc=NumSupprsdSrc+1
          enddo
          close(50)
          write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
          write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
          write(logf,*) "Number of massive sources: ",NumMassiveSrc
          write(logf,*) "Suppressed fraction: ", &
               real(NumSupprsdSrc)/real(NumSupprbleSrc)
       else
          ! Upon restart use the previously calculated suppressed source
          ! list
          open(unit=49,file=sourcelistfilesuppress,status='unknown')
          ! Number of sources
          read(49,*) NumSrc
          close(49)
       endif
       write(logf,*) "Number of sources, with suppression: ",NumSrc
    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif
    
    ! Allocate arrays for this NumSrc
    allocate(srcpos(3,NumSrc))
    allocate(rsrcpos(3,NumSrc))
    allocate(SrcMass(NumSrc))
    allocate(NormFlux(NumSrc))
    allocate(SrcSeries(NumSrc))
    
    if (rank == 0) then
       if (restart == 0 .or. restart == 1) then
          open(unit=50,file=sourcelistfile,status='old')
          ! Number of sources
          read(50,*) NumSrc0
          ! Read in source positions and mass
          ns=0
          do ns0=1,NumSrc0
             read(50,*) srcpos0(1),srcpos0(2),srcpos0(3), &
                  SrcMass00,SrcMass01
             
             if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) then
                ! the cell is still neutral, no suppression
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                ! Source is always at cell centre!!
                rsrcpos(1,ns)=x(srcpos(1,ns))
                rsrcpos(2,ns)=y(srcpos(2,ns))
                rsrcpos(3,ns)=z(srcpos(3,ns))
                SrcMass(ns)=SrcMass00*phot_per_atom1  & !massive sources
                     +SrcMass01*phot_per_atom2      !small sources   
             elseif (SrcMass00 > 0.0d0) then
                !the cell is ionized but source is massive enough to survive
                !and is assumed Pop. II 
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                ! Source is always at cell centre!!
                rsrcpos(1,ns)=x(srcpos(1,ns))
                rsrcpos(2,ns)=y(srcpos(2,ns))
                rsrcpos(3,ns)=z(srcpos(3,ns))
                SrcMass(ns)=SrcMass00*phot_per_atom1
                !else
                !   ! Report
                !   write(logf,*) 'Source dropped: ', &
                !        srcpos0(1),srcpos0(2),srcpos0(3)
                !   write(logf,*) 'mass= ',SrcMass01*(M_grid/M_SOLAR), &
                !        xh(srcpos0(1),srcpos0(2),srcpos0(3),0)
             endif
          enddo
          
          ! Save new source list, without the suppressed ones
          open(unit=49,file=sourcelistfilesuppress,status='unknown')
          write(49,*) NumSrc
          do ns0=1,NumSrc
             write(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
                  SrcMass(ns0)
          enddo
          close(49)
       else ! of restart test
          ! Read source list from file saved previously
          open(unit=49,file=sourcelistfilesuppress,status="old")
          write(logf,*) "Reading ",NumSrc," sources from ", &
               trim(adjustl(sourcelistfilesuppress))
          read(49,*) NumSrc
          do ns0=1,NumSrc
             read(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
                  SrcMass(ns0)
             ! Source is always at cell centre!!
             rsrcpos(1,ns0)=x(srcpos(1,ns0))
             rsrcpos(2,ns0)=y(srcpos(2,ns0))
             rsrcpos(3,ns0)=z(srcpos(3,ns0))
          enddo
          close(49)
       endif ! of restart test
    endif ! of rank 0 test
    
    
#ifdef MPI
    ! Distribute the source parameters to the other nodes
    call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(SrcMass,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Turn masses into luminosities
    do ns=1,NumSrc
       ! Hubble constant h is included in source list (i.e. it contains source 
       ! masses in M_SOLAR	  
       NormFlux(ns)=SrcMass(ns)*M_SOLAR*  &!note that now photons/atom are included in SrcMass
            Omega_B/(Omega0*m_p)/S_star_nominal
       !NormFlux(ns)=NormFlux(ns)/lifetime
       NormFlux(ns)=NormFlux(ns)/lifetime2
    enddo
    
    if (rank == 0) then
       
       write(logf,*) 'Source lifetime=', lifetime2/3.1536e13
       !write(logf,*) 'Source lifetime=', lifetime/3.1536e13
       write(logf,*) 'Total flux= ',sum(NormFlux)
       ! Create array of source numbers for generating random order
       do ns=1,NumSrc
          SrcSeries(ns)=ns
       enddo
       
       ! Make a random order
       call ctrper(SrcSeries(1:NumSrc),1.0)
    endif

#ifdef MPI
    ! Distribute the source series to the other nodes
    call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    return
  end subroutine source_properties

end module sourceprops
