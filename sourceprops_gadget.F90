module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom1, phot_per_atom2, lifetime,S_star_nominal 

  integer :: NumSrc
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: srcMass
  real(kind=dp),dimension(:),allocatable :: NormFlux
  integer,dimension(:),allocatable :: srcSeries

  integer,private :: NumSrc0=0
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01

contains
  
  ! =======================================================================

  subroutine source_properties(zred_now,lifetime2)

    ! Input routine: establish the source properties
    ! Author: Garrelt Mellema
    ! Update: 20-Sep-2006 (3-jan-2005, 15-Apr-2004)

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step

    character(len=180) :: sourcelistfile,sourcelistfilesuppress
    character(len=6) :: zred_str
    integer :: ns,ns0

#ifdef MPI
    integer :: ierror
#endif

    if (NumSrc0 /= 0) then

       write(logf,*) 'Deallocating source arrays'
       deallocate(srcpos)
       deallocate(rsrcpos)
       deallocate(srcMass)
       deallocate(NormFlux)
       deallocate(srcSeries)
    endif

    ! Ask for input
    if (rank == 0) then
       ! Sources are read from file

       ! Construct the file name
       write(zred_str,'(f6.3)') zred_now
       sourcelistfile=trim(adjustl(zred_str))//"_gadget_sources.dat"
      
       open(unit=50,file=sourcelistfile,status='old')
      
       ! Number of sources
       read(50,*) NumSrc0
       NumSrc=NumSrc0 ! no suppression

       ! Report
       write(logf,*) 'Number of sources: ',NumSrc0
       close(50)
    endif
#ifdef MPI
    call MPI_BCAST(NumSrc0,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif

    write(logf,*) 'Number of sources, with suppression: ',NumSrc
    allocate(srcpos(3,NumSrc))
    allocate(rsrcpos(3,NumSrc))
    allocate(SrcMass(NumSrc))
    allocate(NormFlux(NumSrc))
    allocate(SrcSeries(NumSrc))

    if (rank == 0) then
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

          NormFlux(ns)=NormFlux(ns)/S_star_nominal
       enddo
    endif

#ifdef MPI
    ! Distribute the source parameters to the other nodes
    call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(NormFlux,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif

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

    !write(*,*) SrcSeries
    !write(*,*) srcpos
    !write(*,*) rsrcpos
    !write(*,*) NormFlux*S_star_nominal

  end subroutine source_properties

end module sourceprops
