module pmfast

  ! This file contains routine having to do with the PMFAST N-body
  ! simulation.
  ! The routine in here is
  ! - pmfast_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.
  ! It calculates the mass (in particles) of the simulation box
  ! It sets an identifier (id_str) for the resolution (used when
  ! reading in source list files and density fields.

  use sizes
  use astroconstants
  use my_mpi
  use cosmology_parameters

  implicit none

  integer,parameter :: n_box=3248    ! cells/side (in N-body)
  real(kind=8),parameter :: boxsize=100.0  ! Box size in Mpc/h comoving

  ! properties of the box:
  ! M_box      - mass in box
  ! M_particle - mass per particle
  ! M_grid - mean mass per pmfast cell
  real(kind=8) :: M_box,M_particle,M_grid
  
  ! redshift sequence information
  integer :: NumZred               ! number of redshifts
  real(kind=8),dimension(:),allocatable :: zred_array ! array of redshifts
  character(len=8) :: id_str       ! resolution dependent string

  character(len=180),parameter :: dir_dens="/disk/sn-12/garrelt/Simulations/"// &
       "Reionization/100Mpc_WMAP3/203/"
  integer,parameter :: tot_nfiles=7 ! number of files at 812^3

#ifdef MPI
  integer :: ierror
#endif

contains

  subroutine pmfast_ini ()
    
    
    character(len=180) :: redshift_file ! name of file with list of redshifts
    real(kind=8) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter
    integer :: nz ! loop counter

    ! Ask for redshift file
    if (rank == 0) then
       write(*,'(A,$)') 'File with redshifts: '
       read(*,*) redshift_file
       
       ! Open and read redshift file
       open(unit=60,file=redshift_file,form='formatted',status='old')
       read(60,*) NumZred
       allocate(zred_array(NumZred))
       do nz=1,NumZred
          read(60,*) zred_array(nz)
       enddo
       ! print*,nz,zred_array
       close(20)
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
#endif

    ! Parameters of simulations boxes
    ! critical density in M0/Mpc^3
    rho_crit_0_MpcM_sun=rho_crit_0*Mpc**3/M_SOLAR
    rho_crit_0_MpcM_sun=rho_crit_0*Mpc**3/M_SOLAR

    !rho_bar = rho_crit_0_MpcM_sun*Omega0 ! mean matter density in M0/Mpc^3
    rho_matter = rho_crit_0*Omega0        ! mean matter density (g/cm^3)
    M_box      = rho_matter*(boxsize*Mpc/h)**3   ! mass in box (g, not M0) 
    M_particle = 8.0*M_box/(real(n_box)**3) ! mass per particle (g, not M0)
    M_grid = M_particle/8.

    ! Set identifying string (resolution-dependent)
    ! Construct the file name
    if (mesh(1) == 203) id_str="coarsest"
    if (mesh(1) == 406) id_str="coarser"
    if (mesh(1) == 812) id_str="coarse"
    print*,id_str

    return
  end subroutine pmfast_ini

end module pmfast
