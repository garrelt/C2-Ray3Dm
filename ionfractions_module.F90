module ionfractions_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use sm3d, only: read_sm3d_dp_file_routine
  use c2ray_parameters, only: epsilon
  use my_mpi
  
  implicit none

  type ionstates    
     real(kind=dp) :: h(0:1)          !< H  ionization fractions        
     real(kind=dp) :: h_av(0:1)       !< average H  ionization fractions        
     real(kind=dp) :: h_old(0:1)      !< H  ionization fractions from last time step
  end type ionstates

  ! xh - ionization fractions for one cell
#ifdef ALLFRAC
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh
#else
  real(kind=dp),dimension(:,:,:),allocatable :: xh
#endif

#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ===========================================================================
  
  subroutine xfrac_array_init ( )

       ! Allocate ionization fraction array
#ifdef ALLFRAC
       allocate(xh(mesh(1),mesh(2),mesh(3),0:1))
#else
       allocate(xh(mesh(1),mesh(2),mesh(3)))
#endif
       ! Assign ionization fractions. For z = 40 - 20 RECFAST gives 
       ! an average ionization fraction of about 2e-4. We use this
       ! here.
       ! In case of a restart this will be overwritten in xfrac_ini
#ifdef ALLFRAC
       xh(:,:,:,1)=2e-4
       xh(:,:,:,0)=1.0_dp-xh(:,:,:,1)
#else
       xh(:,:,:)=2e-4
#endif

  end subroutine xfrac_array_init

  ! ===========================================================================

  subroutine xfrac_restart_init (zred_now)

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
    real(kind=dp),dimension(:,:,:),target,allocatable :: xh1_real
    !real(kind=si),dimension(:,:,:),allocatable :: xh1_real
    real(kind=dp),dimension(:,:,:),pointer :: ion_fraction

    if (rank == 0) then

       ! Construct file names
       write(zred_str,"(f6.3)") zred_now
       xfrac_file= trim(adjustl(results_dir))// &
            "xfrac3d_"//trim(adjustl(zred_str))//".bin"

       ! Report
       write(unit=logf,fmt="(2A)") "Reading ionization fractions from ", &
            trim(xfrac_file)

       ! GM/110308: If we only work with ionized fractions, do not use 
       ! this anymore (ionized fractions are saved as dp). If we are
       ! working with both neutral and ionized fraction, we need this
       ! since one cannot read in partial arrays. Perhaps some pointer
       ! would help?
       allocate(xh1_real(mesh(1),mesh(2),mesh(3)))
       ion_fraction => xh1_real
       
       call read_sm3d_dp_file_routine(xfrac_file,ion_fraction)

#ifdef ALLFRAC
       xh(:,:,:,1)=ion_fraction(:,:,:)
       xh(:,:,:,0)=1.0_dp-xh(:,:,:,1)
#else
       xh(:,:,:)=ion_fraction(:,:,:)
#endif

       deallocate(xh1_real)

    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
#ifdef ALLFRAC
    call MPI_BCAST(xh,mesh(1)*mesh(2)*mesh(3)*2,MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#else
    call MPI_BCAST(xh,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine xfrac_restart_init

end module ionfractions_module
