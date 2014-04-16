module ionfractions_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
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
    real(kind=dp),dimension(:,:,:),allocatable :: xh1_real
    !real(kind=si),dimension(:,:,:),allocatable :: xh1_real

    if (rank == 0) then
       ! GM/110308: If we only work with ionized fractions, do not use 
       ! this anymore (ionized fractions are saved as dp). If we are
       ! working with both neutral and ionized fraction, we need this
       ! since one cannot read in partial arrays. Perhaps some pointer
       ! would help?
       allocate(xh1_real(mesh(1),mesh(2),mesh(3)))
       write(zred_str,"(f6.3)") zred_now
!       xfrac_file= "./xfrac3d_"//trim(adjustl(zred_str))//".bin"
       xfrac_file= trim(adjustl(results_dir))// &
            !"Ifront3_"//trim(adjustl(zred_str))//".bin"
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
          !read(20) xh
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          !xh(:,:,:,1)=real(xh1_real(:,:,:),dp)
#ifdef ALLFRAC
          xh(:,:,:,1)=xh1_real(:,:,:)
          xh(:,:,:,0)=1.0_dp-xh(:,:,:,1)
#else
          xh(:,:,:)=xh1_real(:,:,:)
#endif
       endif

       ! close file
       close(20)
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
