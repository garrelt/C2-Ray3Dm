module clumping_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use c2ray_parameters, only: type_of_clumping,clumping_factor
  use nbody, only: nbody_type, id_str, dir_dens, NumZred, Zred_array, dir_clump, dir_LLS
  use nbody, only: clumpingformat, clumpingaccess, clumpingheader

  implicit none

  ! Clumping data
  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  real(kind=dp) :: avg_dens !< average density
  character(len=512) :: clumping_fit_file
  public :: set_clumping, clumping_point

#ifdef MPI
  integer,private :: mympierror
#endif

contains
  
  ! ===========================================================================
  
  subroutine set_clumping(z)

    !! 1: constant clumping (with clumping_factor)\n
    !! 2: 3.5Mpc PM, WMAP1 clumping\n
    !! 3: 3.5Mpc PM, WMAP3 clumping\n
    !! 4: 1 Mpc P3M\n
    !! 5: position dependent clumping

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
    real(kind=dp) :: z_read
    real(kind=dp) :: a0
    real(kind=dp) :: a1
    real(kind=dp) :: a2
    real(kind=dp) :: error
    real(kind=dp) :: avg_dens !< average density
    integer :: io_status

    ! clumping in file is in 4B reals, read in via this array
    !real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(clumping_grid))) &
         allocate(clumping_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       write(30,*) "Reading ",id_str," input"
       clump_file=trim(adjustl(dir_clump))// &
            trim(adjustl(zred_str))// &
            "n_all.dat"
            !"c_all.dat"
            !"clump_"//trim(adjustl(id_str))//".dat"
       ! Clumping file should be density with halos included.
       ! This most likely has an identical file name to density
       ! with halos excluded, but resides in a different directory.
       ! GM, 2009-12-09

       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " clumping input from ",trim(clump_file)
       ! Open clumping file: note that the format is determined
       ! by the values of clumpingformat and clumping access,
       ! both set in the nbody module.
       open(unit=20,file=clump_file,form=clumpingformat, &
            access=clumpingaccess,status="old")
       
       ! Read in data
       ! Read in header if there is one
       if (clumpingheader) then
          read(20) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with clumping factors unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif
       ! Read in data and store it in clumping_grid
       read(20) clumping_grid
       write(logf,*) 'Clumping data read'
       !clumping_grid(:,:,:)=clumping_real(:,:,:)
       
       ! close file
       close(20)
       
       ! Normalize to average density
       avg_dens=sum(clumping_grid)/ &
            (real(mesh(1),dp)*real(mesh(2),dp)*real(mesh(3),dp))
       clumping_grid=clumping_grid/avg_dens
       clumping_grid=log10(clumping_grid*clumping_grid)
       ! Report on data: min, max, total
       ! Disabled (GM 2009-12-09): if we read in a density field
       ! there is no sense in reporting these.
       ! assign mean to clumping for reporting in set_clumping
       !clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
       !write(logf,*) "Statistics BEFORE applying clumping fit"
       !write(logf,*) "minimum: ",minval(clumping_grid)
       !write(logf,*) "maximum: ",maxval(clumping_grid)
       !write(logf,*) "average clumping: ",clumping

       ! Read clumping fit (for clumping based on smaller scales
       ! This file contains the parameters for the fit
       ! y=a0+a1*x+a2*x²
       ! where x is the alog10(<n²>_int) and y is alog10(<n²>_Jeans).
       open(unit=21,file=trim(adjustl(dir_clump))//clumping_fit_file, &
          status="old",form="formatted")
       z_read=0.0
       io_status=0
       ! Read in lines until the correct redshift is found.
       do while(z_read-zred_now > 1e-2 .and. io_status == 0)
          read(unit=21,iostat=io_status) z_read,a0,a1,a2,error
       enddo
       close(unit=21)
       if (io_status /= 0) then
          ! Something went wrong, assume no clumping and warn
          a0=0.0
          a1=1.0
          a2=0.0
          write(logf,*) "WARNING: no clumping data found in ", & 
               trim(adjustl(clumping_fit_file))," for redshift ",zred_now
       endif
       clumping_grid(:,:,:)=10.0**(a0+ &
            (a1-1.0)*clumping_grid(:,:,:)+ &
            a2*clumping_grid(:,:,:)*clumping_grid(:,:,:))
    endif

#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_REAL,0,MPI_COMM_NEW,mympierror)
#endif
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(logf,*) "Statistics AFTER applying clumping fit"
       write(logf,*) "clumping fit for redshift ",z_read
       write(logf,*) "minimum: ",minval(clumping_grid)
       write(logf,*) "maximum: ",maxval(clumping_grid)
       write(logf,*) "average clumping: ",clumping
    endif

  end subroutine clumping_init

end module clumping_module
