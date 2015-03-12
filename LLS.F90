!>
!! \brief This module contains data and routines for handling Lyman Limit Systems (LLS) (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 01-April-2014 (27-Jun-2013 (20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))
!!
!! \b Version: 

module lls_module

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: logf
  use my_mpi
  use grid, only: dr, vol, sim_volume
  use cgsconstants, only: c
  use cgsphotoconstants, only: sigma_HI_at_ion_freq
  use astroconstants, only: Mpc
  use cosmology_parameters, only: Omega0, H0
  use nbody, only: id_str,dir_LLS
  use nbody, only: LLSformat, LLSaccess, LLSheader
  use c2ray_parameters, only: type_of_LLS

  implicit none

  ! LLS data
  real(kind=dp),parameter :: opdepth_LL = 2.0 !< typical optical depth of LLS
  real(kind=dp),parameter :: N_1 = opdepth_LL / sigma_HI_at_ion_freq !< typical column density of LLS
  real(kind=dp),public :: n_LLS
  real(kind=dp),public :: coldensh_LLS = 0.0_dp ! Column density of LLSs per cell
  real(kind=dp),public :: mfp_LLS_pMpc
  real,dimension(:,:,:),allocatable :: LLS_grid

  ! LLS parameters
  !> Do not use the LLSs if the mfp is smaller than this.
  real,parameter :: limit_mfp_LLS_pMpc=0.2
  
  ! Different models for LLS redshift evolution
  ! a) Model Prochaska et al. (2010)
  !real(kind=dp),parameter :: C_LLS = 1.9
  !real(kind=dp),parameter :: z_x = 3.7
  !real(kind=dp),parameter,public :: y_LLS = 5.1
  !real(kind=dp),parameter :: beta=1.28 ! not clear what to use here.
  ! b) Model Songaila & Cowie (2010)
  real(kind=dp),parameter :: C_LLS = 2.84
  real(kind=dp),parameter :: z_x = 3.5
  real(kind=dp),parameter,public :: y_LLS = 2.04
  real(kind=dp),parameter :: beta=1.28
  ! c) Model McQuinn et al. (2011)
  !real(kind=dp),parameter :: C_LLS = 2.34
  !real(kind=dp),parameter :: z_x = 3.5
  !real(kind=dp),parameter,public :: y_LLS = 2.85
  !real(kind=dp),parameter :: beta=1.3
  
  public :: set_LLS, LLS_point

  
#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ============================================================================

  subroutine LLS_init ()
    
#ifdef IFORT
    !For gamma function
    use ISO_C_BINDING
#endif
    
#ifdef IFORT
    !For gamma function
    interface
       real(c_double) function tgamma (y) bind(c)
         use iso_c_binding
         real(c_double), value :: y
       end function tgamma
    end interface
#endif

    if (type_of_LLS == 1) then
       ! 1/distance between LLSs expressed in grid cells (z=0)
       n_LLS = C_LLS * (1.0/(1.0 + z_x)) ** y_LLS * dr(1) * H0*sqrt(Omega0) / c

       !n_LLS=n_LLS * ((1.0 + zred)/(1.0+zred_prev))** (y_LLS+1.5)
       !mfp=c/((1.0+z) * Hz * C_LLS * ((1.0 + z)/(1.0 + z_x)) ** y_LLS )
       ! Add the beta correction as explained in Songaila & Cowie (2010).
       ! This corrects for the fact that not all LLS have the same
       ! column density. beta is the slope of the distribution function
       ! of LLS over column densities.
       ! This expression needs the gamma function. For the intel compiler
       ! this is tgamma(). For other compilers (Fortran 2008 standard) this
       ! is gamma().
#ifdef IFORT    
       n_LLS=n_LLS*tgamma(2.0-beta)/(opdepth_LL**(1.0-beta))
#else
       n_LLS=n_LLS*gamma(2.0-beta)/(opdepth_LL**(1.0-beta))
#endif
    else
       ! Set distance between LLS (and mean free path) to infinity
       ! If type_of_LLS = 2 this will be overwritten by cell specific
       ! values.
       n_LLS=0.0d0
    endif

  end subroutine LLS_init

  ! ===========================================================================

  subroutine set_LLS (z)

    !! Two cases:\n
    !! 1: constant LLS optical depth per cell\n
    !! 2: position dependent LLS optical depth

    real(kind=dp),intent(in) :: z

    select case (type_of_LLS)
    case(1)
       ! Calculate the mean free path in pMpc
       mfp_LLS_pMpc=dr(1)/n_LLS/Mpc
       ! Column density per cell due to LLSs
       if (mfp_LLS_pMpc > limit_mfp_LLS_pMpc) then
          coldensh_LLS = N_1 * n_LLS
       else
          coldensh_LLS = 0.0
       endif
    case(2) 
       call read_lls_grid (z)
    end select

    if (rank == 0) then
       write(logf,*) "Average optical depth per cell due to LLSs: ", &
            coldensh_LLS*sigma_HI_at_ion_freq,"(type ", type_of_LLS,")"
       write(logf,*) "Mean free path (pMpc): ", mfp_LLS_pMpc
    endif
    
  end subroutine set_LLS

  ! ===========================================================================

  subroutine LLS_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_LLS /= 2) then
       write(logf,*) &
         "Error: calling position dependent LLS, but array is not initialized."
    else
       coldensh_LLS = LLS_grid(i,j,k)
    endif

  end subroutine LLS_point

  ! ===========================================================================

  subroutine read_LLS_grid (zred_now)

    ! Initializes position dependent LLS optical depth (at redshift zred_now)

    ! Author: Garrelt Mellema

    ! Date: 16-Mar-2011

    ! Version: 

    real(kind=dp),intent(in) :: zred_now
    
    integer :: m1,m2,m3 ! size of mesh in cross section file (header)
    character(len=512):: LLS_file
    character(len=6) :: zred_str

    ! clumping in file is in 4B reals, read in via this array
    !real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(LLS_grid))) &
         allocate(LLS_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       LLS_file=trim(adjustl(dir_LLS))// &
            trim(adjustl(zred_str))// &
            "cross_section.bin"

       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " clumping input from ",trim(LLS_file)

       ! Open clumping file: note that the format is determined
       ! by the values of clumpingformat and clumping access,
       ! both set in the nbody module.
       open(unit=22,file=LLS_file,form=LLSformat, &
            access=LLSaccess,status="old")
       
       ! Read in data
       ! Read in header if there is one
       if (LLSheader) then
          read(22) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with LLS cross sections unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif
       ! Read in data and store it in clumping_grid
       read(22) LLS_grid
       write(logf,*) 'LLS data read'
       
       ! close file
       close(unit=22)
       
       ! Calculate mean free path
       ! Make sure sim_volume is in proper length units
       mfp_LLS_pMpc=sim_volume/(sum(LLS_grid)*Mpc*(1.0+zred_now)**3)

       ! Convert to column density
       LLS_grid=LLS_grid/vol ! 1/(mean free path) = n_LLS
       LLS_grid=N_1 * LLS_grid
    endif

#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(LLS_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_REAL,0,MPI_COMM_NEW,mympierror)
#endif
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    coldensh_LLS=sum(LLS_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(logf,*) "Statistics on LLS column density"
       write(logf,*) "minimum: ",minval(LLS_grid)
       write(logf,*) "maximum: ",maxval(LLS_grid)
       write(logf,*) "average column density: ",coldensh_LLS
       write(logf,*) "mean free path: ",mfp_LLS_pMpc
    endif

  end subroutine read_LLS_grid

end module lls_module
