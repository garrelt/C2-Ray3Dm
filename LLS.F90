!>
!! \brief This module contains data and routines for handling Lyman Limit Systems (LLS) (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 18-May-2023 (25-Jan-2021 (01-April-2014 (27-Jun-2013 ((derives from mat_ini_cosmo2.f))
!!
!! \b Version:
!! The effect of the presence of LLS can be included in different ways. Which
!! way is chosen is determined by the parameter type_of_LLS, set in the module
!! c2ray_parameters:
!! type 1: homogeneous additional absorption
!! type 2: position-dependent additional absorption
!  type 3: finite distance for radiative transfer ("R_max")
!! In type 1 and 2 extra HI column density is added to the grid. This is
!! added at the interface between cells and so it does not contribute to the
!! calculation for a given cell. If a cell has an outgoing column density of
!! N_i then without the LLS model, this will be the ingoing column density of
!! the next cell. With the LLS model, the ingoing column density for the next
!! cell will instead be N_i+N_LLS.
!!
!! If type_of_LLS = 1 then N_LLS is calculated from a mean free path model
!! chosen according to the value of the parameter LLS_model. For LLS_model
!! 1, 2 or 3, we use a redshift-dependent fit given by Worseck et al. (2014). This model, based 
!! on observations, specifies the mean free path of Lyman continuum photons 
!! at a given redshift. The three different models span the uncertainty in
!! this relation (low, standard and high values). For LLS_model 4 and 5 we
!! use a non-evolving proper, respectively comoving value.
!! We convert this quantity to N_LLS by making sure that after a distance
!! correspond to this mfp, a column density corresponding to opdepth_LL
!! (usually 1) is reached.
!! If type_of_LLS = 2 then the programme reads in a grid of N_LLS values
!! which have been calculated using an external programme. The idea is
!! to give higher column densities to higher density cells or cells with
!! collapsed halos in them.

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
  use cosmology_parameters, only: Omega0, H0, h
  use nbody, only: id_str,dir_LLS
  use nbody, only: LLSformat, LLSaccess, LLSheader
  use c2ray_parameters, only: use_LLS, type_of_LLS, LLS_model, R_max_cMpc

  implicit none

  type LLS_model_old
     character(len=25) :: reference
     real(kind=dp) :: C_LLS
     real(kind=dp) :: z_x
     real(kind=dp) :: y_LLS
     real(kind=dp) :: beta
  end type LLS_model_old

  type LLS_model_new
     character(len=25) :: reference
     real(kind=dp) :: A_LLS 
     real(kind=dp) :: z_ref
     real(kind=dp) :: yz_LLS
  end type LLS_model_new

  ! LLS data
  real(kind=dp),parameter :: opdepth_LL = 1.0 !< typical optical depth associated with 1 mean free path
  real(kind=dp),parameter :: N_1 = opdepth_LL / sigma_HI_at_ion_freq !< typical column density of LLS
  real(kind=dp),public :: n_LLS
  real(kind=dp),public :: coldensh_LLS = 0.0_dp ! Column density of LLSs per cell
  real(kind=dp),public :: mfp_LLS_pMpc
  real,dimension(:,:,:),allocatable :: LLS_grid

  ! LLS parameters
  !> Do not use the LLSs if the mfp is smaller than this.
  real,parameter :: limit_mfp_LLS_pMpc=0.2
  real,parameter :: limit_mfp_LLS_cMpc=1.0

  ! LLS models
  type(LLS_model_new),parameter :: low_mfp_LLS=LLS_model_new( &
       reference="W14 mfp low", &
       A_LLS=(37.0-2.0)/(h/0.7), z_ref=4.0, yz_LLS=-5.4-0.4)
  type(LLS_model_new),parameter :: std_mfp_LLS=LLS_model_new( &
       reference="W14 mfp std", &
       A_LLS=37.0/(h/0.7), z_ref=4.0, yz_LLS=-5.4)
  type(LLS_model_new),parameter :: high_mfp_LLS=LLS_model_new(&
       reference="W14 mfp high", &
       A_LLS=(37.0+2.0)/(h/0.7), z_ref=4.0, yz_LLS=-5.4+0.4)
  type(LLS_model_new),parameter :: const_pmfp_LLS=LLS_model_new(&
       reference="constant proper mfp", &
       A_LLS=1.0, z_ref=4.0, yz_LLS=0.0)
  type(LLS_model_new),parameter :: const_cmfp_LLS=LLS_model_new(&
       reference="constant comoving mfp", &
       A_LLS=10.0, z_ref=0.0, yz_LLS=-1.0)

  type(LLS_model_new) :: mfpLLS
  
  real(kind=dp),public :: R_max_LLS
  real(kind=dp),public :: R_max
  

  public :: set_LLS, LLS_point

  
#ifdef MPI
  integer,private :: mympierror
#endif

contains

  ! ============================================================================

  subroutine LLS_init ()

    if (use_LLS) then
       select case (type_of_LLS)
       case(1,2)
          
          ! Choose the appropriate LLS model
          select case (LLS_model)
          case(1)
             mfpLLS=std_mfp_LLS
          case(2)
             mfpLLS=low_mfp_LLS
          case(3)
             mfpLLS=high_mfp_LLS
          case(4)
             mfpLLS=const_pmfp_LLS
          case(5)
             mfpLLS=const_cmfp_LLS
          end select
          
       ! Report
       write(logf,"(A,A)") "Using mean free path model ",mfpLLS%reference

       case(3)
          ! Use a maximum distance from the source, similar to the R_max
          ! implementation in 21cmFAST
          R_max=R_max_cMpc*Mpc
       end select

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
       mfp_LLS_pMpc=mfpLLS%A_LLS*((1.0+z)/(1.0+mfpLLS%z_ref))**mfpLLS%yz_LLS
       mfp_LLS_pMpc=max(mfp_LLS_pMpc,limit_mfp_LLS_cMpc/(1.0+z))
       ! Column density per cell due to LLSs
       n_LLS=dr(1)/(mfp_LLS_pMpc*Mpc)
       coldensh_LLS = N_1 * n_LLS
    case(2) 
       call read_lls_grid (z)
    case(3)
       R_max_LLS=R_max/(1.0+z)
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
    real(kind=dp) :: sim_volume_pMpc

    ! clumping in file is in 4B reals, read in via this array
    !real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(LLS_grid))) &
         allocate(LLS_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       LLS_file=trim(adjustl(dir_LLS))// &
            trim(adjustl(zred_str))// &
            "cross_section_normalized.bin"

       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " LLS input from ",trim(LLS_file)

       ! Open LLS file: note that the format is determined
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
       
       ! Calculate mean free path, this should be 1 pMpc
       ! Make sure sim_volume is in proper length units
       sim_volume_pMpc=sim_volume/(Mpc*(1.0+zred_now))**3
       mfp_LLS_pMpc=sim_volume_pMpc/(sum(LLS_grid)*Mpc*Mpc*Mpc)
       write(logf, *) "Mean free path in file is ",mfp_LLS_pMpc," pMpc"
       if (abs(mfp_LLS_pMpc - 1.0) > 1e-2) write(logf,*) "WARNING: &
            The LLS file contains incorrect values."

       ! Calculate the desired mean free path in pMpc
       mfp_LLS_pMpc=mfpLLS%A_LLS* &
            ((1.0+zred_now)/(1.0+mfpLLS%z_ref))**mfpLLS%yz_LLS

       ! Scale the read-in values to use this mean free path
       LLS_grid=LLS_grid/mfp_LLS_pMpc
       
       ! Do not set LLS column densities if the mfp is too small.
       ! Typically we start using LLS when the mfp is several cells (e.g. 5)
       if (mfp_LLS_pMpc < limit_mfp_LLS_pMpc) then
          LLS_grid(:,:,:)=0.0
       else
          ! Convert units to cgs
          ! fraction of cell covered by LLS (equivalent to n_LLS in 
          ! the routines above)
          LLS_grid(:,:,:)=LLS_grid(:,:,:)*(Mpc/dr(1))*(Mpc/dr(1))
          ! Convert to column density
          LLS_grid(:,:,:)=N_1 * LLS_grid(:,:,:)
       endif
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
