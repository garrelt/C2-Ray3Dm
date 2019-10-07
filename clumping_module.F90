module clumping_module

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use c2ray_parameters, only: type_of_clumping, clumping_factor
  use nbody, only: nbody_type, id_str, dir_dens, NumZred, Zred_array, dir_clump, dir_LLS
  use nbody, only: clumpingformat, clumpingaccess, clumpingheader
  use grid, only: dr
  use density_module, only: ndens
  use mathconstants, only: pi

  implicit none

  ! Clumping data
  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  real(kind=dp) :: avg_dens !< average density
  character(len=512) :: clumping_fit_file
  public :: set_clumping, clumping_point
  real(kind=dp), allocatable, public :: params_scm(:,:,:), params_dcm(:,:), params_gcm(:)
  real(kind=dp), dimension(:), allocatable, public :: redshift_clumping

#ifdef MPI
  integer,private :: mympierror
#endif

contains
  
  ! ===========================================================================

  subroutine load_clumping_model(r)
    !> Set up the clumping model parameters according to the type of clumping specified:
    !! 1: Constant clumping (with clumping_factor)
    !! 2: Globally Averaged Clumping Model GCM
    !! 3: Deterministic Clumping Model DCM (Mao et al. 2019)
    !! 4: Stochastic Clumping Model SCM (Bianco et al. 2020)
    !! 5: Pre-computed grid of clumping

    real(kind=dp), intent(in) :: r  ! resolution of the simulation

    select case (type_of_clumping)
    case(1)
       ! write to logf to inform that you are using this type of clumping
    case(2) 
       call load_GCM_parameters(r)
    case(3)
       call load_DCM_parameters(r)
    case(4)
       call load_SCM_parameters(r)
    case(5)
       ! write to logf to inform that you are using this type of clumping
       ! and that you will use a pre-computed clumping factor grid
    end select

    if(rank == 0) write(logf,*) "Loaded sub-grid clumping factor model, type (", type_of_clumping,")"
      
  end subroutine load_clumping_model
 
  ! ===========================================================================

  subroutine set_clumping(zred)
    !! 1: Constant clumping (with clumping_factor)
    !! 2: Globally Averaged Clumping Model GCM
    !! 3: Deterministic Clumping Model DCM (Mao et al. 2019))
    !! 4: Stochastic Clumping Model SCM (Bianco et al. 2020)
    !! 5: Pre-computed grid of clumping

    real(kind=dp), intent(in) :: zred

    select case (type_of_clumping)
    case(1)
       clumping = clumping_factor
    case(2) 
       clumping = params_gcm(1)*exp(params_gcm(2)*zred+params_gcm(3)*zred*zred)
    case(3)
       if(rank == 0) call deterministic_clumping(zred)
    case(4)
       if(rank == 0) call stochastic_clumping(zred)
    case(5)
       if(rank == 0) call read_clumping_file(zred)
    end select

    ! Report on data: min, max, total assign mean to clumping for reporting in set_clumping, etc...
    if(rank == 0 .and. type_of_clumping==1) then
       write(logf,*) "Setting constant clumping factor to ",clumping," (type ", type_of_clumping,")"
    else if(rank == 0 .and. type_of_clumping==2) then
       write(logf,*) "Setting (mean) global clumping factor to ",clumping," (type ", type_of_clumping,")"
    else if (rank == 0 .and. (type_of_clumping==3 .or. type_of_clumping==4)) then
       write(logf,*) "Setting inhomogeneity dependent clumping factor (type ", type_of_clumping,")"
       clumping = sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) "Statistics AFTER applying such clumping fit:"
       write(logf,*) "clumping fit for redshift ", zred
       write(logf,*) "minimum: ", minval(clumping_grid)
       write(logf,*) "maximum: ", maxval(clumping_grid)
       write(logf,*) "average clumping: ", clumping
    else if(rank == 0 .and. type_of_clumping==5) then
       write(logf,*) "Setting pre-computed clumping factor grid (type ", type_of_clumping,")"
    end if

  end subroutine set_clumping

  ! ===========================================================================

  subroutine clumping_point(i,j,k)
    !> Read clumping grid for position dependent model and allcate value to 
    !! the clumping variable (case = 3, 4 or 5)
    !! 
    integer,intent(in) :: i,j,k

    if (type_of_clumping == 1 .or. type_of_clumping == 2) then
       write(logf,*) &
         "Error: calling position dependent, but array is not initialized."
    else
       clumping = clumping_grid(i,j,k)
    endif
  end subroutine clumping_point

  ! ================ LOAD CLUMPING PARAMETER FILES SUBROUTINES ================

  subroutine load_GCM_parameters(res)
    !> Load precomputed parameters for the Globally Averaged Clumping Model
    !! parameters files should be stored in dir_clump, files should follow the naming system:
    !!   >>  "paramsGCM_[resolution]Mpc.dat"
    !! Table are then stored to a 1dim array with for the following fit (in order): 
    !!   >>  C(z) = C0 * exp(c1*zred + c2*z²) + 1
    !! (in addition, last three variable are error_C0, error_c1, error_c2)
    !!
    !! Author:  Michele Bianco, modified from code by Garrelt Mellema
    !! Date:    04-Oct-2019
    !! 
    implicit none
    real(kind=dp), intent(in) :: res  ! resolution of the simulation

    ! write file name for DCM parameters file
    write(clumping_fit_file, "(a,a,f5.3,a)") trim(dir_clump), "paramsGCM_", res, "Mpc.dat"

    ! open parameter files
    open(unit=13, file=trim(clumping_fit_file), form="unformatted")

    ! read shape of the parameter files and allocate arrays
    allocate(params_gcm(6))

    ! read parameters and save to arrays
    read(13) params_gcm
    close(13)

    write(logf, "(a)") "set parameters for Globally Averaged Clumping Model, fitting: C(zred) = C0 * exp(c1*zred + c2*zred^2) + 1"
    write(logf, "(a,f6.2,a,f5.3)") "parameters are: C0=", params_gcm(1), ", c1=", params_gcm(2)
    write(logf, "(a,f6.5)") " and c2=", params_gcm(3)
  end subroutine load_GCM_parameters

  subroutine load_DCM_parameters(res)
    !> Load precomputed parameters for the Deterministic Clumping Model (Y. Mao et al. 2019)
    !! parameters files should be stored in dir_params, files should follow the naming system:
    !!   >>  "paramsDCM_[resolution]Mpc.dat"
    !! Table are then stored to a 2dim array with the follwing structure:
    !!   >>  params_scm(redshift, quadratic variable)
    !! Quadratic variables should be in order for the quadratic fit: y = a0*x² + a1*x + a2
    !! (first column is redshift and in addition, last three variable are error_a, error_b, error_c)
    !!
    !! Author:  Michele Bianco, modified from code by Garrelt Mellema
    !! Date:    04-Oct-2019
    !! 
    implicit none
    real(kind=dp), intent(in) :: res  ! resolution of the simulation
    integer :: q1, q2

    ! write file name for DCM parameters file
    write(clumping_fit_file, "(a,a,f5.3,a)") trim(dir_clump), "paramsDCM_", res, "Mpc.dat"

    ! open parameter files
    open(unit=13, file=trim(clumping_fit_file), form="unformatted")

    ! read shape of the parameter files and allocate arrays
    read(13) q1, q2
    allocate(params_dcm(q1,q2))

    ! read parameters and save to arrays
    read(13) params_dcm
    close(13)

    write(logf, "(a,i2,a)") "set parameters for Deterministic Clumping Model: ", q1, " redshifts"
    write(logf, "(a,i2,a)") "and ", q2, " parameters"
  end subroutine load_DCM_parameters

  subroutine load_SCM_parameters(res)
    !> Load precomputed parameters for the Stochastic Clumping Model (Bianco et al. 2020)
    !! parameters files should be stored in dir_params, files should follow the naming system:
    !!   >>  "paramsSCM_[resolution]Mpc.dat"
    !! table are then stored to a 3dim array with the follwing structure:
    !!   >>  params_scm(redshift, bin_nr, lognorm variable) 
    !! Lognorm variables should be employed for the fitting to a Lognormal distribution and
    !! in following order: mu, sigma, low, upper and middle value of the density bin
    !!
    !! Author:  Michele Bianco, modified from code by Garrelt Mellema
    !! Date:    04-Oct-2019
    !! 
    implicit none
    real(kind=dp), intent(in) :: res  ! resolution of the simulation
    integer :: l1, l2, l3

    ! load DCM for clumping above or below SCM bins limit, see (Bianco et al. 2020) for more infos
    call load_DCM_parameters(res)

    ! write file name for SCM parameters file
    write(clumping_fit_file, "(a,a,f5.3,a)") trim(dir_clump), "paramsSCM_", res, "Mpc.dat"

    ! open parameter files
    open(unit=12, file=trim(clumping_fit_file), form="unformatted")

    ! read shape of the parameter files and allocate arrays
    read(12) l1, l2, l3
    allocate(params_scm(l1,l2,l3))

    ! read parameters and save to arrays
    read(12) params_scm
    close(12)

    write(logf, "(a,i2,a)") "set parameters for Stochastic Clumping Model for: ", l1, " redshifts, "
    write(logf, "(i2,a,i2,a)") l2, " density bins and ", l3, " lognormal parameters."
  end subroutine load_SCM_parameters

  ! ==================== TOOLS FOR CLUMPING FACTOR METHOD =====================

  function find_nearest(Ndim, arr, val) result(idx)
      !> Look for the closest element in array to a value of choice
      !! Variable:
      !!     > Ndim (int)        : array dimension
      !!     > val (real)        : value to look for
      !!     > arr (1D array)    : array to look
      !! Return:
      !!     > idx (int)     : index of the array element closest to the value
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    04-Oct-2019
      !! 
      implicit none
      integer, intent(in) :: Ndim
      real(kind=dp), intent(in) :: val
      integer :: idx
      real(kind=dp), dimension(Ndim) :: diff_arr, arr, val_arr
      val_arr(:) = val
      diff_arr = abs(arr - val_arr(:))
      idx = minloc(array=diff_arr, dim=1)
  end function find_nearest

  subroutine weight_function(Ndim, arr, val, weight1, weight2, idx)
      !> Weight value for an linear interpolation of parameters.
      !! Here weight1 is the weight assosiacted to the value at position (idx).
      !! While weight2 is associated to the value in position (idx+1) in arr.
      !! i.e: weight1 --> idx, weight2 --> idx+1
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    04-Oct-2019
      !! 
      implicit none
      integer, intent(in) :: Ndim
      real(kind=dp), intent(in) :: val, arr(Ndim)
      real(kind=dp), intent(out) :: weight1, weight2
      integer, intent(out) :: idx
      real(kind=dp) :: dum_var, copy_arr(Ndim)
      integer :: i
      
      copy_arr = arr
      if(arr(1) < arr(Ndim)) then
         ! array must to be in decreasing order to work    
         do i = 1, Ndim
            copy_arr(Ndim+1-i) = arr(i)
         end do
      end if

      ! find element in array closest in value to val
      idx = find_nearest(Ndim, copy_arr, val)

      ! assign weight depending on val 
      if(abs(val - copy_arr(idx)) < 1e-6) then
         ! true when the value correspond to an element of the array
         if(idx == Ndim) then
            weight1 = 0.0
            weight2 = 1.0
            idx = idx-1
         else
            weight1 = 1.0
            weight2 = 0.0
         end if
      else if(val < copy_arr(idx)) then
         ! true when the (idx+1)th is the second closer value to val in arr
         weight2 = abs(val - copy_arr(idx))/abs(copy_arr(idx+1) - copy_arr(idx))
         weight1 = 1.0-weight2
         !print "(i2, f10.6)", idx, arr(idx)
      else if(val > copy_arr(idx)) then
         ! true when the (idx-1)th is the second closer value to val in arr
         weight1 = abs(val - copy_arr(idx))/abs(copy_arr(idx-1) - copy_arr(idx))
         weight2 = 1.0-weight1
         idx = idx-1     ! scale the index of one to respect (w1 -> idx and w2 -> idx+1)
      end if

      if(arr(1) < arr(Ndim)) then
         ! idx and weights need to be adapted to the 
         idx = Ndim-idx
         dum_var = weight2
         weight2 = weight1
         weight1 = dum_var
      end if
  end subroutine weight_function

  function BoxMuller(mu, sigma) result(x0)
      !> Modivied version of the Box-Müller method to generate
      !! indipendent, pseudo-random lognormally distributed numbers
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    04-Oct-2019
      !! 
      real(kind=dp) :: mu, sigma, u0, u1, z0, x0
      call random_seed()
      call random_number(u0)
      call random_number(u1)
      
      z0 = sqrt(-2.0*log(u0))*cos(2.0*pi*u1)
      x0 = exp(z0*sigma + mu)
  end function BoxMuller

  ! ====================== INHOMOGENEITY DEPENDENT MODELS ======================

  subroutine deterministic_clumping(zred)
      !> Fit density field (must be IGM density) to a quadratic equation to infer clumping factor 
      !! for large coarse box. 
      !! The precomputed table are callulated on the simulation:
      !!  >> 6.3 Mpc/h, mesh=3456^3, Npart=1728^3, mpart=5x10^8Msun and RT-mesh=1200^3
      !! This approach is based on Mao et al. 2019 (see for more details).
      !! The method adopted here is the same as presented by Bianco et al. 2020
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    04-Oct-2019
      !! 
      implicit none
      real(kind=dp), intent(in) :: zred
      real(kind=dp) :: w1, w2, paramsdcm(3)
      integer :: i, j, k, idx_z, mesh(3)
      
      ! mesh size of the density and clumping grid
      mesh = shape(ndens)

      ! get weights for redshift
      call weight_function(size(redshift_clumping), redshift_clumping, zred, w1, w2, idx_z)

      ! redshift weighted parameters
      paramsdcm(:) = params_dcm(idx_z,2:4)*w1 + params_dcm(idx_z+1,2:4)*w2

      ! loop over density grid
      do i = 1, mesh(1)
         do j = 1, mesh(2)
            do k = 1, mesh(3)
                  clumping_grid(i,j,k) = paramsdcm(1) * ndens(i,j,k)/avg_dens * ndens(i,j,k)/avg_dens &
                     + paramsdcm(2) * ndens(i,j,k)/avg_dens + paramsdcm(3)
            end do
         end do
      end do

#ifdef MPI 
      ! Distribute the clumping to the other nodes
      call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3),MPI_REAL,0,MPI_COMM_NEW,mympierror)
      call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
  end subroutine deterministic_clumping

  subroutine stochastic_clumping(zred)
      !> Stochastic process to account for scatter when infering clumping factor from the IGM density. 
      !! The precomputed table are callulated on the simulation:
      !!  >> 6.3 Mpc/h, mesh=3456^3, Npart=1728^3, mpart=5x10^8Msun and RT-mesh=1200^3
      !! The method adopted here is firstly presented and explained by Bianco et al. 2020
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    04-Oct-2019
      !! 
      implicit none
      real(kind=dp), intent(in) :: zred
      real(kind=dp) :: w1, w2, paramsscm(size(params_scm(1,:,1)), 5), pardcm(3)
      real(kind=dp) :: wb1, wb2, mu, sigm
      integer :: i, j, k, idx_z, idx_bin, mesh(3)
      
      ! mesh size of the density and clumping grid
      mesh = shape(ndens)

      ! get weights for redshift
      call weight_function(size(redshift_clumping), redshift_clumping, zred, w1, w2, idx_z)

      ! redshift weighted parameters
      paramsscm(:,:) = params_scm(idx_z,:,:)*w1 + params_scm(idx_z+1,:,:)*w2
      print "(5f9.5)", paramsscm
      print *,

      ! loop over density grid
      do i = 1, mesh(1)
         do j = 1, mesh(2)
            do k = 1, mesh(3)
                  if(ndens(i,j,k)/avg_dens < minval(paramsscm(:,5))) then
                     ! redshift weighted parameters for DCM
                     pardcm(:) = params_dcm(idx_z,2:4)*w1 + params_dcm(idx_z+1,2:4)*w2

                     mu = log(pardcm(1) * ndens(i,j,k)/avg_dens * ndens(i,j,k)/avg_dens &
                        + pardcm(2) * ndens(i,j,k)/avg_dens + pardcm(3))

                     sigm = paramsscm(minloc(array=paramsscm(:,5), dim=1),2)

                  else if(ndens(i,j,k)/avg_dens > maxval(paramsscm(:,5))) then
                     ! redshift weighted parameters for DCM
                     pardcm(:) = params_dcm(idx_z,2:4)*w1 + params_dcm(idx_z+1,2:4)*w2

                     ! bin weight mu and sigma for lognormal distribution (in this case )
                     mu = log(pardcm(1) * ndens(i,j,k)/avg_dens * ndens(i,j,k)/avg_dens &
                        + pardcm(2) * ndens(i,j,k)/avg_dens + pardcm(3))
                     
                     sigm = paramsscm(maxloc(array=paramsscm(:,5), dim=1),2)
                  
                  else
                     ! get weights for density bin
                     call weight_function(size(paramsscm(:,5)), paramsscm(:,5), real(ndens(i,j,k)/avg_dens, 8), wb1, wb2, idx_bin)
                  
                     ! bin weight mu and sigma for lognormal distribution
                     mu = paramsscm(idx_bin,1) * wb1 + paramsscm(idx_bin+1,1) * wb2
                     sigm = paramsscm(idx_bin,2) * wb1 + paramsscm(idx_bin+1,2) * wb2
                  end if

                  ! generate lognormal random clumping value
                  clumping_grid(i,j,k) = BoxMuller(mu, sigm)
            end do
         end do
      end do

#ifdef MPI       
      ! Distribute the clumping to the other nodes
      call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3),MPI_REAL,0,MPI_COMM_NEW,mympierror)
      call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
  end subroutine stochastic_clumping

  ! ===========================================================================

  subroutine read_clumping_file(zred)
      !> Initializes position dependent clumping for redshift zred from precomputed clumping grid.
      !! Note: nomenclature has the format from the sub-grid code (Bianco et al. 2020)
      !! i.e: "[zred]_scat.dat"   
      !!
      !! Author:  Michele Bianco, modified from code by Garrelt Mellema
      !! Date:    26-Sep-2019
      !! 
      implicit none
      real(kind=dp), intent(in) :: zred
      character(len=512) :: filename, zstring
      integer :: m1, m2, m3

      ! write file name depending on the redshift value
      write(zstring, "(f6.3)") zred 
      filename = trim(trim(adjustl(dir_clump))//trim(adjustl(zstring))//"_scat.dat")

      !> Open precomputed clumping files
      !! For files calculated with the sub-grid code (Bianco et al. 2020)
      !! variables for open are: form="unformatted", access="stream"
      open(unit=20, file=trim(filename), form=clumpingformat, access=clumpingaccess, status="old")
      ! Read in data header
      if (clumpingheader) then
         read(20) m1, m2, m3
         if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
            write(logf,*) "Warning: file with clumpiness unusable"
            write(logf,*) "mesh found in file: ", m1, m2, m3
            stop
         endif
      endif

      ! Allocate clumping array
      allocate(clumping_grid(m1, m2, m3))

      ! Read in data and store it in clumping_grid
      read(20) clumping_grid
      close(20)

      write(unit=logf,fmt="(2A)") "Reading clumping input from ", trim(filename)

#ifdef MPI       
      ! Distribute the clumping to the other nodes
      call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3),MPI_REAL,0,MPI_COMM_NEW,mympierror)
      call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
  end subroutine read_clumping_file

  ! ===========================================================================

end module clumping_module
