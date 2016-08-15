!>
!! \brief This module contains data and routines for handling subgrid sources
!! of the minihalo (MH) category
!! 
!! \b Author: Kyungjin Ahn, Garrelt Mellema
!!
!! \b Date: 21-June-2016 (imported from Kyungjin's version)
!!
!! \b Version: Derived from Kyungjin's version, modified to be similar
!!  to sourceprops module

module source_sub

  use precision, only: dp, si, i8b
  use my_mpi
  use file_admin, only: logf, results_dir
  use cgsconstants, only: m_p
  use sizes, only: mesh
  use astroconstants, only: M_solar
  use cosmology_parameters, only: Omega_B, Omega0, h
  use nbody, only: id_str, M_grid, dir_src, zred_array
  use material, only: xh, ndens, avg_dens
#ifdef MH  
  use material, only: ndens_previous, avg_dens_previous
#endif
  use grid, only: x,y,z, vol_cMpc3, sim_volume_cMpc3
  use c2ray_parameters, only: phot_per_atom, fstar, &
       lifetime, S_star_nominal, StillNeutral
  use c2ray_parameters, only: emiss, QH_M_real, Ni
#ifdef MH
  use radiation, only: jLW
#endif

  implicit none

#ifdef MH
  ! Model parameters for MH model
  integer, public, parameter :: MHflag = 2 !< Flag indicating type of MH sources

  !> mass of PopIII star in solar mass unit.
  real   (kind=dp), parameter, public :: M_PIIIstar_msun = 300d0 
  !> Name of file with fit parameters
  character(len=100), parameter :: MHfit_file="zred_halodelta1_nMHMpc3_WMAP5"

  integer, public :: NumAGrid !< Number of Active Grid cells with MH sources

  !> Grid with flags, if true there is an active source there, if false not.
  logical,dimension(mesh(1),mesh(2),mesh(3)) :: AGflag

  !> critical nondimensional density over which cell can be minihalo populated,
  !! current and previous redshift
  real(kind=dp) :: densNDcrit, densNDcrit_prev 

  integer,dimension(:,:),allocatable     :: sub_srcpos !< mesh position of subgrid sources
  real(kind=dp),dimension(:),allocatable :: subNormFlux !< normalized ionizing flux of subgrid sources
  real(kind=dp),dimension(:),allocatable :: subNormFlux_LW !< normalized ionizing flux of subgrid sources
  integer,dimension(:),allocatable       :: subSrcSeries !< a randomized list of sources
  
  ! Definitions for reading precalculated # density of minihalos
  integer(kind=i8b),  public :: Nzdata, Ndelta1data  !< LGnMH data size integer header
  real   (kind=dp),  public :: zmin, zmax, LGdelta1min, LGdelta1max  !< LGnMH data domain real*8 header
  real   (kind=dp), public :: dzred, dLGdelta1  !< Needed to search LGnMH data
  real   (kind=dp), dimension(:,:), allocatable, public :: LGnMH_Mpc3 !< log10 of number of MHs per comoving (Mpc)^3, at given redshfit and delta+1, where delta is overdensity.
  real   (kind=dp), public :: tot_subsrcM_msun !< total STELLAR mass in all active (survived after suppression) minihalos

    real(kind=dp)            :: jLWcrit_now, jLWcrit_min

    character(len=512) :: sourcelistfile_sub

contains
  
  ! =======================================================================

  subroutine AGrid_properties(nz,AGlifetime,restart)

    ! Input routine: establish the subgrid source properties
    ! Author: Kyungjin Ahn
    ! Update: 5-Sep-2009, 15-Jul-2010

    ! For random permutation of sources
    !use perm  ! ctrper seems to break if source number too large.

    integer,      intent(in) :: nz
    real(kind=dp),intent(in) :: AGlifetime ! time step
    integer,intent(in) :: restart

    real(kind=dp)            :: zred_prev,zred_now ! previous, current redshift

    character(len=6)   :: z_str !< string value of redshift

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check subgrid props: ",zred_now,nz,AGlifetime
#endif 

    ! Ask for input
    if (rank == 0) then
       
       ! Set values of current and previous redshift
       ! The latter is needed to differentiate mass in time; at the start
       ! of the simulation (nz=1) there is no need to differentiate.
       zred_now  = zred_array(nz)
       if (nz > 1) zred_prev = zred_array(nz-1)

       ! Construct file name for output file with subgrid sources
       write(z_str,'(f6.3)') zred_now
       sourcelistfile_sub=trim(adjustl(results_dir))//&
            trim(adjustl(z_str))//"-"//&
            trim(adjustl(id_str))//"_SUBsources.dat"

       ! Set the critical density (previous and current redshift)
       ! When nz=1 we do not need the previous value
       if (nz > 1) densNDcrit_prev = get_critical_density(zred_prev) 
       densNDcrit = get_critical_density(zred_now)  
       write(logf,*) 'zred, densNDcrit', zred_now, densndcrit
       
       ! Set the current value of the critical LW background
       jLWcrit_now = jLWcrit(zred_now)
       
       ! set marginal value of the LW BG, below this no sources are suppressed
       jLWcrit_min = 0.1d0 * jLWcrit_now 
       
       ! Count the number of active source locations and record
       ! their grid positions
       call establish_number_of_active_MH_sources (zred_now, zred_prev, restart)

    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumAGrid,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    ! Now that we know the number of active sources, allocate source list 
    ! arrays
    if (allocated(sub_srcpos ))  deallocate(sub_srcpos )
    if (allocated(subNormFlux))  deallocate(subNormFlux)
    if (allocated(subNormFlux_LW))  deallocate(subNormFlux_LW)
    allocate(sub_srcpos  (3, NumAGrid))
    allocate(subNormFlux (NumAGrid   ))
    allocate(subNormFlux_LW (NumAGrid   ))
    
    if (rank == 0) then
       ! Fill the source list with values (for the mass)
       call set_MH_sources(zred_now, zred_prev, nz, restart)
       
       ! Turn source list masses into luminosities; photons/atom included in 
       ! subsrcMass
       call assign_uv_luminosities (AGlifetime,nz)
       
       ! Turn source list masses into luminosities; photons/atom included in 
       ! subsrcMass
       call assign_lw_luminosities (AGlifetime,nz)
       
    endif ! end of rank 0 test
    
#ifdef MPI
    ! Distribute source list to the other nodes
    call MPI_BCAST(sub_srcpos, 3*NumAGrid, MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(subNormFlux,  NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(subNormFlux_LW,  NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of Active Grids: ", NumAGrid
#endif
          
    if (rank == 0) then
       
       write(logf,*) 'Subgrid Source lifetime=', AGlifetime/3.1536e13
       write(logf,*) 'Subgrid Total flux= ',sum(subNormFlux)

    endif

    ! Randomize source list
    !call randomize_source_list ()

  end subroutine AGrid_properties

  ! =======================================================================

  function get_critical_density (zred)

    ! Find out the critical dimensionless density, over which cells
    ! can populate minihalos. This is set by matching
    ! the global number density of all minihalos to high-res, small-box
    ! simulation result which resolves minihalos.

    real(kind=dp), intent(in) :: zred

    real(kind=dp) :: get_critical_density

    ! For total number of minihalos in the small box
    integer       :: Ntable, itable
    real(kind=dp) :: NMH_smallbox, NMH_simbox
    real(kind=dp) :: size_smallbox     ! In h^-1 cMpc 
    real(kind=dp) :: volratio    ! (small box volume / simulation box volume)
    real(kind=dp), allocatable :: ztable(:), numMHtable(:)
    real(kind=dp) :: dz_table

    ! For iteration
    integer, parameter :: Nmaxiter = 14
    integer            :: iiter, i, j, k
    real(kind=dp)      :: LGdensND_L, LGdensND_R, LGdensND_N
    real(kind=dp)      :: frac

    ! For fitting
    integer            :: idx_zred, idx_LGdelta1
    real(kind=dp)      :: LGdelta1

    if (MHflag == 1) then
       write(6,*)  'Do something!!!!!!!!!!!!!!!!!!!!!!!'
       stop
    elseif (MHflag == 2) then
       ! Read and find the total number of minihalos in small box
       open(unit=55, file="z_numMH_6.3Mpc_full", status="old")
       read(55,*) size_smallbox  ! in unit of h^-1 cMpc
       read(55,*) Ntable

       allocate(ztable(Ntable))
       allocate(numMHtable(Ntable))

       do itable = 1, Ntable
          read(55,*) ztable(itable), numMHtable(itable)
       enddo
       close(55)

       volratio = (size_smallbox/h)**3/ sim_volume_cMpc3
            !/ (vol_cMpc3*real(mesh(1)*mesh(2)*mesh(3),dp))

       ! Interpolate by finding the nearest neighbor index. Data table given 
       ! should be in the ascending order in z, and has to be in very fine, 
       ! equal size bin. This is the sheer total number of minihalos in the
       ! small box, NOT the number density of minihalos.
       dz_table     = (ztable(Ntable)-ztable(1))/real(Ntable-1, dp)
       NMH_smallbox = numMHtable(int((zred-ztable(1))/dz_table) + 1)
       
       deallocate(ztable)
       deallocate(numMHtable)

       ! Now iterate until we find the right number for the bigbox
       idx_zred     = int((zred - zmin) /dzred) + 1
       ! Very poor extrapolation. Just limit your simulation
       ! within given redshift range.
       if (idx_zred     <= 0         ) idx_zred = 1
       if (idx_zred     > Nzdata     ) idx_zred = Nzdata

       ! bisecting(?) interpolation
       LGdensND_L = -1d0 ! starting L value, should be small enough.
       LGdensND_R = 1d0  ! starting R value, should be large enough.
       iiter      = 0
       do
          iiter = iiter + 1
          LGdensND_N = (LGdensND_R + LGdensND_L)*0.5d0

          NMH_simbox = 0d0
          do k=1,mesh(3)
             do j=1,mesh(2)
                do i=1,mesh(1)
                   if (log10(ndens(i,j,k)/avg_dens) >= LGdensND_N) then
                      idx_LGdelta1 = int((log10(ndens(i,j,k)/avg_dens) - &
                           LGdelta1min) / dLGdelta1) + 1
                      ! Extrapolation, very poor, but only if needed. 
                      if (idx_LGdelta1 <= 0         ) idx_LGdelta1=1
                      if (idx_LGdelta1 > Ndelta1data) idx_LGdelta1=Ndelta1data

                      NMH_simbox = NMH_simbox &
                           + 10d0**LGnMH_Mpc3(idx_zred,idx_LGdelta1)
                   endif
                enddo
             enddo
          enddo

          ! Normalize to the simulation box to get the total number of
          ! minihalos in the simulation box.
          !!!!!!!!! Careful
          ! Because fitting function uses h=0.7, we need to scale further with used h.
          ! SHOULD BE FIXED LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!
          print *, '*************************************'
          print *, 'one box scan -> iiter=', iiter, 'nMH_simbox=', nmh_simbox
          NMH_simbox = NMH_simbox * vol_cMpc3 * (h/0.7d0)**3
          print *, 'one box scan -> iiter=', iiter, 'NMH_simbox=', nmh_simbox
          print *, 'volratio=', volratio
          
          ! fractional difference
          frac = (NMH_simbox*volratio - NMH_smallbox)/NMH_smallbox
          print *, 'frac=', frac
          print *, 'LGdensND_N', LGdensND_N
          print *, '*************************************'
          ! convergence is 1%, chosen just arbitrarily.
          if (abs(frac) <= 0.001d0 .or. iiter == Nmaxiter) exit
          if (NMH_simbox*volratio < NMH_smallbox) then
             LGdensND_R = LGdensND_N
          else
             LGdensND_L = LGdensND_N
          endif
       enddo
       get_critical_density = 10d0**LGdensND_N
    endif   ! end of MHflag=2 case

  end function get_critical_density

  ! =======================================================================

  function jLWcrit(zred)

    ! Get critical jLW at zred
    ! This allows jLWcrit to be redshift dependent but in this version
    ! it is a constant.

    real(kind=dp)             :: jLWcrit
    real(kind=dp), intent(in) :: zred

    jLWcrit = 0.1d0 * 1d-21

  end function jLWcrit

  ! =======================================================================

  subroutine establish_number_of_active_MH_sources (zred_now, zred_prev, restart)

    real(kind=dp),intent(in) :: zred_prev,zred_now ! previous, current redshift
    integer,intent(in) :: restart

    integer :: i, j, k

    ! diff_subsrcMsun is the potential source masses. Identical to removing
    ! previous step contribution to minihalo-sources, due to regulation.
    ! Once a minihalo hosts a Pop III star, it should wait about 100
    ! Myrs (Yoshida et al. 2007, ApJ 663, 687). So better differentiate
    ! like this.
    real(kind=dp)            :: diff_subsrcMsun

    ! Flag grids to activate subgrid physics (or not)
    ! Count the number of cells with MH sources.
    NumAGrid = 0
    AGflag(:,:,:) = .false.
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             diff_subsrcMsun=mass_difference_MH(zred_now, zred_prev, &
                  ndens(i,j,k),ndens_previous(i,j,k))
             ! subgrid is active when neutral & containing minihalos
             ! which are fresh at current redshift!!!!
             ! NOTE THE CRITERION!! XH(x,y,z,1) < StillNeutral
#ifdef ALLFRAC
             if (xh(i,j,k,1) < StillNeutral &
#else
             if (xh(i,j,k) < StillNeutral &
#endif                  
                  .AND. jLW(i,j,k) < jLWcrit_now &
                  .AND. diff_subsrcMsun > 0d0) then
                AGflag(i,j,k) = .true.
                NumAGrid      = NumAGrid + 1
             endif
          enddo
       enddo
    enddo
    
  end subroutine establish_number_of_active_MH_sources

  ! =======================================================================

  subroutine set_MH_sources (zred_now,zred_prev,nz,restart)

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: zred_prev ! previous redshift
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    integer :: i, j, k
    integer :: ns
    integer :: source_counter = 0

    ! diff_subsrcMsun is the potential source masses. Identical to removing
    ! previous step contribution to minihalo-sources, due to regulation.
    ! Once a minihalo hosts a Pop III star, it should wait about 100
    ! Myrs (Yoshida et al. 2007, ApJ 663, 687). So better differentiate
    ! like this.
    real(kind=dp)            :: diff_subsrcMsun

    source_counter = 0
    tot_subsrcM_msun = 0d0
    
    ! Now fill the source list.
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)

             ! If cell is active fill a position in the source list
             if (AGflag(i,j,k)) then
                source_counter = source_counter + 1
                if (source_counter > NumAGrid .and. rank == 0) then
                   write(logf,*) "Error in filling MH source list: ", &
                        "actual source count exceeds previous count"
                   flush(logf)
                endif

                ! Record source position
                sub_srcpos (1, source_counter) = i
                sub_srcpos (2, source_counter) = j
                sub_srcpos (3, source_counter) = k

                ! Determine mass difference between zred_now and zred_prev
                diff_subsrcMsun=mass_difference_MH(zred_now, zred_prev, &
                     ndens(i,j,k),ndens_previous(i,j,k))

                ! Determine source mass using suppression recipe
                if (jLW(i,j,k) <= jLWcrit_min) then
                   ! No suppression, differential mass only.
                   ! Stored in NormFlux variable for now
                   subNormFlux(source_counter) = diff_subsrcMsun
                else 
                   ! suppression by jLW. Use differential mass.
                   ! Stored in NormFlux variable for now
                   subNormFlux(source_counter) = diff_subsrcMsun &
                        * (jLWcrit_now - jLW(i,j,k)) &
                        / (jLWcrit_now - jLWcrit_min)
                endif
             endif
          enddo
       enddo
    enddo

    ! Total source mass of MH sources (stored in NormFlux)
    tot_subsrcM_msun = sum(subNormFlux)

    ! Save source masses in subNormFlux_LW before modifying subNormFlux
    subNormFlux_LW(:) = subNormFlux(:)

    ! Convert subgrid source masses into grid mass unit, and then
    ! multiply phot_per_atom(3). Result is similar to srcMass under
    ! source_properties (but see below for difference when MHflag=2, 
    ! need to divide by fstar(3), since stellar mass is used then.). 
    ! Use it just internally.
    ! We divide by M_grid here for consistency with the processed source
    ! list produced by sourceprops.
    ! See below in assign_uv_luminosities for how these quantities
    ! are further used to produce EUV photon rates.
    if (MHflag == 1) then
       subNormFlux(:) = subNormFlux(:) * (M_solar/M_grid) * phot_per_atom(3)
    elseif (MHflag == 2) then
       subNormFlux(:) = subNormFlux(:) * (M_solar/M_grid) * phot_per_atom(3) &
            / fstar(3)
    endif
    
    ! Save new subgrid source list, just as in sourceprops.
    open(unit=49,file=sourcelistfile_sub,status='unknown')
    write(49,*) NumAGrid
    ! Write MHflag to distinguish physical meaning of
    ! subNormFlux and subNormFlux_LW. When MHflag=2, subNormFlux_LW here
    ! contains the STELLAR BARYON MASS, not BARYON+DM HALO MASS.
    write(49,*) MHflag
    do ns = 1,NumAgrid
       write(49,"(3I5,2e13.4)") sub_srcpos(1,ns), sub_srcpos(2,ns), &
            sub_srcpos(3,ns), subNormFlux(ns), subNormFlux_LW(ns)
    enddo
    close(49)
    

    write(logf,*) "Number of active MH grid cells: ", NumAGrid
       
  end subroutine set_MH_sources

  ! =======================================================================

  function mass_difference_MH (zred_now, zred_prev, nd ,nd_prev)

    real(kind=dp),intent(in) :: zred_prev,zred_now ! previous, current redshift
    real(kind=si),intent(in) :: nd, nd_prev ! current, previous density

    real(kind=dp) :: mass_difference_MH
 
    ! When nz=1, no need to differentiate mass in time!!!!
    if ( zred_now < zred_array(1) ) then
       mass_difference_MH = subsrcM_msun(zred_now, nd/avg_dens, densNDcrit) &
            - subsrcM_msun(zred_prev, nd_prev/avg_dens_previous, &
            densNDcrit_prev)
    else
       mass_difference_MH = subsrcM_msun(zred_now, nd/avg_dens, densNDcrit)
    endif

  end function mass_difference_MH

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime,nz)
    
    real(kind=dp),intent(in) :: lifetime
    integer,intent(in) :: nz

    ! The quantity subNormFlux is the mass (in grid masses) multiplied
    ! with the efficiency factor. Here it is converted into ionizing
    ! photon rate for the source.
    if (MHflag == 1) then
       ! The MH source is a photon efficiency x the mass of the MH:
       ! Mass * phot_per_atom(3) / m_p * Omega_B/Omega0 / lifetime
       ! subNormFlux was previously set to Mass * phot_per_atom(3)
       ! in grid masses.
       subNormFlux(:) = subNormFlux(:) * M_grid * Omega_B/(Omega0*m_p) &
            /S_star_nominal /lifetime
    elseif (MHflag == 2) then
       ! The MH source is a single PopIII star. Its EUV photon rate is
       ! given by
       ! Mass * phot_per_atom(3) / fstar(3) / m_p / lifetime
       ! equivalent to
       ! Mass * Ni(3) * fesc(3) / m_p / lifetime
       ! subNormFlux was previously set to Mass * phot_per_atom(3) / fstar(3)
       ! in grid masses.
       subNormFlux(:) = subNormFlux(:) * M_grid          /        m_p  &
            /S_star_nominal /lifetime
    endif
    
  end subroutine assign_uv_luminosities

  ! =======================================================================

  subroutine assign_lw_luminosities (lifetime,nz)
    
    real(kind=dp),intent(in) :: lifetime
    integer,intent(in) :: nz

    ! Calculate the LW photon rate for MH sources.

    ! There are two options here for different source models. MHflag = 2
    ! appears to be the standard one.
    ! Different subgrid source population schemes reflected by MHflag
    ! MHflag = 1: uniform star formation out of given baryon content
    !          2: one Pop III star/ one minihalo
    ! The quantity subNormFlux_LW contains the total mass of the MH (MHflag=1)
    ! or the stellar mass inside the MH (MHflag=2). Both are in solar masses.
    if (MHflag == 1) then
       ! For MHflag=1 the quantity
       ! subNormFlux_LW(:) * M_solar / m_p * Omega_B/Omega0 * fstar(3)
       ! is the total number of baryons inside the halo in stars
       subNormFlux_LW(:) = subNormFlux_LW(:) * M_solar / m_p * &
            Omega_B/Omega0 * fstar(3)
    elseif (MHflag == 2) then
       ! For MHflag=1 the quantity
       ! subNormFlux_LW(:) * M_solar / m_p 
       ! is the total number of baryons inside the single star in the MH
       subNormFlux_LW(:) = subNormFlux_LW(:) * M_solar / m_p
    endif

    ! Next we convert this number of baryons into a LW luminosity.
    ! - Ni(3) is the number of ionizing photons produced per baryon.
    ! - emiss is the LW emissivity (erg s^-1 Hz^-1 Msun^-1)--- THESE
    !    UNITS DO NOT MAKE SENSE!
    ! - QH_M_real(3) is the number of LW photons produced by a solar
    !    mass (????) in the life time of the population.
    subNormFlux_LW(:) = subNormFlux_LW(:) * Ni(3) * emiss(3) / &
         QH_M_real(3) / lifetime
    
  end subroutine assign_lw_luminosities

  ! =======================================================================

  function subsrcM_msun(zred, densND, densNDcrit)
    
    ! Get subgrid src mass (those not resolved in N-body) in solar mass unit
    ! Cell characteristic given by zred & dimensionless density

    real(kind=dp), intent(in) :: zred
    real(kind=si), intent(in) :: densND
    real(kind=dp), intent(in) :: densNDcrit

    real(kind=dp)             :: subsrcM_msun

    integer                   :: idx_zred, idx_LGdelta1
    real(kind=dp)             :: LGdelta1
    real(kind=dp)             :: N_MH


    if (MHflag == 1) then
       write(6,*)  'Do something!!!!!!!!!!!!!!!!!!!!!!!'
       subsrcM_msun = 0d0
    elseif (MHflag == 2) then
       LGdelta1 = log10(dble(densND))
       if (densND <= 0d0) LGdelta1 = -5  ! for numerical safety.

       idx_zred     = int((zred     - zmin       ) /dzred    ) + 1
       ! Very poor extrapolation. Just limit your simulation
       ! within given redshift range.
       if (idx_zred     <= 0         ) idx_zred = 1
       if (idx_zred     > Nzdata     ) idx_zred = Nzdata

       idx_LGdelta1 = int((LGdelta1 - LGdelta1min) /dLGdelta1) + 1
       ! Extrapolation. LGdelta1 covers -1 to 1, which is not too bad.
       if (idx_LGdelta1 <= 0         ) idx_LGdelta1 = 1
       if (idx_LGdelta1 > Ndelta1data) idx_LGdelta1 = Ndelta1data

       ! densNDcrit is CHOSEN to match the total number of minihalos in
       ! small-box, high-res simulation. So N_MH can be sometimes less than 1,
       ! and more generally fractional numbers.
       if (densND >= densNDcrit) then
          ! number of minihalos contained in the cell, with volume vol_cMpc3 
          ! in comoving Mpc^3 unit.
          ! number of minihalos on the cell with given volume
          ! !!!!!!!! Careful
          ! Because fitting function uses h=0.7, we need to scale further with used h.
          ! SHOULD BE FIXED LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!
          N_MH  = 10d0**LGnMH_Mpc3(idx_zred,idx_LGdelta1) * vol_cMpc3 &
               * (h/0.7d0)**3

          ! Total mass of stellar baryon in the cell.
          subsrcM_msun = N_MH * M_PIIIstar_msun
       else
          subsrcM_msun = 0d0
       endif
    endif

  end function subsrcM_msun

  ! =======================================================================

  subroutine read_LGnMH_Mpc3

#ifdef IFORT
  ! ifort standard for "binary"
  character(len=*),parameter :: tableformat="binary"
  character(len=*),parameter :: tableaccess="sequential"
#else
  ! Fortran2003 standard for "binary"
  character(len=*),parameter :: tableformat="unformatted"
  character(len=*),parameter :: tableaccess="stream"
#endif

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then
       open(unit=20, file=trim(adjustl(MHfit_file)), form=tableformat, &
            access=tableaccess, status="old")
       read(20) Nzdata, Ndelta1data  
       read(20) zmin, zmax, LGdelta1min, LGdelta1max 
       if (rank==0) write(6,*)  '.................'
       if (rank==0) write(6,*)  nzdata, ndelta1data, zmin, zmax, LGdelta1min, LGdelta1max
       
       if (allocated(LGnMH_Mpc3)) deallocate(LGnMH_Mpc3)
       allocate(LGnMH_Mpc3(Nzdata, Ndelta1data))

       read(20) LGnMH_Mpc3
       close(20)

       dzred     = (zmax       -zmin       ) /real(Nzdata     -1, 8)
       dLGdelta1 = (LGdelta1max-LGdelta1min) /real(Ndelta1data-1, 8)
    endif

#ifdef MPI
    call MPI_BCAST(Nzdata,     1,                  MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(Ndelta1data,1,                  MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) then
       if (allocated(LGnMH_Mpc3)) deallocate(LGnMH_Mpc3)
       allocate(LGnMH_Mpc3(Nzdata, Ndelta1data))
    endif
    call MPI_BCAST(zmin,       1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(zmax,       1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGdelta1min,1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGdelta1max,1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGnMH_Mpc3, Nzdata*Ndelta1data, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine read_LGnMH_Mpc3

  ! =======================================================================

  subroutine randomize_source_list ()

    ! For random permutation of sources
    use  m_ctrper

    integer :: ns

    ! Create array for containing randomized source list
    if (allocated(subSrcSeries)) deallocate(subSrcSeries)
    allocate(subSrcSeries( NumAGrid  ))

    if (rank == 0) then

       ! Create array of subgrid source numbers for generating random order
       do ns=1,NumAGrid
          subSrcSeries(ns)=ns
       enddo

       ! Make a random order
       call ctrper(subSrcSeries(1:NumAGrid),1.0)
       !call permi(NumAGrid, subSrcSeries)
       
    endif

#ifdef MPI
    ! Distribute the subgrid source series to the other nodes
    call MPI_BCAST(subSrcSeries,NumAGrid,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine randomize_source_list

#endif

end module source_sub
