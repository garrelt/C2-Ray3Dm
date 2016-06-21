module source_sub

  use precision, only: dp, si, i8b
  use my_mpi
  use file_admin, only: logf, results_dir
  use cgsconstants, only: m_p
  use sizes, only: mesh
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0, h
  use nbody, only: id_str, M_grid, dir_src, zred_array
  use material, only: xh, ndens, ndens_previous, avg_dens, avg_dens_previous
  use grid, only: x,y,z, vol_cMpc3
  use c2ray_parameters, only: phot_per_atom, fstar, &
       lifetime, S_star_nominal, StillNeutral

  implicit none

  integer, public :: NumAGrid !< Number of Active Grids having subgrid sources
  integer, public, parameter :: MHflag = 2
!!$  integer, public, parameter :: MHflag = 0

  integer,dimension(mesh(1),mesh(2),mesh(3)) :: AGflag

  real(kind=dp) :: densNDcrit, densNDcrit_prev !< critical nondimensional density over which cell can be minihalo populated

  integer,dimension(:,:),allocatable     :: sub_srcpos !< mesh position of subgrid sources
  real(kind=dp),dimension(:,:),allocatable :: rsub_srcpos !< grid position of subgrid sources
  real(kind=dp),dimension(:),allocatable :: ssM_msun !< array of mass of subgrid sources (one source) in solar mass unit. Similar to sM00_msun.
  real(kind=dp),dimension(:),allocatable :: subsrcMass !< masses of subgrid sources. Similar to srcMass.
  real(kind=dp),dimension(:),allocatable :: subNormFlux !< normalized ionizing flux of subgrid sources
  integer,dimension(:),allocatable       :: subSrcSeries !< a randomized list of sources
  
  ! Definitions for reading precalculated # density of minihalos
  integer(kind=i8b),  public :: Nzdata, Ndelta1data  !< LGnMH data size integer header
  real   (kind=dp),  public :: zmin, zmax, LGdelta1min, LGdelta1max  !< LGnMH data domain real*8 header
  real   (kind=dp), public :: dzred, dLGdelta1  !< Needed to search LGnMH data
  real   (kind=dp), dimension(:,:), allocatable, public :: LGnMH_Mpc3 !< log10 of number of MHs per comoving (Mpc)^3, at given redshfit and delta+1, where delta is overdensity.
  real   (kind=dp), parameter, public :: M_PIIIstar_msun = 300d0 !< mass of PopIII star in solar mass unit.
  real   (kind=dp), public :: tot_subsrcM_msun !< total STELLAR mass in all active (survived after suppression) minihalos
  character(len=100), parameter :: MHfit_file="zred_halodelta1_nMHMpc3_WMAP5"

contains
  
  ! =======================================================================

  subroutine AGrid_properties(nz,AGlifetime,jLWgrid)

    ! Input routine: establish the subgrid source properties
    ! Author: Kyungjin Ahn
    ! Update: 5-Sep-2009, 15-Jul-2010

    ! For random permutation of sources
    use perm  ! ctrper seems to break if source number too large.

    integer,      intent(in) :: nz
    real(kind=dp),intent(in) :: AGlifetime ! time step
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)),intent(in) :: jLWgrid

    real(kind=dp)            :: zred_prev,zred_now ! previous, current redshift

    real(kind=dp)            :: jLWcrit_now, jLWcrit_min

    ! diff_subsrcMsun is the potentail source masses. Identical to removing
    ! previous step contribution to minihalo-sources, due to regulation.
    ! Once a minihalo hosts a Pop III star, it should wait about 100
    ! Myrs (Yoshida et al. 2007, ApJ 663, 687). So better differentiate
    ! like this.
    real(kind=dp)            :: diff_subsrcMsun

    character(len=512) :: sourcelistfile_sub
    character(len=6)   :: z_str !< string value of redshift

    integer :: ns

    integer :: i, j, k

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check subgrid props: ",zred_now,nz,AGlifetime
#endif 

    ! Ask for input
    if (rank == 0) then
       
       ! When nz=1, no need to differentiate mass in time!!!!
       if (nz > 1) zred_prev = zred_array(nz-1)
       zred_now  = zred_array(nz)
       ! output file for subgrid sources
       write(z_str,'(f6.3)') zred_now
       sourcelistfile_sub=trim(adjustl(results_dir))//&
            trim(adjustl(z_str))//"-"//&
            trim(adjustl(id_str))//"_SUBsources.dat"

       ! When nz=1, no need to differentiate mass in time!!!!
       if (nz > 1) call get_denscrit(zred_prev, densNDcrit_prev)
       call get_denscrit(zred_now,  densNDcrit)
       write(logf,*) 'zred, densNDcrit', zred_now, densndcrit
       write(6,*) 'zred, densNDcrit', zred_now, densndcrit
       
       jLWcrit_now = jLWcrit(zred_now)
       
       !set marginal value under which no sources are suppressed
       jLWcrit_min = 0.1d0 * jLWcrit_now 
       
       ! Flag grids to activate subgrid physics (or not)
       NumAGrid         = 0
       AGflag           = 0
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                ! When nz=1, no need to differentiate mass in time!!!!
                if (nz > 1) then
                   diff_subsrcMsun = &
                        subsrcM_msun(zred_now, &
                        ndens(i,j,k)/avg_dens,densNDcrit) &
                        - subsrcM_msun(zred_prev, &
                        ndens_previous(i,j,k)/avg_dens_previous,    &
                        densNDcrit_prev)
                else
                   diff_subsrcMsun= &
                        subsrcM_msun(zred_now,ndens(i,j,k)/avg_dens,densNDcrit)
                endif
!!$                ! subgrid is active when neutral & containing minihalos
!!$                ! NOTE THE CHANGE OF CRITERION!! XH(x,y,z,0) > StillNeutral
!!$                ! instead of XH(x,y,z,1) < StillNeutral
!!$                if (xh(i,j,k,0) > StillNeutral &
!!$                     .AND. jLWgrid(i,j,k) < jLWcrit_now &
!!$                     .AND. subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit) > 0d0) then
                ! subgrid is active when neutral & containing minihalos
                ! which are fresh at current redshift!!!!
                ! NOTE THE CRITERION!! XH(x,y,z,1) < StillNeutral
                if (xh(i,j,k,1) < StillNeutral &
                     .AND. jLWgrid(i,j,k) < jLWcrit_now &
                     .AND. diff_subsrcMsun > 0d0) then
                   AGflag(i,j,k) = 1
                   NumAGrid      = NumAGrid + 1
                endif
             enddo
          enddo
       enddo
       
       if (allocated(ssM_msun   ))  deallocate(ssM_msun   )
       if (allocated(subsrcMass ))  deallocate(subsrcMass )
       if (allocated(subNormFlux))  deallocate(subNormFlux)
       if (allocated(subSrcSeries)) deallocate(subSrcSeries)
       if (allocated(sub_srcpos ))  deallocate(sub_srcpos )
       if (allocated(rsub_srcpos))  deallocate(rsub_srcpos)
       allocate(ssM_msun    (NumAGrid   ))
       allocate(subsrcMass  (NumAGrid   ))
       allocate(subNormFlux (NumAGrid   ))
       allocate(subSrcSeries( NumAGrid  ))
       allocate(sub_srcpos  (3, NumAGrid))
       allocate(rsub_srcpos (3, NumAGrid))
       
       ssM_msun         = 0d0
       NumAGrid         = 0
       tot_subsrcM_msun = 0d0
       
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                ! When nz=1, no need to differentiate mass in time!!!!
                if (nz > 1) then
                   diff_subsrcMsun = &
                        subsrcM_msun(zred_now, &
                        ndens(i,j,k)/avg_dens,densNDcrit) &
                        - subsrcM_msun(zred_prev, &
                        ndens_previous(i,j,k)/avg_dens_previous,    &
                        densNDcrit_prev)
                else
                   diff_subsrcMsun= &
                        subsrcM_msun(zred_now,ndens(i,j,k)/avg_dens,densNDcrit)
                endif
                ! subgrid is active when neutral & containing minihalos
                if (AGflag(i,j,k) == 1) then
                   NumAGrid                = NumAGrid + 1
                   sub_srcpos (1, NumAGrid) = i
                   sub_srcpos (2, NumAGrid) = j
                   sub_srcpos (3, NumAGrid) = k
                   ! Source is always at cell centre!!
                   rsub_srcpos(1, NumAGrid) = x(i)
                   rsub_srcpos(2, NumAGrid) = y(j)
                   rsub_srcpos(3, NumAGrid) = z(k)
                   if (jLWgrid(i,j,k) <= jLWcrit_min) then
                      ! dens_ND is the dimensionless density.
                      ! subgrid source mass in solar mass unit
!!$                      ! partially suppress by ionized fraction.
!!$                      ssM_msun(NumAgrid) = &
!!$                           subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit) &
!!$                           * xh(i,j,k,0)
                      ! NO PARTIAL SUPPRESSION. differential mass only.
                      ssM_msun(NumAgrid) = &
                           diff_subsrcMsun
                   else ! do interpolation to partially suppress sources
!!$                      ! suppression both by jLW and ionized fraction
!!$                      ssM_msun(NumAgrid) = &
!!$                           subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit) &
!!$                           * (jLWcrit_now - jLWgrid(i,j,k)) &
!!$                           / (jLWcrit_now - jLWcrit_min) &
!!$                           * xh(i,j,k,0)
                      ! suppression only by jLW. Use differential mass.
                      ssM_msun(NumAgrid) = &
                           diff_subsrcMsun &
                           * (jLWcrit_now - jLWgrid(i,j,k)) &
                           / (jLWcrit_now - jLWcrit_min)
                   endif
                endif
             enddo
          enddo
       enddo
       
       tot_subsrcM_msun = sum(ssM_msun)
       
       ! Convert subgrid source masses into grid mass unit, and then
       ! multiply phot_per_atom(3). Result is similar to srcMass under
       ! source_properties(but see below for difference when MHflag=2, 
       ! need to divide by fstar(3), since stellar mass is used then.). 
       ! Use it just internally.
       ! Also see below for the subNormFlux for the difference.
       if (MHflag == 1) then
          subsrcMass = ssM_msun * (M_SOLAR/M_grid) * phot_per_atom(3)
       elseif (MHflag == 2) then
          subsrcMass = ssM_msun * (M_SOLAR/M_grid) * phot_per_atom(3) &
               / fstar(3)
       endif
       
       ! Save new subgrid source list, just 
       open(unit=49,file=sourcelistfile_sub,status='unknown')
       write(49,*) NumAGrid
       ! Write MHflag to distinguish physical meaning of
       ! subsrcMass & ssM_msun. When MHflag=2, ssM_msun is the STELLAR
       ! BARYON MASS, not BARYON+DM HALO MASS.
       write(49,*) MHflag
       do ns = 1,NumAgrid
          write(49,333) sub_srcpos(1,ns),sub_srcpos(2,ns),sub_srcpos(3,ns),&
               subsrcMass(ns), ssM_msun(ns)
333       format(3I5,2e13.4)
       enddo
       close(49)
       

       write(logf,*) "Number of Active Grids: ", NumAGrid
       
       ! Turn masses into luminosities; photons/atom included in subsrcMass
       if (MHflag == 1) then
          subNormFlux = subsrcMass * M_grid * Omega_B/(Omega0*m_p) &
               /S_star_nominal /AGlifetime
       elseif (MHflag == 2) then
          subNormFlux = subsrcMass * M_grid          /        m_p  &
               /S_star_nominal /AGlifetime
       endif

    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumAGrid,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) then
       if (allocated(ssM_msun   ))  deallocate(ssM_msun   )
       if (allocated(subsrcMass ))  deallocate(subsrcMass )
       if (allocated(subNormFlux))  deallocate(subNormFlux)
       if (allocated(subSrcSeries)) deallocate(subSrcSeries)
       if (allocated(sub_srcpos ))  deallocate(sub_srcpos )
       if (allocated(rsub_srcpos))  deallocate(rsub_srcpos)
       allocate(ssM_msun    (NumAGrid   ))
       allocate(subsrcMass  (NumAGrid   ))
       allocate(subNormFlux (NumAGrid   ))
       allocate(subSrcSeries( NumAGrid  ))
       allocate(sub_srcpos  (3, NumAGrid))
       allocate(rsub_srcpos (3, NumAGrid))
    endif
    call MPI_BCAST(ssM_msun,     NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(subsrcMass,   NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(subNormFlux,  NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(sub_srcpos, 3*NumAGrid, MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(rsub_srcpos,3*NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of Active Grids: ", NumAGrid
#endif
          
    if (rank == 0) then
       
       write(logf,*) 'Subgrid Source lifetime=', AGlifetime/3.1536e13
       write(logf,*) 'Subgrid Total flux= ',sum(subNormFlux)

       write(6,*)  'almost ready with MH sources'
       ! Create array of subgrid source numbers for generating random order
       do ns=1,NumAGrid
          subSrcSeries(ns)=ns
       enddo

       ! Make a random order
       call permi(NumAGrid, subSrcSeries)

    endif

#ifdef MPI
    ! Distribute the subgrid source series to the other nodes
    call MPI_BCAST(subSrcSeries,NumAGrid,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    return
  end subroutine AGrid_properties

  ! =======================================================================

  subroutine get_denscrit(zred, densNDcrit)

    ! Find out the critical dimensionless density, over which cells
    ! can populate minihalos. This is set by matching
    ! the global number density of all minihalos to high-res, small-box
    ! simulation result which resolves minihalos.

    real(kind=dp), intent(in) :: zred

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

    ! what we want
    real(kind=dp), intent(out) :: densNDcrit


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

       volratio = (size_smallbox/h)**3 &
            / (vol_cMpc3*real(mesh(1)*mesh(2)*mesh(3),dp))

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
                      idx_LGdelta1 = int((log10(ndens(i,j,k)/avg_dens) - LGdelta1min)&
                           /dLGdelta1) + 1
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
       densNDcrit = 10d0**LGdensND_N
    endif   ! end of MHflag=2 case

  end subroutine get_denscrit

  ! =======================================================================

  function subsrcM_msun(zred, densND, densNDcrit)
    
    ! Get subgrid src mass (those not resolved in N-body) in solar mass unit
    ! Cell characteristic given by zred & dimensionless density

    integer                   :: idx_zred, idx_LGdelta1

    real(kind=dp)             :: subsrcM_msun
    real(kind=dp), intent(in) :: zred
    real(kind=si), intent(in) :: densND
    real(kind=dp), intent(in) :: densNDcrit
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

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then
       open(unit=20, file=trim(adjustl(MHfit_file)), form="binary", status="old")
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

  function jLWcrit(zred)

    ! Get critical jLW at zred

    real(kind=dp)             :: jLWcrit
    real(kind=dp), intent(in) :: zred

    jLWcrit = 0.1d0 * 1d-21

  end function jLWcrit


end module source_sub
