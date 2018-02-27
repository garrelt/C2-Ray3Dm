!>
!! \brief This module contains data and routines for handling the sources
!! of reionization.
!!
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 08-Jul-2015 (30-Jan-2008
!!
!! \b Version: CUBEP3M Nbody simulation, with source suppression, 
!! and different source models; recently extended with partial and
!! gradual suppression.

module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, logf, file_input
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: id_str, M_grid, dir_src, NumZred, n_box, zred_array
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, zeta, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes

  implicit none

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="_sources.dat"
       !sourcelistfile_base="_wsubgrid_sources.dat"
       !sourcelistfile_base="_wsubgrid_gradual_sources.dat"
  character(len=100),parameter,private :: &
       sourcelistfilesuppress_base="_sources_used_wfgamma.dat"

  !> number of columns in source list
  integer,parameter,private :: ncolumns_srcfile=5
  real,dimension(ncolumns_srcfile),private :: srclist

  !> Definitions to help accessing the read in source list array
  integer,parameter :: HMACH=4
  integer,parameter :: LMACH=5
  integer,parameter :: LMACH_SUPPR=6

  !> maximum increase in uv to use up cumulated photons
  real(kind=dp),parameter,private :: cumfrac_max=0.15 

  integer :: NumSrc=0 !< Number of sources
  integer :: Prev_NumSrc !< Previous number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  !real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
  !integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources
  real(kind=dp),dimension(:),allocatable :: uv_array  !< list of UV flux evolution (for some sources models)
  !> The cumulative number of uv photons. We save this number so we can add it
  !! to the uv luminosity the first time sources appear.
  real(kind=dp) :: cumulative_uv=0.0

  character(len=30) :: UV_Model !< type of UV model
  integer :: NumZred_uv !< Number of redshift points in UV model
  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: total_SrcMass
  character(len=6),private :: z_str !< string value of redshift
  character(len=6),private :: z_previous_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources
  real(kind=dp),private :: MassHMACH !< total mass of high mass sources
  real(kind=dp),private :: MassLMACH !< total mass of suppressible sources
  real(kind=dp),private :: MassLMACHsupprsd !< total mass of suppressed sources
  real(kind=dp),private :: dMassHMACH !< change in HMACH mass
  real(kind=dp),private :: dMassLMACHactive !< change in active LMACH mass
  real(kind=dp),private :: MassHMACH_previous !< previous mass in HMACH
  real(kind=dp),private :: MassLMACHactive_previous !< previous mass in active LMACH
  real(kind=dp),private :: cumfrac
  real(kind=dp),private :: xh_ion !< keeps local ionization fraction
  character(len=512),private :: sourcelistfile,sourcelistfilesuppress, &
       sourcelistfile_previous

contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  !! Authors: Garrelt Mellema, Ilian Iliev
  !! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))

  subroutine source_properties(zred_now,nz,lifetime2,restart)

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    integer :: ns,ns0

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    !if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
    !if (allocated(srcSeries)) deallocate(srcSeries)
    
    Prev_NumSrc=NumSrc

    ! Rank 0 reads in sources
    if (rank == 0) then
       
       ! Sources are read from files with redshift in the file name:
       ! construct redshift string
       write(z_str,'(f6.3)') zred_now
       if (nz > 1) write(z_previous_str,'(f6.3)') zred_array(nz-1)
       
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfile_base))
       sourcelistfile_previous=trim(adjustl(dir_src))//&
            trim(adjustl(z_previous_str))//"-"//trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfile_base))
       sourcelistfilesuppress=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfilesuppress_base))

       ! Find the total masses of sources for the previous time step.
       ! These will be used to calculate the change of mass and from
       ! this the luminosity.
       if (UV_model == "Collapsed fraction growth") then
          if (nz>1) then
             call establish_number_of_active_sources (-1,sourcelistfile_previous)
             MassHMACH_previous=MassHMACH
             MassLMACHactive_previous=MassLMACH-MassLMACHsupprsd
          else
             MassHMACH_previous=0.0
             MassLMACHactive_previous=0.0
          endif
       endif

       call establish_number_of_active_sources (restart,sourcelistfile)

       ! Calculate change in mass, this will be used to calculate the luminosity
       if (UV_model == "Collapsed fraction growth") then
          
          dMassHMACH=MassHMACH-MassHMACH_previous
          write(logf,*) "Change in HMACH mass: ",dMassHMACH
          if (dMassHMACH < 0.0) then
             dMassHMACH=0.0
             write(logf,*) "Change in HMACH mass set to zero"
          endif
          dMassLMACHactive=MassLMACH-MassLMACHsupprsd-MassLMACHactive_previous
          write(logf,*) "Change in active LMACH mass: ",dMassLMACHactive
          if (dMassLMACHactive < 0.0) then
             dMassLMACHactive=0.0
             write(logf,*) "Change in activeLMACH mass set to zero"
          endif
       endif
       
    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif
    
    ! Allocate arrays for this NumSrc
    if (NumSrc > 0) then
       allocate(srcpos(3,NumSrc))
       !allocate(SrcMass(NumSrc,0:Number_Sourcetypes))
       allocate(NormFlux(0:NumSrc)) ! 0 will hold lost photons
       !allocate(SrcSeries(NumSrc))

       ! Fill in the source arrays
       if (rank == 0) then
          call read_in_sources (restart)
    
          ! Set cumulative number of uv photons to zero if this is not the
          ! first redshift for which sources are active (this way cumulative_uv
          ! can be used in all cases).
          ! New version: cumulative_uv is slowly reduced
          !if (Prev_NumSrc /= 0) cumulative_uv=0.0
          call assign_uv_luminosities (lifetime2,nz)

          write(logf,*) 'Source lifetime=', lifetime2/(1e6*YEAR),' Myr'
          write(logf,*) 'Total flux= ',sum(NormFlux(1:NumSrc))*S_star_nominal,' s^-1'
          ! Create array of source numbers for generating random order
          !do ns=1,NumSrc
          !   SrcSeries(ns)=ns
          !enddo
          
          ! Make a random order
          !call ctrper(SrcSeries(1:NumSrc),1.0)
       endif

#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       !call MPI_BCAST(SrcMass,(1+Number_Sourcetypes)*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
    
#ifdef MPI
       ! Distribute the source series to the other nodes
       !call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
      
    else

       ! For the prescribed uv evolution model, accumulate the photons
       ! (to be used once the first sources appear)
       if (UV_Model == "Fixed N_gamma") &
            cumulative_uv=cumulative_uv + uv_array(nz)

    endif

  end subroutine source_properties

  ! =======================================================================

  subroutine establish_number_of_active_sources (restart,file_to_use)

    integer,intent(in) :: restart
  character(len=512),intent(in) :: file_to_use
    real :: odens	

    integer :: ns0
    integer :: ns
    real(kind=dp) :: xh_ion
    
    if (restart == 0 .or. restart == 1 .or. restart == -1) then

       ! Initialise counters to zero
       NumSrc = 0
       NumMassiveSrc = 0
       NumSupprbleSrc = 0
       NumSupprsdSrc = 0
       MassHMACH=0.0
       MassLMACH=0.0
       MassLMACHsupprsd=0.0
       open(unit=50,file=file_to_use,status='old')
       ! Number of sources
       read(50,*) NumSrc0
       
       ! Report
       write(logf,*) "Total number of source locations, no suppression: ", &
            NumSrc0
       
       ! Read in source positions and mass to count the number
       ! of non-suppressed sources
       do ns0=1,NumSrc0
          ! If you change the following lines, also change it below in
          ! read_in_sources
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))

          ! Retrieve ionized fraction in source cell
#ifdef ALLFRAC
          xh_ion=xh(srcpos0(1),srcpos0(2),srcpos0(3),1)
#else
          xh_ion=xh(srcpos0(1),srcpos0(2),srcpos0(3))
#endif

          ! Massive sources are never suppressed.
          if (srclist(HMACH) /= 0.0 .or. &
               UV_Model == "Iliev et al partial supp." .or. &
               (UV_model == "Gradual supp." .and. srclist(LMACH_SUPPR) /= 0.0))&
               then
             NumSrc=NumSrc+1
             ! if the cell is still neutral, no suppression (if we use the Iliev
             ! et al source model)   
          elseif (xh_ion < StillNeutral .and. &
               UV_Model == "Iliev et al") then
             NumSrc=NumSrc+1
          endif

          ! Count different types of sources
          if (srclist(HMACH) /= 0.0) then
             NumMassiveSrc=NumMassiveSrc+1
             MassHMACH=MassHMACH+srclist(HMACH)
          endif
          if (srclist(LMACH) /= 0.0) then
             NumSupprbleSrc=NumSupprbleSrc+1
             MassLMACH=MassLMACH+srclist(LMACH)
          endif
          ! How many suppressed?
          if (srclist(LMACH) /= 0.0) then
             if (xh_ion > StillNeutral .or. &
                  (UV_Model /= "Iliev et al" .and. &
                  UV_Model /= "Iliev et al partial supp." .and. &
                  UV_model /= "Gradual supp." .and. &
                  UV_Model /= "Collapsed fraction growth")) then
                NumSupprsdSrc=NumSupprsdSrc+1
                MassLMACHsupprsd=MassLMACHsupprsd+srclist(LMACH)
             endif
          endif
       enddo
       close(50)
       write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
       write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
       write(logf,*) "Number of massive sources: ",NumMassiveSrc
       if (NumSupprbleSrc > 0) write(logf,*) "Suppressed fraction: ", &
            real(NumSupprsdSrc)/real(NumSupprbleSrc)
       ! Collapsed fraction calculation
       ! f_coll = M_halo / M_vol
       ! However, M_halo is given in units of M_grid which is M_vol/n_box^3
       ! so f_coll = M_halo / n_box^3
       write(logf,*) "Collapsed fraction in massive sources: ", &
            MassHMACH/(real(n_box)**3)
       write(logf,*) "Collapsed fraction in suppressable sources: ", &
            MassLMACH/(real(n_box)**3)
       write(logf,*) "Collapsed fraction in suppressed sources: ", &
            MassLMACHsupprsd/(real(n_box)**3)
       if (NumSupprbleSrc > 0) write(logf,*) "Suppressed mass fraction: ", &
            MassLMACHsupprsd/MassLMACH
       
    else
       
       ! Upon restart from intermediate redshift use the previously 
       ! calculated suppressed source list
       open(unit=49,file=sourcelistfilesuppress,status='unknown')
       ! Number of sources
       read(49,*) NumSrc
       close(49)
    endif
    write(logf,*) "Number of sources, with suppression: ",NumSrc

  end subroutine establish_number_of_active_sources

  ! =======================================================================

  subroutine read_in_sources (restart)

    integer,intent(in) :: restart

    integer :: ns
    integer :: ns0

    if (restart == 0 .or. restart == 1) then
       open(unit=50,file=sourcelistfile,status='old')
       ! Number of sources
       read(50,*) NumSrc0
       ! Read in source positions and mass
       ns=0
       do ns0=1,NumSrc0
          ! If you change the following lines, also change it above in
          ! establish_number_of_active_sources
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))
          ! Retrieve ionized fraction in source cell
#ifdef ALLFRAC
          xh_ion=xh(srcpos0(1),srcpos0(2),srcpos0(3),1)
#else
          xh_ion=xh(srcpos0(1),srcpos0(2),srcpos0(3))
#endif

          if (xh_ion < StillNeutral) then
             if (UV_Model == "Iliev et al" .or. &
                  UV_Model == "Iliev et al partial supp." .or. &
                  UV_model == "Gradual supp." .or. &
		  srclist(HMACH) > 0.0d0) then
                ! the cell is still neutral, no suppression
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                ! Collect total source mass (weigthed with efficiency factor
                ! in case of the Iliev et al source model).
                if (UV_Model == "Iliev et al" .or. &
                     UV_Model == "Iliev et al partial supp." .or. &
		     UV_model == "Gradual supp.") then
                   NormFlux(ns)=srclist(HMACH)*phot_per_atom(1) & !massive sources
                        + srclist(LMACH)*phot_per_atom(2)      !small sources  
                elseif (UV_Model == "Collapsed fraction growth") then
                   NormFlux(ns)=srclist(HMACH)*zeta(1)* &
                        dMassHMACH/MassHMACH  & ! high mass sources
                        ! low mass sources  
                        + srclist(LMACH)*zeta(2)* &
                        dMassLMACHactive/(MassLMACH-MassLMACHsupprsd) 
                else
                   NormFlux(ns)=srclist(HMACH)!+srclist(LMACH)
                endif
             endif
          elseif (srclist(HMACH) > 0.0d0 .or. &
               (UV_Model == "Iliev et al partial supp." &
               .and. srclist(LMACH) > 0.0d0) &
               .or. (UV_model == "Gradual supp." &
               .and. srclist(LMACH_SUPPR) > 0.0d0)) then
             !the cell is ionized but source is massive enough to survive
             !and is assumed Pop. II, or source is low-mass, but we assume
             !partial suppression, tuning down its efficiency  
             ns=ns+1
             ! Source positions in file start at 1!
             srcpos(1,ns)=srcpos0(1)
             srcpos(2,ns)=srcpos0(2)
             srcpos(3,ns)=srcpos0(3)
             ! Collect total source mass (weigthed with efficiency factor
             ! in case of the Iliev et al source model), used in 
             ! assign_uv_luminosities to calculate ionizing photon rates
             if (UV_Model == "Iliev et al") then
                NormFlux(ns)=srclist(HMACH)*phot_per_atom(1)  !massive sources
                
             elseif (UV_Model == "Iliev et al partial supp.") then
                !make low-mass sources as efficient as HMACHs
                ! (so typically less efficient)
                NormFlux(ns)=srclist(HMACH)*phot_per_atom(1)  & !massive sources
                     + srclist(LMACH)*phot_per_atom(1)      !low-mass sources
                
	     elseif (UV_Model == "Gradual supp.") then
                !make low-mass sources less efficient
                NormFlux(ns)=srclist(HMACH)*phot_per_atom(1)  & 
                     !low-mass sources with effective mass KLD
                     + srclist(LMACH_SUPPR)*phot_per_atom(2)      
                
             elseif (UV_Model == "Collapsed fraction growth") then
                NormFlux(ns)=srclist(HMACH)*zeta(1)*dMassHMACH/MassHMACH 
             else
                NormFlux(ns)=srclist(HMACH)
             endif
          endif
       enddo
       close(50)
       
       ! Save new source list, without the suppressed ones
       open(unit=49,file=sourcelistfilesuppress,status='unknown')
       write(49,*) NumSrc
       do ns0=1,NumSrc
          write(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
               NormFlux(ns0)
       enddo
       close(49)
    else ! of restart test
       ! Read source list from file saved previously
       open(unit=49,file=sourcelistfilesuppress,status="old")
       write(logf,*) "Reading ",NumSrc," sources from ", &
            trim(adjustl(sourcelistfilesuppress))
       read(49,*) NumSrc
       do ns0=1,NumSrc
          read(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
                NormFlux(ns0)
       enddo
       close(49)
    endif ! of restart test
    
  end subroutine read_in_sources

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)
    
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    integer :: ns

    ! Turn masses into luminosities
    select case (UV_Model)
    case ("Iliev et al", "Iliev et al partial supp.","Gradual supp.")
       do ns=1,NumSrc
          NormFlux(ns)=NormFlux(ns)*M_grid*  &
               Omega_B/(Omega0*m_p)/S_star_nominal
          !note that now photons/atom are already included in NormFlux
          !NormFlux(ns)=NormFlux(ns)/lifetime
          NormFlux(ns)=NormFlux(ns)/lifetime2
       enddo
    case ("Collapsed fraction growth")
       do ns=1,NumSrc
          NormFlux(ns)=NormFlux(ns)*M_grid*  &
               Omega_B/(Omega0*m_p)/S_star_nominal
          NormFlux(ns)=NormFlux(ns)/lifetime2
       enddo
    case ("Fixed N_gamma")
       if (nz <= NumZred_uv) then
          cumfrac=min(cumfrac_max,cumulative_uv/uv_array(nz))
          if (rank == 0) then 
             write(logf,*) 'Cumulative versus current photons: ', &
                  cumulative_uv,uv_array(nz),cumulative_uv/uv_array(nz)
             write(logf,*) 'Cumulative fraction used: ', cumfrac
          endif
          total_SrcMass=sum(NormFlux(1:NumSrc))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             NormFlux(ns)=(1.0+cumfrac)*uv_array(nz)/lifetime2*NormFlux(ns)/total_SrcMass/S_star_nominal
             !NormFlux(ns)=(cumulative_uv+uv_array(nz))/lifetime2*SrcMass(ns,0)/total_SrcMass/S_star_nominal
          enddo
          ! Subtract extra photons from cumulated photons
          cumulative_uv=max(0.0_dp,cumulative_uv-cumfrac*uv_array(nz))
          !write(logf,*) uv_array(nz),SrcMass(:,0),uv_array(nz)/lifetime2
       else
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) "No UV model available, setting fluxes to zero."
       endif
    case ("Fixed Ndot_gamma")
       if (nz <= NumZred_uv) then
          total_SrcMass=sum(NormFlux(1:NumSrc))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             NormFlux(ns)=uv_array(nz)*NormFlux(ns)/total_SrcMass/S_star_nominal
          enddo
       else
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) "No UV model available, setting fluxes to zero."
       endif
    end select
    
  end subroutine assign_uv_luminosities
  
  ! =======================================================================

  !> Initialization routine: determine the source model and optionally read 
  !! in source properties
  !! Author: Garrelt Mellema
  
  
  !! This accomodates different source models
  !! 0: Iliev et al source, Ndot_gamma= f*M_halo/timestep
  !! 1: Fixed total N_gamma, still need to divide by time step
  !! 2: Fixed total Ndot_gamma.
  !! 3: Iliev et al source as above, but partial suppression of LMACHs
  !!    by tuning them down to lower efficiency 
  !! 4: Gradual suppression based on mass
  !! 5: Ndot_gamma proportional to mass but scaled with the total mass
  !!    growth of halos during the time step.
  subroutine source_properties_ini ()
    

    integer :: uv_answer
    real(kind=dp) :: z_in, N_source_nosupp, N_source_supp, N_gamma_nosupp
    character(len=180) :: uv_file ! name of file with uv model for redshifts
    integer :: nz

#ifdef MPI
    integer :: mympierror
#endif

    ! Ask for redshift file
    if (rank == 0) then
       if (.not.file_input) write(*,"(A,$)") "UV Luminosity recipe (0,1,2,3,4): "
       read(stdinput,*) uv_answer
       select case (uv_answer)
       case(0)
          UV_Model = "Iliev et al"
       case(1)
          UV_Model = "Fixed N_gamma"
       case(2)
          UV_Model = "Fixed Ndot_gamma"
       case(3)
          UV_Model = "Iliev et al partial supp."
       case(4)
	  UV_Model = "Gradual supp."
       case(5)
	  UV_Model = "Collapsed fraction growth"
       end select

       if (uv_answer == 1 .or. uv_answer == 2) then
          if (.not.file_input) write(*,"(A,$)") "File with UV data: "
          read(stdinput,*) uv_file
          
          ! Open and read redshift file
          open(unit=60,file=uv_file,form="formatted",status="old")
          read(unit=60,fmt=*) NumZred_uv
          if (NumZred_uv /= NumZred) then
             write(logf,*) "WARNING: Number of redshifts in UV luminosity file (", &
                  NumZred_uv,") does not match number of redshifts in ", &
                  "redshift file (",NumZred,")."
          endif
          allocate(uv_array(NumZred_uv))
          if (uv_answer == 1) then
             do nz=1,NumZred_uv
                read(unit=60,fmt=*) z_in, N_source_nosupp, N_source_supp, & 
                     N_gamma_nosupp, uv_array(nz)
             enddo
          else
             do nz=1,NumZred_uv
                read(unit=60,fmt=*) z_in, uv_array(nz)
             enddo
          endif
          close(60)
       endif
    endif
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(uv_answer,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(UV_Model,30,MPI_CHARACTER,0,MPI_COMM_NEW,mympierror)
    if (uv_answer  == 1 .or. uv_answer == 2) then
       call MPI_BCAST(NumZred_uv,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       if (rank /= 0) allocate(uv_array(NumZred_uv))
       call MPI_BCAST(uv_array,NumZred_uv,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
            mympierror)
    endif
#endif
    
  end subroutine source_properties_ini
  
end module sourceprops
