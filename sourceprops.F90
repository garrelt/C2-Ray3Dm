!>
!! \brief This module contains data and routines for handling the sources
!! of reionization.
!!
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: Aug-2014 (30-Jan-2008)
!!
!! \b Version: CUBEP3M Nbody simulation, LG Nbody simulation, Test simulation; with source suppression, and different source models
!!
!! \b Notes: The cubep3m and LG cases are almost identical.
!! The two differences are
!! 1) For LG the file names are labelled with an integer number and for cubep3m
!!  with redshift
!! 2) For LG the mass units are solar mass (MSOLAR) and for cubep3m grid masses

module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, logf, file_input
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: nbody_type, id_str, M_grid, dir_src, NumZred
  use ionfractions_module, only: xh
  !use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes

  implicit none

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="_sources.dat"
       !sourcelistfile_base="_wsubgrid_sources.dat"
  character(len=100),parameter,private :: &
       sourcelistfilesuppress_base="_sources_used_wfgamma.dat"

  !> number of columns in source list
  integer,parameter,private :: ncolumns_srcfile=5
  real,dimension(ncolumns_srcfile),private :: srclist

  !> maximum increase in uv to use up cumulated photons
  real(kind=dp),parameter,private :: cumfrac_max=0.15 

  integer :: NumSrc=0 !< Number of sources
  integer :: Prev_NumSrc !< Previous number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  !real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
  real(kind=dp),dimension(:),allocatable :: NormFluxPL !< normalized ionizing flux of PL sources
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources
  real(kind=dp),dimension(:),allocatable :: uv_array  !< list of UV flux evolution (for some sources models)
  !> The cumulative number of uv photons. We save this number so we can add it
  !! to the uv luminosity the first time sources appear.
  real(kind=dp) :: cumulative_uv=0.0

  character(len=30) :: UV_Model !< type of UV model
  integer :: NumZred_uv !< Number of redshift points in UV model
  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01,total_SrcMass
  character(len=6),private :: z_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources
  real(kind=dp),private :: cumfrac

  character(len=512),private :: sourcelistfile,sourcelistfilesuppress

contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  !! Authors: Garrelt Mellema, Ilian Iliev
  !! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))

  subroutine source_properties(zred_now,nz,lifetime2,restart)

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
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(NormFluxPL)) deallocate(NormFluxPL)
    
    ! Keep previous number of sources
    Prev_NumSrc=NumSrc

    ! Rank 0 counts the sources in the files
    if (rank == 0) then

       ! Construct file names
       sourcelistfile=construct_sourcefilename(zred_now, nz, &
            sourcelistfile_base)
       sourcelistfilesuppress=construct_sourcefilename(zred_now, nz, &
            sourcelistfilesuppress_base)

       ! Count the active sources (using the source model). This
       ! will set NumSrc to the number of active sources, to be
       ! used in defining the source arrays.
       call count_or_read_in_sources(zred_now, nz, restart, "count")

    endif

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    ! Report
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif
    
    ! Allocate arrays for the NumSrc value found above
    if (NumSrc > 0) then

       allocate(srcpos(3,NumSrc))
       allocate(NormFlux(0:NumSrc)) ! position 0 will hold lost photons

       ! Rank 0 reads in the sources
       if (rank == 0) then
          ! Read in the sources and apply source model
          call count_or_read_in_sources (zred_now, nz, restart, "read")
    
          ! Convert NormFlux from mass units to uv photon production rates
          call assign_uv_luminosities (lifetime2,nz)

          ! Report
          write(logf,*) 'Source lifetime =', lifetime2/(1e6*YEAR),' Myr'
          write(logf,*) 'Total rate of ionizing photons = ', &
               sum(NormFlux(1:NumSrc))*S_star_nominal,' s^-1'

          ! Make a randomized list of source labels.
          !call make_random_srclist ()

       endif

#ifdef MPI
       ! Distribute the source properties to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1,MPI_DOUBLE_PRECISION,0, &
            MPI_COMM_NEW,mympierror)
       ! Distribute the source series to the other nodes
       !call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
    else

       ! If there are no sources and we are using 
       ! the prescribed uv evolution model, accumulate the photons
       ! (to be used once the first sources appear)
       if (UV_Model == "Fixed N_gamma") &
            cumulative_uv=cumulative_uv + uv_array(nz)

    endif

  end subroutine source_properties

  ! ============================================================================

  function construct_sourcefilename(redshift,number_of_redshift,basename)
    
    character(len=512) :: construct_sourcefilename
    real(kind=dp),intent(in) :: redshift
    integer,intent(in) :: number_of_redshift
    character(len=100),intent(in) :: basename

    character(len=6) :: redshift_str
    character(len=3) :: number_of_redshift_str
    character(len=512) :: source_file

    select case (nbody_type)
    case("cubep3m")
       ! Sources are read from files with redshift in the file name:
       ! construct redshift string
       write(redshift_str,'(f6.3)') redshift
       
       ! Construct the file names
       source_file=trim(adjustl(dir_src))//&
            trim(adjustl(redshift_str))//"-"// &
            trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfile_base))

    case("LG")
       ! Sources are read from files with redshift counter nz in the file name:
       ! construct corresponding string
       write(number_of_redshift_str,'(i3.3)') number_of_redshift
       
       ! Construct the file names
       source_file=trim(adjustl(dir_src))//&
            trim(adjustl(number_of_redshift_str))//"-"// &
            trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfile_base))

    case("test")
       source_file=trim(adjustl(dir_src))//"test_sources.dat"

    end select

    ! Report
    write(unit=logf,fmt="(4A)") "Reading ",id_str, &
         " source list from ",trim(source_file)

    ! Copy result to function variable
    construct_sourcefilename=source_file
    
  end function construct_sourcefilename

  ! =======================================================================

  subroutine count_or_read_in_sources (redshift, number_of_redshift, &
       restart, label)

    real(kind=dp),intent(in) :: redshift
    integer,intent(in) :: number_of_redshift
    integer,intent(in) :: restart
    character(len=5),intent(in) :: label ! "count" or "read"

    integer :: ns
    integer :: ns0
    logical :: suppressed

    real(kind=8) :: summed_weighted_mass
    real(kind=8) :: xHII

    ! Initialize the source masses to be read

    if (restart == 0 .or. restart == 1) then
       ! Read in original source list and apply source model
       
       ! Open original source list file
       open(unit=50,file=sourcelistfile,status='old')
       ! Read total number of sources in file
       read(50,*) NumSrc0

       ! Initialize counters and masses to zero
       select case(label)
       case("count")
          ! Report
          write(logf,*) & 
               "Total number of source locations, no suppression: ", &
               NumSrc0
          NumSrc = 0
          NumMassiveSrc = 0
          NumSupprbleSrc = 0
          NumSupprsdSrc = 0
       case("read")
          ns=0
          SrcMass00=0.0
          SrcMass01=0.0
       end select

       ! For every line in the source file
       do ns0=1,NumSrc0
          ! Read source file. The number of columns differs for
          ! different cases
          read(50,*) srclist(1:ncolumns_srcfile)
          ! Extract integer position from scrlist structure
          srcpos0(1:3)=int(srclist(1:3))
          ! First type of sources (if two types these are the 
          ! high mass sources or HMACHs)
          srcMass00=srclist(4)
          ! Second type of sources (if two types these are the 
          ! low mass sources or LMACHs)
          if (ncolumns_srcfile > 4) srcMass01=srclist(5)
          
          ! Count different types of sources
          if (label == "count") then
             if (SrcMass00 > 0.0) NumMassiveSrc=NumMassiveSrc+1
             if (SrcMass01 > 0.0) NumSupprbleSrc=NumSupprbleSrc+1
          endif

          ! Now count or store sources, using the source recipe.
#ifdef ALLFRAC
          ! Extract the ionization fraction in cell (used in some recipes)
          xhii=xh(srcpos0(1),srcpos0(2),srcpos0(3),1)
#else
          xHII=xh(srcpos0(1),srcpos0(2),srcpos0(3))
#endif
          ! Find the total mass of active sources in this cell, multiplied
          ! with their efficiency factor in some cases
          summed_weighted_mass=mass_from_source_models(srcMass00,SrcMass01, &
               xHII,suppressed)

          ! Count the number of suppressed sources
          if (label == "count" .and. SrcMass01 > 0.0 .and. suppressed) &
               NumSupprbleSrc=NumSupprbleSrc+1

          ! Count or store the active sources
          if (summed_weighted_mass > 0.0) then
             select case(label)
             case("count") 
                NumSrc=NumSrc+1
             case("read")
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                NormFlux(ns)=summed_weighted_mass
             end select
          endif
       enddo
       
       select case(label)
       case("count")
          ! Report source counts & statistics on suppression
          write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
          write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
          write(logf,*) "Number of massive sources: ",NumMassiveSrc
          if (NumSupprbleSrc > 0) write(logf,*) "Suppressed fraction: ", &
               real(NumSupprsdSrc)/real(NumSupprbleSrc)
       case("read")
          ! Record source list in sourcelistfilesuppress: first
          ! construct file name
          sourcelistfilesuppress= &
               construct_sourcefilename(redshift,number_of_redshift, &
               sourcelistfilesuppress_base)
          ! Save new source list, without the suppressed ones
          open(unit=49,file=sourcelistfilesuppress,status='unknown')
          write(49,*) NumSrc
          do ns0=1,NumSrc
             write(49,"(3i4,f10.3)") srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
                  NormFlux(ns0)
          enddo
          close(49)
       end select

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
    
  end subroutine count_or_read_in_sources

  ! =======================================================================

  function mass_from_source_models (mass_hmach, mass_lmach, xHII, suppress)

    ! This function takes input read from source files and applies the
    ! source recipe (UV_Model). The result is the value of the source. For most
    ! source recipes this is the mass, possibly multiplied with an
    ! efficiency factor, but it could also be photon production rate.

    real(kind=dp) :: mass_from_source_models
    real(kind=dp),intent(in) :: mass_hmach
    real(kind=dp),intent(in) :: mass_lmach
    real(kind=dp),intent(in) :: xHII
    logical,intent(out) :: suppress

    real(kind=8) :: f_hmach, f_lmach

    select case (UV_Model)
       case("Iliev et al")
          f_hmach=phot_per_atom(1)
          if (xHII > StillNeutral) then
             suppress=.true.
             f_lmach=0.0
          else
             suppress=.false.
             f_lmach=phot_per_atom(2)
          endif
          mass_from_source_models=mass_hmach*f_hmach + mass_lmach*f_lmach
       case("Iliev et al partial supp.")
          f_hmach=phot_per_atom(1)
          if (xHII > StillNeutral) then
             f_lmach=phot_per_atom(1)
             suppress=.true.
          else
             f_lmach=phot_per_atom(2)
             suppress=.false.
          endif
          mass_from_source_models=mass_hmach*f_hmach + mass_lmach*f_lmach
       case default
          mass_from_source_models=mass_hmach
       end select

     end function mass_from_source_models

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)
    
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    integer :: ns

    ! Turn masses into luminosities
    select case (UV_Model)
    case ("Iliev et al", "Iliev et al partial supp.")
       do ns=1,NumSrc
          !note that now photons/atom are already included in NormFlux
          NormFlux(ns)=Luminosity_from_mass(NormFlux(ns),lifetime2)
                       !NormFlux(ns)*M_grid*  &
               !Omega_B/(Omega0*m_p)/S_star_nominal
          !NormFlux(ns)=NormFlux(ns)/lifetime
          !NormFlux(ns)=NormFlux(ns)/lifetime2
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
             NormFlux(ns)=(1.0+cumfrac)*uv_array(nz)/lifetime2* &
                  NormFlux(ns)/(total_SrcMass*S_star_nominal)
          enddo
          ! Subtract extra photons from cumulated photons
          cumulative_uv=max(0.0_dp,cumulative_uv-cumfrac*uv_array(nz))
          !write(logf,*) uv_array(nz),SrcMass(:,0),uv_array(nz)/lifetime2
       else
          ! For high redshifts there may not be a uv model.
          ! Set fluxes to zero.
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) &
               "No UV model available, setting fluxes to zero."
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
          if (rank == 0) write(logf,*) &
               "No UV model available, setting fluxes to zero."
       endif
       case("Test")
          ! Do nothing, the photon production rate is already set from
          ! the source file
          !NormFlux(ns)=NormFlux(ns)
    end select
    
  end subroutine assign_uv_luminosities
  
  ! =======================================================================
  
  function Luminosity_from_mass (Mass, timeperiod)

    ! The normalized flux (luminosity) for normal sources is expressed
    ! in terms of a standard ionizing photon rate, called S_star_nominal.
    ! In radiation the tables have been calculated for a spectrum
    ! with an ionizing photon rate of S_star_nominal.

    ! Mass is supposed to be total mass of the source MULTIPLIED with
    ! the efficiency factor (f) which is the product of the star formation
    ! fraction, the escape fraction and the number of photons produced
    ! per baryon. Because of the latter factor we need to convert the
    ! mass to number of baryons by dividing my the proton mass m_p.
    ! NOTE: number of baryons is the total number of nucleons, which
    ! is why we divide my m_p and NOT by mu*m_p (where mu is mean mass
    ! of atoms/ions).

    real(kind=dp),intent(in) :: Mass !< mass in units of grid masses
    real(kind=dp),intent(in) :: timeperiod !< timeperiod in seconds
    real(kind=dp) :: Luminosity_from_mass !< photon rate in S_star_nominal

    Luminosity_from_mass = Mass*M_grid*Omega_B/(Omega0*m_p)/ &
         (timeperiod*S_star_nominal)

  end function Luminosity_from_mass

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
  !! 4: Test sources (for running the test problems)
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
       if (.not.file_input) write(*,"(A,$)") "UV Luminosity recipe (0 - 4): "
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
          UV_Model = "Test"
       end select

       if (uv_answer == 1 .or. uv_answer == 2) then
          if (.not.file_input) write(*,"(A,$)") "File with UV data: "
          read(stdinput,*) uv_file
          
          ! Open and read redshift file
          open(unit=60,file=uv_file,form="formatted",status="old")
          read(unit=60,fmt=*) NumZred_uv
          if (NumZred_uv /= NumZred) then
             write(logf,*) & 
                  "WARNING: Number of redshifts in UV luminosity file (", &
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

  ! =======================================================================

  subroutine make_random_srclist

    ! Create a random list of source lables in array SrcSeries.
    ! To use a random list of sources, step through the source labels
    ! as given in this array

    ! Module containing the ctrper function
    use m_ctrper

    integer :: ns

    ! Create array of source numbers for generating random order
    if (allocated(srcSeries)) deallocate(srcSeries)
    allocate(SrcSeries(NumSrc))
    
    ! Fill array with sorted values
    do ns=1,NumSrc
       SrcSeries(ns)=ns
    enddo
          
    ! Make a random order
    call ctrper(SrcSeries(1:NumSrc),1.0)
    
  end subroutine make_random_srclist
  
end module sourceprops
