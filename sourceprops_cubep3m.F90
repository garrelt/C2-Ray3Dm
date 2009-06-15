!>
!! \brief This module contains data and routines for handling the sources
!! of reionization.
!!
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 30-Jan-2008
!!
!! \b Version: CUBEP3M Nbody simulation, with source suppression, and different source models

module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, logf, file_input
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: id_str, M_grid, dir_src, NumZred
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes

  integer :: NumSrc !< Number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos !< grid position of sources
  real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources
  real(kind=dp),dimension(:),allocatable :: uv_array  !< list of UV flux evolution (for some sources models)

  character(len=30) :: UV_Model !< type of UV model
  integer :: NumZred_uv !< Number of redshift points in UV model
  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01,total_SrcMass
  character(len=6) :: z_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources

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

    character(len=512) :: sourcelistfile,sourcelistfilesuppress
    integer :: ns,ns0

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(rsrcpos)) deallocate(rsrcpos)
    if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(srcSeries)) deallocate(srcSeries)
    
    ! Ask for input
    if (rank == 0) then
       
       ! Sources are read from file
       write(z_str,'(f6.3)') zred_now
       
       ! Construct the file name
       sourcelistfile=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))//"_sources.dat"
       sourcelistfilesuppress=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            "_sources_used_wfgamma.dat"
       
       if (restart == 0 .or. restart == 1) then
          open(unit=50,file=sourcelistfile,status='old')
          ! Number of sources
          read(50,*) NumSrc0
          
          ! Report
          write(logf,*) "Total number of source locations, no suppression: ", &
               NumSrc0
          
          ! Read in source positions and mass to establish number
          ! of non-suppressed sources
          NumSrc = 0
          NumMassiveSrc = 0
          NumSupprbleSrc = 0
          NumSupprsdSrc = 0
          do ns0=1,NumSrc0
             read(50,*) srcpos0(1),srcpos0(2),srcpos0(3),SrcMass00,SrcMass01
             ! the cell is still neutral, no suppression
             if (SrcMass00 /= 0.0 .or. &
                  xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) &
                  NumSrc=NumSrc+1
             ! Count different types of sources
             if (SrcMass00 /= 0.0) NumMassiveSrc=NumMassiveSrc+1
             if (SrcMass01 /= 0.0) NumSupprbleSrc=NumSupprbleSrc+1
             ! How many suppressed?
             if (SrcMass01 /= 0.0 .and. &
                  xh(srcpos0(1),srcpos0(2),srcpos0(3),1) > StillNeutral) &
                  NumSupprsdSrc=NumSupprsdSrc+1
          enddo
          close(50)
          write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
          write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
          write(logf,*) "Number of massive sources: ",NumMassiveSrc
          if (NumSupprbleSrc > 0) write(logf,*) "Suppressed fraction: ", &
               real(NumSupprsdSrc)/real(NumSupprbleSrc)
       else
          ! Upon restart from intermediate redshift use the previously 
          ! calculated suppressed source list
          open(unit=49,file=sourcelistfilesuppress,status='unknown')
          ! Number of sources
          read(49,*) NumSrc
          close(49)
       endif
       write(logf,*) "Number of sources, with suppression: ",NumSrc
    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif
    
    ! Allocate arrays for this NumSrc
    allocate(srcpos(3,NumSrc))
    allocate(rsrcpos(3,NumSrc))
    allocate(SrcMass(NumSrc,0:Number_Sourcetypes))
    allocate(NormFlux(NumSrc))
    allocate(SrcSeries(NumSrc))
    
    if (rank == 0) then
       if (restart == 0 .or. restart == 1) then
          open(unit=50,file=sourcelistfile,status='old')
          ! Number of sources
          read(50,*) NumSrc0
          ! Read in source positions and mass
          ns=0
          do ns0=1,NumSrc0
             read(50,*) srcpos0(1),srcpos0(2),srcpos0(3), &
                  SrcMass00,SrcMass01
             
             if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) then
                ! the cell is still neutral, no suppression
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                ! Source is always at cell centre!!
                rsrcpos(1,ns)=x(srcpos(1,ns))
                rsrcpos(2,ns)=y(srcpos(2,ns))
                rsrcpos(3,ns)=z(srcpos(3,ns))
                SrcMass(ns,1)=SrcMass00
                SrcMass(ns,2)=SrcMass01
             elseif (SrcMass00 > 0.0d0) then
                !the cell is ionized but source is massive enough to survive
                !and is assumed Pop. II 
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                ! Source is always at cell centre!!
                rsrcpos(1,ns)=x(srcpos(1,ns))
                rsrcpos(2,ns)=y(srcpos(2,ns))
                rsrcpos(3,ns)=z(srcpos(3,ns))
                SrcMass(ns,1)=SrcMass00
                SrcMass(ns,2)=0.0
                !else
                !   ! Report
                !   write(logf,*) 'Source dropped: ', &
                !        srcpos0(1),srcpos0(2),srcpos0(3)
                !   write(logf,*) 'mass= ',SrcMass01*(M_grid/M_SOLAR), &
                !        xh(srcpos0(1),srcpos0(2),srcpos0(3),0)
             endif
          enddo
          
          ! Collect total source mass 
          SrcMass(:,0)=SrcMass(:,1)*phot_per_atom(1)  & !massive sources
               +SrcMass(:,2)*phot_per_atom(2)      !small sources   
          ! Save new source list, without the suppressed ones
          open(unit=49,file=sourcelistfilesuppress,status='unknown')
          write(49,*) NumSrc
          do ns0=1,NumSrc
             write(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
                  SrcMass(ns0,0)
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
                  SrcMass(ns0,0)
             ! Source is always at cell centre!!
             rsrcpos(1,ns0)=x(srcpos(1,ns0))
             rsrcpos(2,ns0)=y(srcpos(2,ns0))
             rsrcpos(3,ns0)=z(srcpos(3,ns0))
          enddo
          close(49)
       endif ! of restart test
    endif ! of rank 0 test
    
    
#ifdef MPI
    ! Distribute the source parameters to the other nodes
    call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(SrcMass,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Turn masses into luminosities
    select case (UV_Model)
    case ("Iliev et al")
       do ns=1,NumSrc
          NormFlux(ns)=SrcMass(ns,0)*M_grid*  &!note that now photons/atom are included in SrcMass
               Omega_B/(Omega0*m_p)/S_star_nominal
          !NormFlux(ns)=NormFlux(ns)/lifetime
          NormFlux(ns)=NormFlux(ns)/lifetime2
       enddo
    case ("Fixed N_gamma")
       if (nz <= NumZred_uv) then
          total_SrcMass=sum(SrcMass(:,0))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             SrcMass(ns,0)=sum(SrcMass(ns,1:Number_Sourcetypes))
             NormFlux(ns)=uv_array(nz)/lifetime2*SrcMass(ns,0)/total_SrcMass/S_star_nominal
          enddo
          write(logf,*) uv_array(nz),SrcMass(:,0),uv_array(nz)/lifetime2
       else
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) "No UV model available, setting fluxes to zero."
       endif
    case ("Fixed Ndot_gamma")
       if (nz <= NumZred_uv) then
          total_SrcMass=sum(SrcMass(:,0))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             SrcMass(ns,0)=sum(SrcMass(ns,1:Number_Sourcetypes))
             NormFlux(ns)=uv_array(nz)*SrcMass(ns,0)/total_SrcMass/S_star_nominal
          enddo
       else
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) "No UV model available, setting fluxes to zero."
       endif
    end select

    if (rank == 0) then
       write(logf,*) 'Source lifetime=', lifetime2/(1e6*YEAR),' Myr'
       !write(logf,*) 'Source lifetime=', lifetime/3.1536e13
       write(logf,*) 'Total flux= ',sum(NormFlux)*S_star_nominal,' s^-1'
       ! Create array of source numbers for generating random order
       do ns=1,NumSrc
          SrcSeries(ns)=ns
       enddo
       
       ! Make a random order
       call ctrper(SrcSeries(1:NumSrc),1.0)
    endif

#ifdef MPI
    ! Distribute the source series to the other nodes
    call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine source_properties

  ! =======================================================================

  !> Initialization routine: determine the source model and optionally read 
  !! in source properties
  !! Author: Garrelt Mellema
  
  
  !! This accomodates different source models
  !! 0: Iliev et al source, Ndot_gamma= f*M_halo/timestep
  !! 1: Fixed total N_gamma, still need to divide by time step
  !! 2: Fixed total Ndot_gamma.
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
       if (.not.file_input) write(*,"(A,$)") "UV Luminosity recipe (0,1,2): "
       read(stdinput,*) uv_answer
       select case (uv_answer)
       case(0)
          UV_Model = "Iliev et al"
       case(1)
          UV_Model = "Fixed N_gamma"
       case(2)
          UV_Model = "Fixed Ndot_gamma"
       end select
       
       if (uv_answer > 0) then
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
    call MPI_BCAST(UV_Model,30,MPI_CHARACTER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(NumZred_uv,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(uv_array(NumZred_uv))
    call MPI_BCAST(uv_array,NumZred_uv,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif
    
  end subroutine source_properties_ini
  
end module sourceprops
