!>
!! \brief This module data and routines which deal with radiative
!!     effects. 
!!
!!  Its main part deal with photo-ionizing radiation, but it
!!     also initializes other radiative properties, such as cooling (which
!!     are contained in different modules).
!!     It can be used in hydrodynamic or stand-alone radiative transfer 
!!     calculations.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 1D version similar to the 3D version.

module radiation
  
  !     This file contains routines having to do with the initialization
  !     of the radiative transport and atomic calculations for Yguazu-a.
  !     It can used for Yguazu-a or non-hydro photo-ionization calculations.
  !
  !     - rad_ini : master routine
  !     - spectrum_parms : Input routine: establish the ionizing spectrum 
  !     - spec_diag : Calculates properties of spectrum
  !     - spec_integr_cores: Calculates spectral integration cores
  !     - spec_integr : Calculates photo ionization integrals
  !     - rad_boundary : Set the radiative boundary condition
  !
  !     Needs following `modules':
  !     romberg : romberg integrators
  !
  !     Author: Garrelt Mellema
  ! 
  !     Date: 31-Jan-2008 (02-Jun-2004 (04-Mar-2004)
  
  ! Version
  ! Simplified version
  ! - Only hydrogen
  ! - Option for Grey photo-ionization cross section
  ! - MPI enabled (broadcasts of radiative parameters to all nodes).

  ! Notes:
  ! - the initialization of the radiative cooling does not really belong
  !   here.
  ! - isothermal is sometimes an input parameter, and sometimes a compile
  !   time parameter. This needs to be streamlined. Probably along similar
  !   lines as the stellar parameters are dealt with.
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigmasb, hplanck, kb, tpic2
  use cgsphotoconstants, only: frth0, frtop1, frtop2, sh0, betah0, sigh
  use astroconstants, only: R_SOLAR, L_SOLAR
  use romberg, only: scalar_romberg,vector_romberg,romberg_initialisation
  use c2ray_parameters, only: teff_nominal, S_star_nominal, isothermal

  implicit none

  !-----------------------------------------------------------------------
  !     NumFreq - Number of integration points in one of the three 
  !               frequency interval.
  !     NumTau -  Number of table points for the optical depth.
  !     NumFreqBnd - Number of frequency bands (1 for hydrogen only)
  !-----------------------------------------------------------------------

  !> Number of integration points in one of the three frequency interval.
  integer,parameter :: NumFreq=128
  !> Number of table points for the optical depth.
  integer,parameter :: NumTau=2000
  !> Number of frequency bands (1 for hydrogen only)
  integer,parameter :: NumFreqBnd=1
  
  !> This parameter sets the optical depth at the entrance of the grid.
  !> It can be used if radiation enters the simulation volume from the
  !> outside.
  real(kind=dp) :: tauHI=0.0

  ! Parameters defining the optical depth entries in the table.
  ! minlogtau is log10(lowest optical depth) (table position 1)
  ! maxlogtau is log10(highest optical depth) (table position NumTau)
  ! dlogtau is the step size in log10(tau) between table entries
  real(kind=dp),parameter :: minlogtau=-20.0 !< log10(lowest optical depth) 
  real(kind=dp),parameter :: maxlogtau=4.0 !< log10(highest optical depth)
  !> step size in log10(tau) between table entries
  real(kind=dp),parameter :: dlogtau=(maxlogtau-minlogtau)/real(NumTau)

  !> Logical that determines the use of grey opacities
  logical,parameter :: grey=.false. ! use grey opacities?

  ! stellar properties
  real(kind=dp) :: teff !< Black body effective temperature
  real(kind=dp) :: rstar !< Black body radius
  real(kind=dp) :: lstar !< Black body luminosity
  real(kind=dp) :: S_star !< Black body ionizing photons rate
  
  ! Photo-ionization integral cores
  real(kind=dp),dimension(NumFreqBnd) :: steph0 !< frequency steps in table 
  !> photo-ionization integral core for H0 (optically thick case)
  real(kind=dp),dimension(:,:,:),allocatable  :: h0int 
  !> photo-ionization heating integral core for H0 (optically thick case)
  real(kind=dp),dimension(:,:,:),allocatable  :: hh0int
  !> photo-ionization integral core for H0 (optically thin case)
  real(kind=dp),dimension(:,:,:),allocatable  :: h0int1
  !> photo-ionization heating integral core for H0 (optically thin case)
  real(kind=dp),dimension(:,:,:),allocatable  :: hh0int1

  ! Photo-ionization integrals (rates)
  !> photo-ionization integral for H0 (optically thick case)
  real(kind=dp),dimension(:,:),allocatable  :: hphot
  !> photo-ionization heating integral for H0 (optically thick case)
  real(kind=dp),dimension(:,:),allocatable  :: hheat
  !> photo-ionization integral for H0 (optically thin case)
  real(kind=dp),dimension(:,:),allocatable  :: hphot1
  !> photo-ionization heating integral for H0 (optically thin case)
  real(kind=dp),dimension(:,:),allocatable  :: hheat1

  !> This type contains all the photo-ionization rates
  !> The in and out rates are used to ensure photon-conservation.
  !> See the C2-Ray paper.
  type photrates
     real(kind=dp) :: h        !< total H ionizing rate
     real(kind=dp) :: hv_h     !< total H heating rate
     real(kind=dp) :: h_in     !< in-rate
     real(kind=dp) :: hv_h_in  !< in-heating rate
     real(kind=dp) :: h_out    !< out-rate
     real(kind=dp) :: hv_h_out !< out-heating rate
  end type photrates

  ! photo-ionization rates (disabled as they are passed as arguments)
  !real(kind=dp),public :: phih,hvphih
  !real(kind=dp),public :: phih_in,phih_out
  !real(kind=dp),public :: hvphih_in,hvphih_out

#ifdef MH
  real(kind=dp),    dimension(mesh(1),mesh(2),mesh(3)), public     :: jLW 
#endif

#ifdef MPI       
    integer,private :: ierror
#endif

contains

!=======================================================================

  !> initializes constants and tables for radiation processes (heating, cooling and ionization)
  subroutine rad_ini ()

    ! initializes constants and tables for radiation processes
    ! (heating, cooling and ionization)

    use radiative_cooling, only: setup_cool

    ! Initialize integration routines
    call romberg_initialisation(NumFreq)

    ! Ask for the parameters of the spectrum
    call spectrum_parms ()

    ! Determine spectrum diagnostics
    call spec_diag ()

    ! Calculate spectral integral cores
    call spec_integr_cores ()

    ! Find the photo-ionization integrals for this spectrum
    call spec_integr ()

    ! Set the radiative boundary conditions
    !call rad_boundary() ! NO LONGER NEEDED

    ! Set source position
    ! call source_position() CALLED ELSEWHERE

    ! Setup cooling
    if (.not.isothermal) call setup_cool () ! SHOULD BE CALLED ELSEWHERE


  end subroutine rad_ini

  !=======================================================================

  !> Input routine: establish the ionizing spectrum
  subroutine spectrum_parms

    ! Input routine: establish the ionizing spectrum
     
    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    use file_admin, only: stdinput
    
    integer :: nchoice
    real(kind=dp) :: totflux

    ! Ask for input

    ! a) Effective temperature

    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (rank == 0 .and. teff_nominal == 0.0) then
       write(*,'(A)') ' '
       teff=0.0
       do while (teff.lt.2000.0.or.teff.gt.200000.) 
          write(*,'(A,$)') 'Give black body effective temperature: '
          read(stdinput,*) teff
          write(*,*)
          if (teff.lt.2000.0.or.teff.gt.200000.) then
             write(*,*) 'Error: Effective temperature out of range. Try again'
             write(*,*) 'Valid range: 2000 to 200,000'
          endif
       enddo

       ! Find total flux (Stefan-Boltzmann law)
       totflux=sigmasb*teff**4
       
       ! b) Luminosity, radius, or ionizing photon rate?
       write(*,'(A)') ' '
       write(*,'(A)') 'You can specify' 
       write(*,'(A)') ' 1) a stellar radius'
       write(*,'(A)') ' 2) a luminosity'
       write(*,'(A)') ' 3) Total number of ionizing photons'
       nchoice=0
       do while (nchoice <= 0 .or. nchoice > 3)
          write(*,'(A,$)') 'Preferred option (1, 2 or 3): '
          read(stdinput,*) nchoice
          if (nchoice <= 0 .or. nchoice > 3) then
             write(*,*) 'Error: Choose between 1 2 or 3'
          endif
       enddo
       if (nchoice.eq.1) then
          write(*,'(A,$)') 'Give radius in solar radii: '
          read(stdinput,*) rstar
          rstar=rstar*r_solar
          lstar=rstar*rstar*(4.0d0*pi*totflux)
          ! Number of photo-ionizing photons set to zero
          ! determined in spec_diag routine
          S_star=0.0
       elseif (nchoice .eq. 2) then
          write(*,'(A,$)') 'Give luminosity in solar luminosities: '
          read(stdinput,*) lstar
          lstar=lstar*l_solar
          rstar=dsqrt(lstar/(4.0d0*pi*totflux))
          ! Number of photo-ionizing photons set to zero
          ! determined in spec_diag routine
          S_star=0.0
       else
          write(*,'(A,$)') 'Give S_* (ionizing photons s^-1): '
          read(stdinput,*) S_star
          ! Assign some fiducial values, these are scaled to correspond 
          ! to S_star in routine spec_diag
          rstar=r_solar
          lstar=rstar*rstar*(4.0d0*pi*totflux)
       endif
    else
       ! teff and S_star are assumed to have been set in the c2ray_parameter 
       ! module
       teff=teff_nominal
       S_star=S_star_nominal
       totflux=sigmasb*teff**4
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       rstar=r_solar
       lstar=rstar*rstar*(4.0d0*pi*totflux)
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(teff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(rstar,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(lstar,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
#endif
    
  end subroutine spectrum_parms

  !=======================================================================
  
  !> Calculates properties of the black body spectrum
  subroutine spec_diag ()

    ! Calculates properties of spectrum
    ! This version: number of ionizing photons, S*, which can be
    ! used to calculate the Stromgren radius and other photon-statistics

    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    ! Tested against numbers listed on 
    ! http://nimbus.pa.uky.edu/plasma2000/input_for_nebular_models.htm
    ! (19 Feb 2004)

    integer :: i
    real(kind=dp) :: rfr,frmax,stepfl,flux
    real(kind=dp) :: fr(0:NumFreq),weight(0:NumFreq),bb(0:NumFreq)
    real(kind=dp) :: S_star_unscaled,scaling
    
    ! This is h/kT (unit 1/Hz, or sec)
    rfr=hplanck/(kb*teff)

    ! Upper limit of frequency integration
    frmax=min(frtop1,10.0*frtop2)

    ! Frequency step
    stepfl=(frmax-frth0)/real(NumFreq)

    ! Fill the arrays (frequency, weight, spectrum)
    do i=0,NumFreq
       fr(i)=frth0+stepfl*real(i)
       weight(i)=stepfl
       bb(i)=tpic2*fr(i)*fr(i)/(exp(fr(i)*rfr)-1.0)
    enddo

    ! Find flux by integrating
    flux=scalar_romberg(bb,weight,NumFreq,NumFreq,0)

    ! Find out what is the S_star for the radius
    ! supplied.
    S_star_unscaled=4.0*pi*rstar*rstar*flux

    ! If S_star is zero, it is set here.
    if (S_star .eq. 0.0) then
       S_star=S_star_unscaled
    else
       ! Find out the factor by which to change the radius
       ! and luminosity to get the required S_star.
       scaling=S_star/S_star_unscaled
       rstar=sqrt(scaling)*rstar
       lstar=scaling*lstar
    endif
    
    ! Report back
    if (rank == 0) then
       write(logf,'(/a)')           'Using a black body with'
       write(logf,'(a,1pe10.3,a)')   ' Teff=       ',teff,' K'
       write(logf,'(a,1pe10.3,a)')   ' Radius=     ',rstar/r_solar, &
            ' R_solar'
       write(logf,'(a,1pe10.3,a)')   ' Luminosity= ',lstar/l_solar, &
            ' L_solar'
       write(logf,'(A,1PE10.3,A//)') ' Number of H ionizing photons: ', &
            S_star,' s^-1'
    endif

  end subroutine spec_diag
  
  !=======================================================================

  !> Calculates spectral integration cores
  subroutine spec_integr_cores ()

    ! Calculates spectral integration cores

    ! Author: Garrelt Mellema
    ! Date: 19-Feb-2004
    ! Version: Simplified version from Coral.

    ! Note: the calculation of the photo-ionization integrals is split
    ! into two parts. The cores (calculated in this routine) are the parts 
    ! that do not change if the effective temperature and luminosity evolve.
    ! For evolving sources, these parts do not need to be recalculated.
    ! In spec_integr the effective temperature part is added, and the
    ! integration over frequency is performed.

    ! Note 2: the cpu time gain of not recalculating these cores should
    !  really be tested.

    ! Note 3: we calculate two integrals over each rate: one for optically
    ! thick cells (ensuring photon-conservation for those cells), and one 
    ! for optically thin cells. The latter are marked with 1.

    integer :: i,n
    real(kind=dp) :: frmax
    real(kind=dp) :: tau(0:NumTau)
    real(kind=dp) :: fr(0:NumFreq)
    real(kind=dp) :: h0ffr(0:NumFreq)
    
    ! Allocate the spectral integral cores
    allocate(h0int(0:NumFreq,0:NumTau,NumFreqBnd))
    allocate(h0int1(0:NumFreq,0:NumTau,NumFreqBnd))
    if (.not.isothermal) then
       allocate(hh0int(0:NumFreq,0:NumTau,NumFreqBnd))
       allocate(hh0int1(0:NumFreq,0:NumTau,NumFreqBnd))
    endif

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do n=1,NumTau
       tau(n)=10.0**(minlogtau+dlogtau*real(n-1))
    enddo
    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    ! frequency band 1
    ! (there is space for NumFreqBnd frequency bands, only
    ! one is used here).
    if (frth0.lt.frtop1) then
       
       ! Upper limit of frequency integration
       frmax=min(frtop1,10.0*frtop2)
       
       ! Step size in frequency 
       steph0(1)=(frmax-frth0)/real(NumFreq)

       do i=0,NumFreq
          fr(i)=frth0+steph0(1)*real(i)
          
          ! Frequency dependence of the absorption
          ! cross section:
          if (grey) then
             h0ffr(i)=1.0
          else
             h0ffr(i)=(betah0*(fr(i)/frth0)**(-sh0)+ &
                  (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
          endif

          do n=0,NumTau
             ! Protect against floating point errors
             ! This needs to be checked. I remember that
             ! -700 is the minimum exponent allowed for
             ! doubleprecision...
             if (tau(n)*h0ffr(i) < 700.0) then
                h0int(i,n,1)=tpic2*fr(i)*fr(i)* &
                     exp(-tau(n)*h0ffr(i))
                h0int1(i,n,1)=tpic2*fr(i)*fr(i)*h0ffr(i)* &
                     exp(-tau(n)*h0ffr(i))
             else
                h0int(i,n,1)=0.0
             endif
             if (.not.isothermal) then
                hh0int(i,n,1)=hplanck*(fr(i)-frth0)*h0int(i,n,1)
                hh0int1(i,n,1)=hplanck*(fr(i)-frth0)*h0int1(i,n,1)
             endif
          enddo
       enddo
    endif
    
  end subroutine spec_integr_cores
  
  ! =======================================================================
  
  !> Calculates photo ionization integrals
  subroutine spec_integr ()

    ! Calculates photo ionization integrals

    ! Author: Garrelt Mellema
    ! Date: 19-Feb-2004
    ! Version: Simplified from Coral

    ! Two types of integrals are evaluated: one for optically thick cells
    ! (hphot, hheat) and one for optically thin cells (hphot1, hheat1).

    
    integer :: i,n,nfrq
    real(kind=dp) :: rstar2,rfr
    real(kind=dp) :: fr(0:NumFreq),func1(0:NumFreq,0:NumTau)
    real(kind=dp) :: func2(0:NumFreq,0:NumTau)
    real(kind=dp) :: weight(0:NumFreq,0:NumTau),phot(0:NumTau)

    ! Allocate photo-ionization tables
    allocate(hphot(0:NumTau,NumFreqBnd))
    allocate(hphot1(0:NumTau,NumFreqBnd))
    if (.not.isothermal) then
       allocate(hheat(0:NumTau,NumFreqBnd))
       allocate(hheat1(0:NumTau,NumFreqBnd))
    endif

    ! This is h/kT
    rfr=hplanck/(kb*teff)

    ! frequency interval 1
    do i=0,NumFreq
       fr(i)=frth0+steph0(1)*real(i)
       do n=0,NumTau
          weight(i,n)=steph0(1)
          func1(i,n)=h0int(i,n,1)/(exp(fr(i)*rfr)-1.0)
          if (.not.isothermal) func2(i,n)=hh0int(i,n,1)/(exp(fr(i)*rfr)-1.0)
       enddo
    enddo
    
    call vector_romberg (func1,weight,NumFreq,NumFreq,NumTau,phot)
    do n=0,NumTau
       hphot(n,1)=phot(n)
    enddo
    
    if (.not.isothermal) then
       call vector_romberg (func2,weight,NumFreq,NumFreq,NumTau,phot)
       do n=0,NumTau
          hheat(n,1)=phot(n)
       enddo
    endif

    ! frequency interval 1
    do i=0,NumFreq
       fr(i)=frth0+steph0(1)*real(i)
       do n=0,NumTau
          weight(i,n)=steph0(1)
          func1(i,n)=h0int1(i,n,1)/(exp(fr(i)*rfr)-1.0)
          if (.not.isothermal) func2(i,n)=hh0int1(i,n,1)/(exp(fr(i)*rfr)-1.0)
       enddo
    enddo
    
    call vector_romberg (func1,weight,NumFreq,NumFreq,NumTau,phot)
    do n=0,NumTau
       hphot1(n,1)=phot(n)
    enddo
    
    if (.not.isothermal) then
       call vector_romberg (func2,weight,NumFreq,NumFreq,NumTau,phot)
       do n=0,NumTau
          hheat1(n,1)=phot(n)
       enddo
    endif

    ! Multiply with 4*pi*r^2 to make it a luminosity
    rstar2=rstar*rstar
    do nfrq=1,NumFreqBnd
       do n=0,NumTau
          hphot(n,nfrq)=4.0*pi*rstar2*hphot(n,nfrq)
          hphot1(n,nfrq)=4.0*pi*rstar2*hphot1(n,nfrq)
        enddo
    enddo
    if (.not.isothermal) then
       do nfrq=1,NumFreqBnd
          do n=0,NumTau
             hheat(n,nfrq)=4.0*pi*rstar2*hheat(n,nfrq)
             hheat1(n,nfrq)=4.0*pi*rstar2*hheat1(n,nfrq)
          enddo
       enddo
    endif

    ! Deallocate the cores
    deallocate(h0int)
    deallocate(h0int1)
    if (.not.isothermal) then
       deallocate(hh0int)
       deallocate(hh0int1)
    endif

  end subroutine spec_integr

  ! =======================================================================
  
  ! Calculates photo-ionization rates
  subroutine photoion (phi,hcolum_in,hcolum_out,vol,nsrc)
    
    ! Calculates photo-ionization rates
    
    ! Author: Garrelt Mellema
    ! Date: 11-May-2005 (f90) (18 feb 2004
    
    ! Version:
    ! Simplified version derived from Coral version, for testing
    ! photon conservation. Only hydrogen is dealt with, and
    ! one frequency band is used.

    use sourceprops, only: NormFlux

    type(photrates),intent(out) :: phi !< result of the routine
    real(kind=dp),intent(in) :: hcolum_in !< H0 column density at front side
    real(kind=dp),intent(in) :: hcolum_out !< H0 column density at back side
    real(kind=dp),intent(in) :: vol !< volume of shell cell is part of
    integer,intent(in) :: nsrc !< number of the source

    real(kind=dp) :: tauh_in,tauh_out
    real(kind=dp) ::  tau1,odpos1,dodpo1
    integer :: iodpo1,iodp11
    
    ! find the optical depths (in and outgoing)
    tauh_in=sigh*hcolum_in
    tauh_out=sigh*hcolum_out
    
    ! find the table positions for the optical depth (ingoing)
    tau1=log10(max(1.0e-20_dp,tauh_in))
    ! odpos1=min(1.0d0*NumTau,max(0.0d0,1.0d0+(tau1-minlogtau)/
    odpos1=min(real(NumTau,dp),max(0.0_dp,1.0+(tau1-minlogtau)/dlogtau))
    iodpo1=int(odpos1)
    dodpo1=odpos1-real(iodpo1,dp)
    iodp11=min(NumTau,iodpo1+1)
    
    ! Find the hydrogen photo-ionization rate (ingoing)
    ! Since all optical depths are hydrogen, we can use
    ! tau1 for all.
    phi%h_in=NormFlux(nsrc)*(hphot(iodpo1,1)+ &
         (hphot(iodp11,1)-hphot(iodpo1,1))*dodpo1)
    if (.not.isothermal) phi%hv_h_in=NormFlux(nsrc)* &
         (hheat(iodpo1,1)+(hheat(iodp11,1)-hheat(iodpo1,1))*dodpo1)
    
    ! Test for optically thick/thin case
    if (abs(tauh_out-tauh_in) > 1e-2) then 
       
       ! find the table positions for the optical depth (outgoing)
       tau1=log10(max(1.0e-20_dp,tauh_out))
       ! odpos1=min(1.0d0*NumTau,max(0.0d0,1.0d0+(tau1-minlogtau)/
       odpos1=min(real(NumTau,dp),max(0.0_dp,1.0+(tau1-minlogtau)/dlogtau))
       iodpo1=int(odpos1)
       dodpo1=odpos1-real(iodpo1)
       iodp11=min(NumTau,iodpo1+1)
        
       ! find the hydrogen photo-ionization rate (outgoing)
       phi%h_out=NormFlux(nsrc)*(hphot(iodpo1,1)+ &
            (hphot(iodp11,1)-hphot(iodpo1,1))*dodpo1)
       if (.not.isothermal) phi%hv_h_out=NormFlux(nsrc)* &
            (hheat(iodpo1,1)+(hheat(iodp11,1)-hheat(iodpo1,1))*dodpo1)
       
       ! The photon conserving photo-ionization rate is the difference between
       ! the one coming in, and the one going out.
       phi%h=(phi%h_in-phi%h_out)/vol
       if (.not.isothermal) phi%hv_h=(phi%hv_h_in-phi%hv_h_out)/vol
       
    else
       
       ! Find the hydrogen photo-ionization rate for the optically thin
       ! case, and from this derive the outgoing rate.
       ! Since all optical depths are hydrogen, we can use
       ! tau1 for all.
       phi%h=NormFlux(nsrc)*(tauh_out-tauh_in)*( &
            hphot1(iodpo1,1)+(hphot1(iodp11,1)-hphot1(iodpo1,1))*dodpo1)/vol
       phi%h_out=phi%h_in-phi%h*vol

       if (.not.isothermal) then
          phi%hv_h=NormFlux(nsrc)*(tauh_out-tauh_in)*( &
               hheat1(iodpo1,1)+(hheat1(iodp11,1)-hheat1(iodpo1,1))*dodpo1)/vol
          phi%hv_h_out=phi%hv_h_in-phi%hv_h*vol
       endif

    endif

  end subroutine photoion

end module radiation

