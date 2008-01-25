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
  !     - source_position : Input routine: establish the source position
  !
  !     Needs following `modules':
  !     romberg : romberg integrators
  !
  !     Author: Garrelt Mellema
  ! 
  !     Date: 02-Jun-2004 (04-Mar-2004)
  
  ! Version
  ! Simplified version from Coral, for testing
  ! - Only hydrogen
  ! - Grey photo-ionization cross section
  
  use precision
  use my_mpi
  use romberg
  use c2ray_parameters

  implicit none

  !-----------------------------------------------------------------------
  !     NumFreq - Number of integration points in one of the three 
  !               frequency interval.
  !     NumTau -  Number of table points for the optical depth.
  !-----------------------------------------------------------------------

  integer,parameter :: NumFreq=128
  integer,parameter :: NumTau=2000
  integer,parameter :: NumFreqBnd=1
  
  real(kind=dp),parameter :: minlogtau=-20.0
  real(kind=dp),parameter :: maxlogtau=4.0
  real(kind=dp),parameter :: dlogtau=(maxlogtau-minlogtau)/real(NumTau)

  ! stellar properties
  real(kind=dp) :: teff,rstar,lstar,S_star
  
  real(kind=dp),dimension(NumFreqBnd) :: steph0
  !real(kind=dp),dimension(0:NumFreq,0:NumTau,NumFreqBnd) :: h0int
  real(kind=dp),dimension(:,:,:),allocatable :: h0int
  !real(kind=dp),dimension(0:NumFreq,0:NumTau,NumFreqBnd) :: hh0int
  real(kind=dp),dimension(:,:,:),allocatable :: hh0int

  real(kind=dp),dimension(0:NumTau,NumFreqBnd) :: hphot
  !real(kind=dp),dimension(0:NumTau,NumFreqBnd) :: hheat
  real(kind=dp),dimension(:,:),allocatable :: hheat
  
  real(kind=dp) :: tauHI=0.0

  type photrates
     real(kind=dp) :: h
     real(kind=dp) :: hv_h
     real(kind=dp) :: h_in
     real(kind=dp) :: hv_h_in
     real(kind=dp) :: h_out
     real(kind=dp) :: hv_h_out
  end type photrates

  ! photo-ionization rates
  !real(kind=dp),public :: phih,hvphih
  !real(kind=dp),public :: phih_in,phih_out
  !real(kind=dp),public :: hvphih_in,hvphih_out

#ifdef MPI       
    integer,private :: ierror
#endif

contains

!=======================================================================

  subroutine rad_ini ()

    ! initializes constants and tables for radiation processes
    ! (heating, cooling and ionization)

    use radiative_cooling

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

  subroutine spectrum_parms

    ! Input routine: establish the ionizing spectrum
     
    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    use cgsconstants
    use astroconstants
    use file_admin, only: stdinput
    
    integer :: nchoice
    real(kind=dp) :: totflux

    ! Ask for input

    ! a) Effective temperature

    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
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
       ! teff and S_star were set in the parameter module
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
  
  subroutine spec_diag ()

    ! Calculates properties of spectrum
    ! This version: number of ionizing photons, S*, which can be
    ! used to calculate the Stromgren radius.

    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    ! Tested against numbers listed on 
    ! http://nimbus.pa.uky.edu/plasma2000/input_for_nebular_models.htm
    ! (19 Feb 2004)

    use cgsconstants
    use cgsphotoconstants
    use astroconstants

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
       write(30,'(/a)')           'Using a black body with'
       write(30,'(a,1pe10.3,a)')   ' Teff=       ',teff,' K'
       write(30,'(a,1pe10.3,a)')   ' Radius=     ',rstar/r_solar, &
            ' R_solar'
       write(30,'(a,1pe10.3,a)')   ' Luminosity= ',lstar/l_solar, &
            ' L_solar'
       write(30,'(A,1PE10.3,A//)') ' Number of H ionizing photons: ', &
            S_star,' s^-1'
    endif

  end subroutine spec_diag
  
  !=======================================================================

  subroutine spec_integr_cores ()

    ! Calculates spectral integration cores

    ! Author: Garrelt Mellema
    ! Date: 19-Feb-2004
    ! Version: Simplified version from Coral.

    use cgsconstants
    use cgsphotoconstants
    
    logical,parameter :: grey=.false. ! use grey opacities?

    integer :: i,n
    real(kind=dp) :: frmax
    real(kind=dp) :: tau(0:NumTau)
    real(kind=dp) :: fr(0:NumFreq)
    real(kind=dp) :: h0ffr(0:NumFreq)
    
    ! Allocate the spectral integral cores
    allocate(h0int(0:NumFreq,0:NumTau,NumFreqBnd))
    if (.not.isothermal) allocate(hh0int(0:NumFreq,0:NumTau,NumFreqBnd))

    ! fill the optica depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    !dlogtau=(maxlogtau-minlogtau)/real(NumTau)
    do n=1,NumTau
       tau(n)=10.0**(minlogtau+dlogtau*real(n-1))
    enddo
    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey) write(30,*) 'WARNING: Using grey opacities'

    ! frequency band 1
    ! (there is space for NumFreqBnd frequency bands, only
    ! one is used here).
    if (frth0.lt.frtop1) then
       
       ! Upper limit of frequency integration
       frmax=min(frtop1,10.0*frtop2)
       
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
             if (tau(n)*h0ffr(i) .lt. 700.0) then
                h0int(i,n,1)=tpic2*fr(i)*fr(i)* &
                     exp(-tau(n)*h0ffr(i))
             else
                h0int(i,n,1)=0.0
             endif
             if (.not.isothermal) hh0int(i,n,1)= &
                  hplanck*(fr(i)-frth0)*h0int(i,n,1)
          enddo
       enddo
    endif
    
  end subroutine spec_integr_cores
  
  ! =======================================================================
  
  subroutine spec_integr ()

    ! Calculates photo ionization integrals

    ! Author: Garrelt Mellema
    ! Date: 19-Feb-2004
    ! Version: Simplified from Coral
    ! It still uses THREE frequency intervals (needed
    ! for H + He), but the actually only the first
    ! is used.
    
    use cgsconstants
    use cgsphotoconstants
    
    integer :: i,n,nfrq
    real(kind=dp) :: rstar2,rfr
    real(kind=dp) :: fr(0:NumFreq),func1(0:NumFreq,0:NumTau)
    real(kind=dp) :: func2(0:NumFreq,0:NumTau)
    real(kind=dp) :: weight(0:NumFreq,0:NumTau),phot(0:NumTau)

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
       allocate(hheat(0:NumTau,NumFreqBnd))
       call vector_romberg (func2,weight,NumFreq,NumFreq,NumTau,phot)
       do n=0,NumTau
          hheat(n,1)=phot(n)
       enddo
    endif

    ! Multiply with 4*pi*r^2 to make it a luminosity
    rstar2=rstar*rstar
    do nfrq=1,NumFreqBnd
       do n=0,NumTau
          hphot(n,nfrq)=4.0*pi*rstar2*hphot(n,nfrq)
          if (.not.isothermal) hheat(n,nfrq)=4.0*pi*rstar2*hheat(n,nfrq)
       enddo
    enddo

    ! Deallocate the cores
    deallocate(h0int)
    if (.not.isothermal) deallocate(hh0int)
  end subroutine spec_integr

  ! =======================================================================
  
  subroutine photoion (phi,hcolum_in,hcolum_out,vol,nsrc)
    
    ! Calculates photo-ionization rates
    
    ! Author: Garrelt Mellema
    ! Date: 11-May-2005 (f90) (18 feb 2004
    
    ! Version:
    ! Simplified version derived from Coral version, for testing
    ! photon conservation. Only hydrogen is dealt with, and
    ! one frequency band is used.

    ! Notes:
    ! Need to add optically thin limit (GM, 060919)

    use cgsphotoconstants
    use sourceprops, only: NormFlux

    type(photrates),intent(out) :: phi
    real(kind=dp),intent(in) :: hcolum_in,hcolum_out,vol
    integer,intent(in) :: nsrc

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
    
    ! find the table positions for the optical depth (outgoing)
    tau1=log10(max(1.0e-20_dp,tauh_out))
    ! odpos1=min(1.0d0*NumTau,max(0.0d0,1.0d0+(tau1-minlogtau)/
    odpos1=min(real(NumTau,dp),max(0.0_dp,1.0+(tau1-minlogtau)/dlogtau))
    iodpo1=int(odpos1)
    dodpo1=odpos1-real(iodpo1,dp)
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

  end subroutine photoion

end module radiation

