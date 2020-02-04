!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_sed_parameters
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigma_SB, &                    ! Stefan-Boltzmann constant
                          hplanck, &                     ! Planck constant
                          k_B, &                         ! Boltzmann constant
                          two_pi_over_c_square, &        ! two times pi over c aquare
                          ev2fr                          ! eV to Hz conversion
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               sigma_HI_at_ion_freq      ! HI cross section at its ionzing frequency
  use astroconstants, only: R_SOLAR, &                   ! Solar radius
                            L_SOLAR                      ! Solar luminosity
  use romberg, only: scalar_romberg, &                   ! 1D integration function
                     vector_romberg, &                   ! 1D integration subroutine
                     romberg_initialisation              ! Romberg initialisation procedure
  use c2ray_parameters, only: T_eff_nominal => teff_nominal,&            ! Black body  effective temperature for for nominal BB SED
                              S_star_nominal, &          ! Ionizing photon rate for for nominal BB SED
                              pl_index_nominal,&         ! Power law index for for nominal PL SED
                              EddLeff_nominal,&          ! Eddington efficiency for for nominal PL SED
                              EddLum, &                  ! Eddington luminosity for for nominal PL SED
                              pl_S_star_nominal, &       ! Ionizing photon rate for nominal PL SED
                              pl_MinFreq_nominal, &      ! Lowest frequency for nominal PL SED
                              pl_MaxFreq_nominal         ! Highest frequency for nominal PL SED
  !use material, only: isothermal
  use c2ray_parameters, only: isothermal

  use radiation_sizes, only: NumFreq,NumFreqBnd
  use radiation_sizes, only: freq_min,freq_max

  implicit none

  ! Number of photon of the power law source
  real(kind=dp) :: pl_input_flux = 0.0

  ! Type of source, B=black body, P=power law source
  Character :: sourcetype = " "

  ! Stellar properties
  real(kind=dp) :: T_eff  = 0.0      ! Black body effective temperature
  real(kind=dp) :: R_star = 0.0      ! Black body radius
  real(kind=dp) :: L_star = 0.0      ! Black body luminosity
  real(kind=dp) :: L_star_ion = 0.0 ! Black body ionizing luminosity
  real(kind=dp) :: S_star = 0.0      ! Black body ionizing photons rate
  real(kind=dp) :: R_star2 = 0.0     ! Square of R_star
  real(kind=dp) :: h_over_kT    ! Planck constant over k_B * T_eff
  real(kind=dp) :: bb_MaxFreq           ! Maximum frequency for integration of total power of BB
  !real(kind=dp) :: freq_max_src ! Maximum frequency for integration of total power

  ! The lowest and highest frequency subbands used for the bb and pl source
  integer :: bb_FreqBnd_UpperLimit=NumFreqBnd
  integer :: pl_FreqBnd_LowerLimit
  integer :: pl_FreqBnd_UpperLimit

  ! Power law source properties
  real(kind=dp) :: pl_index = 1.0            ! Power law index
  real(kind=dp) :: pl_MinFreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_MaxFreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: Edd_Efficiency = 0.0 ! Eddinton efficieny
  real(kind=dp) :: pl_S_star = 1.0

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms (freq_max_src)
    
    use file_admin, only: stdinput, file_input
    
    real(kind=dp),intent(out) :: freq_max_src

    integer :: i_choice = 0                 ! option number
    real(kind=dp) :: bb_luminosity_unscaled ! black body total flux
    
    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (T_eff_nominal == 0.0) then
       
       if (rank == 0) then
          do while ((sourcetype  /= 'B') .and. (sourcetype /= 'P'))
             
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), power law source (P), ", &
                  "or both (A)"
             read(stdinput,*) sourcetype ! Read source type, either blackbody or power law source
             
          enddo
          
          !In blackbody case, ask for effective temperature of the source. 
          !The temperature should be bounded below by 2000 and above by 1000000.
          !And then ask for some parameters of the blackbody.
          if (sourcetype == 'B' .or. sourcetype == "A") then
             
             write(logf,*) 'Black body source'
             
             do while (T_eff < 2000.0 .or. T_eff > 1000000.) 
                if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature (2000 <= T <= 1000000): '
                read(stdinput,*) T_eff      ! Read temperature of black body
                write(logf,*) 'Temperature of the black body : ', T_eff
                if (T_eff < 2000.0 .or. T_eff > 1000000.) then
                   write(*,*) 'Error: Effective temperature out of range. Try again'
                endif
             enddo
             
             ! Find total flux of blackbody (Stefan-Boltzmann law)
             bb_luminosity_unscaled = sigma_SB*T_eff*T_eff*T_eff*T_eff
          
             ! Ask for radius, luminosity, ionizing luminosity or ionizing photon rate?
             if (.not.file_input) then
                write(*,'(A)') 'You can specify' 
                write(*,'(A)') ' 1) a stellar radius'
                write(*,'(A)') ' 2) a total luminosity'
                write(*,'(A)') ' 3) Total ionizing luminosity'
                write(*,'(A)') ' 4) Total number of ionizing photons'
             endif
             
             ! Report error if options are not 1, 2, 3 and 4
             do while (i_choice <= 0 .or. i_choice > 4)
                if (.not.file_input) write(*,'(A,$)') 'Preferred option (1, 2, 3 or 4): '
                read(stdinput,*) i_choice       ! Read option from the input, 1 to 4
                if (i_choice <= 0 .or. i_choice > 4) then
                   write(*,*) 'Error: Choose between 1 2 3 or 4'
                endif
             enddo
             
             select case (i_choice)
                
             case (1)
                write(logf,*) 'A stellar radius is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give radius in solar radius: '
                read(stdinput,*) R_star      ! Read radius of the black body
                write(logf,*) 'The radius is ', R_star, ' solar radius'
                R_star=R_star*r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                S_star=0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (2)
                write(logf,*) 'A total luminosity is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give total luminosity in solar luminosity: '
                read(stdinput,*) L_star      ! Read luminosity of the black body
                write(logf,*) 'The luminosity is ', L_star, ' solar luminosity'
                L_star=L_star*l_solar
                R_star=dsqrt(L_star/(4.0d0*pi*bb_luminosity_unscaled))
                S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (3)
                write(logf,*) 'Total ionizing luminosity is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give total ionizing luminosity in solar luminosity: '
                read(stdinput,*) L_star_ion   ! Read ionizing luminosity of the black body
                write(logf,*) 'Total ionizing luminosity is ', L_star_ion, ' solar luminosty'
                L_star_ion=L_star_ion*l_solar
                ! Assign some fiducial values, these are overwritten in routine spec_diag
                R_star=r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (4)
                write(logf,*) 'Rate of ionzing photons (S_star) is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give the number of ionizing photons per second: '
                read(stdinput,*) S_star
                write(logf,*) 'The number of photons per second is ', S_star
                ! Assign some fiducial values for R_star and L_star, 
                ! these are scaled to correspond to S_star in routine spec_diag
                R_star=r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                
             end select

          endif
          ! In power-law case, we ask for number of ionizing photons per second or Eddinton luminosity efficiency 
          if (sourcetype == 'P' .or. sourcetype == 'A') then
             
             write(logf,*) 'Power law source'
             if (.not.file_input) then
                write(*,'(A)') 'You can specify'
                write(*,'(A)') ' 1)  Number of ionizing photons per second '
                write(*,'(A)') ' 2)  Efficiency parameter assuming a 1e6 solar mass BH' 
             endif
             
             ! Report error if options are not 1 and 2
             do while (i_choice <= 0 .or. i_choice > 2)
                if (.not.file_input) write(*,'(A,$)') 'Preferred option (1 or 2): '
                read(stdinput,*) i_choice
                if (i_choice <= 0 .or. i_choice > 2) then
                   write(*,*) 'Error: Choose between 1 or 2'
                endif
             enddo
             
             select case (i_choice)
                
             case (1)
                write(logf,*) 'Rate of ionizing photons is specified'
                if (.not.file_input) write(*,'(A,$)') 'give number of ionizing photons per second '
                read(stdinput,*) pl_S_star          ! Read ionizing photons per second
                write(logf,*) 'The rate is ', pl_S_star
                
                ! Set the Eddinton luminosity efficiency to the nominal value 
                Edd_Efficiency=EddLeff_nominal
                pl_scaling=1.0 ! fiducial value, to be updated in spec_diag
                
             case (2)
                write(logf,*) 'Efficiency parameter is specified'
                if (.not.file_input) write(*,'(A,$)') 'give efficiency parameter '
                read(stdinput,*) Edd_Efficiency         ! Read Eddington efficiency
                write(logf,*) 'The efficiency parameter is ', Edd_Efficiency
                ! Set some fiducial value, to be updated in spec_diag
                pl_S_star=0.0
                pl_scaling=1.0
             end select
             
             if (.not.file_input) write(*,'(A,$)') 'Specify power law index (for number of photons, not energy) '
             read(stdinput,*) pl_index      ! Read power law index, this number equal to one plus that of energy 
             write(logf,*) 'Power law index is ', pl_index
             if (.not.file_input) write(*,'(A,$)') 'give lower and upper frequency limits in eV '
             read(stdinput,*) pl_MinFreq,pl_MaxFreq     ! Read lower and upper frequency limits in eV	
             write(logf,*) 'The lower energy limit is ', pl_MinFreq, ' eV'
             write(logf,*) 'The upper energy limit is ', pl_MaxFreq, ' eV'

             ! Convert eVs to Hz for the frequency limits
             pl_MinFreq = pl_MinFreq * ev2fr
             pl_MaxFreq = pl_MaxFreq * ev2fr

             ! set some fiducial values for the BB source here, though they are not useful
             R_star=r_solar
             S_star=0.0
             L_star=0.0
             T_eff=1.0e5
             
          endif
       endif
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(T_eff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(R_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(L_star,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(L_star_ion,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(Edd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
#endif
    
       ! In case source properties are taken from c2ray_parameters
    else
       ! T_eff and S_star are assumed to have been set in the c2ray_parameter module
       T_eff=T_eff_nominal
       S_star=S_star_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       R_star=r_solar
       L_star=R_star*R_star*4.0d0*pi*sigma_SB*T_eff**4
       
       ! Power law source properties set to nominal values from c2ray_parameter module
       pl_index=pl_index_nominal
       pl_MinFreq=pl_MinFreq_nominal
       pl_MaxFreq=pl_MaxFreq_nominal
       pl_S_star=pl_S_star_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       pl_scaling=1.0
       Edd_Efficiency=0.0
    endif
    
    ! This is h/kT
    h_over_kT = hplanck/(k_B*T_eff)

    ! Set maximum frequency for integration. For the BB source this is calculated
    ! from the exponent in the BB equation, for the PL source this has been specified.
    ! We take the maximum of the two
    bb_MaxFreq = 25.0/h_over_kT
    freq_max_src=max(pl_MaxFreq,bb_MaxFreq)

  end subroutine spectrum_parms
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine spectrum diagnostics
  ! This routine integrates the SEDs in order to achieve proper scaling
  ! of them. The spectrum_parms subroutine has set some of the SED properties
  ! but they are not necessarily consistent. Consistency is set here
  subroutine spec_diag ()
    
    real(kind=dp) :: bb_ionizing_luminosity_unscaled
    real(kind=dp) :: S_star_unscaled, S_scaling
    real(kind=dp) :: pl_S_star_unscaled
    real(kind=dp) :: pl_S_star_wanted
    real(kind=dp) :: pl_ionizing_luminosity_unscaled
    real(kind=dp) :: pl_ionizing_luminosity
    real(kind=dp) :: pl_ionizing_luminosity_wanted
    
    ! Case L_star_ion is provided. Find R_star, L_star and blackblody flux
    if (L_star_ion /= 0.0d0) then
       
       ! Black-body ionizing luminosity (energy sense)
       bb_ionizing_luminosity_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","L")
       ! Find radius from the scaled and specified ionizing luminosities
       R_star = sqrt(L_star_ion/bb_ionizing_luminosity_unscaled)*R_star
       ! Find total luminosity from Stefan-Boltzmann law
       L_star = R_star*R_star*4.0_dp*pi*sigma_SB*T_eff**4
       ! Find the S_star
       S_star = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","S")
       
    else
       
       ! Case L_star_ion is not provided.  
       ! Now we know R_star and L_star (except in the case where S_star is provided).
       ! So we can continue to find blackbody flux 
       
       ! Black-body flux (photon sense)
       S_star_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","S")
       ! If S_star is zero, it is set here.
       if (S_star == 0.0) then
          S_star=S_star_unscaled
       else
          ! If S_star is specified, then R_star and L_star have to be tuned accordingly.
          S_scaling = S_star/S_star_unscaled
          R_star = sqrt(S_scaling)*R_star
          L_star = S_scaling*L_star
          L_star_ion = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","L")
       endif
    endif
       
    ! This is R_star^2
    R_star2=R_star*R_star

    ! Report back to the log file
    if (rank == 0) then
       if (sourcetype == 'B' .or. sourcetype == " " .or. sourcetype == "A") then
          write(logf,'(/a)')           'Using a black body with'
          write(logf,'(a,es10.3,a)')   ' Teff =       ', T_eff, ' K'
          write(logf,'(a,es10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
          write(logf,'(a,es10.3,a)')   ' Luminosity = ', L_star/l_solar, ' L_solar'
          write(logf,'(a,es10.3,a)')   ' Ionzing photon rate = ', S_star, ' s^-1'
       endif
    endif

    ! Determine the scaling factor pl_scaling
    ! needed to achieve either the specified ionizing photon rate or ionizing luminosity
    if (pl_S_star > 0.0) then
       ! Total power-law ionizing photon rate is specified (photon sense)
       pl_S_star_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
       !pl_S_star_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","S")
       pl_S_star_wanted = pl_S_star
       pl_scaling = pl_S_star_wanted/pl_S_star_unscaled
       pl_ionizing_luminosity = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       !pl_ionizing_luminosity = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","L")
       Edd_Efficiency = pl_ionizing_luminosity / EddLum
    else
       ! The power-law ionizing luminosity is specified (energy sense). 
       pl_ionizing_luminosity = EddLum*Edd_Efficiency
       pl_ionizing_luminosity_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       !pl_ionizing_luminosity_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","L")
       pl_ionizing_luminosity_wanted = pl_ionizing_luminosity
       pl_scaling = pl_ionizing_luminosity_wanted/pl_ionizing_luminosity_unscaled
       pl_S_star = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
       !pl_S_star = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","S")
    endif
    
    ! Report back to the log file
    if (rank == 0) then
       if (sourcetype == 'P' .or. sourcetype == " " .or. sourcetype == "A") then
          write(logf,'(/a)')           'Using a power law source with'
          write(logf,'(a,es10.3)')   ' Power law index = ', pl_index
          write(logf,'(a,es10.3)')   ' Efficiency parameter = ', Edd_Efficiency
          write(logf,'(a,es10.3)')   ' Ionizing photon rate = ', pl_S_star
          write(logf,'(a,es10.3)')   ' Ionizing luminosity = ', pl_ionizing_luminosity
          write(logf,'(a,2f8.3,a)')   ' between energies ', &
               pl_MinFreq/(1e3*ev2fr),pl_MaxFreq/(1e3*ev2fr),' kEv'
       endif
    endif
    
  end subroutine spec_diag

  function integrate_sed(freq_min, freq_max, sourcetype, sedtype)
    
    ! function type
    real(kind=dp) :: integrate_sed

    ! arguments
    real(kind=dp),intent(in) :: freq_min
    real(kind=dp),intent(in) :: freq_max
    character(len=1),intent(in) :: sourcetype ! P or B
    character(len=1),intent(in) :: sedtype ! L or S

    integer :: i_freq
    real(kind=dp) :: freq_step
    real(kind=dp),dimension(0:NumFreq) :: frequency
    real(kind=dp),dimension(0:NumFreq) :: weight
    real(kind=dp),dimension(0:NumFreq) :: integrand

    ! Set default answer
    integrate_sed=0.0

    ! Set the frequency step for the integration
    freq_step=(freq_max-freq_min)/real(NumFreq)
    
    ! Fill the frequency array and the weight array
    do i_freq=0,NumFreq
       frequency(i_freq)=freq_min+freq_step*real(i_freq)
       weight(i_freq)=freq_step
    enddo

    select case (sourcetype)

    case("B")
       do i_freq=0,NumFreq
          if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
             ! this blackbody is in number of photon sense
             integrand(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                  (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
             ! when the argument of the exponential function gets too high
          else
             integrand(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                  (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                  (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
          endif
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = 4.0*pi*R_star*R_star*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)

    case("P")
       do i_freq=0,NumFreq
          ! this power-law is in number of photon sense
          integrand(i_freq)=frequency(i_freq)**(-pl_index)
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = pl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)

    end select

  end function integrate_sed

end module radiation_sed_parameters
