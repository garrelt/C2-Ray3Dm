!     This module contains data and routines which deal with radiative
!     effects. Its main part deals with photo-ionizing radiation, but it
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
  use sed_parameters, only: stellar_SED_type, & ! Type of stellar SED
       bb_Teff,&            ! Black body  effective temperature for for nominal BB SED
       bb_S_star, &          ! Ionizing photon rate for for nominal BB SED
       bb_MinFreq, &      ! Lowest frequency for nominal PL SED
       bb_MaxFreq, &         ! Highest frequency for nominal PL SED
       pl_index,&         ! Power law index for for nominal PL SED
       EddLeff_nominal,&          ! Eddington efficiency for for nominal PL SED
       EddLum, &                  ! Eddington luminosity for for nominal PL SED
       pl_S_star, &       ! Ionizing photon rate for nominal PL SED
       pl_MinFreq, &      ! Lowest frequency for nominal PL SED
       pl_MaxFreq         ! Highest frequency for nominal PL SED
  !use material, only: isothermal
  use c2ray_parameters, only: isothermal

  use radiation_sizes, only: NumFreq,NumFreqBnd
  use radiation_sizes, only: freq_min,freq_max

  implicit none

  ! Number of photon of the power law source
  real(kind=dp) :: pl_input_flux = 0.0

  ! Stellar SED properties
  character :: sourcetype = " " ! Type of source, B=black body, &
  ! P=power law source
  character,dimension(2) :: sourcetype_string=(/ "B", "P" /)
  real(kind=dp) :: L_star_ion = 0.0  ! SED Ionizing luminosity
  real(kind=dp) :: S_star = 0.0      ! SED ionizing photons rate
  real(kind=dp) :: MinFreq = 0.0     ! SED minimum frequency
  real(kind=dp) :: MaxFreq = 0.0     ! SED maximum frequency

  ! X-ray SED properties
  character :: xray_sourcetype = " " ! Type of source, H=HMXB, &
  ! Q=quasar
  character,dimension(2) :: xray_sourcetype_string=(/ "H", "Q" /)
  real(kind=dp) :: L_xray = 0.0  ! X-ray Ionizing luminosity
  real(kind=dp) :: S_star_xray = 0.0  ! X-ray ionizing photons rate

  ! Stellar properties
  real(kind=dp) :: T_eff  = 0.0      ! Black body effective temperature
  real(kind=dp) :: R_star = 0.0      ! Black body radius
  real(kind=dp) :: R_star2 = 0.0     ! Square of R_star
  real(kind=dp) :: L_star = 0.0      ! Black body luminosity
  real(kind=dp) :: h_over_kT         ! Planck constant over k_B * T_eff

  ! Power law source properties
  real(kind=dp) :: S_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: Edd_Efficiency = 0.0 ! Eddinton efficieny

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms (freq_min_src, freq_max_src)
    
    use file_admin, only: stdinput, file_input
    
    real(kind=dp),intent(out) :: freq_min_src, freq_max_src

    integer :: i_choice = 0                 ! option number
    real(kind=dp) :: bb_luminosity_unscaled ! black body total flux
    
    ! All SED parameters are expected to be set in sed_parameters.f90
    ! Previous versions of the code allowed interactive specification
    ! of the SED but this has been dropped in this version.

    ! Set source type (from sed_parameters)
    sourcetype=sourcetype_string(stellar_SED_type)

    select case (sourcetype)

       ! Black body case
    case ("B")

       ! Report to log file
       write(logf,*) 'Black body source'
       
       ! In blackbody case, limit the effective temperate. 
       ! It should be higher than 2000 and lower than 1000000.
       T_eff=max(min(bb_teff,1e6),2000.)
       write(logf,*) 'Temperature of the black body : ', T_eff
             
       ! Set rate of ionzing photons (S_star)
       S_star=bb_S_star

       ! Set minimum and maximum frequency
       MinFreq=bb_MinFreq
       MaxFreq=bb_MaxFreq

       ! Find total flux of blackbody (Stefan-Boltzmann law)
       bb_luminosity_unscaled = sigma_SB*T_eff*T_eff*T_eff*T_eff

       ! Assign some fiducial values for R_star and L_star, 
       ! these are scaled to correspond to S_star in routine spec_diag
       R_star=r_solar
       L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
       
       ! This is h/kT
       h_over_kT = hplanck/(k_B*T_eff)
       
    case("P")
    
       write(logf,*) 'Power law source'
       write(logf,*) 'Power law index is ', pl_index

       ! Set rate of ionzing photons (S_star)
       S_star=pl_S_star

       ! Set minimum and maximum frequency
       MinFreq=pl_MinFreq
       MaxFreq=pl_MaxFreq

    end select

    ! Report
    write(logf,*) 'Rate of ionzing photons (S_star) is specified'
    write(logf,*) 'The number of ionizing photons per second is ', S_star
    write(logf,*) 'The lower energy limit is ', MinFreq/ev2fr, ' eV'
    write(logf,*) 'The upper energy limit is ', MaxFreq/ev2fr, ' eV'
    
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(T_eff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
#endif
    
    ! Set maximum frequency for integration.
    freq_max_src=MaxFreq

  end subroutine spectrum_parms
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine spectrum diagnostics
  ! This routine integrates the SEDs in order to achieve proper scaling
  ! of them. The spectrum_parms subroutine has set some of the SED properties
  ! but they are not necessarily consistent. Consistency is set here
  subroutine spec_diag ()
    
    real(kind=dp) :: bb_ionizing_luminosity_unscaled
    real(kind=dp) :: S_star_unscaled
    real(kind=dp) :: Ionizing_Luminosity

    if (sourcetype == "B") then

       ! Now we know R_star and L_star (except in the case where S_star is provided).
       ! So we can continue to find blackbody flux 
          
       ! Black-body flux (photon sense)
       S_star_unscaled = integrate_sed(MinFreq,MaxFreq,"B","S")
       ! If S_star is specified, then R_star and L_star have to be tuned accordingly.
       S_scaling = S_star/S_star_unscaled
       R_star = sqrt(S_scaling)*R_star
       L_star = S_scaling*L_star
       Ionizing_Luminosity = integrate_sed(MinFreq,MaxFreq,"B","L")
 
       ! This is R_star^2
       R_star2=R_star*R_star
       
       ! Report back to the log file
       if (rank == 0) then
          write(logf,'(/a)')           'Using a black body with'
          write(logf,'(a,es10.3,a)')   ' Teff =       ', T_eff, ' K'
          write(logf,'(a,es10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
          write(logf,'(a,es10.3,a)')   ' Ionizing Luminosity = ', Ionizing_Luminosity/l_solar, ' L_solar'
          write(logf,'(a,es10.3,a)')   ' Ionzing photon rate = ', S_star, ' s^-1'
       endif
    endif
    
    if (sourcetype == "P") then
       ! Determine the scaling factor S_scaling
       ! needed to achieve either the specified ionizing photon rate or ionizing luminosity
       ! Total power-law ionizing photon rate is specified (photon sense)
       S_star_unscaled = integrate_sed(MinFreq,MaxFreq,"P","S")
       S_scaling = S_star/S_star_unscaled
       write(logf,*) S_star,S_star_unscaled,S_scaling
       Ionizing_Luminosity = integrate_sed(MinFreq,MaxFreq,"P","L")
       
       ! Report back to the log file
       if (rank == 0) then
          write(logf,'(/a)')           'Using a power law source with'
          write(logf,'(a,es10.3)')   ' Power law index = ', pl_index
          write(logf,'(a,es10.3)')   ' Ionizing photon rate = ', S_star
          write(logf,'(a,es10.3)')   ' Ionizing luminosity = ', Ionizing_Luminosity/L_solar
          write(logf,'(a,2f8.3,a)')   ' between energies ', &
               pl_MinFreq/(ev2fr),pl_MaxFreq/(ev2fr),' Ev'
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
       integrate_sed = S_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)

    end select

  end function integrate_sed

end module radiation_sed_parameters
