!     This module contains data and routines for tables with photoionization
!     and heating rates 
!
!     These can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_tables
  
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
                               sigma_HI_at_ion_freq      ! HI cross section at i
  use astroconstants, only: R_SOLAR, &                   ! Solar radius
                            L_SOLAR                      ! Solar luminosity
  use romberg, only: scalar_romberg, &                   ! 1D integration function
                     vector_romberg, &                   ! 1D integration subroutine
                     romberg_initialisation              ! Romberg initialisation procedure

  !use material, only: isothermal
  use c2ray_parameters, only: isothermal

  use radiation_sizes, only: NumFreq,NumFreqBnd,NumBndin1
  use radiation_sizes, only: Numheatbin
  use radiation_sizes, only: freq_min,freq_max,delta_freq
  use radiation_sizes, only: NumTau
  use radiation_sizes, only: pl_index_cross_section_HI
  use radiation_sizes, only: setup_scalingfactors

  use radiation_sed_parameters, only: sourcetype
  use radiation_sed_parameters, only: T_eff, R_star, R_star2, h_over_kT
  use radiation_sed_parameters, only: pl_index, pl_scaling
  use radiation_sed_parameters, only: pl_minfreq, pl_maxfreq
  use radiation_sed_parameters, only: spectrum_parms, spec_diag
  
  implicit none

  ! Parameters defining the optical depth entries in the table.
  real(kind=dp),parameter :: minlogtau = -20.0                             ! Table position starts at log10(minlogtau) 
  real(kind=dp),parameter :: maxlogtau = 4.0                               ! Table position ends at log10(maxlogtau) 
  real(kind=dp),parameter :: dlogtau = (maxlogtau-minlogtau)/real(NumTau)  ! dlogtau is the step size in log10(tau)

  ! Logical that determines the use of grey opacities
  logical,parameter :: grey = .false. 

  ! The lowest and highest frequency subbands used for the bb and pl source
  integer :: bb_FreqBnd_LowerLimit=1
  integer :: bb_FreqBnd_UpperLimit=NumFreqBnd
  integer :: pl_FreqBnd_LowerLimit
  integer :: pl_FreqBnd_UpperLimit

  ! Integrands ( frequency, optical depth )
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HI

  ! Frequency dependence of cross section (subband dependent)
  real(kind=dp), dimension(0:NumFreq) :: cross_section_freq_dependence

  ! Weights
  real(kind=dp), dimension(0:NumFreq, 0:NumTau) :: vector_weight

  ! Frequency array
  real(kind=dp), dimension(0:NumFreq) :: frequency

  ! Optical depth array
  real(kind=dp), dimension(0:NumTau) :: tau

  ! Integration table ( optical depth, sub-bin )
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thick_table 
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thin_table

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                             !
  subroutine rad_ini ()                                                                      !

    ! Contains the maximum frequency, determined by source input
    ! passed from spectrum_parms to setup_scalingfactors
    real(kind=dp) :: freq_max_src

    !
    ! Ask for the parameters of the spectrum                                                 !
    call spectrum_parms (freq_max_src)                                                                   ! 
                                                                                             ! 
    ! Initializes constants and tables for radiation processes                               ! 
    ! (heating, cooling and ionization)                                                      !
    call setup_scalingfactors (freq_max_src)                                                             !
                                                                                             !
    ! Initialize integration routines                                                        !
    call romberg_initialisation (NumFreq)                                                    !
                                                                                             !
    ! Determine spectrum diagnostics                                                         !
    call spec_diag ()                                                                        !

#ifdef MPILOG
    write(logf,*) 'about to integrate'
#endif
                                                                                           !                                                                                             !
    ! Generate photoionization tables and heating tables                                     !
    call spec_integration ()                                                                 !                            

#ifdef MPILOG
    write(logf,*) 'end of radiation'
#endif
                                                                                           !
  end subroutine rad_ini                                                                     !                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Generate photoionization tables and heating tables   
  subroutine spec_integration ()

    integer :: i_tau
    integer :: i_subband

    ! Allocate integrands and table arrays
    call allocate_spec_integration_arrays

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do i_tau = 1,NumTau
       tau(i_tau) = 10.0**(minlogtau+dlogtau*real(i_tau-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    if (rank == 0) then
       write(logf,"(A,I3)") "Using BB up to frequency ", &
            freq_max(NumBndin1)
       write(logf,"(A,F10.2,A)") "  this is energy ", &
            freq_max(NumBndin1)/ev2fr," eV"
    endif

    ! In frequency band 1, fill in integrands and make tables

    ! Go through all the sub-bin in band 1
    do i_subband=1,NumBndin1  
       
       ! Fill frequency array
       call set_frequency_array(i_subband)
       
       ! Set frequency dependence of cross section
       call set_cross_section_freq_dependence(i_subband, &
            pl_index_cross_section_HI(i_subband),grey)
       
       ! Make photo integrands
       call fill_photo_integrands(i_subband)
       
       ! Assign values to the heating integrands
       if (.not.isothermal) then
          call fill_heating_integrands_HI
       endif
       
       ! Set integration weights
       call set_integration_weights(i_subband)
       
       ! Make photo tables
       call make_photo_tables(i_subband)
       
       ! Make heating tables
       if (.not.isothermal) then
          call make_heat_tables_HI(1)
       endif
       
    enddo
    
#ifdef MPILOG
    write(logf,*) 'Tables made, cleaning up'
#endif
    
    ! deallocate the useless photo integrand
    call deallocate_integrand_arrays

#ifdef MPILOG
    write(logf,*) 'Clean!'
#endif

    ! report a table value
    if (rank == 0) then
       write(logf,*) "bb_photo_thick_table: ",sum(bb_photo_thick_table(0,:))
       write(logf,*) "bb_photo_thin_table: ",sum(bb_photo_thin_table(0,:))
       write(logf,*) "bb_heat_thick_table: ",(bb_heat_thick_table(0,1))
       write(logf,*) "bb_heat_thin_table: ",(bb_heat_thin_table(0,1))
    endif
    
  end subroutine spec_integration
  
!---------------------------------------------------------------------------
  
  subroutine allocate_spec_integration_arrays
    
    call allocate_integrand_arrays
    call allocate_table_arrays

  end subroutine allocate_spec_integration_arrays

!---------------------------------------------------------------------------

  subroutine allocate_integrand_arrays

    ! Photoionization integrand as a function of frequency and tau
    allocate(bb_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(bb_photo_thin_integrand(0:NumFreq, 0:NumTau)) 
    allocate(pl_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(pl_photo_thin_integrand(0:NumFreq, 0:NumTau)) 

    ! Heating integrand as a function of frequency and tau
    if (.not.isothermal) then
       allocate(bb_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))  
       allocate(bb_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
    endif

  end subroutine allocate_integrand_arrays

!---------------------------------------------------------------------------

  subroutine deallocate_integrand_arrays
    
    deallocate(bb_photo_thick_integrand)
    deallocate(bb_photo_thin_integrand)
    deallocate(pl_photo_thick_integrand)
    deallocate(pl_photo_thin_integrand)
    
    ! deallocate the useless heating integrand
    if (.not.isothermal) then
       deallocate(bb_heat_thick_integrand_HI)
       deallocate(bb_heat_thin_integrand_HI)
       deallocate(pl_heat_thick_integrand_HI)
       deallocate(pl_heat_thin_integrand_HI)
    endif
    
  end subroutine deallocate_integrand_arrays

!---------------------------------------------------------------------------

  subroutine allocate_table_arrays

    ! Allocate the table arrays

    ! Photoionization table as a function of photo sub-bin and tau
    allocate(bb_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(bb_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    
    ! Heating table as a function of heating sub-bin and tau
    if (.not.isothermal) then
       allocate(bb_heat_thick_table(0:NumTau, 1:NumheatBin))
       allocate(bb_heat_thin_table(0:NumTau, 1:NumheatBin))
       allocate(pl_heat_thick_table(0:NumTau, 1:NumheatBin))
       allocate(pl_heat_thin_table(0:NumTau, 1:NumheatBin))
    endif
    
  end subroutine allocate_table_arrays

!---------------------------------------------------------------------------

  subroutine set_frequency_array(i_subband)
    
    integer,intent(in) :: i_subband
    
    integer :: i_freq
    
    do i_freq=0,NumFreq
       frequency(i_freq) = freq_min(i_subband)+ &
            delta_freq(i_subband)*real(i_freq)
    enddo
    
  end subroutine set_frequency_array

!---------------------------------------------------------------------------

  subroutine set_cross_section_freq_dependence(i_subband,pl_index,grey)

    integer,intent(in) :: i_subband
    real(kind=dp),intent(in) :: pl_index
    logical,intent(in) :: grey
    
    integer :: i_freq

    if (grey) then
       do i_freq=0,NumFreq
          cross_section_freq_dependence(i_freq) = 1.0
       enddo
    else
       do i_freq=0,NumFreq
          cross_section_freq_dependence(i_freq) = &
               (frequency(i_freq)/freq_min(i_subband))**(-pl_index)        
       enddo
    endif
    
  end subroutine set_cross_section_freq_dependence
    
!---------------------------------------------------------------------------

  subroutine fill_photo_integrands(i_subband)

    integer,intent(in) :: i_subband

    integer :: i_tau
    integer :: i_freq
    
    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition 
       do i_freq=0,NumFreq
          
          ! Assign values to the photo integrands
          if (tau(i_tau)*cross_section_freq_dependence(i_freq) < 700.0) then
             ! GM/130729 For these high frequencies this
             ! BB exponential term can overflow. Test for this.
             if (frequency(i_freq)*h_over_kT < 700.0) then  
                bb_photo_thick_integrand(i_freq,i_tau) = &
                     4.0_dp*pi*R_star2*two_pi_over_c_square* &
                     frequency(i_freq)*frequency(i_freq)* &
                     exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))/ &
                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
                bb_photo_thin_integrand(i_freq,i_tau) = &
                     4.0_dp*pi*R_star2*two_pi_over_c_square* &
                     frequency(i_freq)*frequency(i_freq)* &
                     cross_section_freq_dependence(i_freq)* &
                     exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))/ &
                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
             else
                bb_photo_thick_integrand(i_freq,i_tau) = 0.0
                bb_photo_thin_integrand(i_freq,i_tau) = 0.0
             endif
             pl_photo_thick_integrand(i_freq,i_tau) = &
                  pl_scaling*frequency(i_freq)**(-pl_index)* &
                  exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))
             pl_photo_thin_integrand(i_freq,i_tau) = &
                  pl_scaling*frequency(i_freq)**(-pl_index)* &
                  cross_section_freq_dependence(i_freq)* &
                  exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))
          else
             bb_photo_thick_integrand(i_freq,i_tau) = 0.0
             bb_photo_thin_integrand(i_freq,i_tau) = 0.0
             pl_photo_thick_integrand(i_freq,i_tau) = 0.0
             pl_photo_thin_integrand(i_freq,i_tau) = 0.0
          endif
          
       enddo
    enddo

  end subroutine fill_photo_integrands

!---------------------------------------------------------------------------

  subroutine fill_heating_integrands_HI

    integer :: i_tau
    integer :: i_freq
    
    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition
       do i_freq=0,NumFreq
          
          bb_heat_thick_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               bb_photo_thick_integrand(i_freq,i_tau)
          bb_heat_thin_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               bb_photo_thin_integrand(i_freq,i_tau)
          pl_heat_thick_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               pl_photo_thick_integrand(i_freq,i_tau)
          pl_heat_thin_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               pl_photo_thin_integrand(i_freq,i_tau)
          
       enddo
       
    enddo

  end subroutine fill_heating_integrands_HI

!---------------------------------------------------------------------------

  subroutine set_integration_weights(i_subband)

    ! Sets the weights for the call to vector_romberg
    integer,intent(in) :: i_subband

    vector_weight(:,:) = delta_freq(i_subband)  

  end subroutine set_integration_weights

!---------------------------------------------------------------------------

  subroutine make_photo_tables(i_subband)

    integer,intent(in) :: i_subband

    real(kind=dp), dimension(0:NumTau) :: answer

    ! Make photo tables
    call vector_romberg (bb_photo_thick_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_photo_thick_table(:,i_subband) = answer
    call vector_romberg (bb_photo_thin_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_photo_thin_table(:,i_subband) = answer
    call vector_romberg (pl_photo_thick_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_photo_thick_table(:,i_subband) = answer
    call vector_romberg (pl_photo_thin_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_photo_thin_table(:,i_subband) = answer  

  end subroutine make_photo_tables

!---------------------------------------------------------------------------

  subroutine make_heat_tables_HI(table_position)
    
    integer,intent(in) :: table_position

    real(kind=dp), dimension(0:NumTau) :: answer

    call vector_romberg (bb_heat_thick_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thick_table(:,table_position) = answer
    call vector_romberg (bb_heat_thin_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thin_table(:,table_position) = answer
    call vector_romberg (pl_heat_thick_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thick_table(:,table_position) = answer
    call vector_romberg (pl_heat_thin_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thin_table(:,table_position) = answer
    
  end subroutine make_heat_tables_HI
  
!---------------------------------------------------------------------------

end module radiation_tables
