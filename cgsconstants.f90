!>
!! \brief This module contains definitions of physical constants 
!! and conversion factors, as well as routines to set recombination
!! and collisional ionization coefficients.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema, Martina M. Friedrich
!!
!! \b Date: 2015-06-30 (previously undated)
!!
!! \b Version: cgs units
!!

module cgsconstants

  use precision, only: dp
  use mathconstants, only: pi

  implicit none

  ! A collection of physical constants and conversion factors
  ! Units: cgs

  !> proton mass
  real(kind=dp), parameter :: m_p=1.672661e-24_dp
  !> speed of light
  real(kind=dp), parameter :: c=2.997925e+10_dp
  !> Planck constant
  real(kind=dp), parameter :: hplanck=6.6260755e-27_dp
  !> Stefan-Boltzmann constant
  real(kind=dp), parameter :: sigma_SB=5.670e-5_dp
  !> Boltzmann constant
  real(kind=dp), parameter :: k_B=1.381e-16_dp
  !> Gravitational constant
  real(kind=dp), parameter :: G_grav=6.6732e-8_dp

  !> ev2k   - conversion factor between evs and kelvins
  real(kind=dp),parameter  :: ev2k=1.0/8.617e-05
  !> ev2erg  - conversion factor between evs and ergs
  real(kind=dp),parameter  :: ev2erg=1.602e-12
  !> ev2j   - conversion factor between ergs and Joules
  real(kind=dp),parameter  :: erg2j=1e-7
  
  ! The following are scaled to frequency scaling

  !> Frequency scaling factor,
  !! this scaling parameter is independent of any main program scaling
  !! (see scaling.f90), and may only be used in the radiation physics 
  !! subroutines (currently switched off)
  real(kind=dp),parameter  :: sclfre=1.0e15
  !> conversion between evs and frequency
  real(kind=dp),parameter  :: ev2fr=0.241838e15!/sclfre      

  ! h/k, Planck/Boltzm
  ! Check this number, seems scaled
  !real(kind=dp),parameter  :: hoverk=47979.72484
  !> Planck constant scaled
  real(kind=dp),parameter  :: hscl=hplanck!*sclfre 
  !> tpic2  - 2*pi/c^2 times scaling factors needed for the integral cores
  real(kind=dp),parameter :: two_pi_over_c_square=2.0*pi/(c*c)

  !> Hydrogen recombination parameter (power law index)
  real(kind=dp),parameter :: albpow=-0.7_dp  !in the old code -0.79
  !> Hydrogen recombination parameter (value at 10^4 K)
  real(kind=dp),parameter :: bh00=2.59e-13_dp ! OTS value, alpha_B
  !> Helium0 recombination parameter (power law index)
  real(kind=dp), parameter :: alcpow=-0.672_dp
  !> Helium0 recombination parameter (value at 10^4 K)
  real(kind=dp), parameter :: bhe00=4.26e-13_dp !alpha_b+alpha_1
  !> Helium1 recombination parameter (value at 10^4 K)
  real(kind=dp), parameter :: bhe10=1.53e-12_dp !different in the book! 
  !here it was 5.91e-12 I replace with book value of 2*1.1e-12

  !> Hydrogen ionization energy (in eV)
  real(kind=dp), parameter :: eth0=13.598
  !> Hydrogen ionization energy (in erg)
  real(kind=dp),parameter :: hionen=eth0*ev2erg
  !> Hydrogen ionization energy expressed in K
  real(kind=dp),parameter :: temph0=eth0*ev2k
  !> Hydrogen collisional ionization parameter 1
  real(kind=dp),parameter :: xih0=1.0
  !> Hydrogen collisional ionization parameter 2
  real(kind=dp),parameter :: fh0=0.83
  !> Hydrogen collisional ionization parameter
  real(kind=dp),parameter :: colh0=1.3e-8*fh0*xih0/(eth0*eth0)
  ! critical electron density at which colisions become important for y-fraction ionization (Osterbrock)
  real(kind=dp),parameter :: n_el_crit=4.0e3

  !> Helium ionization energy (in eV)
  real(kind=dp), dimension(0:1), parameter :: ethe=(/24.587,54.416/)
  !> Helium ionization energy (in erg)
  real(kind=dp), dimension(0:1), parameter :: heionen=(/ethe(0)*ev2erg,ethe(1)*ev2erg/)
  !> Helium ionization energy expressed in K
  real(kind=dp), dimension(0:1), parameter :: temphe=(/ethe(0)*ev2k,ethe(1)*ev2k/)
  !> Helium collisional ionization parameter 1
  real(kind=dp),dimension(0:1),parameter :: xihe=(/2.0,1.0/)
  !> Helium collisional ionization parameter 2
  real(kind=dp),dimension(0:1),parameter :: fhe=(/0.63,1.30/)
  !> Helium collisional ionization parameter
  real(kind=dp),dimension(0:1),parameter :: colhe= &
       (/1.3e-8*fhe(0)*xihe(0)/(ethe(0)*ethe(0)), &
       1.3e-8*fhe(1)*xihe(1)/(ethe(1)*ethe(1))/)

  !> Hydrogen 0 A recombination parameter
  real(kind=dp) :: arech0
  !> Hydrogen 0 B recombination parameter
  real(kind=dp) :: brech0

  !> Helium   0 A recombination parameter
  real(kind=dp) :: areche0
  !> Helium   0 B recombination parameter
  real(kind=dp) :: breche0
  !> Helium   0 1 recombination parameter
  real(kind=dp) :: oreche0
  
  !> Helium   + A recombination parameter
  real(kind=dp) :: areche1
  !> Helium   + B recombination parameter
  real(kind=dp) :: breche1
  !> Helium   + 2 recombination parameter
  real(kind=dp) :: treche1
  
  
  !> H0 collisional ionization parameter at T=temp0
  real(kind=dp) :: colli_HI
  !> He0 collisional ionization parameter at T=temp0
  real(kind=dp) :: colli_HeI
  !> He1 collisional ionization parameter at T=temp0
  real(kind=dp) :: colli_HeII
  !> Fraction fo He++ -> He+ recombination photons that goes into 2photon decay
  real(kind=dp) :: v

contains

  !> This subroutine contains all temperature dependent recombination and 
  !! collisional ionization rates. It should be called every time after the
  !! temperature is recalculated.
  subroutine ini_rec_colion_factors (temperature)
    
    real(kind=dp),intent(in):: temperature

    call ini_hydrogen_recombination (temperature)
    
    call ini_helium0_recombination (temperature)
    
    call ini_helium1_recombination (temperature)

    call ini_hydrogen_helium_collisional_ionization (temperature)

  end subroutine ini_rec_colion_factors

  ! ---------------------------------------------------------------------------

  subroutine  ini_hydrogen_recombination (temperature)
    
    ! Recombination from H^+ to H^0 coefficients.

    ! MMF: Implement the new fit as mentioned in email from 2008-06-08
    ! 
    ! fit1=bh00*1.4*(T/1e4)^albpow*0.9 
    ! fit2=bh00*5*(T/1e4)^albpow*1.95
    ! brech0=1/[(1/fit1)+(1/fit2)]
    ! 
    real(kind=dp),intent(in):: temperature
    !real(kind=dp) :: fit1,fit2
    real(kind=dp) :: lambda

    ! Hydrogen A & B recombination rates, fits from Hui & Gnedin (1997)
    lambda =  2.0_dp*(temph0/temperature)
    arech0=1.269e-13*lambda**1.503_dp/(1.0_dp+(lambda/0.522)**0.470)**1.923
    brech0=2.753e-14*lambda**1.500_dp/(1.0_dp+(lambda/2.740)**0.407)**2.242

  end subroutine ini_hydrogen_recombination

  !----------------------------------------------------------------------------

  subroutine ini_helium0_recombination (temperature)

    ! Recombination from He^+ to He^0 coefficients.
    ! MMF: See email from 2008-08-11
    ! For now, use old fit.

    real(kind=dp),intent(in):: temperature
    real(kind=dp) :: lambda
    real(kind=dp) :: dielectronic

    ! He0 recombination rates including the dielectronic recombination
    if (temperature < 9.e3_dp) then   ! was 1.5e4
       lambda=2.0_dp*(temph0/temperature)
       areche0=1.269e-13_dp*lambda**1.503_dp/(1.0_dp+(lambda/0.522)**0.470)**1.923
       breche0=2.753e-14_dp*lambda**1.500_dp/(1.0_dp+(lambda/2.740)**0.407)**2.242
    else
       lambda=   2.0_dp*(temphe(0)/temperature)
       dielectronic = 1.9e-3_dp*temperature**(-1.5_dp)*exp(-4.7e5_dp/temperature)*(1.0_dp+0.3_dp*exp(-9.4e4_dp/temperature))
       areche0= 3.000e-14_dp*lambda**0.654_dp + dielectronic
       breche0= 1.260e-14_dp*lambda**0.750_dp + dielectronic
    endif
    oreche0=areche0-breche0
    ! above 1.5e4 I follow the fit from Hui & Gnedin (1997) for He0 for radiation recombination, 
    ! this fits the data from Hummer&Storey (1998) good.
    ! below that, the data from Hummer&Storey is fitted by the fit for H from Hui%Gnedin1997
    ! comparing the He0 and H data from Hummer&Storey (1998) and Hummer (1994) shows, that 
    ! the two never diviate more than 4% above T=1.5e4. 

    ! without dielectronic recombinations
    !areche0 = bhe00*(temperature/1.0e4_dp)**alcpow
    !oreche0 =1.54e-13_dp*(temperature/1.0e4_dp)**(-0.5_dp)  ! values from Osterbrock interpolated
    !breche0 = areche0-oreche0                         !B recombination
    ! NOT USED ANYMORE

  end subroutine ini_helium0_recombination

  !----------------------------------------------------------------------------

  subroutine ini_helium1_recombination (temperature)

    ! Recombination from He^2+ to He^+ coefficients.
    ! MMF: Implement the new fit as mentioned in email from 2008-06-08

    real(kind=dp),intent(in):: temperature
    real(kind=dp) :: lambda

    ! fits to He1 recombination rates to the data from Hummer (1994)
    ! using the fits from Hui & Gnedin (1997), they are excellent fits 
    ! to the table data from Hummer (1994)
    ! The same is true for the A rates
    lambda=2.0_dp*(temphe(1)/temperature)
    breche1=5.5060e-14_dp*lambda**1.5_dp/(1.0_dp+(lambda/2.740_dp)**0.407_dp)**2.242_dp
    areche1=2.538e-13*lambda**1.503_dp/(1.0_dp+(lambda/0.522_dp)**0.470_dp)**1.923_dp

    ! For the B-C values, that is recombination to n=2, there is no data in
    ! Hummer (1994), neither a fit in Hui&Gnedin (1997)
    treche1 = 3.4e-13_dp*(temperature/1.0e4_dp)**(-0.6_dp)  ! extrapolate Osterbrok B value minus C value, p 38

    ! fraction of He++ recombination photons that go into two-photon decay.
    v=0.285_dp*(temperature/1.0e4_dp)**0.119_dp ! fit to Hummer and Seaton, 1964 table five

  end subroutine ini_helium1_recombination

  !----------------------------------------------------------------------------

  subroutine ini_hydrogen_helium_collisional_ionization (temperature)

    ! Collisional ionization coefficients for H^0 to H^+ (temperature dependent)

    real(kind=dp),intent(in):: temperature
    real(kind=dp) :: sqrtt0

    ! Fits from Cox (1970)
    sqrtt0 =sqrt(temperature)
    colli_HI =colh0*sqrtt0*exp(-temph0/temperature)
    colli_HeI=colhe(0)*sqrtt0*exp(-temphe(0)/temperature)
    colli_HeII=colhe(1)*sqrtt0*exp(-temphe(1)/temperature)

    ! I might consider using the Hui&Gnedin fits instead. However, the difference is not so big
    !sqrtt0 = temperature**(-1.5_dp)*exp(-temph0/temperature)
    !acolh0 =21.11_dp*sqrtt0*exp(-temph0/temperature)* &
    !(2.0_dp*temph0/temperature)**-1.089_dp/(1.0_dp+((2.0_dp*temph0/temperature)/0.354_dp)**0.874)**1.101_dp
    !acolhe0=32.38_dp*sqrtt0*exp(-temphe(0)/temperature)* &
    ! (2.0_dp*temphe(0)/temperature)**-1.146_dp/(1.0_dp+((2.0_dp*temphe(0)/temperature)/0.416_dp)**0.987)**1.056;
    !acolhe1=19.95_dp*sqrtt0*exp(-temphe(1)/temperature)* &
    !(2.0_dp*temphe(1)/temperature)**-1.089_dp/(1.0_dp+((2.0_dp*temphe(1)/temperature)/0.553_dp)**0.735)**1.275;

  end subroutine ini_hydrogen_helium_collisional_ionization

  !----------------------------------------------------------------------------

  subroutine ini_recombinations_collisional_ionizations_constant

    ! SET THEM HERE EXPLIXITLY FOR T=10 000 K TO THE ONES FROM GABRIEL
    brech0  = 2.59182e-13_dp     ! H0 b-recombination coefficient
    breche0 = 2.61613e-13_dp     ! He0 b-recombination coefficient
    breche1 = 1.54528e-12_dp     ! He1 b-recombination coefficient

    areche0 = 4.22471e-13_dp     ! He0 a-recombination coefficient
    areche1 = 2.22561e-12_dp     ! He1 a-recombination coefficient
    arech0  = 4.29695e-13_dp     ! H0 a-recombination coefficient

    oreche0 =areche0-breche0     ! he0 1 recombination coefficient
    treche1 = 3.46e-13_dp        !Osterbrok B value minus C value, p 38    
    v=0.285_dp

    colli_HI  = 8.96396e-16_dp  ! H collisional ionization coefficient
    colli_HeI = 7.46415e-22_dp  ! He0 collisional ionization coefficient
    colli_HeII = 2.28059e-37_dp  ! He1 collisional ionization coefficient

  end subroutine ini_recombinations_collisional_ionizations_constant

end module cgsconstants


