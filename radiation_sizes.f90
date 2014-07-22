!     This module sets the sizes of radiation quantities

module radiation_sizes

  use precision, only: dp
  use cgsphotoconstants, only: ion_freq_HI => frth0,&             ! HI ionization energy in frequency
                               sigma_HI_at_ion_freq => sigh    ! HI cross section at its ionzing frequency
  !use material, only: isothermal
  use c2ray_parameters, only: isothermal

  
  implicit none

  integer,parameter :: NumFreq = 512      ! Number of integration points in each of the frequency bins
  integer,parameter :: NumTau = 2000      ! Number of table points for the optical depth
  integer,parameter :: NumBndin1 = 1      ! Number of frequency sub-bins in interval 1 
  integer,parameter :: NumFreqBnd=NumBndin1       ! Total number of frequency bins
  integer,parameter :: NumheatBin=NumBndin1   ! Total number of heating bins 

  ! Optical depths at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the outside.
  real(kind=dp) :: boundary_tauHI = 0.0

  real(kind=dp), dimension(:), allocatable :: delta_freq      ! Frequency width of integration 
  real(kind=dp), dimension(:), allocatable :: freq_max        ! Maximum freqeucny of integration 
  real(kind=dp), dimension(:), allocatable :: freq_min        ! Minimum freqeucny of integration

  ! Power law fit parameter for frequency range 1:3
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HI    ! Power law index of cross section of HI

  ! Cross section of atoms
  real(kind=dp), dimension(:), allocatable :: sigma_HI       ! Cross section of HI

contains

  ! set up some scaling factors arrays
  subroutine setup_scalingfactors (freq_max_src)

    ! Contains the maximum frequency, determined by source input
    real(kind=dp),intent(in) :: freq_max_src

    integer :: i_subband

    ! Allocate size of arrays
    allocate(delta_freq(1:NumFreqBnd))  
    allocate(freq_max(1:NumFreqBnd))  
    allocate(freq_min(1:NumFreqBnd))  
    allocate(pl_index_cross_section_HI(1:NumFreqBnd))
    allocate(sigma_HI(1:NumFreqBnd))

    ! Assignment of maximum frequency in the sub-bin partition.
    select case (NumBndin1)

    case (1)

       freq_max(NumBndin1) = freq_max_src

    end select

    ! Assignment of minimum frequency in the sub-bin partition.

    freq_min(NumBndin1) = ion_freq_HI

    ! calculate the width of frequency sub-bin
    do i_subband=1,NumFreqBnd 
       delta_freq(i_subband) = (freq_max(i_subband)-freq_min(i_subband))/real(NumFreq)
    enddo

    ! Assign value sigma of HI, HeI and HeII at different frequencies 
    select case (NumBndin1)

    case (1)

       sigma_HI(1) = sigma_HI_at_ion_freq

    end select

    ! Assign power-law index of HI, HeI and HeII at different frequencies (about absorption)
    select case (NumBndin1)

    case (1)

       pl_index_cross_section_HI(1) = 2.761_dp

    end select

  end subroutine setup_scalingfactors

end module radiation_sizes
