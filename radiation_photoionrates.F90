!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_photoionrates
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: hplanck                     ! Planck constant
  use cgsphotoconstants, only: ion_freq_HI            ! HI ionization energy in frequency

  !use material, only: isothermal
  use c2ray_parameters, only: isothermal

  use radiation_sizes, only: NumFreqBnd
  use radiation_sizes, only: sigma_HI
  use radiation_sizes, only: NumTau, NumheatBin
  
  use radiation_tables, only: minlogtau, dlogtau
  use radiation_tables, only: bb_photo_thick_table, bb_photo_thin_table 
  use radiation_tables, only: pl_photo_thick_table, pl_photo_thin_table 
  use radiation_tables, only: bb_heat_thick_table, bb_heat_thin_table 
  use radiation_tables, only: pl_heat_thick_table, pl_heat_thin_table 
  use radiation_tables, only: bb_FreqBnd_UpperLimit, bb_FreqBnd_LowerLimit
  use radiation_tables, only: pl_FreqBnd_UpperLimit, pl_FreqBnd_LowerLimit

  implicit none

  ! photrates contains all the photo-ionization rates and heating rates
  type photrates    
     real(kind=dp) :: photo_cell_HI          ! HI photoionization rate of the cell    
     real(kind=dp) :: heat_cell_HI           ! HI heating rate of the cell       
     real(kind=dp) :: photo_in_HI            ! HI photoionization rate incoming to the cell    
     real(kind=dp) :: heat_in_HI             ! HI heating rate incoming to the cell
     real(kind=dp) :: photo_out_HI           ! HI photoionization rate outgoing from the cell
     real(kind=dp) :: heat_out_HI            ! HI heating rate outgoing from the cell
     real(kind=dp) :: heat                   ! Total heating rate of the cell
     real(kind=dp) :: photo_in               ! Total photoionization rate incoming to the cell
     real(kind=dp) :: photo_out               ! Total photoionization rate incoming to the cell
  end type photrates

  ! This definition allows adding two variables of type photrates using the 
  ! + sign.
  ! The function photrates_add is defined below.
  interface operator (+)
     module procedure photrates_add
  end interface operator (+)

  ! tablepos helps to locate correct position of the photoionization and heating tables
  type tablepos
    real(kind=dp), dimension(NumFreqBnd) :: tau            
    real(kind=dp), dimension(NumFreqBnd) :: odpos          
    real(kind=dp), dimension(NumFreqBnd) :: residual       
    integer, dimension(NumFreqBnd)       :: ipos           
    integer, dimension(NumFreqBnd)       :: ipos_p1        
  end type tablepos 

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  function photoion_rates (colum_in_HI,colum_out_HI, &
       vol,nsrc,i_state)

    use sourceprops, only: NormFlux,NormFluxPL
    !use cgsphotoconstants

    ! Function type
    type(photrates) :: photoion_rates

    ! Incoming and outgoing HI column density
    real(kind=dp), intent(in) :: colum_in_HI, colum_out_HI

    ! Volume of shell cell
    real(kind=dp), intent(in) :: vol

    ! Ionization state of cell
    real(kind=dp), intent(in) :: i_state

    ! Number of the source
    integer, intent(in) :: nsrc 

    integer :: i_subband
    real(kind=dp) :: colum_cell_HI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_in_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_out_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HI

    type(tablepos) :: tau_pos_in, tau_pos_out
    type(photrates) :: phi

    ! New source position, set local photo-ionization and heating rates 
    ! to zero. The structure phi is ultimately copied to the result of this function
    call set_photrates_to_zero (phi)

    ! Set the column densities (HI, HeI, HeII) of the current cell
    colum_cell_HI = colum_out_HI-colum_in_HI

    ! Calculate the optical depth (incoming, HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
       tau_in_all(i_subband) = colum_in_HI*sigma_HI(i_subband)
    enddo

    ! total tau_out (including HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
       tau_out_all(i_subband) = colum_out_HI*sigma_HI(i_subband)
    enddo

    ! find the table positions for the optical depth (ingoing and outgoing)
    tau_pos_in = set_tau_table_positions(tau_in_all)
    tau_pos_out = set_tau_table_positions(tau_out_all)

    ! Find the photo-ionization rates by looking up the values in
    ! the (appropriate) photo-ionization tables and add to the
    ! rates
    if (NormFlux(nsrc) > 0.0) &  
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFlux(nsrc),"B",vol)
    !if (colum_in_HI == 0.0) write(logf,*) "After photolookup: ", &
    !     phi%photo_cell_HI, phi%photo_cell_HeI, &
    !              phi%photo_cell_HeII, phi%heat
    if (allocated(NormFluxPL)) then
       if (NormFluxPL(nsrc) > 0.0) &  
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFluxPL(nsrc),"P",vol)
    endif
    ! Find the heating rates rates by looking up the values in
    ! the (appropriate) photo-ionization tables and using the
    ! secondary ionization. Add them to the rates.
    if (.not.isothermal) then
       
       ! The optical depths (HI, HeI, HeII) at current cell
       ! These are only needed in heat_lookuptable
       do i_subband=1,NumFreqBnd
          tau_cell_HI(i_subband) = colum_cell_HI*sigma_HI(i_subband)
       enddo
              
       !if (colum_in_HI == 0.0) then
       !   write(logf,*) "Before heatlookup: ", &
       !        phi%photo_cell_HI, phi%photo_cell_HeI, &
       !        phi%photo_cell_HeII, phi%heat
       !endif
       if (allocated(NormFlux)) then
          if (NormFlux(nsrc) > 0.0) & 
               phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
               tau_in_all,tau_out_all, &
               tau_cell_HI,NormFlux(nsrc),"B", &
               vol,i_state)
       endif
       !if (colum_in_HI == 0.0) write(logf,*) "After heatlookup: ", &
       !     phi%photo_cell_HI, phi%photo_cell_HeI, &
       !     phi%photo_cell_HeII, phi%heat

        if (allocated(NormFluxPL)) then
           if (NormFluxPL(nsrc) > 0.0) &  
                phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
                tau_in_all,tau_out_all, &
                tau_cell_HI,NormFluxPL(nsrc),"P", &
                vol,i_state)
        endif
    endif

    ! Assign result of function
    photoion_rates = phi

  end function photoion_rates
 
!---------------------------------------------------------------------------

  ! Calculates the table position data for an optical depth tau
  function set_tau_table_positions (tau)

    real(kind=dp), dimension(1:NumFreqBnd),intent(in) :: tau
    type(tablepos) :: set_tau_table_positions
    type(tablepos) :: tau_position

    integer :: i_subband
    
    
    ! fill the table positions structure for the optical depth tau
    do i_subband=1,NumFreqBnd  
       tau_position%tau(i_subband) = log10(max(1.0e-20_dp,tau(i_subband)))
       tau_position%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
            (tau_position%tau(i_subband)-minlogtau)/dlogtau))
       tau_position%ipos(i_subband) = int(tau_position%odpos(i_subband))
       tau_position%residual(i_subband) = tau_position%odpos(i_subband)- &
            real(tau_position%ipos(i_subband),dp)
       tau_position%ipos_p1(i_subband) = min(NumTau, &
            tau_position%ipos(i_subband)+1)
    enddo
    
    ! Set the return value
    set_tau_table_positions=tau_position

  end function set_tau_table_positions

!---------------------------------------------------------------------------

  function read_table(table,tablesize,table_position,i_subband,i_subband2)
    
    integer,intent(in) :: tablesize
    real(kind=dp), pointer, dimension(:,:),intent(in) :: table
    !real(kind=dp),dimension(0:NumTau, 1:tablesize),intent(in) :: table
    type(tablepos),intent(in) :: table_position
    integer,intent(in) :: i_subband
    integer,intent(in) :: i_subband2
    
    real(kind=dp) :: read_table
    
    read_table = table(table_position%ipos(i_subband),i_subband2)+ &
         ( table(table_position%ipos_p1(i_subband),i_subband2)- &
         table(table_position%ipos(i_subband),i_subband2) ) * &
         table_position%residual(i_subband)
    
  end function read_table

  !---------------------------------------------------------------------------

  ! find out the correct position in the photo and heating tables
  function photo_lookuptable(tau_pos_in,tau_pos_out, &
       tau_in_all,tau_out_all, &
       NFlux,table_type, &
       vol)
    
    !use cgsphotoconstants
    
    ! Function type
    type(photrates) :: photo_lookuptable
    
    ! Optical depth below which we should use the optically thin tables
    real(kind=dp),parameter :: tau_photo_limit = 1.0e-7 

    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol
    real(kind=dp), dimension(NumFreqBnd), intent(in) :: tau_in_all, tau_out_all
    character,intent(in) :: table_type
    
    integer ::  i_subband
    real(kind=dp) :: phi_photo_in_all, phi_photo_out_all, phi_photo_all
    real(kind=dp), pointer, dimension(:,:) :: photo_thick_table, photo_thin_table
    integer :: Minimum_FreqBnd
    integer :: Maximum_FreqBnd
        
    ! New source. Set all the rates to zero to initialize them.
    call set_photrates_to_zero (photo_lookuptable)

    ! pointers point to the correct tables to use, BB or PL source
    ! Set the maximum frequency band to consider (and limit the
    ! loop over the subbands below)
    if (table_type == "B") then 
       photo_thick_table => bb_photo_thick_table
       photo_thin_table => bb_photo_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=1
    elseif (table_type == "P") then
       photo_thick_table => pl_photo_thick_table
       photo_thin_table => pl_photo_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=1
    endif
    
    ! loop through the relevant frequency bands
    do i_subband=Minimum_FreqBnd, Maximum_FreqBnd
       
       ! Incoming total photoionization rate
       phi_photo_in_all = NFlux* &
            read_table(photo_thick_table,NumFreqBnd,tau_pos_in, &
            i_subband,i_subband)
       photo_lookuptable%photo_in = &
            photo_lookuptable%photo_in + phi_photo_in_all

       ! Total cell photo-ionization rate, calculated differently
       ! for optically thick and thin cells
       if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
            tau_photo_limit) then
          
          ! When current cell is optically thick
          phi_photo_out_all = NFlux* &
               read_table(photo_thick_table,NumFreqBnd, &
               tau_pos_out,i_subband,i_subband)
          phi_photo_all = phi_photo_in_all-phi_photo_out_all
          
       else
          
          ! When current cell is optically thin
          phi_photo_all = NFlux* &
               (tau_out_all(i_subband)-tau_in_all(i_subband))* &
               read_table(photo_thin_table,NumFreqBnd, &
               tau_pos_in,i_subband,i_subband)
          phi_photo_out_all = phi_photo_in_all-phi_photo_all

       endif

       ! Collect all outgoing photons
       photo_lookuptable%photo_out = &
            photo_lookuptable%photo_out + phi_photo_out_all

       ! Assign to the HI photo-ionization rate
       photo_lookuptable%photo_cell_HI = photo_lookuptable%photo_cell_HI + &
            phi_photo_all/vol
       
    enddo
    
  end function photo_lookuptable
  
 !---------------------------------------------------------------------------
 
 ! find out the correct position in the photo and heating tables.
 ! it updates phi
  function heat_lookuptable (tau_pos_in,tau_pos_out, &
       tau_in_all,tau_out_all, &
       tau_cell_HI, &
       NFlux,table_type, &
       vol,i_state)

    ! Function type
    type(photrates) :: heat_lookuptable

    ! Optical depth below which we should use the optically thin tables
    real(kind=dp),parameter :: tau_heat_limit = 1.0e-4

    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol, i_state
    real(kind=dp), dimension(1:NumFreqBnd), intent(in) :: tau_in_all, &
         tau_out_all
    character,intent(in) :: table_type
    real(kind=dp), dimension(1:NumFreqBnd),intent(in) :: tau_cell_HI

    integer ::  i_subband, i
    real(kind=dp) :: phi_heat_HI
    real(kind=dp) :: phi_heat_in_HI
    real(kind=dp) :: phi_heat_out_HI
    real(kind=dp) :: f_heat
    real(kind=dp) :: df_heat
    real(kind=dp), pointer, dimension(:,:) :: heat_thick_table, heat_thin_table
    integer :: Minimum_FreqBnd
    integer :: Maximum_FreqBnd
    ! Related to secondary ionizations

    ! New source. Set all the rates to zero to initialize them.
    call set_photrates_to_zero (heat_lookuptable)

    ! pointers point to the correct tables to use, BB or PL source
    ! Set the maximum frequency band to consider (and limit the
    ! loop over the subbands below)
    if (table_type == "B") then 
       heat_thick_table => bb_heat_thick_table
       heat_thin_table => bb_heat_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=1
    elseif (table_type == "P") then
       heat_thick_table => pl_heat_thick_table
       heat_thin_table => pl_heat_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=1
    endif

    ! Current cell individual heating rates of HI, HeI, HeII
    ! loop through the frequency bands
    do i_subband=Minimum_FreqBnd,Maximum_FreqBnd
       
       ! For every subband these will contain the heating due to
       ! the different species by photons in that subband.
       ! The sum over all subbands is collected in f_heat.
       ! All variables starting will phi_ are rates for one subband
       ! and local variables. The heat_lookuptable% is final answer for
       ! all subbands and sources.
       phi_heat_HI = 0.0_dp

       ! Incoming current cell HI heating rate at band 1
       phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
            tau_pos_in,i_subband,i_subband)
          
       ! When current cell is HI optically thick (in total tau)
       if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
            tau_heat_limit) then
          phi_heat_out_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_out,i_subband,i_subband)
          phi_heat_HI = (phi_heat_in_HI-phi_heat_out_HI)/vol
          
          ! When current cell is HI optically thin
       else
          phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
               read_table(heat_thin_table,NumheatBin, &
               tau_pos_in,i_subband,i_subband)
          !phi_heat_out_HI = phi_heat_in_HI-phi_heat_HI
          phi_heat_HI = phi_heat_HI/vol
       endif
          
       ! Save the heating in f_heat variable
       df_heat = phi_heat_HI
       
       ! Add the subband contribution to the total rates
       f_heat = f_heat + df_heat

    enddo
       
    ! Total heating rate on current cell
    ! Needs to be cumulative because one cell may contain different
    ! types of sources.
    heat_lookuptable%heat = f_heat 
    
  end function heat_lookuptable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function photrates_add (rate1, rate2)
    
    type(photrates) :: photrates_add
    type(photrates),intent(in) :: rate1, rate2
    
    photrates_add%photo_cell_HI = rate1%photo_cell_HI + rate2%photo_cell_HI
    photrates_add%heat_cell_HI = rate1%heat_cell_HI + rate2%heat_cell_HI
    photrates_add%photo_in_HI = rate1%photo_in_HI + rate2%photo_in_HI
    photrates_add%heat_in_HI = rate1%heat_in_HI + rate2%heat_in_HI
    photrates_add%photo_out_HI = rate1%photo_out_HI + rate2%photo_out_HI
    photrates_add%heat_out_HI = rate1%heat_out_HI + rate2%heat_out_HI
    photrates_add%heat = rate1%heat + rate2%heat
    photrates_add%photo_in = rate1%photo_in + rate2%photo_in
    photrates_add%photo_out = rate1%photo_out + rate2%photo_out
    
  end function photrates_add
  
  subroutine set_photrates_to_zero (rate1)
    
    type(photrates),intent(out) :: rate1
    
    rate1%photo_cell_HI = 0.0
    rate1%heat_cell_HI = 0.0
    rate1%photo_in_HI = 0.0
    rate1%heat_in_HI = 0.0
    rate1%photo_out_HI = 0.0
    rate1%heat_out_HI = 0.0
    rate1%heat = 0.0
    rate1%photo_in = 0.0
    rate1%photo_out = 0.0
    
  end subroutine set_photrates_to_zero
  
end module radiation_photoionrates
