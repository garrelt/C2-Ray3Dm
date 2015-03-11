! This module contains routines having to do with the calculation of
! the thermal evolution of a single point/cell. 

module thermalevolution

  use precision, only: dp
  use c2ray_parameters, only: minitemp, relative_denergy, cosmological
  use radiative_cooling, only: coolin
  use tped, only: temper2pressr, pressr2temper, electrondens
  !use cgsconstants
  use atomic, only: gamma1
  use cosmology, only: cosmo_cool
  !use radiation, only: photrates
  use radiation_photoionrates, only: photrates
  use ionfractions_module, only: ionstates

  implicit none

contains

  ! calculates the thermal evolution of one grid point
  subroutine thermal (dt,end_temper,avg_temper,ndens_electron,ndens_atom,ion,phi)!,pos)

    ! The time step
    real(kind=dp), intent(in) :: dt
    ! end time temperature of the cell
    real(kind=dp), intent(inout) :: end_temper
    ! average temperature of the cell
    real(kind=dp), intent(out) :: avg_temper
    ! Electron density of the cell
    real(kind=dp), intent(in) :: ndens_electron
    ! Number density of atoms of the cell
    real(kind=dp), intent(in) :: ndens_atom
    ! Photoionization rate and heating rate
    type(photrates), intent(in) :: phi
    ! Ionized fraction of the cell
    type(ionstates), intent(in) :: ion
    ! mesh position of cell 
    !integer, intent(in) :: pos
 
    ! initial temperature
    real(kind=dp) :: initial_temp
    ! timestep taken to solve the ODE
    real(kind=dp) :: dt_ODE
    ! timestep related to thermal timescale
    real(kind=dp) :: dt_thermal
    ! record the time elapsed
    real(kind=dp) :: cumulative_time
    ! internal energy of the cell
    real(kind=dp) :: internal_energy
    ! thermal timescale, used to calculate the thermal timestep
    real(kind=dp) :: thermal_timescale
    ! heating rate
    real(kind=dp) :: heating
    ! cooling rate
    real(kind=dp) :: cooling
    ! difference of heating and cooling rate
    real(kind=dp) :: thermal_rate
    ! cosmological cooling rate
    real(kind=dp) :: cosmo_cool_rate
    ! Counter of number of thermal timesteps taken
    integer :: i_heating

    ! heating rate
    heating = phi%heat

    ! Find initial internal energy
    internal_energy = temper2pressr(end_temper,ndens_atom, &
         electrondens(ndens_atom,ion%h_old))/(gamma1)
    !internal_energy = temper2pressr(end_temper,ndens_atom,electrondens(ndens_atom,ion%begin_HII, &
    !                                ion%begin_HeII,ion%begin_HeIII))/(gamma1)

    ! Set the cosmological cooling rate
    if (cosmological) then
       ! Disabled for testing
       cosmo_cool_rate=cosmo_cool(internal_energy)
    else
       cosmo_cool_rate=0.0
    endif

    ! Thermal process is only done if the temperature of the cell 
    ! is larger than the minimum temperature requirement
    if (end_temper.gt.minitemp) then

       ! stores the time elapsed is done
       cumulative_time = 0.0 
   
       ! initialize the counter
       i_heating = 0

       ! initialize time averaged temperature
       avg_temper = 0.0 

       ! initial temperature
       initial_temp = end_temper  

       ! thermal process begins
       do

          ! update counter              
          i_heating = i_heating+1 
         
          ! update cooling rate from cooling tables
          cooling = coolin(ndens_atom,ndens_electron,ion%h_av, &
               end_temper)+cosmo_cool_rate
          !cooling = coolin(ndens_atom,ndens_electron,ion%avg_HI,ion%avg_HII,ion%avg_HeI,ion%avg_HeII,&
          !                 ion%avg_HeIII, end_temper)+cosmo_cool_rate

          ! Find total energy change rate
          thermal_rate = max(1d-50,abs(cooling-heating))

          ! Calculate thermal time scale
          thermal_timescale = internal_energy/abs(thermal_rate)

          ! Calculate time step needed to limit energy change
          ! to a fraction relative_denergy
          dt_thermal = relative_denergy*thermal_timescale

          ! Time step to large, change it to dt_thermal
          ! Make sure we do not integrate for longer than the
          ! total time step
          dt_ODE = min(dt_thermal,dt-cumulative_time)

          ! Find new internal energy density
          internal_energy = internal_energy+dt_ODE*(heating-cooling)

          ! Update avg_temper sum (first part of dt_thermal sub time step)
          avg_temper = avg_temper+0.5*end_temper*dt_ODE

          ! Find new temperature from the internal energy density
          end_temper = pressr2temper(internal_energy*gamma1,ndens_atom, &
               electrondens(ndens_atom,ion%h_av))
          !end_temper = pressr2temper(internal_energy*gamma1,ndens_atom,electrondens(ndens_atom,&
          !                           ion%avg_HII,ion%avg_HeII,ion%avg_HeIII))

          ! Update avg_temper sum (second part of dt_thermal sub time step)
          avg_temper = avg_temper+0.5*end_temper*dt_ODE
                    
          ! Take measures if temperature drops below minitemp
          if (end_temper.lt.minitemp) then
             internal_energy = temper2pressr(minitemp,ndens_atom, &
                  electrondens(ndens_atom,ion%h_av))
             !internal_energy = temper2pressr(minitemp,ndens_atom,electrondens(ndens_atom,ion%avg_HII,&
             !                                ion%avg_HeII,ion%avg_HeIII))
             end_temper = minitemp
          endif
                    
          ! Update fractional cumulative_time
          cumulative_time = cumulative_time+dt_ODE
  
          ! Exit if we reach dt
          if (cumulative_time.ge.dt.or.abs(cumulative_time-dt).lt.1e-6*dt) exit

          ! In case we spend too much time here, we exit
          if (i_heating.gt.10000) exit
        	     	
       enddo
              
       ! Calculate the averaged temperature
       if (dt.gt.0.0) then
          avg_temper = avg_temper/dt
       else
          avg_temper = initial_temp
       endif
       
       ! Calculate the final temperature 
       end_temper = pressr2temper(internal_energy*gamma1,ndens_atom, &
            electrondens(ndens_atom,ion%h))
       !end_temper = pressr2temper(internal_energy*gamma1,ndens_atom,electrondens(ndens_atom,&
       !                           ion%end_HII,ion%end_HeII,ion%end_HeIII))
       
    endif
    
  end subroutine thermal
  
end module thermalevolution
