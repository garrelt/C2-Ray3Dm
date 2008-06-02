module cosmology

  ! This file contains routines having to do with the cosmological 
  ! evolution
  
  ! - cosmology_init: initializes cosmological time and sets lengths
  !            and volumes from comoving to proper scaling.
  ! - time2zred: convert time to redshift
  ! - zred2time: convert redshift to time
  ! - redshift_evol: calculate redshift from time, and zfactor between the
  !    current and previous time
  ! - cosmo_evol: cosmological evolution of space, density
  ! - cosmo_cool: cosmological adiabatic cooling rate
  ! - compton_cool: Compton cooling wrt the CMB.

  use precision, only: dp
  use cosmology_parameters

  implicit none

  real(kind=dp) :: zred_t0 ! initial redshift
  real(kind=dp) :: t0      ! time of initial redshift
  real(kind=dp) :: zred    ! current redshift
  real(kind=dp),private :: zfactor ! scaling factor between two redshifts

contains
  ! =======================================================================

  subroutine cosmology_init (zred0,time)
    
    real(kind=dp),intent(in) :: zred0
    real(kind=dp),intent(in) :: time
    
    ! Cosmological time corresponding to (initial) redshift zred0
    ! NOTE: Good only for high-z!!!
    t0 = 2.*(1.+zred0)**(-1.5)/(3.*H0*sqrt(Omega0))
    
    ! Initialize redshift
    zred_t0=zred0 ! keep initial redshift
    zred=0.0 ! needs to be zero, so comoving will be changed to proper
    
    ! Initialize lengths and volumes to proper units
    call redshift_evol(time)
    call cosmo_evol( )
    
  end subroutine cosmology_init

  ! =======================================================================

  function time2zred (time)

    ! Calculates the cosmological redshift for a given time

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 21-May-2005)
    
    ! Version: f90

    ! History: - 20-Aug-2006, conversion to f90

    real(kind=dp) :: time2zred
    real(kind=dp),intent(in) :: time

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    time2zred = -1+(1.+zred_t0)*(t0/(t0+time))**(2./3.)

    return
  end function time2zred

  ! =======================================================================

  function zred2time (zred1)

    ! Calculates the time for a given cosmological redshift

    ! Author: Garrelt Mellema

    ! Date: 30-Sep-2006
    
    ! Version: f90

    ! History: 

    real(kind=dp) :: zred2time
    real(kind=dp),intent(in) :: zred1

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    zred2time = t0*( ((1.0+zred_t0)/(1.0+zred1))**1.5 - 1.0 )

    return
  end function zred2time

  ! =======================================================================

  subroutine redshift_evol (time)

    ! Calculates the cosmological redshift from time
    ! and the scale factor zfactor for use in cosmo_evol
    
    ! Author: Garrelt Mellema
    ! Date: 10-Mar-2006
    ! Version: F90

    ! History:
    ! - 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: time

    real(kind=dp) :: zred_prev

    ! Calculate the change since the previous redshift.
    ! Note: the initial redshift should be ZERO since
    ! the variables are initialized as comoving!
    ! NOTE: Good only for high-z!!!
    zred_prev=zred
    zred = -1+(1.+zred_t0)*((t0+time)/t0)**(-2./3.)

    ! Take the average zfactor between zred_prev and zred
    zfactor=(1.0+zred_prev)/(1.+zred)

  end subroutine redshift_evol

  ! =======================================================================

  subroutine cosmo_evol ()

    ! Calculates the cosmological evolution of space and densities

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: F90 first version

    ! History:
    ! - 19-Nov-2004: first version f77

    use sizes
    use grid
    use sourceprops, only: rsrcpos
    use material
    
    real(kind=dp) :: zfactor3

    zfactor3=zfactor*zfactor*zfactor

    ! Change the grid coordinates
    x(:)=x(:)*zfactor
    y(:)=y(:)*zfactor
    z(:)=z(:)*zfactor

    dr(:)=dr(:)*zfactor
    
    vol=vol*zfactor3

    ! Source positions (multiple source version)
    if (allocated(rsrcpos)) rsrcpos(:,:)=rsrcpos(:,:)*zfactor
    
    ! Change the densities
    ndens(:,:,:)=ndens(:,:,:)/zfactor3

    ! Change the volumes (not needed???)
    ! volx(i,j,k)=volx(i,j,k)*zfactor3
    ! voly(i,j,k)=voly(i,j,k)*zfactor3
    ! volz(i,j,k)=volz(i,j,k)*zfactor3

  end subroutine cosmo_evol

  ! =======================================================================

  real(kind=dp) function cosmo_cool (e_int)

    ! Calculates the cosmological adiabatic cooling

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: f90

    ! History:
    ! 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: e_int

    real(kind=dp) :: dzdt

    ! Cosmological cooling rate:
    ! 2*(da/dt)/a*e
    ! or
    ! 2*(dz/dt)/(1+z)*e
    ! with a the cosmological scale factor

    ! dz/dt (for flat LambdaCDM)
    dzdt=H0*(1.+zred)*sqrt(Omega0*(1.+zred)**3+1.-Omega0)

    !Cooling rate
    cosmo_cool=e_int*2.0/(1.0+zred)*dzdt

  end function cosmo_cool

  ! =======================================================================

  real(kind=dp) function compton_cool (temper,eldens)
    
    ! Calculates the (cosmological) Compton cooling rate
    
    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: first version
    
    ! History:
    ! 16-May-2005: first version f77
    

    ! parameter reference?

    real(kind=dp),intent(in) :: temper ! temperature
    real(kind=dp),intent(in) :: eldens ! electron density
    
    !Cooling rate
    compton_cool=5.65e-36*eldens*(1.0+zred)**4* &
         (temper-cmbtemp*(1.0+zred))
    
  end function compton_cool
  
end module cosmology
