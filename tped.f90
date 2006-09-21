module tped

  ! This file contains routines having to do with the calculation of
  ! temperature and electron density.

  ! - temper2pressr : find temperature from pressure
  ! - pressr2temper : find pressure from temperature
  ! - electrondens:   find electron density

  use precision
  use cgsconstants
  use abundances

  implicit none

contains

  !=======================================================================

  function temper2pressr (temper,ndens,eldens) result(pressr)
    
    real(kind=dp),intent(in) :: ndens
    real(kind=dp),intent(in) :: temper
    real(kind=dp),intent(in) :: eldens
    
    real(kind=dp) :: pressr

    pressr=(ndens+eldens)*kb*temper

    !temper2pressr=pressr

  end function temper2pressr

  ! =======================================================================

  function pressr2temper (pressr,ndens,eldens) result(temper)

    real(kind=dp),intent(in) :: pressr
    real(kind=dp),intent(in) :: ndens
    real(kind=dp),intent(in) :: eldens
    
    real(kind=dp) :: temper

    temper=pressr/(kb*(ndens+eldens))
    
    !pressr2temper=temper
    
  end function pressr2temper
      
  ! =======================================================================

  function electrondens(ndens,xh)
      
    real(kind=dp) :: electrondens
    real(kind=dp),intent(in) :: ndens
    real(kind=dp),intent(in) :: xh(0:1)

    electrondens=ndens*(xh(1)+abu_c)

  end function electrondens

  ! =======================================================================

  function rho2n(rho)
      
    ! Calculates number density from mass density
    
    real(kind=dp) :: rho2n
    real(kind=dp),intent(in) :: rho

    rho2n=rho/(mu*m_p)

  end function rho2n

  ! =======================================================================

  function n2rho(ndens)
      
    ! Calculates number density from mass density

    real(kind=dp) :: n2rho
    real(kind=dp),intent(in) :: ndens
    
    n2rho=ndens*m_p*mu

  end function n2rho

end module tped
