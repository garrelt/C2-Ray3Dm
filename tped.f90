!>
!! \brief This module contains routines having to do with the calculation of
!!  temperature and electron density.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-Mar-08 (but older)
!!
!! \b Doxygen \b note: This module contains "elemental" functions (i.e.
!! functions without side effects working on scalars). Doxygen does not
!! recognize these as functions, and therefore does not list them.
!!
!! The elemental functions in this module are:
!! - temper2pressr (calculates pressure from temperature, density and electron density).
!! - pressr2temper (calculates temperature from pressure, density and electron density).
!! - rho2n (calculates number density from mass density)
!! - n2rho (calculates mass density from number density)
!!
module tped

  ! This file contains routines having to do with the calculation of
  ! temperature and electron density.

  ! - temper2pressr : find temperature from pressure
  ! - pressr2temper : find pressure from temperature
  ! - electrondens:   find electron density

  use precision, only: dp
  use cgsconstants, only: kb, m_p
  use abundances, only: abu_c, mu

  implicit none

contains

  !=======================================================================

  !> find temperature from pressure
  elemental function temper2pressr (temper,ndens,eldens) result(pressr)
    
    real(kind=dp),intent(in) :: ndens !< number density
    real(kind=dp),intent(in) :: temper !< temperature
    real(kind=dp),intent(in) :: eldens !< electron density
    
    real(kind=dp) :: pressr !< pressure

    pressr=(ndens+eldens)*kb*temper

    !temper2pressr=pressr

  end function temper2pressr

  ! =======================================================================

  !> find pressure from temperature
  elemental function pressr2temper (pressr,ndens,eldens) result(temper)

    real(kind=dp),intent(in) :: pressr !< pressure
    real(kind=dp),intent(in) :: ndens !< number density
    real(kind=dp),intent(in) :: eldens !< electron density
    
    real(kind=dp) :: temper !< temperature

    temper=pressr/(kb*(ndens+eldens))
    
    !pressr2temper=temper
    
  end function pressr2temper
      
  ! =======================================================================

  !> find electron density
  function electrondens(ndens,xh)
      
    real(kind=dp) :: electrondens 
    real(kind=dp),intent(in) :: ndens !< number density
    real(kind=dp),intent(in) :: xh(0:1) !< H ionization fractions

    electrondens=ndens*(xh(1)+abu_c)

  end function electrondens

  ! =======================================================================

  !> Find number density from mass density
  elemental function rho2n(rho)
      
    ! Calculates number density from mass density
    
    real(kind=dp) :: rho2n
    real(kind=dp),intent(in) :: rho !< mass density

    rho2n=rho/(mu*m_p)

  end function rho2n

  ! =======================================================================

  !> Find mass density from number density
  elemental function n2rho(ndens)
      
    ! Calculates number density from mass density

    real(kind=dp) :: n2rho
    real(kind=dp),intent(in) :: ndens !< number density
    
    n2rho=ndens*m_p*mu

  end function n2rho


end module tped
