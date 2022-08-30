!>
!! \brief This module contains routines for calculating the ionization evolution of a single grid point.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!!

module doric_module
  
  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of a single grid point.
  ! It can used for Yguazu-a or non-hydro photo-ionization calculations.

  ! - doric : time dependent solution of the ionization equation
  ! - coldens : column density of a single cell.
  ! - coldensh_bndry : Boundary condition for H column density

  use precision, only: dp

  implicit none

  ! No public data associated with this module

contains
  
  !=======================================================================

  !> calculates time dependent ionization state for hydrogen
  subroutine doric (dt,temp0,rhe,rhh,xfh,xfh_av,phih)
    
    ! calculates time dependent ionization state for hydrogen
    
    ! Author: Garrelt Mellema
    ! Date: 11-May-2005 (f90) (15 June 2004
    
    ! Version:
    ! Simplified version of the original doric as used in Coral,
    ! - H only
    
    ! Modifications:
    ! 15-Jun-2004 (GM): test for (1.0d0-ee)/deltht to avoid fp errors.
    ! 11-May-2005 (GM): f90
    
    use cgsconstants, only: bh00,albpow,colh0,temph0
    use clumping_module, only: clumping
    !use material, only: clumping
    !use radiation
    use tped, only: electrondens ! should this really be used inside doric?
    use c2ray_parameters, only: epsilon
    
    real(kind=dp),parameter :: sqrtt_isothermal=sqrt(1e4)
    real(kind=dp),parameter :: acolh0_isothermal=colh0* &
         sqrtt_isothermal*exp(-temph0/1e4)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: temp0 !< local temperature
    real(kind=dp),intent(inout) :: rhe !< electron density
    real(kind=dp),intent(in) :: rhh !< H density (or total density?)
    real(kind=dp),dimension(0:1),intent(inout) :: xfh !< H ionization fractions
    real(kind=dp),dimension(0:1),intent(inout) :: xfh_av !< H ionization fractions (time-averaged)
    real(kind=dp),intent(in) :: phih !< photo-ionization rate
    
    real(kind=dp) :: brech0,sqrtt0,acolh0
    real(kind=dp) :: rhe0,xfh1old,xfh0old,aih0,delth,eqxfh0,eqxfh1
    real(kind=dp) :: aphoth0
    real(kind=dp) :: deltht,ee
    real(kind=dp) :: avg_factor
    
    ! find the hydrogen recombination rate at the local temperature
    brech0=clumping*bh00*(temp0/1e4)**albpow
    
    ! find the hydrogen collisional ionization rate at the local 
    ! temperature
    sqrtt0=sqrt(temp0)
    acolh0=colh0*sqrtt0*exp(-temph0/temp0)
    
    ! Find the true photo-ionization rate
    aphoth0=phih
    
    ! determine the hydrogen and helium ionization states and 
    ! electron density
    ! (schmidt-voigt & koeppen 1987)
    ! The electron density is the time-averaged value.
    ! This needs to be iterated, in this version the iteration
    ! is assumed to take place outside the doric routine.
    
    ! Save old values
    !rhe0=rhe
    xfh1old=xfh(1)
    xfh0old=xfh(0)

    ! Calculate equilibrium values
    aih0=aphoth0+rhe*acolh0
    delth=aih0+rhe*brech0
    eqxfh1=aih0/delth
    eqxfh0=rhe*brech0/delth
    
    ! Assign equilibrium values
    xfh(1)=eqxfh1
    xfh(0)=eqxfh0
    
    ! The question here is whether it would be better to calculate
    ! xfh_av(0) first, and xfh_av(1) from it.
    xfh_av(1)=eqxfh1
    xfh_av(0)=eqxfh0
    
  end subroutine doric
  
  ! =======================================================================
  
  !> Calculates the column density (of hydrogen)
  !! for a cell of ionization fraction xh, length path, 
  !! and density ndenstime dependent ionization state for hydrogen
  function coldens(path,xh0,ndens)
    
    ! Calculates the column density (of hydrogen)
    ! for a cell of ionization fraction xh, length dr, 
    ! and density ndens
    
    real(kind=dp) :: coldens
    real(kind=dp),intent(in) :: path !< path length through cell
    real(kind=dp),intent(in) :: xh0  !< neutral H fraction
    real(kind=dp),intent(in) :: ndens !< number density of H
    
    ! Column density over a distance dr (cell)
    coldens=xh0*ndens*path
    
  end function coldens
  
  !=======================================================================

  !> Sets the boundary condition for the hydrogen column density
  function coldensh_bndry()
    
    ! Sets the boundary condition for the hydrogen column density
    
    use cgsphotoconstants, only: sigh => sigma_HI_at_ion_freq
    use radiation_sizes, only: boundary_tauHI
    
    real(kind=dp):: coldensh_bndry
    
    ! Column density at the left of the first cell
    coldensh_bndry=boundary_tauHI/sigh
    
  end function coldensh_bndry
  
end module doric_module
  
