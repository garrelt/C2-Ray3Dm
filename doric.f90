module doric_module
  
  ! This file contains routines having to do with the calculation of
  ! the ionization evolution of a single point.
  ! It can used for Yguazu-a or non-hydro photo-ionization calculations.

  ! - doric : time dependent solution of the ionization equation
  ! - photoion : photo-ionization rates for given optical depths (moved
  !              to radiation module).
  ! - coldens : column density of a single cell.
  ! - coldensh_bndry : Boundary condition for H column density

  use precision

  implicit none

contains

    !=======================================================================
    subroutine doric (dt,temp0,rhe,rhh,xfh,xfh_av,phih,finalpass)
      
      ! calculates time dependent ionization state for hydrogen
      
      ! Author: Garrelt Mellema
      ! Date: 11-May-2005 (f90) (15 June 2004
      
      ! Version:
      ! Simplified version of the original doric as used in Coral,
      ! for testing photon conservation
      ! - H only
      
      ! Modifications:
      ! 15-Jun-2004 (GM): test for (1.0d0-ee)/deltht to avoid fp errors.
      ! 11-May-2005 (GM): f90
      
      use cgsconstants
      use subgrid_clumping
      use radiation
      use tped
      use c2ray_parameters, only: epsilon

      real(kind=dp),intent(in) :: dt,rhh
      real(kind=dp),intent(inout) :: rhe
      real(kind=dp),dimension(0:1),intent(out) :: xfh
      real(kind=dp),dimension(0:1),intent(inout) :: xfh_av
      
      integer :: itt
      real(kind=dp) :: temp0,brech0,sqrtt0,acolh0
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
      ! (divide by (time averaged) density of neutrals)
      aphoth0=phih/((xfh_av(0)+epsilon)*rhh)
      
      ! determine the hydrogen and helium ionization states and 
      ! electron density
      ! (schmidt-voigt & koeppen 1987)
      ! The electron density is the time-averaged value.
      ! This needs to be iterated, in this version the iteration
      ! is assumed to take place outside the doric routine.
      
      ! Save old values
      rhe0=rhe
      xfh1old=xfh(1)
      xfh0old=xfh(0)
      
      aih0=aphoth0+rhe*acolh0
      delth=aih0+rhe*brech0
      eqxfh1=aih0/delth
      eqxfh0=rhe*brech0/delth
      deltht=delth*dt
      ee=exp(-deltht)
      xfh(1)=(xfh1old-eqxfh1)*ee+eqxfh1
      xfh(0)=(xfh0old-eqxfh0)*ee+eqxfh0
      rhe=electrondens(rhh,xfh)
      
      ! determine neutral densities (take care of precision fluctuations)
      if (xfh(0).lt.epsilon.and.abs(xfh(0)).lt.1.0e-10) xfh(0)=epsilon
      
      ! Determine average ionization fraction over the time step
      ! Mind fp fluctuations. (1.0-ee)/deltht should go to 1.0 for
      ! small deltht, but finite precision leads to values slightly
      ! above 1.0 and for very small values even to 0.0.
      if (deltht.lt.1.0e-8) then
         avg_factor=1.0
      else
         avg_factor=(1.0-ee)/deltht
      endif
      !     The question here is whether it would be better to calculate
      ! xfh_av(0) first, and xfh_av(1) from it.
      xfh_av(1)=eqxfh1+(xfh1old-eqxfh1)*avg_factor
      xfh_av(0)=1.0-xfh_av(1)
      
      ! Take care of precision
      if (xfh_av(0).lt.epsilon.and.abs(xfh_av(0)).lt.1.0e-10) xfh_av(0)=epsilon
      
      return
    end subroutine doric
    
    ! =======================================================================
    
    function coldens(path,xh0,ndens)
      
      
      ! Calculates the column density (of hydrogen)
      ! for a cell of ionization fraction xh, length dr, 
      ! and density ndens

      real(kind=dp) :: coldens
      real(kind=dp),intent(in) :: path
      real(kind=dp),intent(in) :: xh0,ndens
      
      ! Column density over a distance dr (cell)
      
      coldens=xh0*ndens*path
      
      return
    end function coldens
    
    !=======================================================================
    
    function coldensh_bndry()
      
      ! Sets the boundary condition for the hydrogen column density
      
      use cgsconstants
      use cgsphotoconstants
      use radiation
      
      real(kind=dp):: coldensh_bndry

      ! Column density at the left of the first cell
      coldensh_bndry=tauHI/sigh
      
      return
    end function coldensh_bndry
    
  end module doric_module
  
