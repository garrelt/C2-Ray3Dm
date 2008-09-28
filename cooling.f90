!>
!! \brief This module contains data and subroutines for radiative cooling
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! \b Version: Collisional Ionization Equilibrium cooling

module radiative_cooling

  ! This module file contains subroutines for radiative cooling
  ! - coolin     - calculate cooling rate
  ! - setup_cool - setup cooling table

  ! Version: CIE cooling

  use precision, only: dp
  
  implicit none

  private

  integer,parameter :: temppoints=61 !< number of points in cooling table
  real(kind=dp),dimension(temppoints) :: cie_cool !< CIE cooling curve
  real(kind=dp) :: mintemp !< lowest log10(temperature) in table
  real(kind=dp) :: dtemp !< steps in log10(temperature) in table

  public :: coolin, setup_cool

contains
  
  !===========================================================================

  !> Calculate the cooling rate
  function coolin(nucldens,eldens,xh,temp0)
    
    real(kind=dp) :: coolin
    
    real(kind=dp),intent(in) :: nucldens !< number density
    real(kind=dp),intent(in) :: eldens !< electron density
    real(kind=dp),dimension(0:1),intent(in) :: xh !< H ionization fractions
    real(kind=dp),intent(in) :: temp0 !< temperature
    
    real(kind=dp) :: tpos, dtpos
    integer :: itpos,itpos1
    
    tpos=(log10(temp0)-mintemp)/dtemp+1.0d0
    itpos=min(temppoints-1,max(1,int(tpos)))
    dtpos=tpos-real(itpos)
    itpos1=min(temppoints,itpos+1)
    
    ! Cooling curve
    coolin=nucldens*eldens* &
         (cie_cool(itpos)+(cie_cool(itpos1)-cie_cool(itpos))*dtpos)
    
    return
  end function coolin
  
  !===========================================================================

  !> Read in and set up the cooling table(s)
  subroutine setup_cool ()

    real(kind=dp),dimension(temppoints) :: temp
    integer :: itemp
    integer :: element,ion,nchck

    ! Open cooling table (CIE)
    open(unit=22,file='tables/corocool.tab',status='old')
    ! Read the cooling data
    do itemp=1,temppoints
       read(22,*) temp(itemp),cie_cool(itemp)
    enddo
    close(22)
    
    mintemp=temp(1)
    dtemp=temp(2)-temp(1)
    ! not needed: maxtemp=temp(temppoints)
    
    ! Convert cooling to from log to linear 
    do itemp=1,temppoints
       cie_cool(itemp)=10.0d0**cie_cool(itemp)
    enddo

    return
  end subroutine setup_cool
  
end module radiative_cooling
