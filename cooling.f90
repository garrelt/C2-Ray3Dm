module radiative_cooling

  ! This module file contains subroutines for radiative cooling
  ! - coolin     - calculate cooling rate
  ! - setup_cool - setup cooling table

  ! Version: CIE cooling

  use precision, only: dp
  
  implicit none

  integer,parameter,private :: temppoints=61
  real(kind=dp),dimension(temppoints),private :: cie_cool
  real(kind=dp),private :: mintemp, dtemp

contains
  
  !===========================================================================

  function coolin(nucldens,eldens,xh,temp0)
    
    real(kind=dp) :: coolin
    
    real(kind=dp),intent(in) :: nucldens,eldens,xh(0:1),temp0
    
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

  subroutine setup_cool ()

    real(kind=dp),dimension(temppoints) :: temp
    integer :: itemp
    integer :: element,ion,nchck

    ! Open cooling table (CIE)
    open(unit=22,file='cool.tab',status='old')
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
