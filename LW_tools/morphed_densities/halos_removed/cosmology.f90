module cosmosimple

  ! straightforward zred2time, Kyungjin Ahn

  use precision, only: dp
  use cosmology_parameters, only: H0,Omega0

  implicit none

contains
!---------------------------------------------------------------
  real(kind=dp) function zred2time (zred)
    real(kind=dp), intent(in) :: zred
    
    ! Cosmological time corresponding to redshift zred
    ! NOTE: Good only for high-z!!!
    zred2time = 2d0*(1d0+zred)**(-1.5d0)/(3d0*H0*sqrt(Omega0))
    
  end function zred2time



!---------------------------------------------------------------
  real(kind=dp) function time2zred (time)
    real(kind=dp), intent(in) :: time
    ! NOTE: Good only for high-z!!!
    time2zred = (2d0/(3d0*H0*sqrt(Omega0)*time))**(2d0/3d0) -1d0
    
  end function time2zred

  
end module cosmosimple
