module astroconstants

  use precision, only: dp

  ! A collection of astronomical units and conversion factors
  ! Units: cgs
  
  real(kind=dp),parameter :: R_SOLAR=6.9599e10
  real(kind=dp),parameter :: L_SOLAR=3.826e33
  real(kind=dp),parameter :: M_SOLAR=1.98892d33
  
  real(kind=dp),parameter :: YEAR=3.15576E+07 ! Julian year

  real(kind=dp),parameter :: pc=3.086e18
  real(kind=dp),parameter :: kpc=1e3*pc
  real(kind=dp),parameter :: Mpc=1e6*pc
  real(kind=dp),parameter :: lightyear=9.463e17
  real(kind=dp),parameter :: AU=1.49597870E+13

end module astroconstants
