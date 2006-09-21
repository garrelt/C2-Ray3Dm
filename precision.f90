module precision

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2003-012-09
  ! This module is also accepted by the F compiler (Dec 9, 2003)

  ! This module is used by the MPI_Hydrodynamics programme.
  ! It contains the precision definitions

  ! This is single precision
  integer,public,parameter :: si=selected_real_kind(6,37) 
  integer,public,parameter :: dp=selected_real_kind(15,307) 
  integer,public,parameter :: li=selected_int_kind(9)
  real(kind=dp),public,parameter :: tiny_dp=tiny(1.0_dp) ! smallest dp

end module precision
