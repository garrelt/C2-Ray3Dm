module subgrid_clumping

  ! This module specifies the clumping behaviour of the matter

  use cosmology, only: zred
  
  implicit none

  integer,parameter,private :: type_of_clumping=1
  real,parameter,private :: clumping_factor=1.0
  real,public :: clumping

contains
  subroutine set_clumping(a)

    real,intent(in) :: a
    real :: z

    select case (type_of_clumping)
    case(1)
       clumping = clumping_factor
    case(2) 
       z=1.0/a-1.0
       clumping = 27.466*exp(-0.114*z+0.001328*z*z)
    case(3)
       z=1.0/a-1.0
       clumping=26.2917*exp(-0.1822*zred+0.003505*zred*zred)
    case(4)
       z=1.0/a-1.0
       clumping=17.57*exp(-0.101*zred+0.0011*zred*zred)
    end select
    
  end subroutine set_clumping

end module subgrid_clumping
