module subgrid_clumping

  ! This module specifies the clumping behaviour of the matter

  use c2ray_parameters, only: type_of_clumping,clumping_factor
  use file_admin, only: log

  implicit none

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
       clumping=26.2917*exp(-0.1822*z+0.003505*z*z)
    case(4)
       z=1.0/a-1.0
       clumping=17.57*exp(-0.101*z+0.0011*z*z)
    end select
    write(log,*) 'Setting global clumping factor to ',clumping,'(type ', &
         type_of_clumping,')'
    
  end subroutine set_clumping

end module subgrid_clumping
