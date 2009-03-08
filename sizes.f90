!>
!! \brief This module contains basic size parameter definitions
!!
!! Module for C2-Ray
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2006-08-20
module sizes

  ! Module for C2-Ray
  ! Author: Garrelt Mellema
  ! Date: 2006-08-20
  ! This module is also accepted by the F compiler

  ! This module contains basic parameter definitions
  
  ! use precision

  implicit none

  private
  !> Number of spatial dimensions
  integer,parameter,public :: Ndim=3

  !> Size of the mesh for spatial coordinate.


  !for test problem
  integer,dimension(Ndim),parameter,public :: mesh=(/ 100, 100, 100 /)

  !for 3072^3 particles, 6144^3 grid, 114/h Mpc box
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 256, 256, 256 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 384, 384, 384 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 512, 512, 512 /)

  !for 1728^3 particles, 3456^3 grid, 64/h Mpc box
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 216, 216, 216 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 432, 432, 432 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 864, 864, 864 /)

  !for LG simulation: 1024^3 particles, 1024^3 grid, 64/h Mpc box
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 512, 512, 512 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 256, 256, 256 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 128, 128, 128 /)

  !for old-style PMFAST simulation, 1624^3 particles, 3248^3 grid, 100/h or 35/h Mpc box
  !integer,dimension(Ndim),parameter,public :: mesh=(/ 203, 203, 203 /)
  !integer,dimension(Ndim),parameter,public :: mesh=(/406, 406, 406/)
  !integer,dimension(Ndim),parameter,public :: mesh=(/812, 812, 812/)

end module sizes
