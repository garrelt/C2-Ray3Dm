module sizes

  ! Module for C2-Ray
  ! Author: Garrelt Mellema
  ! Date: 2006-08-20
  ! This module is also accepted by the F compiler

  ! This module contains basic parameter definitions
  
  ! use precision
  implicit none

  private
  integer,parameter,public :: Ndim=3

  !Size of the mesh for spatial coordinate.
  integer,dimension(Ndim),parameter,public :: mesh=(/203, 203, 203/)
  !integer,dimension(Ndim),parameter :: mesh=(/406, 406, 406/)
  !integer,dimension(Ndim),parameter :: mesh=(/812, 812, 812/)

end module sizes
