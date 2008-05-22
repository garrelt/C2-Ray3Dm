module file_admin

  ! This module contains parameters having to do
  ! with file in/ouput

  implicit none

  character(len=50),parameter :: results_dir="../results/"
  integer,public,parameter :: stdinput=5 ! set to not 5 if input via file (not implemented in all codes, check main program!)
  integer,public,parameter :: logf=30     ! unit number of log file(s)
  integer,public,parameter :: ah3=40     ! unit number of ah3 files (if applicable)
  logical,public :: file_input=.false.

contains
  
  subroutine flag_for_file_input(flag)

    logical,intent(in) :: flag

    ! Sets flag for inpit via file (and not via std input)
    file_input=flag
    
  end subroutine flag_for_file_input

end module file_admin
