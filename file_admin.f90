module file_admin

  ! This module contains parameters having to do
  ! with file in/ouput

  integer,public,parameter :: stdinput=5 ! set to not 5 if input via file (not implemented in all codes, check main program!)
  integer,public,parameter :: log=30     ! unit number of log file(s)
  integer,public,parameter :: ah3=40     ! unit number of ah3 files (if applicable)

end module file_admin
