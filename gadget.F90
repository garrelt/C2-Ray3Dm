module gadget

  ! This file contains routine having to do with the GADGET N-body
  ! simulation.
  ! The routine in here is
  ! - pmfast_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput
  use my_mpi

  implicit none

  ! redshift sequence information
  integer, public :: NumZred               ! number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array ! array of redshifts

#ifdef MPI
  integer,private :: ierror
#endif

contains

  subroutine gadget_ini ()
    
    character(len=180) :: redshift_file ! name of file with list of redshifts
    integer :: nz ! loop counter

    ! In some cases a special file system is used, and its name is
    ! found from an environment variable.
       
    ! Ask for redshift file
    if (rank == 0) then
       NumZred=1
       allocate(zred_array(NumZred))
       write(*,'(A,$)') 'Initial redshift: '
       read(stdinput,*) zred_array(1)
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
#endif

    return
  end subroutine gadget_ini

end module gadget
