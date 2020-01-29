module readinput

  use precision, only: dp,si
  use file_admin, only: file_input
  use my_mpi



  ! redshift sequence information
  integer, public                               :: NumZred               !< number of redshifts
  real(kind=dp),dimension(:),allocatable,public :: zred_array !< array of redshifts 
  real(kind=dp), public                         :: boxsize               !< boxsize in Mpc/h


contains
  subroutine read_input(nz0)
    integer,intent(out) :: nz0    ! nz0 is the number of the starting slice
    integer             :: nz
    character(len=180)  :: redshift_file ! name of file with list of redshifts
    integer             :: mympierror

    if (rank == 0) then
       if (.not.file_input) write(*,"(A,$)") "boxsize in Mpc/h"
       read(*,*) boxsize
       if (.not.file_input) write(*,"(A,$)") "Number of starting slice: "
       read(*,*) nz0
       if (.not.file_input) write(*,"(A,$)") "File with redshifts: "
       read(*,*) redshift_file
       
       ! Open and read redshift file
       open(unit=60,file=redshift_file,form="formatted",status="old")
       read(unit=60,fmt=*) NumZred
       allocate(zred_array(NumZred))
       do nz=1,NumZred
          read(unit=60,fmt=*) zred_array(nz)
       enddo
       close(60)
    endif

#ifdef MPI       
    call MPI_BCAST(nz0,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(NumZred,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) allocate(zred_array(NumZred))
    call MPI_BCAST(zred_array,NumZred,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         mympierror)
#endif

  end subroutine read_input

end module readinput
