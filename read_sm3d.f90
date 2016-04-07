module sm3d

use precision, only: dp,si

contains
  subroutine read_sm3d_dp_file_routine(filename,data_array)

    character(len=*),intent(in) :: filename
    real(kind=dp),pointer,dimension(:,:,:),intent(inout) :: data_array

    integer :: m1,m2,m3
    integer,dimension(3) :: shape_of_data

    shape_of_data=shape(data_array)

    ! Open sm3d data file
    open(unit=20,file=filename,form="unformatted",status="old")       
    
    ! Read in data
    read(20) m1,m2,m3
    if (m1 /= shape_of_data(1) .or. m2 /= shape_of_data(2) .or. &
         m3 /= shape_of_data(3)) then
       write(logf,*) "Warning: file ",filename," unusable."
       write(logf,*) "since mesh found in file: ",m1,m2,m3
       write(logf,*) "and ",shape_of_data(:)," needed."
    else
       read(20) data_array
    endif
    ! close file
    close(20)

  end subroutine read_sm3d_dp_file_routine

  subroutine read_sm3d_si_file_routine(filename,data_array)

    character(len=*),intent(in) :: filename
    real(kind=si),pointer,dimension(:,:,:),intent(inout) :: data_array

    integer :: m1,m2,m3
    integer,dimension(3) :: shape_of_data

    shape_of_data=shape(data_array)

    ! Open sm3d data file
    open(unit=20,file=filename,form="unformatted",status="old")       
    
    ! Read in data
    read(20) m1,m2,m3
    if (m1 /= shape_of_data(1) .or. m2 /= shape_of_data(2) .or. &
         m3 /= shape_of_data(3)) then
       write(logf,*) "Warning: file ",filename," unusable."
       write(logf,*) "since mesh found in file: ",m1,m2,m3
       write(logf,*) "and ",shape_of_data(:)," needed."
    else
       read(20) data_array
    endif
    ! close file
    close(20)

  end subroutine read_sm3d_si_file_routine

  subroutine write_sm3d_dp_file_routine(filename,data_array)

    character(len=*),intent(in) :: filename
    real(kind=dp),pointer,dimension(:,:,:),intent(inout) :: data_array

    integer,dimension(3) :: shape_of_data

    shape_of_data=shape(data_array)

    ! Open sm3d data file
    open(unit=20,file=filename,form="unformatted",status="unknown")       
    
    ! Write header and data
    write(20) shape_of_data(1:3)
    write(20) data_array

    ! close file
    close(20)

  end subroutine write_sm3d_dp_file_routine

  subroutine write_sm3d_si_file_routine(filename,data_array)

    character(len=*),intent(in) :: filename
    real(kind=si),pointer,dimension(:,:,:),intent(inout) :: data_array

    integer,dimension(3) :: shape_of_data

    shape_of_data=shape(data_array)

    ! Open sm3d data file
    open(unit=20,file=filename,form="unformatted",status="unknown")       
    
    ! Write header and data
    write(20) shape_of_data(1:3)
    write(20) data_array

    ! close file
    close(20)

  end subroutine write_sm3d_si_file_routine

end module sm3d
