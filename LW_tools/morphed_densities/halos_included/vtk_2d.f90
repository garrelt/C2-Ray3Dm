module vtk_2d

  use precision, only: dp, si

contains
!--------------------------------------------------
  subroutine vtk2d(array, N1, N2, unitnum, varname, varlen)
    integer, intent(in)               :: N1, N2, unitnum, varlen
    real(kind=dp), intent(in)         :: array(N1,N2)
    character(len=varlen), intent(in) :: varname
    integer                           :: i, j

222 format(A)
    write(unitnum,222)  "# vtk DataFile Version 2.0"
    write(unitnum,222)  "Volume example"
    write(unitnum,222)  "ASCII"
    write(unitnum,222)  "DATASET STRUCTURED_POINTS"
    write(unitnum,223)  "DIMENSIONS ", N1, N2, 1
223 format(A,3I5)
    write(unitnum,222)  "ORIGIN 0. 0. 0."
    write(unitnum,222)  "SPACING 1. 1. 1."
    write(unitnum,224)  "POINT_DATA ", N1*N2
224 format(A,I10)
    write(unitnum,222)  ""
    write(unitnum,225) "SCALARS ", varname, " float"
225 format(3A)
    write(unitnum,222) "LOOKUP_TABLE default"
    do j=   1, N2
       do i=1, N1
          write(1,'(E14.7)') array(i,j)
       enddo
    enddo
    
  end subroutine vtk2d

end module vtk_2d
