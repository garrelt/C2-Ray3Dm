module grid

  use precision, only: dp
  use sizes, only: Ndim, mesh
  use astroconstants, only: Mpc
  use cosmology_parameters, only: h
  use my_mpi
  use file_admin, only: stdinput
  !use nbody_simulation, only :: boxsize

  implicit none

  ! Contains grid data
  ! dr - cell size
  ! x,y,z - x,y,z coordinates
  ! vol - volume of one cell
  real(kind=dp),dimension(Ndim) :: dr
  real(kind=dp),dimension(mesh(1)) :: x
  real(kind=dp),dimension(mesh(2)) :: y
  real(kind=dp),dimension(mesh(3)) :: z
  real(kind=dp) :: vol
  
contains

  ! =======================================================================

  subroutine grid_ini()
    
    ! Initializes grid properties
    
    ! Author: Garrelt Mellema
    
    ! Date: 20-Aug-2006 (f77 version: 15-Apr-2004)
    
    ! Version:
    ! Three-dimensional cartesian grid
    
    ! dr - cell size
    ! x,y,z - x,y,z coordinates
    ! vol - volume of one cell
    ! contained in common block in grid.h
    
    integer :: i,j,k
    real(kind=dp) :: xgrid,ygrid,zgrid
    
#ifdef MPI
    integer :: ierror
#endif

    ! Ask for grid size
    if (rank == 0) then
       write(*,'(A,$)') 'Enter comoving size of grid in x,y,z (Mpc/h): '
!       write(*,'(A,$)') 'Enter physical (comoving) size of grid in x,y,z (Mpc): '
       read(stdinput,*) xgrid,ygrid,zgrid
    endif
    
#ifdef MPI
    ! Distribute the total grid size over all processors
    call MPI_BCAST(xgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(ygrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(zgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif
    
    xgrid=xgrid*Mpc/h
    ygrid=ygrid*Mpc/h
    zgrid=zgrid*Mpc/h
    
    ! Calculate cell sizes
    dr(1)=xgrid/real(mesh(1))
    dr(2)=ygrid/real(mesh(2))
    dr(3)=zgrid/real(mesh(3))
    
    ! coordinates of a cell
    do k=1,mesh(3)
       z(k)=(real(k)-0.5)*dr(3)
    enddo
    
    do j=1,mesh(2)
       y(j)=(real(j)-0.5)*dr(2)
    enddo
    
    do i=1,mesh(1)
       x(i)=(real(i)-0.5)*dr(1)
    enddo
    
    ! Volume of a cell.
    ! do k=1,mesh(3)
    ! do j=1,mesh(2)
    ! do i=1,mesh(1)
    ! vol(i,j,k)=dr(1)*dr(2)*dr(3)
    ! enddo
    ! enddo
    !      enddo
    ! Scalar version
    vol=dr(1)*dr(2)*dr(3)
    
    return
  end subroutine grid_ini
  
end module grid
