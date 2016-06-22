!>
!! \brief This module contains data and routines for handling the physical grid (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 20-Aug-2006 (f77 version: 15-Apr-2004)
!!
!!
module grid

  use precision, only: dp
  use sizes, only: Ndim, mesh
  use astroconstants, only: Mpc
  use cosmology_parameters, only: h
  use my_mpi
  use file_admin, only: stdinput,logf,file_input
  use nbody, only : boxsize

  implicit none

  ! Contains grid data
  ! dr - cell size
  ! x,y,z - x,y,z coordinates
  ! vol - volume of one cell
  real(kind=dp),dimension(Ndim) :: dr !< cell size
  real(kind=dp),dimension(Ndim) :: dr_cMpc !< cell size in comoving Mpc
  real(kind=dp),dimension(:),allocatable :: x !< spatial coordinate x
  real(kind=dp),dimension(:),allocatable :: y !< spatial coordinate y
  real(kind=dp),dimension(:),allocatable :: z !< spatial coordinate z
  real(kind=dp) :: vol !< volume of grid cell
  real(kind=dp) :: vol_cMpc3 !< volume of grid cell in cMpc^3
  real(kind=dp) :: sim_volume !< volume of entire simulation box
  real(kind=dp) :: sim_volume_cMpc3 !< volume of entire simulation box in cMpc^3
  
contains

  ! =======================================================================

  !> Initializes grid properties
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
    integer :: alloc_status
    real(kind=dp) :: xgrid,ygrid,zgrid
    
#ifdef MPI
    integer :: ierror
#endif

    ! Ask for grid size (if rank 0 and not set in nbody module)
    if (rank == 0) then
       if (boxsize == 0.0) then
          if (.not.file_input) then
             write(*,'(A,$)') 'Enter comoving size of grid in x,y,z (Mpc/h): '
          endif
          read(stdinput,*) xgrid,ygrid,zgrid
       else
          xgrid=boxsize
          ygrid=boxsize
          zgrid=boxsize
       endif
       ! Report
       write(logf,*) "Box size is ",xgrid," Mpc/h (comoving)"
       flush(logf)
    endif
    
#ifdef MPI
    ! Distribute the total grid size over all processors
    call MPI_BCAST(xgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#ifdef MPILOG
    write(logf,*) ierror
#endif
    call MPI_BCAST(ygrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#ifdef MPILOG
    write(logf,*) ierror
#endif
    call MPI_BCAST(zgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#ifdef MPILOG
    write(logf,*) ierror
#endif
#endif
    
    ! Simulation volume (comoving)
    sim_volume_cMpc3=xgrid*ygrid*zgrid/(h**3)
    sim_volume=sim_volume_Mpc3*Mpc**3

    ! Divide out the Hubble factor
    xgrid=xgrid/h
    ygrid=ygrid/h
    zgrid=zgrid/h
    
    ! Calculate cell sizes in cMpc
    dr_cMpc(1)=xgrid/real(mesh(1))
    dr_cMpc(2)=ygrid/real(mesh(2))
    dr_cMpc(3)=zgrid/real(mesh(3))

    ! Calculate cell sizes in cm
    dr(:)=dr_cMpc(:)*Mpc
    
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
    vol_cMpc3=vol/(Mpc**3)
    
    ! allocate coordinates
    allocate(x(mesh(1)),stat=alloc_status)
#ifdef MPILOG
    write(logf,*) alloc_status
#endif
    allocate(y(mesh(2)),stat=alloc_status)
#ifdef MPILOG
    write(logf,*) alloc_status
#endif
    allocate(z(mesh(3)),stat=alloc_status)
#ifdef MPILOG
    write(logf,*) alloc_status
#endif

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
    
#ifdef MPILOG
    write(logf,*) "End of grid_ini"
    flush(logf)
#endif

  end subroutine grid_ini
  
end module grid
