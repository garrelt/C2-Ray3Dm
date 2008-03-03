module subgrid_clumping

  ! This module specifies the clumping behaviour of the matter

  use precision, only: dp, si
  use sizes, only: mesh
  use my_mpi
  use c2ray_parameters, only: type_of_clumping,clumping_factor
  use file_admin, only: log
  use pmfast, only: id_str, dir_dens, tot_nfiles

  implicit none

  private

  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  public :: set_clumping, clumping_point

contains

  ! ===========================================================================

  subroutine set_clumping(z)

    real(kind=dp),intent(in) :: z

    select case (type_of_clumping)
    case(1)
       clumping = clumping_factor
    case(2) 
       clumping = 27.466*exp(-0.114*z+0.001328*z*z)
    case(3)
       clumping=26.2917*exp(-0.1822*z+0.003505*z*z)
    case(4)
       clumping=17.57*exp(-0.101*z+0.0011*z*z)
    case(5)
       call clumping_init (z)
    end select

    if (rank == 0) write(log,*) "Setting (mean) global clumping factor to ", &
         clumping,"(type ", type_of_clumping,")"
    
  end subroutine set_clumping

  ! ===========================================================================

  subroutine clumping_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_clumping /= 5) then
       write(log,*) &
         'Error: calling position dependent, but array is not initialized.'
    else
       clumping = clumping_grid(i,j,k)
    endif

  end subroutine clumping_point

  ! ===========================================================================

  subroutine clumping_init (zred_now)

    ! Initializes position dependent clumping (at redshift zred_now)

    ! Author: Ilian Iliev, modified from code by Garrelt Mellema

    ! Date: 03-Mar-2008 (5-Mar-2007( 20-Aug-2006, 19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004))

    ! Version: 
    ! - Three-dimensional. 
    ! - gas clumping field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    
    integer :: i,j,k,n,nfile ! loop counters
    character(len=512):: clump_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: summed_density
    real(kind=dp) :: avg_dens

    ! clumping in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: clumping_real

#ifdef MPI
    integer :: mympierror
#endif

    if (.not.(allocated(clumping_grid))) &
         allocate(clumping_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,'(f6.3)') zred_now
       if (id_str == "coarsest") then
          ! Allocate array needed to read in data
          allocate(clumping_real(mesh(1),mesh(2),mesh(3)))
          write(30,*) 'Reading ',id_str,' input'
          clump_file=trim(adjustl(dir_dens))//"coarser_densities/"// &
               trim(adjustl(zred_str))// &
               "clump_"//trim(adjustl(id_str))//".dat"
          
          ! Open clumping file: note that it is in `binary' form
          open(unit=20,file=clump_file,form='binary',status='old')
          
          ! Read in data and store it in clumping_grid
          read(20) clumping_real
          clumping_grid(:,:,:)=clumping_real(:,:,:)
          
          ! close file
          close(20)

          ! Deallocate array needed for reading in the data.
          deallocate(clumping_real)

       else if (id_str == "coarser") then
          ! For the highest resolution the density is spread out
          ! over tot_nfiles. Otherwise the same.
          allocate(clumping_real(mesh(1),mesh(2),mesh(3)/tot_nfiles))
          do nfile=0,tot_nfiles-1
             write(nfile_str,'(I1)') nfile
             clump_file=trim(adjustl(dir_dens))// &
                  trim(adjustl(zred_str))// &
                  "clump"//nfile_str//".dat"
             write(30,*) clump_file
             open(unit=20,file=clump_file,form='binary',status='old')
             read(20) clumping_real
             clumping_grid(:,:, &
                  1+nfile*(mesh(3)/tot_nfiles): &
                  (1+nfile)*(mesh(3)/tot_nfiles))=clumping_real(:,:,:)
             close(20)
          enddo
          deallocate(clumping_real)
       else!812^3 run
          print*,"Cannot do that, data unavailable"!for 812^3 runs just use the 406^3 values in each 2x2x2? 
       endif
    endif
    !write(30,*) 'Distributing the clumping'
#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(clumping,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    !write(30,*) 'Density distributed'
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(log,*) 'minimum: ',minval(clumping_grid)
       write(log,*) 'maximum: ',maxval(clumping_grid)
       write(log,*) 'mean clumping: ',clumping
    endif

  end subroutine clumping_init

end module subgrid_clumping
