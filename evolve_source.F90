!>
!! \brief This module contains routines for calculating the effect (ionization 
!! and heating) of one source (3D).
!! The routines in this module handle the stepping through the (3D) grid.
!! The actual work is done per grid point using the routines in the
!! module evolve_point.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-09-05
!!
!! \b Version: 3D, MPI & OpenMP

module evolve_source

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,timefile,iterdump, results_dir, dump_dir
  use abundances, only: abu_he
  use c2ray_parameters, only: subboxsize, max_subbox
  use sizes, only: Ndim, mesh
  use sourceprops, only: NumSrc, srcpos, NormFlux, NormFluxPL !SrcSeries
  use radiation_sed_parameters, only: S_star, pl_S_star
  use photonstatistics, only: photon_loss

  use evolve_data, only: periodic_bc
  use evolve_data, only: coldensh_out
  use evolve_data, only: photon_loss_src_thread
  use evolve_data, only: last_l,last_r
  use evolve_data, only: tn

  use evolve_point, only: evolve0d
  
  implicit none

  save

  private

  ! mesh positions of end points for RT
  integer,dimension(Ndim) :: lastpos_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: lastpos_r !< mesh position of right end point for RT
  integer,public :: sum_nbox !< sum of all nboxes (on one processor)
  integer,public :: sum_nbox_all !< sum of all nboxes (on all processors)

  real(kind=dp) :: photon_loss_src

  public do_source

contains

  ! ===========================================================================
  
  !> Does the ray-tracing over the entire 3D grid for one source.
  !! The number of this source in the current list is ns1.
  subroutine do_source(dt,ns1,niter)

    ! Does the ray-tracing over the entire 3D grid for one source.
    ! The number of this source in the current list is ns1.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer, intent(in) :: ns1 !< number of the source being done
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: naxis,nplane,nquadrant
    integer :: ns
    integer :: k
    integer :: nbox
    integer :: nnt
    integer :: logf1

    ! Total ionizing photon rate of all contributions to the source
    real(kind=dp) :: total_source_flux

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: rtpos
      
    ! Pick up source number from the source list
    !ns=SrcSeries(ns1)
    ns=ns1

    ! Report on source
    !write(logf,*) "Source number: ",ns
    !write(logf,*) NormFlux(ns)
    !write(logf,*) srcpos(:,ns)

    ! reset column densities for new source point
    ! coldensh_out is unique for each source point
    coldensh_out(:,:,:)=0.0

    ! Find the mesh position for the end points of the loop
    ! We trace until we reach max_subbox (set in c2ray_parameters)
    ! or the end of the grid. In the periodic case the end of the
    ! grid is always mesh/2 away from the source. If the grid is
    ! even-sized we trave mesh/2 cells to the left and mesh/2-1
    ! cell to the right. If it is odd, it is mesh/2 in either direction.
    ! The mod(mesh,2) takes care of handling this.
    if (periodic_bc) then
       lastpos_r(:)=srcpos(:,ns)+min(max_subbox,mesh(:)/2-1+mod(mesh(:),2))
       lastpos_l(:)=srcpos(:,ns)-min(max_subbox,mesh(:)/2)
    else
       lastpos_r(:)=min(srcpos(:,ns)+max_subbox,mesh(:))
       lastpos_l(:)=max(srcpos(:,ns)-max_subbox,1)
    endif

    ! Loop through grid in the order required by 
    ! short characteristics
    
    ! Transfer is done in a set of cubes of increasing size.
    ! If the HII region is small we do not waste time calculating
    ! column densities of parts of the grid where no radiation
    ! penetrates. To test whether the current subbox is large
    ! enough we use the photon_loss_src. If this is non-zero,
    ! photons are leaving this subbox and we need to do another
    ! one. We also stop once we have done the whole grid.
    nbox=0 ! subbox counter
    total_source_flux=NormFlux(ns)*S_star ! &
         !! + NormFluxPL(ns)*pl_S_star
    photon_loss_src=total_source_flux !-1.0 ! to pass the first while test
    last_r(:)=srcpos(:,ns) ! to pass the first while test
    last_l(:)=srcpos(:,ns) ! to pass the first while test

    ! Loop through boxes of increasing size
    ! NOTE: make this limit on the photon_loss a fraction of
    ! a source flux loss_fraction*NormFlux(ns)*S_star)
    do while (photon_loss_src > 1e-6*total_source_flux &
    !do while (all(photon_loss_src(:) /= 0.0) &
         .and. last_r(3) < lastpos_r(3) &
         .and. last_l(3) > lastpos_l(3))
       nbox=nbox+1 ! increase subbox counter
       photon_loss_src = 0.0 ! reset photon_loss_src to zero
       photon_loss_src_thread(:) = 0.0 ! reset photon_loss_src to zero
       last_r(:)=min(srcpos(:,ns)+subboxsize*nbox,lastpos_r(:))
       last_l(:)=max(srcpos(:,ns)-subboxsize*nbox,lastpos_l(:))

       ! OpenMP: if we have multiple OpenMP threads (nthreads > 1) we 
       ! parallelize over the threads by doing independent parts of
       ! the mesh.
       if (nthreads > 1) then ! OpenMP parallelization

          ! First do source point (on first pass)
          if (nbox == 1) then
             rtpos(:)=srcpos(:,ns)
             call evolve0D(dt,rtpos,ns,niter)
          endif

          ! do independent areas of the mesh in parallel using OpenMP
          !$omp parallel default(shared) private(tn)
          !!!reduction(+:photon_loss_src)

          ! Find out your thread number
#ifdef MY_OPENMP
          tn=omp_get_thread_num()+1
#else
          tn=1
#endif
          
          ! Then do the the axes
          !$omp do schedule(dynamic,1)
          do naxis=1,6
             call evolve1D_axis(dt,ns,niter,naxis)
          enddo
          !$omp end do

          ! Then the source planes
          !$omp do schedule (dynamic,1)
          do nplane=1,12
             call evolve2D_plane(dt,ns,niter,nplane)
          end do
          !$omp end do

          ! Then the quadrants
          !$omp do schedule (dynamic,1)
          do nquadrant=1,8
             call evolve3D_quadrant(dt,ns,niter,nquadrant)
          end do
          !$omp end do

          !$omp end parallel
          ! Collect photon losses for each thread
          do nnt=1,nthreads
             photon_loss_src=photon_loss_src + &
                  photon_loss_src_thread(nnt)
          enddo

       else ! No OpenMP parallelization

          ! 1. transfer in the upper part of the grid 
          !    (srcpos(3)-plane and above)
          do k=srcpos(3,ns),last_r(3)
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do
          
          ! 2. transfer in the lower part of the grid (below srcpos(3))
          do k=srcpos(3,ns)-1,last_l(3),-1
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do

          ! No OpenMP threads so we use position 1
          ! GM/121127: previous versions of the code did not have
          ! the variable tn set to 1 if we were not running OpenMP.
          ! This led to non-photon-conservations (and should have
          ! led to memory errors...)
          photon_loss_src=photon_loss_src_thread(1)

       endif

    enddo

    ! Record the final photon loss, this is the photon loss that leaves
    ! the grid.
    photon_loss(1)=photon_loss(1) + photon_loss_src

    ! Sum the total number of subboxes used for reporting later
    sum_nbox=sum_nbox+nbox

  end subroutine do_source

  ! ===========================================================================

  !> Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
  !! directions.
  subroutine evolve2D(dt,rtpos,ns,niter)

    ! Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
    ! directions.
    
    real(kind=dp),intent(in) :: dt      !! passed on to evolve0D
    integer,dimension(Ndim),intent(inout) :: rtpos !< mesh position, pos(3) is
                                                 !! intent(in)
    integer,intent(in) :: ns           !< current source
    integer,intent(in) :: niter        !< passed on to evolve0D

    integer :: i,j ! mesh positions
    
    ! sweep in `positive' j direction
    do j=srcpos(2,ns),last_r(2)
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
          
       end do
    end do
    
    ! sweep in `negative' j direction
    do j=srcpos(2,ns)-1,last_l(2),-1
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
       end do
    end do

  end subroutine evolve2D

  ! ===========================================================================

  ! Ray tracing for the axes going through the source point
  ! should be called after having done the source point
  subroutine evolve1D_axis(dt,ns,niter,naxis)

    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: naxis        ! axis to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (naxis)
    case(1)
       ! sweep in +i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)+1,last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `negative' k
       end do
    end select
    
  end subroutine evolve1D_axis

  ! ===========================================================================

  !> Ray tracing for planes containing the source point
  !! should be called after evolve1D_axis
  subroutine evolve2D_plane(dt,ns,niter,nplane)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nplane        ! plane to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nplane)
    case(1)
       ! sweep in +i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
       
    end select
    
  end subroutine evolve2D_plane

  ! ===========================================================================

  !> Ray tracing for the 8 octants 
  !! should be called after evolve2D_plane
  subroutine evolve3D_quadrant(dt,ns,niter,nquadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    
    real(kind=dp),intent(in) :: dt     ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nquadrant    ! which quadrant to do    

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nquadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    end select

  end subroutine evolve3D_quadrant

end module evolve_source
