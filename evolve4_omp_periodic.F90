module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (3D).
     
  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! Version:
  ! - MPI parallelization over the sources
  ! - OpenMP parallelization over the octants

  ! Needs:
  ! doric : ionization calculation for one point + photo-ionization rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI definitions
  use sizes, only: Ndim, mesh
  use grid, only: x,y,z,vol,dr
  use material, only: ndens, xh, temper
  use sourceprops, only: SrcSeries, NumSrc, srcpos
  use photonstatistics, only: state_before, calculate_photon_statistics, &
       photon_loss
  use c2ray_parameters, only: convergence_fraction

  implicit none

  private

  public :: evolve3D, phih_grid

  real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)) :: phih_grid
  real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1) :: xh_av
  real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1) :: xh_intermed
  real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)) :: coldensh_out
  real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)) :: buffer
  real(kind=dp) :: photon_loss_all

contains
  ! =======================================================================
  
  subroutine evolve3D (dt)

    ! Calculates the evolution of the hydrogen ionization state
     
    ! Author: Garrelt Mellema
     
    ! Date: 21-Aug-2006 (f77/OMP: 13-Jun-2005)

    ! Version: Multiple sources / Using average fractions to converge
    ! loop over sources
    
    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz

    ! For random permutation of sources:
    use  m_ctrper, only: ctrper

    ! The time step
    real(kind=dp),intent(in) :: dt

    ! Will contains the integer position of the cell being treated
    integer,dimension(Ndim) :: pos
      
    ! Loop variables
    integer :: i,j,k,l,nx,ns,ns1,niter
    integer :: naxis,nplane,nquadrant

    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag

#ifdef MPI
    integer :: ierror
#endif

    ! End of declarations

    ! Initial state (for photon statistics)
    call state_before ()

    ! initialize average and intermediate results to initial values
    xh_av(:,:,:,:)=xh(:,:,:,:)
    xh_intermed(:,:,:,:)=xh(:,:,:,:)
    
    ! Loop over sources
    niter=0
    do
       ! Loop counter
       niter=niter+1
         
       ! reset global rates to zero for this iteration
       phih_grid(:,:,:)=0.0
       
       ! reset photon loss counter
       photon_loss=0.0
       
       ! Make a random order :: call in serial
       if ( rank == 0 ) call ctrper (SrcSeries(1:NumSrc),1.0)
       
#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
       
       ! Source Loop - distributed for the MPI nodes
       do ns1=1+rank,NumSrc,npr
          ns=SrcSeries(ns1)
          write(30,*) 'Processor ',rank,' doing source at:',srcpos(:,ns)
          
          ! reset column densities for new source point
          ! coldensh_out is unique for each source point
          coldensh_out(:,:,:)=0.0
          
          ! Loop through grid in the order required by short characteristics
          ! and useful for parallelization

          ! First do source point
          pos(:)=srcpos(:,ns)
          call evolve0D(dt,pos,ns,niter)

          ! do independent areas of the mesh in parallel using OpenMP
          !$omp parallel default(shared)

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

          write(30,*) sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
       enddo ! sources loop
       
       ! End of parallelization
       
       ! Report photon losses over grid boundary
       write(30,*) 'photon loss counter: ',photon_loss
       
#ifdef MPI
       ! accumulate threaded photon loss
       call MPI_ALLREDUCE(photon_loss, photon_loss_all, 1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, ierror)
       
       ! accumulate threaded phih_grid
       
       call MPI_ALLREDUCE(phih_grid, buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, ierror)
       
       phih_grid(:,:,:)=buffer(:,:,:)
       
       ! accumulate threaded xh_intermed
       call MPI_ALLREDUCE(xh_av(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, ierror)
       
       xh_av(:,:,:,1) = buffer(:,:,:)
       xh_av(:,:,:,0) = max(0.0_dp,min(1.0_dp,1.0-xh_av(:,:,:,1)))
       
       ! accumulate threaded xh_intermed
       call MPI_ALLREDUCE(xh_intermed(:,:,:,1), buffer, &
            mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_NEW, ierror)
       
       xh_intermed(:,:,:,1)=buffer(:,:,:)
       xh_intermed(:,:,:,0)=max(0.0_dp,min(1.0_dp,1.0-xh_intermed(:,:,:,1)))
#else
       photon_loss_all=photon_loss
#endif
       photon_loss=photon_loss_all/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
       
       ! Apply total photo-ionization rates from all sources (phih_grid)
       write(30,*) 'Applying Rates'
       conv_flag=0 ! will be used to check for convergence
       
       ! Loop through the entire mesh
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                pos=(/ i,j,k /)
                call evolve0D_global(dt,pos,conv_flag)
             enddo
          enddo
       enddo
       write(30,*) 'Number of non-converged points: ',conv_flag
       
       write(30,*) sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))

       ! Update xh if converged and exit
       if (conv_flag <= int(convergence_fraction*mesh(1)*mesh(2)*mesh(3))) then
          xh(:,:,:,:)=xh_intermed(:,:,:,:)
          exit
       else
          if (niter > 50) then
             ! Complain about slow convergence
             write(30,*) 'Multiple sources not converging'
             exit
          endif
       endif
    enddo

    ! Calculate photon statistics
    call calculate_photon_statistics (dt)

    return
  end subroutine evolve3D

  ! ===========================================================================

  subroutine evolve1D_axis(dt,ns,niter,naxis)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: naxis        ! axis to do

    integer :: i,j,k
    integer,dimension(Ndim) :: pos ! mesh position

    select case (naxis)
    case(1)
       ! sweep in +i direction
       pos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
          pos(1)=i
          call evolve0D(dt,pos,ns,niter) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       pos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
          pos(1)=i
          call evolve0D(dt,pos,ns,niter) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       pos(1)=srcpos(1,ns)
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
          pos(2)=j
          call evolve0D(dt,pos,ns,niter) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       pos(1)=srcpos(1,ns)
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
          pos(2)=j
          call evolve0D(dt,pos,ns,niter) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       pos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          call evolve0D(dt,pos,ns,niter) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       pos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          call evolve0D(dt,pos,ns,niter) !# `negative' k
       end do
    end select
    
    return
  end subroutine evolve1D_axis

  ! ===========================================================================

  subroutine evolve2D_plane(dt,ns,niter,nplane)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nplane        ! plane to do

    integer,dimension(Ndim) :: pos ! mesh position

    integer :: i,j,k

    select case (nplane)
    case(1)
       ! sweep in +i,+j direction
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
          pos(2)=j
          do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
          pos(2)=j
          do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
          pos(2)=j
          do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       pos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
          pos(2)=j
          do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       pos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       pos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       pos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       pos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
             pos(1)=i
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       pos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       pos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       pos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       pos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             call evolve0D(dt,pos,ns,niter)
          enddo
       enddo
       
    end select
    
    return
  end subroutine evolve2D_plane

  ! ===========================================================================

  subroutine evolve3D_quadrant(dt,ns,niter,nquadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    
    real(kind=dp),intent(in) :: dt     ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nquadrant    ! which quadrant to do    
    integer :: i,j,k
    integer,dimension(Ndim) :: pos ! mesh position

    select case (nquadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
                pos(1)=i
                call evolve0D(dt,pos,ns,niter)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,ns)+1,srcpos(3,ns)+mesh(3)/2-1+mod(mesh(3),2)
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)+1,srcpos(2,ns)+mesh(2)/2-1+mod(mesh(2),2)
             pos(2)=j
             do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             do i=srcpos(1,ns)+1,srcpos(1,ns)+mesh(1)/2-1+mod(mesh(1),2)
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,ns)-1,srcpos(3,ns)-mesh(3)/2,-1
          pos(3)=k
          do j=srcpos(2,ns)-1,srcpos(2,ns)-mesh(2)/2,-1
             pos(2)=j
             do i=srcpos(1,ns)-1,srcpos(1,ns)-mesh(1)/2,-1
                pos(1)=i
                call evolve0D(dt,pos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    end select
    return
  end subroutine evolve3D_quadrant

  !=======================================================================

  subroutine evolve0D(dt,rtpos,ns,niter)
    
    ! Calculates the photo-ionization rate for one cell and adds it to 
    ! the collective rate.
    
    ! Author: Garrelt Mellema
    
    ! Date: 01-Feb-2008 (21-Aug-2006, 20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: multiple sources, fixed temperature
    
    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated.
    ! For the first pass (niter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence in the
    ! case of isolated sources.

    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use radiation, only: photoion, photrates
    use c2ray_parameters, only: epsilon,convergence1,convergence2
    use mathconstants, only: pi
    
    real(kind=dp),parameter :: max_coldensh=2e19_dp ! column density for stopping chemisty
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: rtpos ! cell position (for RT)
    integer,intent(in)      :: ns ! source number 
    integer,intent(in)      :: niter ! pass number
    
    integer :: nx,nd,nit,idim ! loop counters
    integer,dimension(Ndim) :: pos
    integer,dimension(Ndim) :: srcpos1
    real(kind=dp) :: dist,path,vol_ph
    real(kind=dp) :: xs,ys,zs
    real(kind=dp) :: coldensh_in
    real(kind=dp) :: coldensh_cell
    real(kind=dp) :: ndens_p
    real(kind=dp) :: avg_temper
    real(kind=dp) :: de
    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0
    real(kind=dp) :: yh_av0
    real(kind=dp) :: convergence
    
    type(photrates) :: phi

    ! set convergence tolerance
    convergence=convergence1

    ! Map pos to mesh pos, assuming a periodic mesh
    pos(:)=modulo(rtpos(:)-1,mesh(:))+1

    ! Initialize local ionization states to the global ones
    do nx=0,1
       yh0(nx)=xh(pos(1),pos(2),pos(3),nx)
       yh_av(nx)=xh_av(pos(1),pos(2),pos(3),nx)
    enddo

    ! Initialize local temperature and density
    avg_temper=temper
    ndens_p=ndens(pos(1),pos(2),pos(3))

    ! Find the column density at the entrance point of the cell (short
    ! characteristics)

    if (rtpos(1) == srcpos(1,ns).and.rtpos(2) == srcpos(2,ns).and. &
         rtpos(3) == srcpos(3,ns)) then
       ! Do not call cinterp for the source point.
       ! Set coldensh and path by hand
       coldensh_in=0.0
       path=0.5*dr(1)
       
       ! Find the distance to the source (average?)
       dist=0.5*dr(1)         ! this makes vol=dx*dy*dz
       !vol_ph=4.0/3.0*pi*dist**3
       vol_ph=dr(1)*dr(2)*dr(3)

    else

       ! For all other points call cinterp to find the column density
       srcpos1(:)=srcpos(:,ns)
       call cinterp(rtpos,srcpos1,coldensh_in,path)
       path=path*dr(1)
          
       ! Find the distance to the source
       xs=dr(1)*real(rtpos(1)-srcpos(1,ns))
       ys=dr(2)*real(rtpos(2)-srcpos(2,ns))
       zs=dr(3)*real(rtpos(3)-srcpos(3,ns))
       dist=sqrt(xs*xs+ys*ys+zs*zs)
         
       ! Find the volume of the shell this cell is part of 
       ! (dilution factor).
       vol_ph=4.0*pi*dist*dist*path

    endif

    ! Only do chemistry if this is the first pass over the sources,
    ! and if column density is below the maximum
    ! On the first pass it may be beneficial to assume isolated sources,
    ! but on later passes the effects of multiple sources has to be
    ! taken into account. Therefore no changes to xh, xh_av, etc.
    ! should happen on later passes!
    if (niter == 1 .and. coldensh_in < max_coldensh) then

       ! Iterate to get mean ionization state (column density / optical depth) 
       ! in cell
       nit=0
       do 
          nit=nit+1

          ! Debug write
          if (niter > 1 .and. nit > 2) write(*,*) niter, nit, pos(1:3)

          ! Store the value of yh_av found in the previous iteration
          ! (for convergence test)
          yh_av0=yh_av(0)

          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(path,yh_av(0),ndens_p)

          ! Calculate (photon-conserving) photo-ionization rate
          call photoion(phi,coldensh_in,coldensh_in+coldensh_cell,vol_ph,ns)
          phi%h=phi%h/(yh_av(0)*ndens_p)

          ! Restore yh to initial values (for doric)
          yh(:)=yh0(:)
              
          ! Calculate (mean) electron density
          de=electrondens(ndens_p,yh_av)

          ! Calculate the new and mean ionization states (yh and yh_av)
          call doric(dt,avg_temper,de,ndens_p,yh,yh_av,phi%h)

          ! Test for convergence on ionization fraction
          ! Depending on how the multiple sources are handled this
          ! could be convergence on yh_av or on yh
          if ((abs((yh_av(0)-yh_av0)/yh_av(0)) < convergence .or. &
               (yh_av(0) < 1e-12))) exit

          ! Warn about non-convergence
          if (nit > 5000) then
             write(30,*) 'Convergence failing (source ',ns,')'
             write(30,*) 'xh: ',yh_av(0),yh_av0
             exit
          endif
       enddo ! end of iteration

       ! Copy ionic abundances back if this is the first iteration
       ! and the first source. This will speed up convergence if
       ! the sources are isolated and only ionizing up.
       ! In other cases it does not make a difference.
       xh_intermed(pos(1),pos(2),pos(3),1)=max(yh(1), &
            xh_intermed(pos(1),pos(2),pos(3),1))
       xh_intermed(pos(1),pos(2),pos(3),0)=1.0- &
            xh_intermed(pos(1),pos(2),pos(3),1)
       xh_av(pos(1),pos(2),pos(3),1)=max(yh_av(1), &
            xh_av(pos(1),pos(2),pos(3),1))
       xh_av(pos(1),pos(2),pos(3),0)=1.0- &
            xh_av(pos(1),pos(2),pos(3),1)
       ! if (niter == 1.and.ns == 1) then
       !do nx=0,1
       ! xh_intermed(pos(1),pos(2),pos(3),nx)=yh(nx)
       !   xh_av(pos(1),pos(2),pos(3),nx)=yh_av(nx)
       !enddo
       !endif

    endif


    ! For niter > 1, just ray trace and exit. Do not touch the ionization
    ! fractions. They are updated using phih_grid in evolve0d_global
    
    ! Add the (time averaged) column density of this cell
    ! to the total column density (for this source)
    coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
         coldens(path,yh_av(0),ndens_p)

    ! Calculate (photon-conserving) photo-ionization rate
    if (coldensh_in < max_coldensh) then
       call photoion(phi,coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
            vol_ph,ns)
       phi%h=phi%h/(yh_av(0)*ndens_p)
    else
       phi%h=0.0
       phi%h_out=0.0
    endif

    ! Save photo-ionization rates, they will be applied in evolve0D_global
    phih_grid(pos(1),pos(2),pos(3))= &
         phih_grid(pos(1),pos(2),pos(3))+phi%h

    ! Photon statistics: register number of photons leaving the grid
    if ( &
         (rtpos(1) == srcpos(1,ns)-1-mesh(1)/2).or. &
         (rtpos(1) == srcpos(1,ns)+mesh(1)/2).or. &
         (rtpos(2) == srcpos(2,ns)-1-mesh(2)/2).or. &
         (rtpos(2) == srcpos(2,ns)+mesh(2)/2).or. &
         (rtpos(3) == srcpos(3,ns)-1-mesh(3)/2).or. &
         (rtpos(3) == srcpos(3,ns)+mesh(3)/2)) then
       
       photon_loss=photon_loss + phi%h_out*vol/vol_ph
    endif

    return
  end subroutine evolve0D

  ! =======================================================================

  subroutine evolve0D_global(dt,pos,conv_flag)

    ! Calculates the evolution of the hydrogen ionization state.

    ! Author: Garrelt Mellema

    ! Date: 20-May-2005 (5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources, Global update, no ray tracing

    ! Multiple sources
    ! Final pass: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence
    
    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use c2ray_parameters, only: convergence1,convergence2

    real(kind=dp),intent(in) :: dt
    integer,dimension(Ndim),intent(in) :: pos
    integer,intent(inout) :: conv_flag

    integer :: nx,nit ! loop counter
    real(kind=dp) :: de
    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0
    real(kind=dp) :: avg_temper
    real(kind=dp) :: ndens_p

    real(kind=dp) :: phih

    real(kind=dp) :: yh_av0
    real(kind=dp) :: convergence
    
    ! This routine does the final (whole grid) pass
    ! Also set convergence tolerance
    convergence=convergence2

    ! Initialize local ionization states to global ones
    do nx=0,1
       yh0(nx)=xh(pos(1),pos(2),pos(3),nx)
       yh(nx)=yh0(nx)
       yh_av(nx)=xh_av(pos(1),pos(2),pos(3),nx) ! use calculated xh_av
    enddo

    ! Initialize local scalars for density and temperature
    ndens_p=ndens(pos(1),pos(2),pos(3))
    avg_temper=temper

    ! Use the collected rates
    phih=phih_grid(pos(1),pos(2),pos(3))
    
    ! Add lost photons
    ! (if the cell is ionized, add a fraction of the lost photons)
    !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
    phih=phih + photon_loss/(vol*xh_av(pos(1),pos(2),pos(3),0)*ndens_p)

    nit=0
    do 
       nit=nit+1
       ! Save the values of yh_av found in the previous
       ! iteration
       yh_av0=yh_av(0)

       ! Copy ionic abundances back to initial values for doric
       yh(:)=yh0(:)
              
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,yh_av)

       ! Calculate the new and mean ionization states
       call doric(dt,avg_temper,de,ndens_p,yh,yh_av,phih)

       ! Test for convergence on ionization fraction
       if ((abs((yh_av(0)-yh_av0)/yh_av(0)) < convergence2 &
             .or. (yh_av(0) < 1e-12))) exit
                  
       ! Warn about non-convergence
       if (nit > 5000) then
          write(30,*) 'Convergence failing (global)'
          write(30,*) 'xh: ',yh_av(0),yh_av0
          exit
       endif
    enddo

    ! Test for convergence on ionization fraction
    yh_av0=xh_av(pos(1),pos(2),pos(3),0) ! use previously calculated xh_av
    if (abs((yh_av(0)-yh_av0)) > convergence2 .and. &
         (abs((yh_av(0)-yh_av0)/yh_av(0)) > convergence2 .and. &
         (yh_av(0) > 1e-12))) then
       conv_flag=conv_flag+1
    endif

    ! Copy ionic abundances back to intermediate global arrays.
    do nx=0,1
       xh_intermed(pos(1),pos(2),pos(3),nx)=yh(nx)
       xh_av(pos(1),pos(2),pos(3),nx)=yh_av(nx)
    enddo

    return
  end subroutine evolve0D_global

  ! ===========================================================================

  subroutine cinterp (pos,srcpos,cdensi,path)
    
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point srcpos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    integer,dimension(Ndim) :: pos
    integer srcpos(Ndim)
    real(kind=dp) :: cdensi
    real(kind=dp) :: path

    integer :: i,j,k,i0,j0,k0

    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4
    real(kind=dp) :: di,dj,dk


    ! map to local variables (should be pointers ;)
    i=pos(1)
    j=pos(2)
    k=pos(3)
    i0=srcpos(1)
    j0=srcpos(2)
    k0=srcpos(3)
    
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel=i-i0
    jdel=j-j0
    kdel=k-k0
    idela=abs(i-i0)
    jdela=abs(j-j0)
    kdela=abs(k-k0)
    
    ! Find coordinates of points closer to source
    sgni=sign(1,idel)
!      if (idel.eq.0) sgni=0
    sgnj=sign(1,jdel)
!      if (jdel.eq.0) sgnj=0
    sgnk=sign(1,kdel)
!      if (kdel.eq.0) sgnk=0
    im=i-sgni
    jm=j-sgnj
    km=k-sgnk
    di=real(i-i0)
    dj=real(j-j0)
    dk=real(k-k0)

    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       if(kdel >= 0) then
          ! ray crosses z plane below destination point
          alam=(real(km-k0)+0.5)/dk
       else
          ! ray crosses z plane above destination point
          alam=(real(km-k0)-0.5)/dk
       end if
              
       xc=alam*di+real(i0) ! x of crossing point on z-plane 
       yc=alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx=2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1=(1.-dx)*(1.-dy)    ! interpolation weights of
       s2=(1.-dy)*dx         ! corner points to c-point
       s3=(1.-dx)*dy
       s4=dx*dy
       
       !s1=((1.-dx)*(1.-dy))**2    ! interpolation weights of
       !s2=((1.-dy)*dx)**2         ! corner points to c-point
       !s3=((1.-dx)*dy)**2
       !s4=(dx*dy)**2

       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jp=modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)    !# column densities at the
       c2=coldensh_out(ip,jmp,kmp)     !# four corners
       c3=coldensh_out(imp,jp,kmp)
       c4=coldensh_out(ip,jp,kmp)
       
       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       ! column density at the crossing point
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4) 

       ! Take care of diagonals
       ! if (kdela.eq.idela.or.kdela.eq.jdela) then
       ! if (kdela.eq.idela.and.kdela.eq.jdela) then
       ! cdensi=sqrt(3.0)*cdensi
       !else
       !cdensi=sqrt(2.0)*cdensi
       !endif
       !endif

       if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
          if (idela == 1.and.jdela == 1) then
             cdensi=sqrt(3.0)*cdensi
          else
             cdensi=sqrt(2.0)*cdensi
          endif
       endif
       ! if (kdela.eq.1) then
       ! if ((w3.eq.1.0).or.(w2.eq.1.0)) cdensi=sqrt(2.0)*cdensi
       ! if (w1.eq.1.0) cdensi=sqrt(3.0)*cdensi
       ! write(30,*) idela,jdela,kdela
       !endif

       ! Path length from c through d to other side cell.
       dxp=di/dk
       dyp=dj/dk
       path=sqrt(dxp*dxp+dyp*dyp+1.0) ! pathlength from c to d point  


       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
          
       if(jdel >= 0) then
          alam=(real(jm-j0)+0.5)/dj
       else
          alam=(real(jm-j0)-0.5)/dj
       end if
       zc=alam*dk+real(k0)
       xc=alam*di+real(i0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dx=2.0*abs(xc-(real(im)+0.5*sgni))
       s1=(1.-dx)*(1.-dz)
       s2=(1.-dz)*dx
       s3=(1.-dx)*dz
       s4=dx*dz
       !s1=((1.-dx)*(1.-dz))**2
       !s2=((1.-dz)*dx)**2
       !s3=((1.-dx)*dz)**2
       !s4=(dx*dz)**2
       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)
       c2=coldensh_out(ip,jmp,kmp)
       c3=coldensh_out(imp,jmp,kp)
       c4=coldensh_out(ip,jmp,kp)

       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       ! Take care of diagonals
       if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
          if (idela == 1.and.kdela == 1) then
             !write(30,*) 'error',i,j,k
             cdensi=sqrt(3.0)*cdensi
          else
             !write(30,*) 'diagonal',i,j,k
             cdensi=sqrt(2.0)*cdensi
          endif
       endif

       dxp=di/dj
       dzp=dk/dj
       path=sqrt(dxp*dxp+1.0+dzp*dzp)
       

       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then
       
       if(idel >= 0) then
          alam=(real(im-i0)+0.5)/di
       else
          alam=(real(im-i0)-0.5)/di
       end if
       zc=alam*dk+real(k0)
       yc=alam*dj+real(j0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1=(1.-dz)*(1.-dy)
       s2=(1.-dz)*dy
       s3=(1.-dy)*dz
       s4=dy*dz
       !s1=((1.-dz)*(1.-dy))**2
       !s2=((1.-dz)*dy)**2
       !s3=((1.-dy)*dz)**2
       !s4=(dz*dy)**2

       imp=modulo(im-1,mesh(1))+1
       jp=modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=coldensh_out(imp,jmp,kmp)
       c2=coldensh_out(imp,jp,kmp)
       c3=coldensh_out(imp,jmp,kp)
       c4=coldensh_out(imp,jp,kp)
       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1)
       w2=s2*weightf(c2)
       w3=s3*weightf(c3)
       w4=s4*weightf(c4)
       
       cdensi=(c1*w1+c2*w2+c3*w3+c4*w4)/(w1+w2+w3+w4)
       
       if ( idela == 1 .and. ( jdela == 1 .or. kdela == 1 ) ) then
          if ( jdela == 1 .and. kdela == 1 ) then
             ! WRITE(30,*) 'error',i,j,k
             cdensi=sqrt(3.0)*cdensi
          else
             ! WRITE(30,*) 'diagonal',i,j,k
             cdensi=sqrt(2.0)*cdensi
          endif
       endif
       
       dyp=dj/di
       dzp=dk/di
       path=sqrt(1.0+dyp*dyp+dzp*dzp)
       
    end if
    
    return
  end subroutine cinterp

  ! =========================================================================

  real(kind=dp) function weightf (cd)

    use cgsphotoconstants, only: sigh

    real(kind=dp),intent(in) :: cd

    ! weightf=1.0
    ! weightf=1.0/max(1.0d0,cd**0.54)
    ! weightf=exp(-min(700.0,cd*0.15*6.3d-18))
    weightf=1.0/max(0.6d0,cd*sigh)

    ! weightf=1.0/log(max(e_ln,cd))

    return
  end function weightf

end module evolve
