module photonstatistics
  
  ! This module handles the calculation of the photon statistics
  ! For: C2-Ray

  ! Author: Garrelt Mellema

  ! Date: 26-Feb-2008

  ! Photon statistics
  ! photon_loss is a sum over all sources, summing is done in evolve0d
  ! and evolve3d (in case of parallelization).
 
  use precision, only: dp
  use cgsconstants, only: albpow,bh00,colh0,temph0
  use sizes, only: mesh
  use grid, only: vol
  use material, only: ndens, xh, temper, clumping
  use tped, only: electrondens

  logical,parameter :: do_photonstatistics=.true.
  real(kind=dp) :: totrec
  real(kind=dp) :: totcollisions
  real(kind=dp) :: dh0
  real(kind=dp) :: total_ion
  real(kind=dp) :: grtotal_ion
  real(kind=dp) :: photon_loss

  real(kind=dp),private :: h0_before,h0_after,h1_before,h1_after
  integer,private :: i,j,k

contains
  subroutine initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0

  end subroutine initialize_photonstatistics

  subroutine calculate_photon_statistics (dt)

    real(kind=dp),intent(in) :: dt

    ! Call the individual routines needed for this calculation

    call state_after () ! number of neutrals after integration
    call total_rates (dt) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  subroutine state_before ()

    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_before=h0_before+ndens(i,j,k)*xh(i,j,k,0)
             h1_before=h1_before+ndens(i,j,k)*xh(i,j,k,1)
          enddo
       enddo
    enddo

    h0_before=h0_before*vol
    h1_before=h1_before*vol

  end subroutine state_before

  subroutine total_rates(dt)

    real(kind=dp),intent(in) :: dt

    real(kind=dp),dimension(0:1) :: yh
 
    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             yh(0)=xh(i,j,k,0)
             yh(1)=xh(i,j,k,1)
             totrec=totrec+ndens(i,j,k)*xh(i,j,k,1)*    &
                  electrondens(ndens(i,j,k),yh)*  &
                  clumping*bh00*(temper/1e4)**albpow
             totcollisions=totcollisions+ndens(i,j,k)*   &
                  xh(i,j,k,0)*electrondens(ndens(i,j,k),yh)* &
                  colh0*sqrt(temper)*exp(-temph0/temper)
          enddo
       enddo
    enddo

    totrec=totrec*vol*dt
    totcollisions=totcollisions*vol*dt

  end subroutine total_rates
  
  subroutine state_after()
    
    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_after=h0_after+ndens(i,j,k)*xh(i,j,k,0)
             h1_after=h1_after+ndens(i,j,k)*xh(i,j,k,1)
          enddo
       enddo
    enddo
    h0_after=h0_after*vol
    h1_after=h1_after*vol
    
  end subroutine state_after
  
  subroutine total_ionizations ()
    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    total_ion=totrec+dh0
    
  end subroutine total_ionizations

end module photonstatistics
