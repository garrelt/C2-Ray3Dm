!>
!! \brief This module contains data and routines for calculating the photon statistics
!!
!! Module for C2-Ray
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:  26-Feb-2008
!!
!! \b Version: 3D

module photonstatistics
  
  ! Photon statistics
  ! photon_loss: this is kept here, but calculated in the evolve module.

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: logf
  use cgsconstants, only: albpow,bh00,colh0,temph0
  use sizes, only: mesh
  use grid, only: vol
  use material, only: ndens, temper, clumping, clumping_point
  use tped, only: electrondens
  use sourceprops, only: NormFlux, NumSrc
  use radiation, only: S_star, NumFreqBnd
  use c2ray_parameters, only: type_of_clumping

  implicit none

  !> true if checking photonstatistics
  logical,parameter :: do_photonstatistics=.true.
  !> Total number of recombinations
  real(kind=dp) :: totrec
  !> Total number of collisional ionizations
  real(kind=dp) :: totcollisions
  !> Change in number of neutral H atoms
  real(kind=dp) :: dh0
  !> Total number of ionizing photons used
  real(kind=dp) :: total_ion
  !> Grand total number of ionizing photons used
  real(kind=dp) :: grtotal_ion
  !> Number of photons leaving the grid
  real(kind=dp) :: photon_loss(NumFreqBnd)

  real(kind=dp),private :: h0_before !< number of H atoms at start of time step
  real(kind=dp),private :: h0_after !< number of H atoms at end of time step
  real(kind=dp),private :: h1_before !< number of H ions at start of time step
  real(kind=dp),private :: h1_after !< number of H ions at end of time step

  integer,private :: i !< mesh loop index (x)
  integer,private :: j !< mesh loop index (y)
  integer,private :: k !< mesh loop index (z)

contains

  !----------------------------------------------------------------------------

  !> Initialize the photon statistics
  subroutine initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0

  end subroutine initialize_photonstatistics

  !----------------------------------------------------------------------------

  !> Call the individual routines needed for photon statistics calculation
  subroutine calculate_photon_statistics (dt,xh_l,xh_r)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_r

    ! Call the individual routines needed for this calculation

    call state_after (xh_l) ! number of neutrals after integration
    call total_rates (dt,xh_r) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the start of the time step
  subroutine state_before (xh_l)

    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l

    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_before=h0_before+ndens(i,j,k)*xh_l(i,j,k,0)
             h1_before=h1_before+ndens(i,j,k)*xh_l(i,j,k,1)
          enddo
       enddo
    enddo

    h0_before=h0_before*vol
    h1_before=h1_before*vol

  end subroutine state_before

  !----------------------------------------------------------------------------

  !> Calculates total number of recombinations and collisions
  subroutine total_rates(dt,xh_l)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l

    real(kind=dp),dimension(0:1) :: yh
    real(kind=dp) :: ndens_p ! needed because ndens may be single precision

    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             yh(0)=xh_l(i,j,k,0)
             yh(1)=xh_l(i,j,k,1)
             ndens_p=ndens(i,j,k)
             ! Set clumping to local value if we have a clumping grid
             if (type_of_clumping == 5) &
                  call clumping_point (i,j,k)
             totrec=totrec+ndens_p*xh_l(i,j,k,1)*    &
                  electrondens(ndens_p,yh)*  &
                  clumping*bh00*(temper/1e4)**albpow
             totcollisions=totcollisions+ndens_p*   &
                  xh_l(i,j,k,0)*electrondens(ndens_p,yh)* &
                  colh0*sqrt(temper)*exp(-temph0/temper)
          enddo
       enddo
    enddo

    totrec=totrec*vol*dt
    totcollisions=totcollisions*vol*dt

  end subroutine total_rates
  
  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the end of the time step
  subroutine state_after(xh_l)
    
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l

    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_after=h0_after+ndens(i,j,k)*xh_l(i,j,k,0)
             h1_after=h1_after+ndens(i,j,k)*xh_l(i,j,k,1)
          enddo
       enddo
    enddo
    h0_after=h0_after*vol
    h1_after=h1_after*vol
    
  end subroutine state_after
  
  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine total_ionizations ()
    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    total_ion=totrec+dh0
    
  end subroutine total_ionizations

  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine report_photonstatistics (dt)

    real(kind=dp),intent(in) :: dt

    real(kind=dp) :: totalsrc,photcons,total_photon_loss

    total_photon_loss=sum(photon_loss)*dt* &
         real(mesh(1))*real(mesh(2))*real(mesh(3))
    totalsrc=sum(NormFlux(1:NumSrc))*s_star*dt
    photcons=(total_ion-totcollisions)/totalsrc
    if (rank == 0) then
       write(logf,"(7(1pe10.3))") &
            total_ion, totalsrc, &
            photcons, &
            dh0/total_ion, &
            totrec/total_ion, &
            total_photon_loss/totalsrc, &
            totcollisions/total_ion
    endif
    
  end subroutine report_photonstatistics

end module photonstatistics
