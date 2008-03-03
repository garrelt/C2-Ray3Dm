module output_module
  
  ! This file contains routines having to do with the output
  ! of Ifront programs.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: stdinput

  implicit none
  
  private

  integer,parameter :: max_input_streams=5
  integer,dimension(max_input_streams) :: streams

  real(kind=dp) :: grtotalsrc

  public :: setup_output, output, close_down

contains
  !----------------------------------------------------------------------------

  subroutine setup_output ()
    
    ! Sets up output stream
    
    ! Version: Five streams
    ! Stream1:
    ! Ifront1.out contains a line of constant y and z going through the
    ! source for all timesteps. (formatted)
    
    ! Stream2: 
    ! Ionization fractions for the full data cube (unformatted)
    ! Ifront3',f5.3,'.bin'
    
    ! Stream3: 
    ! Temperature for the full data cube (unformatted)
    ! temper3',f5.3,'.bin'
    
    ! Stream 4:
    ! Ionization fractions in a plane for one time step
    ! Ifront2_xy_',f5.3,'.bin'
    ! Ifront2_xz_',f5.3,'.bin'
    ! Ifront2_yz_',f5.3,'.bin'
    
    ! Stream 5:
    ! Densities in a plane for one time step
    ! ndens_xy_',f5.3,'.bin'
    ! ndens_xz_',f5.3,'.bin'
    ! ndens_yz_',f5.3,'.bin'
    
    ! photon statistics
    use photonstatistics, only: do_photonstatistics, &
         initialize_photonstatistics

    if (rank == 0) then
       write(*,*) 'Which output streams do you want?'
       write(*,*) 'Enter a mask for streams 1 through 5'
       write(*,*) 'E.g. 1,0,1,0,0 means streams 1 and 3, but not 2, 4, 5'
    
       read(stdinput,*) streams(1),streams(2),streams(3),streams(4),streams(5)
       
       ! Open files
       if (do_photonstatistics) then
          open(unit=90,file='PhotonCounts.out', &
               form='formatted',status='unknown')
          write(90,*) 'redshift, photon conservation number, fraction new ionization, fraction recombinations, fraction photon losses, fraction collisional ionization, grand total photon conservation number'
          open(unit=95,file='PhotonCounts2.out', &
               form='formatted',status='unknown')
          write(95,*) 'redshift, total number of ions, grand total ionizing photons, mean ionization fraction (by volume and mass)'
       endif
    endif
    if (do_photonstatistics) then
       call initialize_photonstatistics ()
       grtotalsrc=0.0
    endif

  end subroutine setup_output
  
  !-----------------------------------------------------------------------------
  subroutine close_down ()
    
    ! Closes down
    
    if (rank == 0) then
       if (streams(1).eq.1) close(51)
       if (streams(2).eq.1) close(52)
       if (streams(3).eq.1) close(53)
       close(90)
    endif

  end subroutine close_down
  
  !----------------------------------------------------------------------------

  subroutine output(zred_now,time,dt)

    ! Simple output routine.

    ! Output format:
    ! The five output streams (see setup_output)

    !PhotonCounts: time
    !              Number of (ionizations + recombinations) / photons 
    !                   during time step
    !              Number of ionizations /(ionizations + recombinations)
    !              Number of recombinations /(ionizations + recombinations)
    !              Number of (ionizations + recombinations) / photons 
    !              Number of (ionizations + recombinations) / photons 
    !                   since t=0

    use sizes, only: mesh
    use grid, only: x, vol
    use material, only: xh, temper, ndens
    use evolve, only: phih_grid
    use sourceprops, only: srcpos, NormFlux, NumSrc
    use photonstatistics, only: do_photonstatistics, total_ion, totrec, totcollisions, dh0, grtotal_ion, photon_loss
    use radiation, only: teff,rstar,lstar,S_star

    real(kind=dp),intent(in) :: zred_now,time,dt

    integer :: i,j,k,ns
    character(len=6) :: zred_str
    character(len=40) :: file1,file2,file3,file4,file5,file6
    real(kind=dp) :: totalsrc
    real(kind=dp) :: totions,totphots,volfrac,massfrac
    logical crossing,recording_photonstats

    if (rank == 0) then
       ! Stream 1
       if (streams(1).eq.1) then
          write(file1,'(f6.3)') zred_now
          file1='Ifront1_'//trim(adjustl(file1))//'.dat'
          open(unit=51,file=file1,form='unformatted',status='unknown')
          do i=1,mesh(1)
             write(51,'(5(1pe10.3,1x))') x(i), &
                  xh(i,srcpos(2,1),srcpos(3,1),0), &
                  xh(i,srcpos(2,1),srcpos(3,1),1), &
                  temper, &
                  ndens(i,srcpos(2,1),srcpos(3,1))
          enddo
          close(51)
       endif
       
       ! Stream 2
       if (streams(2).eq.1) then
          write(file1,'(f6.3)') zred_now
          file1='xfrac3d_'//trim(adjustl(file1))//'.bin'
          open(unit=52,file=file1,form='unformatted',status='unknown')
          write(52) mesh(1),mesh(2),mesh(3)
          write(52) (((xh(i,j,k,1),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(52)
       endif
       
       ! Stream 3
!       if (streams(3).eq.1) then
!          write(file1,'(f6.3)') zred_now
!          file1='Temper3_'//trim(adjustl(file1))//'.bin'
!          open(unit=53,file=file1,form='unformatted',status='unknown')
!          write(53) mesh(1),mesh(2),mesh(3)
!          write(53) (((real(temper),i=1,mesh(1)),j=1,mesh(2)), &
!               k=1,mesh(3))
!          close(53)
!       endif
       
       ! Stream 3
       if (streams(3).eq.1) then
          write(file1,'(f6.3)') zred_now
          file1='IonRates3_'//trim(adjustl(file1))//'.bin'
          open(unit=53,file=file1,form='unformatted',status='unknown')
          write(53) mesh(1),mesh(2),mesh(3)
          write(53) (((real(phih_grid(i,j,k)),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(53)
       endif
       
       ! Stream 4
       if (streams(4).eq.1) then
          write(zred_str,'(f6.3)') zred_now
          file1='Ifront2_xy_'//trim(adjustl(zred_str))//'.bin'
          file2='Ifront2_xz_'//trim(adjustl(zred_str))//'.bin'
          file3='Ifront2_yz_'//trim(adjustl(zred_str))//'.bin'
          open(unit=54,file=file1,form='unformatted',status='unknown')
          open(unit=55,file=file2,form='unformatted',status='unknown')
          open(unit=56,file=file3,form='unformatted',status='unknown')
          ! xy cut through source 
          write(54) mesh(1),mesh(2)
          !        write(54) ((real(xh(i,j,srcpos(3,1),1)),i=1,mesh(1)),
          write(54) ((real(xh(i,j,mesh(3)/2,1)),i=1,mesh(1)), &
               j=1,mesh(2))
          ! xz cut through source 
          write(55) mesh(1),mesh(3)
          !        write(55) ((real(xh(i,srcpos(2,1),k,1)),i=1,mesh(1)),
          write(55) ((real(xh(i,mesh(2)/2,k,1)),i=1,mesh(1)), &
               k=1,mesh(3))
          ! yz cut through source 
          write(56) mesh(2),mesh(3)
          !        write(56) ((real(xh(srcpos(1,1),j,k,1)),j=1,mesh(2)),
          write(56) ((real(xh(mesh(1)/2,j,k,1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(54)
          close(55)
          close(56)
       endif
       
       ! Stream 5
       if (streams(5).eq.1) then
          write(zred_str,'(f6.3)') zred_now
          file4='ndens_xy_'//trim(adjustl(zred_str))//'.bin'
          file5='ndens_xz_'//trim(adjustl(zred_str))//'.bin'
          file6='ndens_yz_'//trim(adjustl(zred_str))//'.bin'
          open(unit=57,file=file4,form='unformatted',status='unknown')
          open(unit=58,file=file5,form='unformatted',status='unknown')
          open(unit=59,file=file6,form='unformatted',status='unknown')
          ! xy cut through source 
          write(57) mesh(1),mesh(2)
          !        write(57) ((real(ndens(i,j,srcpos(3,1))),i=1,mesh(1)),
          write(57) ((real(ndens(i,j,mesh(3)/2)),i=1,mesh(1)),j=1,mesh(2))
          ! xz cut through source 
          write(58) mesh(1),mesh(3)
          !        write(58) ((real(ndens(i,srcpos(2,1),k)),i=1,mesh(1)),
          write(58) ((real(ndens(i,mesh(2)/2,k)),i=1,mesh(1)),k=1,mesh(3))
          ! yz cut through source 
          write(59) mesh(2),mesh(3)
          !        write(59) ((real(ndens(srcpos(1,1),j,k)),j=1,mesh(2)),
          write(59) ((real(ndens(mesh(1)/2,j,k)),j=1,mesh(2)),k=1,mesh(3))
          close(57)
          close(58)
          close(59)
       endif
       
       ! Check if we are tracking photon conservation
       if (do_photonstatistics) then
          ! Photon Statistics
          total_ion=total_ion!+photon_loss
          totalsrc=sum(NormFlux(1:NumSrc))*s_star*dt
          grtotal_ion=grtotal_ion+total_ion-totcollisions
          grtotalsrc=grtotalsrc+totalsrc
          if (time.gt.0.0) then
             write(90,'(f6.3,6(1pe10.3))') &
                  zred_now, &
                  (total_ion-totcollisions)/totalsrc, &
                  dh0/total_ion, &
                  totrec/total_ion, &
                  photon_loss/total_ion, &
                  totcollisions/total_ion, &
                  grtotal_ion/grtotalsrc
          endif
          totions=sum(ndens(:,:,:)*xh(:,:,:,1))*vol
          volfrac=sum(xh(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          massfrac=sum(ndens(:,:,:)*xh(:,:,:,1))/sum(ndens)
          write(95,'(f6.3,4(1pe10.3))') zred_now,totions,grtotalsrc,volfrac,massfrac

       endif
    endif
    
    return
  end subroutine output

end module output_module
