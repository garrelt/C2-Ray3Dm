!>
!! \brief This module contains routines for file output

module output_module
  
  ! This file contains routines having to do with the output
  ! of Ifront programs.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, results_dir, file_input, logf

  implicit none
  
  private

  integer,parameter :: max_input_streams=5 !< maximum number of input streams
  integer,dimension(max_input_streams) :: streams !< flag array for input streams

  public :: setup_output, output, close_down

contains
  !----------------------------------------------------------------------------

  !> Initializes output streams
  subroutine setup_output ()
    
    ! Sets up output stream
    
    ! Version: Five streams
    ! Stream1:
    ! Ifront1.out contains a line of constant y and z going through the
    ! source for all timesteps. (formatted)
    
    ! Stream2: 
    ! Ionization fractions for the full data cube (unformatted)
    ! Ifront3",f5.3,".bin"
    
    ! Stream3: 
    ! Temperature for the full data cube (unformatted)
    ! temper3",f5.3,".bin"
    
    ! Stream 4:
    ! Ionization fractions in a plane for one time step
    ! Ifront2_xy_",f5.3,".bin"
    ! Ifront2_xz_",f5.3,".bin"
    ! Ifront2_yz_",f5.3,".bin"
    
    ! Stream 5:
    ! Densities in a plane for one time step
    ! ndens_xy_",f5.3,".bin"
    ! ndens_xz_",f5.3,".bin"
    ! ndens_yz_",f5.3,".bin"
    
    ! photon statistics
    use photonstatistics, only: do_photonstatistics, &
         initialize_photonstatistics

    if (rank == 0) then
       if (.not.file_input) then
          write(*,*) "Which output streams do you want?"
          write(*,*) "Enter a mask for streams 1 through 5"
          write(*,*) "E.g. 1,0,1,0,0 means streams 1 and 3, but not 2, 4, 5"
       endif
       read(stdinput,*) streams(1),streams(2),streams(3),streams(4),streams(5)
       
       ! Open files
       if (do_photonstatistics) then
          !open(unit=90,file=trim(adjustl(results_dir))//"PhotonCounts.out", &
          !     form="formatted",status="unknown",access="append")
          ! Fortran 2003 standard
          open(unit=90,file=trim(adjustl(results_dir))//"PhotonCounts.out", &
               form="formatted",status="unknown",position="append")
          write(90,*) "redshift, photon conservation number, ", &
               "fraction new ionization, fraction recombinations, ", &
               "fraction photon losses, fraction collisional ionization, ", &
               "grand total photon conservation number"
          !open(unit=95,file=trim(adjustl(results_dir))//"PhotonCounts2.out", &
          !     form="formatted",status="unknown",access="append")
          ! Fortran 2003 standard
          open(unit=95,file=trim(adjustl(results_dir))//"PhotonCounts2.out", &
               form="formatted",status="unknown",position="append")
          write(95,*) "redshift, total number of ions, ", &
               "grand total ionizing photons, mean ionization fraction ", &
               "(by volume and mass)"
       endif
#ifdef MPILOG
       write(logf,*) "Making output streams according to: ", &
            streams(1),streams(2),streams(3),streams(4),streams(5)
#endif
    endif
    if (do_photonstatistics) then
       call initialize_photonstatistics ()
    endif

#ifdef MPILOG
    write(logf,*) "End of setup output"
#endif

  end subroutine setup_output
  
  !-----------------------------------------------------------------------------

  !> Closes down output streams
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

  !> Produce output for a time frame
  subroutine output(zred_now,time,dt,photcons_flag)

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
    use photonstatistics, only: do_photonstatistics, total_ion, totrec, totcollisions, dh0, grtotal_ion, photon_loss, LLS_loss, grtotal_src
    use radiation, only: teff,rstar,lstar,S_star

    real(kind=dp),intent(in) :: zred_now,time,dt
    integer,intent(out) :: photcons_flag

    integer :: i,j,k,ns
    character(len=6) :: zred_str
    character(len=40) :: file1,file2,file3,file4,file5,file6
    real(kind=dp) :: totalsrc,photcons,total_photon_loss,total_LLS_loss
    real(kind=dp) :: totions,totphots,volfrac,massfrac
    logical crossing,recording_photonstats

#ifdef MPI
    integer :: mympierror
#endif

    ! Set photon conservation flag to zero on all processors
    photcons_flag=0

    ! Only produce output on rank 0
    if (rank == 0) then
       ! Stream 1
       if (streams(1).eq.1) then
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"Ifront1_"//trim(adjustl(file1))//".dat"
          open(unit=51,file=file1,form="formatted",status="unknown")
          do i=1,mesh(1)
             write(51,"(5(1pe10.3,1x))") x(i), &
                  xh(i,srcpos(2,1),srcpos(3,1),0), &
                  xh(i,srcpos(2,1),srcpos(3,1),1), &
                  !1.0_dp-xh(i,srcpos(2,1),srcpos(3,1)), &
                  !xh(i,srcpos(2,1),srcpos(3,1)), &
                  temper, &
                  ndens(i,srcpos(2,1),srcpos(3,1))
          enddo
          close(51)
       endif
       
       ! Stream 2
       if (streams(2).eq.1) then
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"xfrac3d_"//trim(adjustl(file1))//".bin"
          open(unit=52,file=file1,form="unformatted",status="unknown")
          write(52) mesh(1),mesh(2),mesh(3)
          write(52) (((xh(i,j,k,1),i=1,mesh(1)),j=1,mesh(2)),k=1,mesh(3))
          !write(52) xh
          close(52)
       endif
       
       ! Stream 3
!       if (streams(3).eq.1) then
!          write(file1,"(f6.3)") zred_now
!          file1="Temper3_"//trim(adjustl(file1))//".bin"
!          open(unit=53,file=file1,form="unformatted",status="unknown")
!          write(53) mesh(1),mesh(2),mesh(3)
!          write(53) (((real(temper),i=1,mesh(1)),j=1,mesh(2)), &
!               k=1,mesh(3))
!          close(53)
!       endif
       
       ! Stream 3
       if (streams(3).eq.1) then
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"IonRates3_"//trim(adjustl(file1))//".bin"
          open(unit=53,file=file1,form="unformatted",status="unknown")
          write(53) mesh(1),mesh(2),mesh(3)
          write(53) (((real(phih_grid(i,j,k)),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(53)
       endif
       
       ! Stream 4
       if (streams(4).eq.1) then
          write(zred_str,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"Ifront2_xy_"//trim(adjustl(zred_str))//".bin"
          file2=trim(adjustl(results_dir))//"Ifront2_xz_"//trim(adjustl(zred_str))//".bin"
          file3=trim(adjustl(results_dir))//"Ifront2_yz_"//trim(adjustl(zred_str))//".bin"
          open(unit=54,file=file1,form="unformatted",status="unknown")
          open(unit=55,file=file2,form="unformatted",status="unknown")
          open(unit=56,file=file3,form="unformatted",status="unknown")
          ! xy cut through source 
          write(54) mesh(1),mesh(2)
          !        write(54) ((real(xh(i,j,srcpos(3,1),1)),i=1,mesh(1)),
          write(54) ((real(xh(i,j,mesh(3)/2,1)),i=1,mesh(1)), &
          !write(54) ((real(xh(i,j,mesh(3)/2)),i=1,mesh(1)), &
               j=1,mesh(2))
          ! xz cut through source 
          write(55) mesh(1),mesh(3)
          !        write(55) ((real(xh(i,srcpos(2,1),k,1)),i=1,mesh(1)),
          write(55) ((real(xh(i,mesh(2)/2,k,1)),i=1,mesh(1)), &
          !write(55) ((real(xh(i,mesh(2)/2,k)),i=1,mesh(1)), &
               k=1,mesh(3))
          ! yz cut through source 
          write(56) mesh(2),mesh(3)
          !        write(56) ((real(xh(srcpos(1,1),j,k,1)),j=1,mesh(2)),
          write(56) ((real(xh(mesh(1)/2,j,k,1)),j=1,mesh(2)), &
          !write(56) ((real(xh(mesh(1)/2,j,k)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(54)
          close(55)
          close(56)
       endif
       
       ! Stream 5
       if (streams(5).eq.1) then
          write(zred_str,"(f6.3)") zred_now
          file4=trim(adjustl(results_dir))//"ndens_xy_"//trim(adjustl(zred_str))//".bin"
          file5=trim(adjustl(results_dir))//"ndens_xz_"//trim(adjustl(zred_str))//".bin"
          file6=trim(adjustl(results_dir))//"ndens_yz_"//trim(adjustl(zred_str))//".bin"
          open(unit=57,file=file4,form="unformatted",status="unknown")
          open(unit=58,file=file5,form="unformatted",status="unknown")
          open(unit=59,file=file6,form="unformatted",status="unknown")
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
          ! total_ion is the total number of new ionization, plus
          ! the total number of recombinations, and here is also
          ! added the number of photons lost from the grid. Since
          ! this number was divided by the number of cells, we
          ! multiply by this again.
          total_photon_loss=sum(photon_loss)*dt* &
               real(mesh(1))*real(mesh(2))*real(mesh(3))
          total_LLS_loss = LLS_loss*dt
          !total_ion=total_ion + total_photon_loss
          totalsrc=sum(NormFlux(1:NumSrc))*s_star*dt
          photcons=(total_ion+LLS_loss-totcollisions)/totalsrc
          if (time > 0.0) then
             write(90,"(f6.3,9(1pe10.3))") &
                  zred_now, &
                  total_ion, totalsrc, &
                  photcons, &
                  dh0/total_ion, &
                  totrec/total_ion, &
                  total_LLS_loss/totalsrc, &
                  total_photon_loss/totalsrc, &
                  totcollisions/total_ion, &
                  grtotal_ion/grtotal_src
          endif
          totions=sum(ndens(:,:,:)*xh(:,:,:,1))*vol
          !totions=sum(ndens(:,:,:)*xh(:,:,:))*vol
          volfrac=sum(xh(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          !volfrac=sum(xh(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
          massfrac=sum(ndens(:,:,:)*xh(:,:,:,1))/sum(real(ndens,dp))
          !massfrac=sum(ndens(:,:,:)*xh(:,:,:))/sum(real(ndens,dp))
          write(95,"(f6.3,4(1pe10.3))") zred_now,totions,grtotal_src, &
               volfrac,massfrac

          if (abs(1.0-photcons) > 0.15) then
             if ((1.0-photcons) > 0.15 .and. &
                  total_photon_loss/totalsrc < (1.0-photcons) ) then
                photcons_flag=1
                ! Report photon conservation
                write(logf,"(A,2(1pe10.3,x))") &
                     "Photon conservation problem: ", &
                     photcons, total_photon_loss/totalsrc

             endif
          endif
       endif
    endif
    
#ifdef MPI
     call MPI_BCAST(photcons_flag,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine output

end module output_module
