!>
!! \brief This module contains routines for file output
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2014-07-22 (but older)
!!
!! \b Version: 3D, hydrogen only

module output_module
  
  ! This file contains routines having to do with the output
  ! of the C2-Ray program.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, results_dir, file_input, logf
  use c2ray_parameters, only: isothermal
  use sizes, only: mesh
  use grid, only: x, vol
  use material, only: xh, ndens
  use material, only: temper_val, temperature_grid, temperature_states
  use material, only: temperature_states_dbl, get_temperature_point
  use evolve, only: phih_grid
  use sourceprops, only: srcpos, NormFlux, NumSrc
  use photonstatistics, only: initialize_photonstatistics
  use photonstatistics, only: do_photonstatistics, total_ion, totrec
  use photonstatistics, only: totcollisions, dh0, grtotal_ion, photon_loss
  use photonstatistics, only: LLS_loss, grtotal_src
  use radiation, only: teff,rstar,lstar,S_star
    

  implicit none
  
  private

  ! To controle what kind of output is produced we use an array of flags.
  ! There can be at most max_output_streams types of output.
  integer,parameter :: max_output_streams=5 !< maximum number of output streams
  integer,dimension(max_output_streams) :: streams !< flag array for output streams

  public :: setup_output, output, close_down

contains
  !----------------------------------------------------------------------------

  !> Initializes output streams
  subroutine setup_output ()
    
    ! Sets up output stream
    
    ! Version: Five streams
    ! Stream1:
    ! Ifront1.out contains a line of constant y and z going through the
    ! centre of the grid for all timesteps. (formatted)
    
    ! Stream2: 
    ! Ionization fractions for the full data cube (unformatted)
    ! xfrac3d_",f5.3,".bin"
    
    ! Stream3: 
    ! Ionization rate for the full data cube (unformatted)
    ! "Ionrates3_",f5.3,".bin"
    ! If non-isothermal: Temperature for the full data cube (unformatted)
    ! "Temper3d_",f5.3,".bin"
    
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
    
    if (rank == 0) then
       if (.not.file_input) then
          write(*,*) "Which output streams do you want?"
          write(*,*) "Enter a mask for streams 1 through 5"
          write(*,*) "E.g. 1,0,1,0,0 means streams 1 and 3, but not 2, 4, 5"
       endif
       read(stdinput,*) streams(1),streams(2),streams(3),streams(4),streams(5)
       
#ifdef MPILOG
       write(logf,*) "Making output streams according to: ", &
            streams(1),streams(2),streams(3),streams(4),streams(5)
#endif
    endif

    ! Initialize photon statistics (open files and initialize some quantities)
    if (do_photonstatistics) then
       if (rank == 0) call open_photonstatistics_files()
       call initialize_photonstatistics ()
    endif

#ifdef MPILOG
    write(logf,*) "End of setup output"
#endif

  end subroutine setup_output
  
  !-----------------------------------------------------------------------------

  !> Closes down output streams
  subroutine close_down ()
    
    ! Rank 0 takes care of file i/o
    if (rank == 0) then
       ! Close any open output files
       ! There are none
       if (do_photonstatistics) then
          ! Close the photon statistics files
          close(90)
          close(95)
       endif
    endif

  end subroutine close_down
  
  !----------------------------------------------------------------------------

  subroutine open_photonstatistics_files ()

    ! Open files
    if (do_photonstatistics .and. rank == 0) then

       ! Open file 1
       open(unit=90,file=trim(adjustl(results_dir))//"PhotonCounts.out", &
            form="formatted",status="unknown",position="append")
       ! Write header (content information)
       write(90,*) "Columns: redshift, ", &
            "total number of photons used on the grid, ", &
            "total number of photons produced on the grid, ",&
            "photon conservation number, ", &
            "fraction new ionization, fraction recombinations, ", &
            "fraction LLS losses (seems to be wrong), ", &
            "fraction photon losses, fraction collisional ionization, ", &
            "grand total photon conservation number"

       ! Open file 2
       open(unit=95,file=trim(adjustl(results_dir))//"PhotonCounts2.out", &
            form="formatted",status="unknown",position="append")
       ! Write header (content information)
       write(95,*) "Columns: redshift, total number of ions, ", &
            "grand total ionizing photons, mean ionization fraction ", &
            "(by volume and mass)"
    endif

  end subroutine open_photonstatistics_files

  !----------------------------------------------------------------------------

  !> Produce output for a time frame
  subroutine output(zred_now,time,dt,photcons_flag)

    ! Simple output routine.

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag
#ifdef MPILOG     
     write(logf,*) 'output 1'
#endif 
    if (streams(1) == 1) call write_stream1 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 2'
#endif 
    if (streams(2) == 1) call write_stream2 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 3'
#endif 
    if (streams(3) == 1) call write_stream3 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 4'
#endif 
    if (streams(4) == 1) call write_stream4 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 5'
#endif 
    if (streams(5) == 1) call write_stream5 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 6'
#endif 
    call write_photonstatistics (zred_now,time,dt,photcons_flag)

#ifdef MPILOG     
     write(logf,*) 'output 7'
#endif 

   end subroutine output

  !----------------------------------------------------------------------------

   !> produces output for a time frame. See below for format
   subroutine write_stream1 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift
     
    character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".dat"
    character(len=512) :: file1
    integer :: i,j,k
    type(temperature_states_dbl) :: temperature_point
    real,dimension(mesh(1)) :: temperature_profile

    ! Only produce output on rank 0
    if (rank == 0) then

       ! Stream 1
       if (streams(1) == 1) then
          ! Construct file name
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//basename//trim(adjustl(file1)) &
               //base_extension
          ! Open file
          open(unit=51,file=file1,form="formatted",status="unknown")

          ! Get temperature profile
          do i=1,mesh(1)
             call get_temperature_point(i,srcpos(2,1),srcpos(3,1), &
                  temperature_point)
             temperature_profile(i)=temperature_point%current
          enddo

          do i=1,mesh(1)
             write(51,"(5(es10.3,1x))") x(i), &
#ifdef ALLFRAC
                  xh(i,srcpos(2,1),srcpos(3,1),0), &
                  xh(i,srcpos(2,1),srcpos(3,1),1), &
#else
                  1.0_dp-xh(i,srcpos(2,1),srcpos(3,1)), &
                  xh(i,srcpos(2,1),srcpos(3,1)), &
#endif
                  temperature_profile(i), &
                  ndens(i,srcpos(2,1),srcpos(3,1))
          enddo
          close(51)
       endif
       
    endif
    
  end subroutine write_stream1
  
  !> produces output for a time frame. See below for format
  subroutine write_stream2 (zred_now)
    
    real(kind=dp),intent(in) :: zred_now !< current redshift
    
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=512) :: file1
    integer :: i,j,k
    
    if (rank == 0) then
       ! Stream 2
       if (streams(2) == 1) then
          ! Construct file name
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))// &
               "xfrac3d_"//trim(adjustl(file1))//base_extension
          ! Open file
          open(unit=52,file=file1,form="unformatted",status="unknown")
          ! Write data header
          write(52) mesh(1),mesh(2),mesh(3)
          ! Write data (depending on -DALLFRAC compile option)
#ifdef ALLFRAC
          write(52) (((xh(i,j,k,1),i=1,mesh(1)),j=1,mesh(2)),k=1,mesh(3))
#else
          write(52) xh
#endif
          close(52)
       endif
    endif
       
  end subroutine write_stream2

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream3 (zred_now)
    
    real(kind=dp),intent(in) :: zred_now !< current redshift
    
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=512) :: file1
    integer :: i,j,k
    
    if (rank == 0) then
       ! Stream 3
       if (streams(3) == 1) then
          ! Construct filename
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))// &
               "IonRates3_"//trim(adjustl(file1))//".bin"
          ! Open file
          open(unit=53,file=file1,form="unformatted",status="unknown")
          ! Write data header
          write(53) mesh(1),mesh(2),mesh(3)
          ! Write data
          write(53) (((real(phih_grid(i,j,k)),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(53)
          
          if (.not.isothermal) then
             file1="Temper3d_"//trim(adjustl(file1))//".bin"
             open(unit=53,file=file1,form="unformatted",status="unknown")
             write(53) mesh(1),mesh(2),mesh(3)
             write(53) ((((temperature_grid%current),i=1,mesh(1)),j=1,mesh(2)), &
                  k=1,mesh(3))
             close(53)
             
             ! Should add heating rate to this output stream
          endif
       endif
    endif

  end subroutine write_stream3

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream4 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=512) :: file1,file2,file3
    character(len=6) :: zred_str
    integer :: i,j,k

    if (rank == 0) then
       ! Stream 4
       if (streams(4) == 1) then
          write(zred_str,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))// &
               "Ifront2_xy_"//trim(adjustl(zred_str))//".bin"
          file2=trim(adjustl(results_dir))// &
               "Ifront2_xz_"//trim(adjustl(zred_str))//".bin"
          file3=trim(adjustl(results_dir))// &
               "Ifront2_yz_"//trim(adjustl(zred_str))//".bin"
          open(unit=54,file=file1,form="unformatted",status="unknown")
          open(unit=55,file=file2,form="unformatted",status="unknown")
          open(unit=56,file=file3,form="unformatted",status="unknown")
          ! xy cut through source 
          write(54) mesh(1),mesh(2)
#ifdef ALLFRAC
          write(54) ((real(xh(i,j,mesh(3)/2,1)),i=1,mesh(1)), &
               j=1,mesh(2))
#else
          write(54) ((real(xh(i,j,mesh(3)/2)),i=1,mesh(1)), &
               j=1,mesh(2))
#endif
          close(54)
          
          ! xz cut through source 
          write(55) mesh(1),mesh(3)
#ifdef ALLFRAC
          write(55) ((real(xh(i,mesh(2)/2,k,1)),i=1,mesh(1)), &
               k=1,mesh(3))
#else
          write(55) ((real(xh(i,mesh(2)/2,k)),i=1,mesh(1)), &
               k=1,mesh(3))
#endif
          close(55)
          
          ! yz cut through source 
          write(56) mesh(2),mesh(3)
#ifdef ALLFRAC
          write(56) ((real(xh(mesh(1)/2,j,k,1)),j=1,mesh(2)), &
               k=1,mesh(3))
#else
          write(56) ((real(xh(mesh(1)/2,j,k)),j=1,mesh(2)), &
               k=1,mesh(3))
#endif
          close(56)
       endif
    endif

  end subroutine write_stream4

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream5 (zred_now)
    
    real(kind=dp),intent(in) :: zred_now !< current redshift
    
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=512) :: file1,file2,file3
    character(len=6) :: zred_str
    integer :: i,j,k

    if (rank == 0) then
       ! Stream 5
       if (streams(5) == 1) then
          write(zred_str,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))// &
               "ndens_xy_"//trim(adjustl(zred_str))//".bin"
          file2=trim(adjustl(results_dir))// &
               "ndens_xz_"//trim(adjustl(zred_str))//".bin"
          file3=trim(adjustl(results_dir))// &
               "ndens_yz_"//trim(adjustl(zred_str))//".bin"
          
          open(unit=57,file=file1,form="unformatted",status="unknown")
          ! xy cut through source 
          write(57) mesh(1),mesh(2)
          !        write(57) ((real(ndens(i,j,srcpos(3,1))),i=1,mesh(1)),
          write(57) ((real(ndens(i,j,mesh(3)/2)),i=1,mesh(1)),j=1,mesh(2))
          close(57)
          
          ! xz cut through source 
          open(unit=58,file=file2,form="unformatted",status="unknown")
          write(58) mesh(1),mesh(3)
          !        write(58) ((real(ndens(i,srcpos(2,1),k)),i=1,mesh(1)),
          write(58) ((real(ndens(i,mesh(2)/2,k)),i=1,mesh(1)),k=1,mesh(3))
          close(58)
          
          ! yz cut through source 
          open(unit=59,file=file3,form="unformatted",status="unknown")
          write(59) mesh(2),mesh(3)
          !        write(59) ((real(ndens(srcpos(1,1),j,k)),j=1,mesh(2)),
          write(59) ((real(ndens(mesh(1)/2,j,k)),j=1,mesh(2)),k=1,mesh(3))
          close(59)
       endif
    endif

  end subroutine write_stream5

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_photonstatistics (zred_now,time,dt,photcons_flag)
    
    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag
    
    real(kind=dp) :: totalsrc,photcons,total_photon_loss
    real(kind=dp) :: total_LLS_loss
    real(kind=dp) :: totions,totphots,volfrac,massfrac
    
#ifdef MPI
    integer :: mympierror
#endif
    
    if (rank == 0) then
 
       ! Check if we are tracking photon conservation
       if (do_photonstatistics) then
          ! Photon Statistics
          ! Calculate quantities for PhotonCounts.dat:
          ! total_photon_loss: the number of photons lost from the grid. 
          !    Since this number was divided by the number of cells, we
          !    multiply by this again.
          total_photon_loss=sum(photon_loss)*dt* &
               real(mesh(1))*real(mesh(2))*real(mesh(3))
          ! total_LLS_loss: the number of photons lost due to LLSs. This
          !    number does not appear to be calculated correctly
          total_LLS_loss = LLS_loss*dt
          ! total_src: total number of photons produced by sources
          totalsrc=sum(NormFlux(1:NumSrc))*s_star*dt
          ! photcons: photon conservation number. total_ion is the total 
          !    number of new ionizations plus the total number of 
          !    recombinations. We subtract the total number of ionizations
          !    due to collisions. If photon conservation is perfect
          !    total_ion-totcollisions should equal the number of photons
          !    produced by sources. The number will be <1 if there are
          !    photon losses over the boundary of the grid or due to LLS.
          !    In the multiple sources case even without these losses the
          !    number will be only approximately 1.
          photcons=(total_ion-totcollisions)/totalsrc
          !photcons=(total_ion+LLS_loss-totcollisions)/totalsrc

          write(logf,*) "Change in output: ",total_ion
          write(logf,*) "Rates in output: ",totrec,totcollisions
          
          ! Write PhotonCounts.dat
          if (time > 0.0) then
             write(90,"(f6.3,9(es10.3))") &
                  zred_now, &
                  total_ion, totalsrc, &
                  photcons, &
                  dh0/total_ion, &
                  totrec/total_ion, &
                  total_LLS_loss/totalsrc, &
                  total_photon_loss/totalsrc, &
                  totcollisions/total_ion, &
                  grtotal_ion/grtotal_src
             flush(90) ! force writing of output
          endif

          ! Calculate quantities for PhotonCounts2.dat:
          ! totions: total number of ions in volume
          ! volfrac: Global average ionized fraction (by volume)
          ! massfrac: Global average ionized fraction (by mass)
#ifdef ALLFRAC
          totions=sum(ndens(:,:,:)*xh(:,:,:,1))*vol
          volfrac=sum(xh(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          massfrac=sum(ndens(:,:,:)*xh(:,:,:,1))/sum(real(ndens,dp))
#else
          totions=sum(ndens(:,:,:)*xh(:,:,:))*vol
          volfrac=sum(xh(:,:,:))/real(mesh(1)*mesh(2)*mesh(3))
          massfrac=sum(ndens(:,:,:)*xh(:,:,:))/sum(real(ndens,dp))
#endif
          ! Write PhotonCounts2.dat
          write(95,"(f6.3,4(es10.3))") zred_now,totions,grtotal_src, &
               volfrac,massfrac
          flush(95) ! force writing of output
          
          ! Report and flag for non-conservations of photons.
          ! This flag will only have an effect if stop_on_photon_violation
          ! in c2ray_parameters is .true.
          if (abs(1.0-photcons) > 0.15) then
             if ((1.0-photcons) > 0.15 .and. &
                  total_photon_loss/totalsrc < (1.0-photcons) ) then
                photcons_flag=1
                ! Report photon conservation
                write(logf,"(A,2(es10.3,x))") &
                     "Photon conservation problem: ", &
                     photcons, total_photon_loss/totalsrc

             endif
          endif
       endif
    endif
    
#ifdef MPI
    call MPI_BCAST(photcons_flag,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine write_photonstatistics

end module output_module


