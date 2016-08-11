!>
!! \brief This module contains data and routines for handling the Lyman-Werner
!! background. This is needed to calculate the suppression of mini-halo
!! sources
!! 
!! \b Author: Kyungjin Ahn, Garrelt Mellema
!!
!! \b Date: 23-June-2016 (imported from Kyungjin's version)
!!
!! \b Version: Derived from Kyungjin's version

module getjLW

  use precision, only: dp
  use my_mpi
  use astroconstants, only: M_SOLAR
  use cgsconstants, only: m_p
  use sizes, only: mesh
  use nbody, only: zred_array
  use cosmology, only: zred2time, zred_array_out
  use sourceprops, only: NumSrc, NumMassiveSrc, NumSupprbleSrc, NumSupprsdSrc, srcpos00, srcpos01, sM00_msun, sM01_msun
#ifdef MH
  use source_sub, only: NumAGrid, sub_srcpos, ssM_msun, MHflag
  ! Catch 22: source_sub needs jLW from getjLW and getjLW needs variables
  ! from source_sub
  ! Solution: new data module for MH source lists.
  use jLWgreen, only: greenK, read_greenK, get_rLW, rLW_zobs
  use radiation, only: jLW
#endif
  use c2ray_parameters, only: phot_per_atom, Ni, fstar
  use cosmology_parameters, only: Omega0, Omega_B, HcOm
  use file_admin, only: logf, results_dir

  implicit none
#ifdef IFORT
  ! ifort standard for "binary"
  character(len=*),parameter :: binaryformat="binary"
  character(len=*),parameter :: binaryaccess="sequential"
#else
  ! Fortran2003 standard for "binary"
  character(len=*),parameter :: binaryformat="unformatted"
  character(len=*),parameter :: binaryaccess="stream"
#endif

  ! LW emissivity for each type of source
  real(kind=dp),parameter :: emis00  = 1.67d21  !< LW emissivity (erg s^-1 Hz^-1 Msun^-1) of high-mass source
  real(kind=dp),parameter :: emis01  = 3d21   !< LW emissivity (erg s^-1 Hz^-1 Msun^-1) of low-mass source
  real(kind=dp),parameter :: emissub = 3d21   !< LW emissivity (erg s^-1 Hz^-1 Msun^-1) of subgrid source

  real(kind=dp) :: CC00   !< conversion from difference between C2Ray lifetime and true lifetime of high-mass source
  real(kind=dp) :: QH_M_C2ray00
!  real(kind=dp),parameter :: QH_M_real00 = 10d0**46.8d0
  real(kind=dp),parameter :: QH_M_real00 = 6.309573445d46 ! = 10^46.8

  real(kind=dp) :: CC01   !< conversion from difference between C2Ray lifetime and true lifetime of low-mass source
  real(kind=dp) :: QH_M_C2ray01
  real(kind=dp),parameter :: QH_M_real01 = 1.2d48

  real(kind=dp) :: CCsub   !< conversion from difference between C2Ray lifetime and true lifetime of subgrid source
  real(kind=dp) :: QH_M_C2Raysub
  real(kind=dp),parameter :: QH_M_realsub = 1.2d48

  real(kind=dp),    dimension(:,:,:),allocatable                   :: srclum
  complex(kind=dp), dimension(mesh(1)/2+1,mesh(2),mesh(3)), public :: srclumK

  real(kind=dp),    dimension(:,:,:),allocatable                   :: jLWtemp
  complex(kind=dp), dimension(:,:,:),allocatable                   :: gK_sK

  real(kind=dp) :: jLWaver

#ifdef MPI
    integer :: mympierror
#endif

contains

#ifdef MH
! ======================================================================

  subroutine read_jLW(zred_now)

    ! Read jLW from a file
    ! This is only done in case of a restart.

    real(kind=dp),intent(in) :: zred_now

    character(len=512) :: filej
    integer            :: m1,m2,m3

    ! I/O is done by rank 0
    ! Note: Kyungjin's version seemed to have all nodes reading!
    if (rank == 0) then

       ! Construct file name from zred_now
       write(filej,"(f6.3)") zred_now
       filej=trim(adjustl(results_dir))//"jLW3d_"//trim(adjustl(filej))//".bin"
       
       ! Open file and read data
       open(unit=33, file=filej, form=binaryformat, access=binaryaccess, &
            status="old")
       read(33) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with jLW unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
          stop
       endif
       
       write(logf,*) '**************************************'
       write(logf,"(A,f6.3)") 'reading jLW file at z=', zred_now 
       write(logf,*) '**************************************'
       write(6,*) '**************************************'
       write(6,"(A,f6.3)") 'reading jLW file at z=', zred_now 
       write(6,*) '**************************************'

       read(33) jLW
       close(33)
       
    endif

#ifdef MPI
    ! Share the values of jLW with the other nodes
    CALL MPI_BCAST(jLW, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_NEW, mympierror)
#endif
    
  end subroutine read_jLW

  ! ======================================================================

  subroutine get_jLW (zred,nz_LW)

    ! This routines finds jLW for the current redshift zred. This redshift
    ! should correspond to zred_array_out(nz_LW).

    include 'fftw3.f'

    integer, intent(in) :: nz_LW !< index of redshift in zred_array_out
    real(kind=dp), intent(in) :: zred !< redshift

    integer :: nzz ! loop counter for loop over redshift list, for LW calc.
    integer :: nz_LWbegin ! index of starting redshift for LW calculation.
    integer :: n_past ! loop counter for loop over past slices for LW.
    integer :: n_pastslice ! number of past slices needed for LW calculation.
    integer :: ii, jj, kk
    real(kind=dp)       :: zobs
    real(kind=dp)       :: zstart

    ! for FFT
    integer(kind=8)          :: PLAN

    if (rank == 0) then

       ! Calculate srclumK, the Fourier transform of the LW luminosity field.
       ! We calculate this for the previous redshift as it depends on the
       ! source properties which have not yet been calculated.
       call get_srclumK(zred_array_out(nz_LW-1))
       ! Write this to a file for future use
       call output_srclumK(zred_array_out(nz_LW-1))

       ! Determine lookback redshift for LW field calculation
       zobs = zred ! should be equal to zred_array_out(nz_out)
       ! This routine provides the LW horizon (rLW_zobs) for zobs
       call get_rLW(zobs)
       ! Find the highest redshift which we need to look back to.
       ! This assumes a matter dominated universe.
       zstart = ((1d0+zobs)**(-0.5d0) - rLW_zobs*0.5d0*HcOm)**(-2d0) - 1d0
       ! Find which number this corresponds to
       if (zstart > zred_array_out(1)) then
          nz_LWbegin = 1
       else
          do nzz = 1, nz_LW ! or nz_LW-1 ?????
             if (zred_array_out(nzz) >= zstart .and. zstart > zred_array_out(nzz+1)) then
                nz_LWbegin = nzz
             endif
          enddo
       endif

       ! required number of past slices to consider for LW calculation
       n_pastslice  = nz_LW - nz_LWbegin

       write(logf,*) 'n_pastslice, nz_lwbegin', n_pastslice, nz_lwbegin
       write(6,*) 'n_pastslice, nz_lwbegin', n_pastslice, nz_lwbegin

       ! jLw: sum over contributions from old slices
       jLW     =  0d0
       n_past  = 0
       do nzz = nz_LWbegin, nz_LW-1
          n_past    = n_past + 1
          ! Read in the precalculated k-space Green's function
          call read_greenK(zred_array_out(nzz), zred_array_out(nzz+1), &
               zred_array_out(nz_LW))
          ! Read in the MH source properties
          call read_srclumK(zred_array_out(nzz))

          ! Allocate arrays
          if (allocated(gK_sK)) deallocate(gK_sK)
          allocate(gK_sK(mesh(1)/2+1,mesh(2),mesh(3)))

          if (allocated(jLWtemp)) deallocate(jLWtemp)
          allocate(jLWtemp(mesh(1),mesh(2),mesh(3)))

          ! convolution theorem
          gK_sK = greenK * srclumK

          ! FFT gK_sK to get jLWtemp
          call DFFTW_PLAN_DFT_C2R_3D(PLAN,MESH(1),MESH(2),MESH(3), &
               gK_sK, jLWtemp, FFTW_ESTIMATE)
          call DFFTW_EXECUTE_DFT_C2R(PLAN, gK_sK, jLWtemp)
          call DFFTW_DESTROY_PLAN(PLAN)

          jLW   = jLW + jLWtemp

          deallocate(gK_sK)
          deallocate(jLWtemp)
       enddo

    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute jLW to the other nodes
    CALL MPI_BCAST(jLW, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_NEW, mympierror)
#endif

    ! Record average value to log file
    IF (rank == 0) THEN
       jLWaver = sum(jLW)/dble(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) "jLW volume average", jLWaver
    ENDIF

    return

  end subroutine get_jLW

! ======================================================================

  subroutine get_srclumK(nz_LW)

    include 'fftw3.f'
    
    integer, intent(in) :: nz_LW
    integer             :: n00, n01, nsub
    integer(kind=8)     :: PLAN
    real(kind=dp)       :: C2Ray_lifetime
    real(kind=dp)       :: coeff00, coeff01
    real(kind=dp)       :: coeffsub

    ! for MPI
    integer             :: NPROCS, MYRANK, IERR



    ! Allocate the LW luminosity grid
    if (allocated(srclum)) deallocate(srclum)
    allocate(srclum(mesh(1),mesh(2),mesh(3)))
    srclum(:,:,:) = 0d0

    ! Add LW contribution of LMACHs and HMACHs to LW luminosity grid
    do n00 = 1, NumSrc
       srclum(srcpos00(1, n00), srcpos00(2, n00), srcpos00(3, n00)) = &
            srclum(srcpos00(1, n00), srcpos00(2, n00), srcpos00(3, n00)) + &
            NormFlux_LW(n00)
    enddo

    ! LW contribution of MHs to LW luminosity grid
    do nsub = 1, NumAGrid
       srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub))= &
            srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub)) &
            + ssM_msun(nsub) * coeffsub
    enddo

    ! FFT srclum into srclumK
    call DFFTW_PLAN_DFT_R2C_3D(PLAN,MESH(1),MESH(2),MESH(3), &
         srclum,srclumK, FFTW_ESTIMATE)
    call DFFTW_EXECUTE_DFT_R2C(PLAN,srclum,srclumK)
    call DFFTW_DESTROY_PLAN(PLAN)
    ! srclum is no longer needed
    deallocate(srclum)

    ! normalize k-space component to follow H&E eqtns (6-96), (6-97)
    srclumK = srclumK /dble(mesh(1)*mesh(2)*mesh(3))
    
  end subroutine get_srclumK

! ======================================================================

  subroutine output_srclumK(nz_LW)

    ! Write the FT of the LW luminosity field to a file for future use.

    integer, intent(in) :: nz_LW
    character(len=6)    :: z_str
    character(len=512)  :: dir_sK, fname

    ! Directory for file
    dir_sK = results_dir !"../results/"

    ! Construct file name
    write(z_str,'(f6.3)') zred_array_out(nz_LW)
    fname = trim(adjustl(dir_sK))//trim(adjustl(z_str))//"-srcK"

    ! Open file
    open(unit=10, file=fname, form=binaryformat, access=binaryaccess, &
         status='unknown')
    ! Write data in SM3D format
    write(10) mesh(1), mesh(2), mesh(3)
    write(10) srclumK
    ! Close file
    close(10)

  end subroutine output_srclumK

! ======================================================================

  subroutine read_srclumK(zsrc)

    ! Read the FT of the LW luminosity field at redshift zsrc from a file 

    real(kind=dp), intent(in) :: zsrc !< redshift
    character(len=6)          :: z_str
    character(len=512)        :: dir_sK, fname
    integer                   :: m1,m2,m3

    ! Directory for file
    dir_sK = results_dir !"../results/"

    ! Construct file name
    write(z_str,'(f6.3)') zsrc
    fname = trim(adjustl(dir_sK))//trim(adjustl(z_str))//"-srcK"

    ! Open file
    open(unit=10, file=fname, form=binaryformat, access=binaryaccess, &
         status='old')
    ! Read data (SM3D format)
    read(10) m1, m2, m3
    if (m1 /= mesh(1) .or. m2 /= mesh(2) .or. m3 /= mesh(3)) then
       write(6,*) 'mesh number not matched for k-space source luminosity'
       close(10)
       stop
    endif
    read(10) srclumK
    ! Close file
    close(10)

  end subroutine read_srclumK

#endif
end module getjLW
