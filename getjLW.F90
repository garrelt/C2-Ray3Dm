module getjLW
!
  use precision, only: dp
  use my_mpi
  use astroconstants, only: M_SOLAR
  use cgsconstants, only: m_p
  use sizes, only: mesh
  use nbody, only: zred_array
  use cosmology, only: zred2time 
  use sourceprops, only: NumSrc, NumMassiveSrc, NumSupprbleSrc, NumSupprsdSrc, srcpos00, srcpos01, sM00_msun, sM01_msun
  use source_sub, only: NumAGrid, sub_srcpos, ssM_msun, MHflag
  use jLWgreen, only: greenK, read_greenK, get_rLW, rLW_zobs
  use c2ray_parameters, only: phot_per_atom, Ni, fstar
  use cosmology_parameters, only: Omega0, Omega_B, HcOm
  use file_admin, only: logf, results_dir
  
!  use nbody, only: boxsize  ! in Mpc/h  <--- retrieve later

  implicit none

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

  real(kind=dp),    dimension(mesh(1),mesh(2),mesh(3)), public     :: jLW 
  real(kind=dp),    dimension(:,:,:),allocatable                   :: jLWtemp
  complex(kind=dp), dimension(:,:,:),allocatable                   :: gK_sK

  real(kind=dp) :: jLWaver

contains
! ======================================================================
  subroutine read_jLW(zred_now)

    real(kind=dp),intent(in) :: zred_now

    character(len=512) :: filej
    integer            :: m1,m2,m3

    write(filej,"(f6.3)") zred_now
    filej=trim(adjustl(results_dir))//"jLW3d_"//trim(adjustl(filej))//".bin"
    open(unit=33,file=filej,form="binary",status="old")
    read(33) m1,m2,m3
    if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
       write(logf,*) "Warning: file with jLW unusable"
       write(logf,*) "mesh found in file: ",m1,m2,m3
       stop
    endif

    if (rank == 0) then
       write(logf,*) '**************************************'
       write(logf,"(A,f6.3)") 'reading jLW file at z=', zred_now 
       write(logf,*) '**************************************'
       write(6,*) '**************************************'
       write(6,"(A,f6.3)") 'reading jLW file at z=', zred_now 
       write(6,*) '**************************************'
    endif

    read(33) jLW
    close(33)


  end subroutine read_jLW

! ======================================================================
  subroutine get_jLW(nz0, nz)

    include 'fftw3.f'

    integer, intent(in) :: nz0, nz
    integer :: nzz ! loop counter for loop over redshift list, for LW calc.
    integer :: nz_LWbegin ! index of starting redshift for LW calculation.
    integer :: n_past ! loop counter for loop over past slices for LW.
    integer :: n_pastslice ! number of past slices needed for LW calculation.
    integer :: ii, jj, kk
    real(kind=dp)       :: zobs
    real(kind=dp)       :: zstart

    ! for FFT
    integer(kind=8)          :: PLAN

#ifdef MPI
    ! for MPI
    integer             :: NPROCS, MYRANK, IERR

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    IF (MYRANK == 0) THEN
#endif

       ! Get srclumK, the Fourier component of src-luminosity distribution
       call get_srclumK(nz)
       call output_srclumK(nz)

       ! Determine lookback redshift for LW field calculation
       zobs = zred_array(nz+1)
       call get_rLW(zobs)
       zstart = ((1d0+zobs)**(-0.5d0) - rLW_zobs*0.5d0*HcOm)**(-2d0) - 1d0
       if (zstart > zred_array(1)) then
          nz_LWbegin = 1
       else
          do nzz = 1, nz
             if (zred_array(nzz) >= zstart .and. zstart > zred_array(nzz+1)) then
                nz_LWbegin = nzz
             endif
          enddo
       endif
       ! required number of past slices to consider for LW calculation
       n_pastslice  = nz+1 - nz_LWbegin

       write(logf,*) 'n_pastslice, nz_lwbegin', n_pastslice, nz_lwbegin
       write(6,*) 'n_pastslice, nz_lwbegin', n_pastslice, nz_lwbegin

       ! jLw: sum over contributions from old slices
       jLW     =  0d0
       n_past  = 0
       do nzz = nz_LWbegin, nz
          n_past    = n_past + 1
          call read_greenK(zred_array(nzz), zred_array(nzz+1), zred_array(nz+1))
          call read_srclumK(zred_array(nzz))
          ! convolution theorem
          if (allocated(gK_sK)) deallocate(gK_sK)
          allocate(gK_sK(mesh(1)/2+1,mesh(2),mesh(3)))

          if (allocated(jLWtemp)) deallocate(jLWtemp)
          allocate(jLWtemp(mesh(1),mesh(2),mesh(3)))

          gK_sK = greenK * srclumK

          call DFFTW_PLAN_DFT_C2R_3D(PLAN,MESH(1),MESH(2),MESH(3), &
               gK_sK, jLWtemp, FFTW_ESTIMATE)
          call DFFTW_EXECUTE_DFT_C2R(PLAN, gK_sK, jLWtemp)
          call DFFTW_DESTROY_PLAN(PLAN)

          jLW   = jLW + jLWtemp

          deallocate(gK_sK)
          deallocate(jLWtemp)
       enddo
#ifdef MPI
    ENDIF

    CALL MPI_BCAST(jLW, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

    IF (MYRANK == 0) THEN
       jLWaver = sum(jLW)/dble(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) "jLW volume average", jLWaver
    ENDIF
#endif

    return

  end subroutine get_jLW

! ======================================================================
  subroutine get_srclumK(nz)

    include 'fftw3.f'
    
    integer, intent(in) :: nz
    integer             :: n00, n01, nsub
    integer(kind=8)     :: PLAN
    real(kind=dp)       :: C2Ray_lifetime
    real(kind=dp)       :: coeff00, coeff01
    real(kind=dp)       :: coeffsub

    ! for MPI
    integer             :: NPROCS, MYRANK, IERR


    C2ray_lifetime = abs(zred2time(zred_array(nz)) - zred2time(zred_array(nz+1)))

    QH_M_C2ray00  = Ni(1) * (M_solar/m_p)  /C2ray_lifetime
    QH_M_C2ray01  = Ni(2) * (M_solar/m_p)  /C2ray_lifetime

    CC00          = QH_M_C2ray00 / QH_M_real00
    CC01          = QH_M_C2ray01 / QH_M_real01

    coeff00       = emis00 * CC00 * fstar(1) * Omega_B/Omega0
    coeff01       = emis01 * CC01 * fstar(2) * Omega_B/Omega0

    write(6,*) 'Sanity check: fesc00 = ', phot_per_atom(1) /(Ni(1) *fstar(1) )
    write(6,*) 'Sanity check: fesc01 = ', phot_per_atom(2) /(Ni(2) *fstar(2) )

    if (allocated(srclum)) deallocate(srclum)
    allocate(srclum(mesh(1),mesh(2),mesh(3)))
    srclum(:,:,:) = 0d0

    do n00 = 1, NumMassiveSrc
       srclum(srcpos00(1, n00), srcpos00(2, n00), srcpos00(3, n00)) = &
            srclum(srcpos00(1, n00), srcpos00(2, n00), srcpos00(3, n00)) + &
            sM00_msun(n00) * coeff00
    enddo
    do n01 = 1, NumSupprbleSrc-NumSupprsdSrc
       srclum(srcpos01(1, n01), srcpos01(2, n01), srcpos01(3, n01)) = &
            srclum(srcpos01(1, n01), srcpos01(2, n01), srcpos01(3, n01)) + &
            sM01_msun(n01) * coeff01
    enddo
    ! Different subgrid source population schemes reflected by MHflag
    ! MHflag = 1: uniform star formation out of given baryon content
    !          2: one Pop III star/ one minihalo
    ! Below ssM_msun differs for different MHflags. Also note difference
    ! in fstar(3). ssM_msun is total subgrid stellar mass for MHflag=2.
    if (MHflag == 1) then
       QH_M_C2raysub  = Ni(3) * (M_solar/m_p)  /C2ray_lifetime
       CCsub          = QH_M_C2raysub / QH_M_realsub
       coeffsub       = emissub * CCsub * fstar(3) * Omega_B/Omega0
       write(6,*) 'Sanity check: fesc_sub = ', phot_per_atom(3) /(Ni(3) *fstar(3))
       do nsub = 1, NumAGrid
          srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub))= &
               srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub)) + ssM_msun(nsub) * coeffsub
       enddo
    elseif (MHflag == 2) then
       QH_M_C2raysub  = Ni(3) * (M_solar/m_p)  /C2ray_lifetime
       CCsub          = QH_M_C2raysub / QH_M_realsub
       ! Check out the difference from MHflag=1 case!!!
       coeffsub       = emissub * CCsub 
       write(6,*) 'Sanity check: fesc_sub = ', phot_per_atom(3) /(Ni(3) *fstar(3))

       do nsub = 1, NumAGrid
          srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub))= &
               srclum(sub_srcpos(1,nsub), sub_srcpos(2,nsub), sub_srcpos(3,nsub)) + ssM_msun(nsub) * coeffsub
       enddo
    endif

    call DFFTW_PLAN_DFT_R2C_3D(PLAN,MESH(1),MESH(2),MESH(3), &
         srclum,srclumK, FFTW_ESTIMATE)
    call DFFTW_EXECUTE_DFT_R2C(PLAN,srclum,srclumK)
    call DFFTW_DESTROY_PLAN(PLAN)
    deallocate(srclum)
    ! normalize k-space component to follow H&E eqtns (6-96), (6-97)
    srclumK = srclumK /dble(mesh(1)*mesh(2)*mesh(3))

    return

  end subroutine get_srclumK

! ======================================================================
  subroutine output_srclumK(nz)

    integer, intent(in) :: nz
    character(len=6)    :: z_str
    character(len=512)  :: dir_sK, fname

    dir_sK = "../results/"

    write(z_str,'(f6.3)') zred_array(nz)
    fname = trim(adjustl(dir_sK))//trim(adjustl(z_str))//"-srcK"
    open(unit=10, file=fname, form='binary', status='unknown')

    write(10) mesh(1), mesh(2), mesh(3)
    write(10) srclumK
    close(10)

  end subroutine output_srclumK

! ======================================================================
  subroutine read_srclumK(zsrc)

    real(kind=dp), intent(in) :: zsrc
    character(len=6)          :: z_str
    character(len=512)        :: dir_sK, fname
    integer                   :: m1,m2,m3

    dir_sK = "../results/"

    write(z_str,'(f6.3)') zsrc
    fname = trim(adjustl(dir_sK))//trim(adjustl(z_str))//"-srcK"
    open(unit=10, file=fname, form='binary', status='old')

    read(10) m1, m2, m3
    if (m1 .ne. mesh(1) .or. m2 .ne. mesh(2) .or. m3 .ne. mesh(3)) then
       write(6,*) 'mesh number not matched for k-space source luminosity'
       close(10)
       stop
    endif
    read(10) srclumK
    close(10)

  end subroutine read_srclumK

end module getjLW
