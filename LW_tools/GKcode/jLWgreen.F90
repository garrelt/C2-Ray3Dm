module jLWgreen
!
  use precision, only: dp
  use my_mpi
  use sizes, only: mesh
  use mathconstants, only: pi
  use cgsconstants, only: c
  use cosmology_parameters, only: Omega0, h
  use astroconstants, only: Mpc
  use readinput, only: boxsize  ! in Mpc/h

  implicit none

  real(kind=dp)        :: rLW_zobs !< past LW horizon to be used for green.
  real(kind=dp),public :: HcOm     !< H/c*Omega0**0.5 in Mpc^-1, for rLW_zobs.

  real(kind=dp),    dimension(mesh(1),mesh(2),mesh(3))     :: green !< periodic spatial green function for LW flux
  complex(kind=dp), dimension(mesh(1)/2+1,mesh(2),mesh(3)) :: greenK !< Fourier component of green


contains
! ==========================================================================
  subroutine jLW_green(zsbegin, zsend, zobs)

    !     
    ! This calculates the Fourier component of Green function, G_klm, where
    ! k,l,m are integers that represent wavenumbers in x,y,z, repectively. 
    !    
    ! First the spatial Green ftn is constructed by periodic spacing of a
    ! unit source. First take about 1/8 of the spatial box, get the Green
    ! ftn there, and assign the Green ftn on the rest (7/8) of the box.
    ! Even inside this 1/8 of the box, we use the symmetry to reduce the
    ! number of operations maximally (1/6). Altogether, 1/48 * mesh**3
    ! operations needed from each source.
    !
    ! Then this is FFTed to create Green ftn in the k-space, namely G_klm.
    !
    ! zsbegin corresponds to some of N-body snapshot redshifts. As the source
    ! luminosity is assumed to be constant during each adjacent N-body
    ! redshifts, zsbegin is also identical to the "beginning" redshift of
    ! newly created sources, and zsend the "end". This is observed at zobs.
    ! zsend can be equal to zobs, but not smaller(later in time) than zobs.
    !
    ! WARNING: rLW, dtauMpc_sm, & dtauMpc_bg use approximate expressions, 
    !          valid only when the universe is at high z (matter dominated). 
    !          Pls contact the author for correction, when universe at z~1
    !          should be treated.
    !
    ! Author: Kyungjin Ahn (31-May-2009)
    ! Update: 
    !

    include 'fftw3.f'

    ! beginning source zred; ending source zred; observing zred
    real(kind=dp),intent(in) :: zsbegin, zsend, zobs 
    ! conformal lookback times to zsend and zsbegin from zobs.
    real(kind=dp)            :: dtauMpc_sm, dtauMpc_bg
    ! boxsize in Mpc (without h); cellsize in Mpc (without h)
    real(kind=dp)            :: boxMpc, d_boxMpc ! boxsize in Mpc (without h)
    ! luminosity distance (cm); (source redshift)+1; some constant (Mpc^-1); 
    ! rLW scaling factor; source to observer comoving distance (Mpc)
    real(kind=dp)            :: DL, zs1, scfacomm, scfactor, ros
    integer,dimension(:,:),allocatable :: ijksrc ! periodic src locations
    real(kind=dp)            :: zobs1, zs1_zobs1
    integer                  :: i, j, k, iextend, isrc

    real(kind=dp)            :: greenijk
    integer                  :: iswap, jswap, kswap

    ! for FFT
    integer(kind=8)          :: PLAN

#ifdef MPI
    ! for MPI
    integer                  :: NPROCS, MYRANK, IERR, IJK
#endif

    integer                  :: ISTA, IEND

!    real(kind=dp), dimension(mesh(1)*mesh(2)*mesh(3))     :: green1D
!    real(kind=dp), dimension(mesh(1)*mesh(2)*mesh(3))     :: green1Dsummed
    real(kind=dp), dimension(:),allocatable     :: green1D
    real(kind=dp), dimension(:),allocatable     :: green1Dsummed

#ifdef MPI
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)
#endif
    scfacomm = (h/0.7d0)**(-1d0) * (Omega0/0.27d0)**(-0.5d0) ! 

    boxMpc   = boxsize/h

    d_boxMpc = boxMpc/mesh(1) !works only when the box is a cube.

    zobs1      = 1d0 + zobs
    dtauMpc_sm = 2d0/HcOm * (zobs1**(-0.5d0) - (1d0+zsend)**(-0.5d0))
    dtauMpc_bg = 2d0/HcOm * (zobs1**(-0.5d0) - (1d0+zsbegin)**(-0.5d0))

    ! periodicity
    call get_rLW(zobs)
    iextend  = min(int(rLW_zobs/boxMpc)+1, int(dtauMpc_bg/boxMpc)+1)
    if (allocated(ijksrc)) deallocate(ijksrc)
    allocate(ijksrc((2*iextend+2)**3,3))

    ! src number & its position, of course coming from periodicity.
    isrc = 0
    do k=-iextend, iextend+1
       do j=-iextend, iextend+1
          do i=-iextend, iextend+1
             isrc           = isrc+1
             ijksrc(isrc,1) = i * mesh(1) + 1
             ijksrc(isrc,2) = j * mesh(2) + 1
             ijksrc(isrc,3) = k * mesh(3) + 1
          enddo
       enddo
    enddo

    ! spatial Green function
    green(:,:,:) = 0d0  ! Absolutely necessary, especially for MPI jobs.

#ifdef MPI
    ! parallelize over the periodic sources
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!!$    call PARA_RANGE(1, (2*iextend+2)**3, NPROCS, MYRANK, ISTA, IEND)
    ! Don't waste too much time communicating, setting 64 as maximum.
    call PARA_RANGE(1, (2*iextend+2)**3, min(64,NPROCS), MYRANK, ISTA, IEND)
#else
    ISTA = 1
    IEND = (2*iextend+2)**3
#endif
    ! Thread specific calculation. Value of 'green' varies over threads.
    do isrc = ISTA, IEND
       ! 1/6 of 1/8 of the whole box
       do k=1,mesh(3)/2+1
          do j=1,k
             do i=1,j
                ros  = d_boxMpc * &
                     sqrt(dble( (i-ijksrc(isrc,1))**2+(j-ijksrc(isrc,2))**2 &
                     +(k-ijksrc(isrc,3))**2 ) )
                ! region of effect is ros=[dtauMpc_sm, dtauMpc_bg]
                if (dtauMpc_sm <= ros .and. ros <= dtauMpc_bg) then
                   zs1  = ((1d0+zobs)**(-0.5d0) - ros*0.5d0*HcOm)**(-2d0)
                   zs1_zobs1 = zs1/zobs1
                   DL        = ros/zobs1 * zs1_zobs1
                   scfactor  = scfacomm * (zs1/21d0)**(-0.5d0)
                   if (ros < 0.1d0*d_boxMpc) then
                      green(i,j,k) = green(i,j,k) + 3d0/(4d0*pi*(d_boxMpc/zobs1)**2)
                   else
                      green(i,j,k) = green(i,j,k) + 1d0/(4d0*pi*DL*DL) &
                           * zs1_zobs1 * fmod(ros, scfactor)
                   endif
                endif

             enddo
          enddo
       enddo
    enddo

    ! Now calculate 1/8 of the whole box
    do k=1,mesh(3)/2+1
       do j=1,k
          do i=1,j
             greenijk     = green(i,j,k)
             green(i,k,j) = greenijk
             green(j,i,k) = greenijk
             green(j,k,i) = greenijk
             green(k,i,j) = greenijk
             green(k,j,i) = greenijk
          enddo
       enddo
    enddo

    ! Now calculate the whole box
    ! Dodge those indices which are out of range.
    do k=1,      mesh(3)/2+1
       kswap          = mesh(3)-k+2
       do j=1,   mesh(2)/2+1
          jswap       = mesh(2)-j+2
          do i=1,mesh(1)/2+1
             iswap    = mesh(1)-i+2
             greenijk = green(i,j,k)
             if       (kswap <= mesh(3)) then
                green      (i    ,     j, kswap) = greenijk
                if    (iswap <= mesh(1)) then
                   green   (iswap,     j, kswap) = greenijk
                   if (jswap <= mesh(2)) then
                      green(iswap, jswap, kswap) = greenijk
                   endif
                endif
                if    (jswap <= mesh(2)) then
                   green   (i    , jswap, kswap) = greenijk
                endif
             endif


             if (iswap <= mesh(1)) then
                green   (iswap, j    , k    ) = greenijk
                if (jswap <= mesh(2)) then
                   green(iswap, jswap, k    ) = greenijk
                endif
             endif
             if (jswap <= mesh(2)) then
                green   (i    , jswap, k    ) = greenijk
             endif

          enddo
       enddo
    enddo
             

    ! Normalize the green function: This is the correct normalization
    ! for a green function in the discrete Fourier transform, which 
    ! follows the convention of Hockney & Eastwood (see eqtn 6-96).
    ! H&E convention: X(p) = \sum_{k} X(k) exp(2*pi*i*p*k/N), where
    ! p & k represent spatial & k-space indices, repectively.
    ! This convention enforces the following Green ftn convetion:
    ! phi(p) = 1/N * \sum_{p'} G(p,p') g(p'), and convolution:
    ! phi(k) = G(k) * g(k). 1/N here should be 1/N^3 in 3D case.
    ! After getting phi(k), phi(p) = FFTWbackward(phi(k)), simply.
    ! green is still in Mpc^-2 unit. See below for unit change.
    green = green * dble(mesh(1)*mesh(2)*mesh(3))

    ! Change the unit into cm^-2 ster^-1
    green = green / Mpc**2 /(4d0*pi)

    
#ifdef MPI
    if (allocated(green1D))       deallocate(green1D)
    if (allocated(green1Dsummed)) deallocate(green1Dsummed)

    allocate(green1D      (mesh(1)*mesh(2)*mesh(3)))
    allocate(green1Dsummed(mesh(1)*mesh(2)*mesh(3)))

    IJK        = 0
    green1D(:) = 0d0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             IJK = IJK + 1
             green1D(IJK) = green(i,j,k)
          enddo
       enddo
    enddo

    ! perform MPI sum into green1Dsummed
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    CALL MPI_REDUCE(green1D, green1Dsummed, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)

    deallocate(green1D)

    ! replace green with the summed value
    if (MYRANK == 0) then
       IJK = 0
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                IJK = IJK + 1
                green(i,j,k) = green1Dsummed(IJK)
             enddo
          enddo
       enddo
!!$       call outputvtk_jLWgreen(zobs)
    endif
    
    deallocate(green1Dsummed)

    CALL MPI_BCAST(green, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
#endif

#ifdef MPI
    IF (MYRANK == 0) then
#endif
    call DFFTW_PLAN_DFT_R2C_3D(PLAN,MESH(1),MESH(2),MESH(3), &
         green,greenK, FFTW_ESTIMATE)
    call DFFTW_EXECUTE_DFT_R2C(PLAN,green,greenK)
    call DFFTW_DESTROY_PLAN(PLAN)
    ! normalize k-space component to follow H&E eqtns (6-96), (6-97)
    greenK = greenK /dble(mesh(1)*mesh(2)*mesh(3))

#ifdef MPI
    ENDIF
    CALL MPI_BCAST(greenK, (mesh(1)/2+1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, IERR)
#endif

    return

  end subroutine jLW_green

! =========================================================================
  function fmod(ros, scfactor)

    ! picket-fence modulation factor
    ! 
    ! Author: Kyungjin Ahn (27-May-2009)
    ! Update: 
    !
    real(kind=dp)             :: fmod
    real(kind=dp), intent(in) :: ros, scfactor

    if (ros/scfactor <= 97.39d0) then
       fmod = 1.7d0*exp(-(ros/(116.29d0*scfactor))**0.68d0) - 0.7d0
    else
       fmod = 0d0
    endif

  end function fmod
    
! =========================================================================
  subroutine get_rLW(zobs)

    ! get past LW horizon from redshift zobs
    ! comoving Mpc over which Lyman delta redshifted into Lyman gamma
    ! From eqtn 8 in Ahn et al. 2009 (eqtn 8 is faulty in the paper!!!!!)
    ! 
    ! Author: Kyungjin Ahn (27-May-2009)
    ! Update: 
    !

    real(kind=dp),intent(in) :: zobs 

    rLW_zobs = 2d0/HcOm * (zobs+1d0)**(-0.5d0) * &
         (1d0 - ((15d0/16d0)/(8d0/9d0))**(-0.5d0))

    return

  end subroutine get_rLW

! =========================================================================
  subroutine get_HcOm

    ! Simple constant calculation. In unit of Mpc^-1.

    HcOm = h*100d0*1d5 /c *Omega0**0.5d0

    return

  end subroutine get_HcOm

!===================================================================
  subroutine outputvtk_jLWgreen(zsbegin, zsend, zobs)
    ! plot green function of jLW field

    use file_admin, only : results_dir

    real(kind=dp),intent(in) :: zsbegin, zsend, zobs 
    character(len=6)         :: zb_str, ze_str, zo_str
    character(len=300)       :: fileout
    integer                  :: i, j

    write(zb_str,"(f6.3)") zsbegin
    write(ze_str,"(f6.3)") zsend
    write(zo_str,"(f6.3)") zobs
    fileout=trim(adjustl(results_dir))//"jLWgr_"//trim(adjustl(zb_str))//"-"//trim(adjustl(ze_str))//"-"//trim(adjustl(zo_str))//".vtk"
    
    open(unit=60,file=fileout,status="unknown")
222 format(A)
    write(60,222)  "# vtk DataFile Version 2.0"
    write(60,222)  "Volume example"
    write(60,222)  "ASCII"
    write(60,222)  "DATASET STRUCTURED_POINTS"
    write(60,223)  "DIMENSIONS ", mesh(1), mesh(2), 1
223 format(A,3I5)
    write(60,222)  "ORIGIN 0. 0. 0."
    write(60,222)  "SPACING 1. 1. 1."
    write(60,224)  "POINT_DATA ", mesh(1)*mesh(2)
224 format(A,I10)
    write(60,222)  ""
    write(60,222) "SCALARS jLWgreen float"
    write(60,222) "LOOKUP_TABLE default"
    
    do j=   1,mesh(2)
       do i=1,mesh(1)
          write(60,*) real(green(i,j,1))
       enddo
    enddo
    close(60)
    
  end subroutine outputvtk_jLWgreen


!===================================================================
  subroutine output_jLWgreenK(zsbegin, zsend, zobs)
    ! plot green function of jLW field

    use file_admin, only : results_dir

    real(kind=dp),intent(in) :: zsbegin, zsend, zobs 
    character(len=6)         :: zb_str, ze_str, zo_str
    character(len=300)       :: fileout

    write(zb_str,"(f6.3)") zsbegin
    write(ze_str,"(f6.3)") zsend
    write(zo_str,"(f6.3)") zobs
    fileout=trim(adjustl(results_dir))//"gK_"//trim(adjustl(zb_str))//"-"//trim(adjustl(ze_str))//"-"//trim(adjustl(zo_str))//"_dat"
    
    open(unit=60,file=fileout,form='binary',status="unknown")

    write(60) mesh(1), mesh(2), mesh(3)
    write(60) greenK
    close(60)
    
  end subroutine output_jLWgreenK
  
  
end module jLWgreen
