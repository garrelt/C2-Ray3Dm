module jLWgreen
!
  use precision, only: dp
  use sizes, only: mesh
  use mathconstants, only: pi
  use cgsconstants, only: c
  use cosmology_parameters, only: Omega0, h, HcOm
  use astroconstants, only: Mpc
  use nbody, only: boxsize  ! in Mpc/h

  implicit none

  real(kind=dp)        :: rLW_zobs !< past LW horizon to be used for green.

  complex(kind=dp), dimension(mesh(1)/2+1,mesh(2),mesh(3)) :: greenK !< Fourier component of green

contains
! ==========================================================================
  subroutine read_greenK(zsbegin, zsend, zobs)

    ! This routine reads in precalculated k-space green function in 
    ! mesh1/2+1, mesh2, mesh3 dimension in 8-byte complex variables.
    !
    ! --- How k-space green function is precalculated: --------------------
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
    ! ---------------------------------------------------------------------
    ! Update: KA (12-Mar-2010)
    !

    ! beginning source zred; ending source zred; observing zred
    real(kind=dp),intent(in) :: zsbegin, zsend, zobs 
    ! string corresponding to above redshifts
    character(len=6)         :: zb_str, ze_str, zo_str
    character(len=512)       :: dir_gK, fname
    integer                  :: m1, m2, m3

    dir_gK = "../GKresults/"

    write(zb_str, '(f6.3f)') zsbegin
    write(ze_str, '(f6.3f)') zsend
    write(zo_str, '(f6.3f)') zobs

    fname=trim(adjustl(dir_gK))//"gK_"//trim(adjustl(zb_str))//"-"//trim(adjustl(ze_str))//"-"//trim(adjustl(zo_str))//"_dat"

    open(unit=15, file=fname, form='binary', status='old')
    read(15) m1, m2, m3
    if (m1 /= mesh(1) .or. m2 /= mesh(2) .or. m3 /= mesh(3)) then
       write(6,*) 'mesh number not matched for greenK!! Aborting!'
       close(15)
       stop
    endif
    read(15) greenK
    close(15)

    return
    
  end subroutine read_greenK

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

end module jLWgreen
