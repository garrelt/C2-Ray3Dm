!>
!! \brief Main program for C2Ray-3Dm (3D, multiple sources)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 22-May-2008
!<
Program GreenK

  ! Authors: Kyungjin Ahn

  ! Date: 1-Mar-2010

  ! Precalculates k-space green function and dump it.

  use precision, only: dp
  use file_admin, only: flag_for_file_input
  use my_mpi !, only: mpi_setup, mpi_end, rank
  use readinput
  use jLWgreen, only: jLW_green, get_HcOm, get_rLW, rLW_zobs, HcOm, output_jLWgreenK

  implicit none

  integer       :: nz0, nz, nzz, nz_LWbegin
  real(kind=dp) :: zsbegin, zsend, zobs, zstart

#ifdef MPI
  call mpi_setup()
#endif

  if (iargc() > 0) then
     if (rank == 0) call flag_for_file_input(.true.)
  endif
    
  call read_input(nz0)

  ! set the constant HcOm for LW calculation
  call get_HcOm

  ! Determine lookback redshift for LW field calculation
  do nz=nz0+1,NumZred
     zobs = zred_array(nz)
     call get_rLW(zobs)
     zstart = ((1d0+zobs)**(-0.5d0) - rLW_zobs*0.5d0*HcOm)**(-2d0) - 1d0
     if (rank == 0) then
        write(*,'(A,f6.3,A,f6.3,A,e14.7)'), 'zobs = ', zobs, ' zstart = ', zstart, ' rLW=', rLW_zobs
     endif
     if (zstart > zred_array(nz0)) then
        nz_LWbegin = nz0
     else
        do nzz = nz0, nz
           if (zred_array(nzz) >= zstart .and. zstart > zred_array(nzz+1)) then
              nz_LWbegin = nzz
           endif
        enddo
     endif
     do nzz = nz_LWbegin, nz-1
        zsbegin = zred_array(nzz)
        zsend   = zred_array(nzz+1)
        call jLW_green       (zsbegin, zsend, zobs)
        if (rank == 0) call output_jLWgreenK(zsbegin, zsend, zobs)
     enddo
  enddo

end Program GreenK
