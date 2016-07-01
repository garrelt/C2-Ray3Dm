!>
!! \brief This module contains data and routines for cosmological problems
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! This module keeps track of the current redshift
!!

module cosmology

  ! This file contains routines having to do with the cosmological 
  ! evolution
  
  ! - cosmology_init: initializes cosmological time and sets lengths
  !            and volumes from comoving to proper scaling.
  ! - time2zred: convert time to redshift
  ! - zred2time: convert redshift to time
  ! - redshift_evol: calculate redshift from time, and zfactor between the
  !    current and previous time
  ! - cosmo_evol: cosmological evolution of space, density
  ! - cosmo_cool: cosmological adiabatic cooling rate
  ! - compton_cool: Compton cooling wrt the CMB.

  use precision, only: dp
  use my_mpi
  use cosmology_parameters, only: H0,Omega0,cmbtemp,cosmo_id
  use file_admin, only: logf
  use sizes, only: mesh
  use cgsconstants, only: c

  implicit none

  real(kind=dp) :: zred_t0 !< initial redshift
  real(kind=dp) :: t0      !< time of initial redshift
  real(kind=dp) :: zred    !< current redshift
  real(kind=dp) :: Hz      !< Hubble constant at current redshift
  real(kind=dp),private :: zfactor !< scaling factor between two redshifts
  !> list of output redshifts
  real(kind=dp),dimension(:),allocatable :: zred_array_out
  !> counter of output redshifts
  integer :: nz_out=0

contains
  ! =======================================================================

  !> initializes cosmological time and sets lengths and volumes from 
  !! comoving to proper scaling.
  subroutine cosmology_init (zred0,time)
    
    real(kind=dp),intent(in) :: zred0 !< initial redshift
    real(kind=dp),intent(in) :: time  !< initial time
    
    ! Report cosmological parameter set used
    if (rank == 0) then
       write(logf,*) "Cosmology used: ",trim(adjustl(cosmo_id))
    endif
    
    ! Cosmological time corresponding to (initial) redshift zred0
    ! NOTE: Good only for high-z!!!
    t0 = 2.*(1.+zred0)**(-1.5)/(3.*H0*sqrt(Omega0))
    
    ! Initialize redshift
    zred_t0=zred0 ! keep initial redshift
    zred=0.0 ! needs to be zero, so comoving will be changed to proper
    
    call initialize_zred_array_out()

    ! Initialize lengths and volumes to proper units
    call redshift_evol(time)
    call cosmo_evol( )
    
  end subroutine cosmology_init

  ! =======================================================================

  subroutine initialize_zred_array_out

    integer :: ierror=0 
    integer :: ioflag=0

    ! Allocate array for storing output redshifts (which may be more finely
    ! spaced than the zred_array from the N-body outputs).
    allocate(zred_array_out(number_timesteps*NumZred))

    ! If this is a restart there will be a file with previous
    ! values of output redshifts. This file needs to be read as
    ! we need to know these previous values for the LW lightcone

    if (rank == 0) then
       open(unit=71,file=trim(adjustl(results_dir))//"zred_out.dat", &
            status="old",iostat=ierror)
       if .not.(ierror) then
          do while (ioflag >= 0) ! EOF gives a -1 here. Some of the files have
             ! other errors, which should be ignored.
             ! Read in catalogues and count sources
             read(71,*,iostat=ioflag) zred_array_out(nz_out)
             nz_out=nz_out + 1
          enddo
          nz_out=nz_out-1 ! we counted one too many with the above procedure
       else
          write(logf,'(A)') "Error opening zred_out.dat"
       endif
    endif
    
#IFDEF mpi
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(nz_out,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(zred_array_out,number_timesteps*NumZred, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine initialize_zred_array_out

  ! =======================================================================

  !> Calculates the cosmological redshift for a given time
  real(kind=dp) function time2zred (time)

    ! Calculates the cosmological redshift for a given time

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 21-May-2005)
    
    ! Version: f90

    ! History: - 20-Aug-2006, conversion to f90

    real(kind=dp),intent(in) :: time

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    time2zred = -1+(1.+zred_t0)*(t0/(t0+time))**(2./3.)

  end function time2zred

  ! =======================================================================

  !> Calculates the time for a given cosmological redshift
  real(kind=dp) function zred2time (zred1)

    ! Calculates the time for a given cosmological redshift

    ! Author: Garrelt Mellema

    ! Date: 30-Sep-2006
    
    ! Version: f90

    ! History: 

    real(kind=dp),intent(in) :: zred1

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    zred2time = t0*( ((1.0+zred_t0)/(1.0+zred1))**1.5 - 1.0 )

  end function zred2time

  ! =======================================================================

  !> Calculates the cosmological redshift from time
  !! and the scale factor zfactor for use in cosmo_evol
  subroutine redshift_evol (time)

    ! Calculates the cosmological redshift from time
    ! and the scale factor zfactor for use in cosmo_evol
    
    ! Author: Garrelt Mellema
    ! Date: 10-Mar-2006
    ! Version: F90

    ! History:
    ! - 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: time

    real(kind=dp) :: zred_prev

    integer :: nz_loop

    ! Calculate the change since the previous redshift.
    ! Note: the initial redshift should be ZERO since
    ! the variables are initialized as comoving!
    ! NOTE: Good only for high-z!!!
    zred_prev=zred
    zred = -1+(1.+zred_t0)*((t0+time)/t0)**(-2./3.)

    ! Take the average zfactor between zred_prev and zred
    zfactor=(1.0+zred_prev)/(1.+zred)

    Hz=H0*(1.+zred)**(1.5)*sqrt(Omega0) ! Hubble constant at current redshift (cgs)

    ! Store redshift in list
    nz_out=nz_out+1
    zred_array_out(nz_out)=zred

    ! Write output redshift file
    ! This file should be overwritten every time.
    if (rank == 0) then
       open(unit=71,file=trim(adjustl(results_dir))//"zred_out.dat", &
            status="replace",iostat=ierror)
       do nz_loop=1,nz_out
          write(71,'(f6.3)') zred_array_out(nz_loop)
       enddo
       close(71)
    endif

  end subroutine redshift_evol

  ! =======================================================================

  !> Calculates the cosmological evolution of space and densities\n
  !! Changes variables: ndens, r, dr, vol
  subroutine cosmo_evol ()

    ! Calculates the cosmological evolution of space and densities

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: F90 first version

    ! History:
    ! - 19-Nov-2004: first version f77

    !use sizes
    use grid, only: x,y,z,dr,vol
    use material, only: ndens, n_LLS, y_LLS
    
    real(kind=dp) :: zfactor3

    zfactor3=zfactor*zfactor*zfactor

    ! Change the grid coordinates
    x(:)=x(:)*zfactor
    y(:)=y(:)*zfactor
    z(:)=z(:)*zfactor

    dr(:)=dr(:)*zfactor
    
    vol=vol*zfactor3

    ! Change the densities
    ndens(:,:,:)=ndens(:,:,:)/zfactor3

    ! Change the volumes (not needed???)
    ! volx(i,j,k)=volx(i,j,k)*zfactor3
    ! voly(i,j,k)=voly(i,j,k)*zfactor3
    ! volz(i,j,k)=volz(i,j,k)*zfactor3

    ! Evolution of LLS
    n_LLS=n_LLS * zfactor**(-y_LLS-1.5)

  end subroutine cosmo_evol
  
  ! =======================================================================

  !> Calculates the cosmological adiabatic cooling
  real(kind=dp) function cosmo_cool (e_int)

    ! Calculates the cosmological adiabatic cooling

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: f90

    ! History:
    ! 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: e_int

    real(kind=dp) :: dzdt

    ! Cosmological cooling rate:
    ! 2*(da/dt)/a*e
    ! or
    ! 2*(dz/dt)/(1+z)*e
    ! with a the cosmological scale factor

    ! dz/dt (for flat LambdaCDM)
    dzdt=H0*(1.+zred)*sqrt(Omega0*(1.+zred)**3+1.-Omega0)

    !Cooling rate
    cosmo_cool=e_int*2.0/(1.0+zred)*dzdt

  end function cosmo_cool

  ! =======================================================================

  !> Calculates the (cosmological) Compton cooling rate (against the CMB)
  real(kind=dp) function compton_cool (temper,eldens)
    
    ! Calculates the (cosmological) Compton cooling rate
    
    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: first version
    
    ! History:
    ! 16-May-2005: first version f77
    

    ! parameter reference?

    real(kind=dp),intent(in) :: temper ! temperature
    real(kind=dp),intent(in) :: eldens ! electron density
    
    !Cooling rate
    compton_cool=5.65e-36*eldens*(1.0+zred)**4* &
         (temper-cmbtemp*(1.0+zred))
    
  end function compton_cool
  
end module cosmology
