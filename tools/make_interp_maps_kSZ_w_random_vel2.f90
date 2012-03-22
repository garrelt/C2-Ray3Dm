program make_map

  use sizes, only: mesh
  use nbody
  use cosmology_parameters
  use cgsconstants
  use mathconstants

  implicit none

  integer :: nz, direction, ranx, rany, i,j,k
  real*8 :: zcurrent, z_old, zmiddle, Deltaz, weight, Hz
  real*8 :: meandeltaT,rmsdeltaT
  real*8 :: mean1,mean2,mean3,rms1,rms2,rms3
  character(len=6) :: zred_str
  character(len=180) :: file1, file2, file3
  character(len=180) :: redshift_file ! name of file with list of redshifts
  real(kind=8), dimension(100) :: zred_now
  !     assumes cubical mesh!!
  real,dimension(3,mesh(1),mesh(1)) :: DeltaT_real, DeltaT_old_real
  real,dimension(3,mesh(1),mesh(1)) :: tau_real, tau_old_real
  real*8,dimension(3,mesh(1),mesh(1)) :: DeltaT,DeltaT_old,DeltaTadd 
  real*8,dimension(3,mesh(1),mesh(1)) :: tau,tau_old,tauadd
  real*8,dimension(mesh(1),mesh(1)) :: DeltaTtot, DeltaTadd1 
  real*8,dimension(mesh(1),mesh(1)):: DeltaTadd2, DeltaTadd3
  real*8,dimension(mesh(1),mesh(1)):: tauadd1, tauadd2, tauadd3
  !     DeltaTtot would collect the data map
  !     assumes cubical mesh 
  real*8 :: vrandom, rantheta, vrms_box, vrms,ranampl
  integer :: rant, rana

!  integer :: NumZred
  integer,parameter :: MaxNumZred=200

  file2='partial_integrals_w_random_vel/output_redshifts.dat'
  open(100,file=file2)

  !     Ask for redshift file
  write(*,'(A,$)') 'File with redshifts: '
  read(*,*) redshift_file

  !     Open and read redshift file
  open(unit=60,file=redshift_file,form='formatted',status='old')
  read(60,*) NumZred
  if (NumZred > MaxNumZred) write(*,*)&
       'Your program is about to die in a horrible way'

  do nz=1,NumZred
     read(60,*) zred_now(nz)
  end do
  zcurrent=zred_now(1)      !start at first output redshift
  !      zcurrent=zred_now(2)      !start at second output redshift
  direction=1               !start with x-direction
  !      direction=2               !start with y-direction
  !      direction=3               !start with z-direction
  nz=1
  !     
  do
     !loop until overlap
     Hz=100*sqrt(Omega0*(1.+zcurrent)**3+(1-Omega0)) !in km/s/Mpc*h
     Deltaz=Hz*boxsize/(c/1e5) !c is in cm/s, boxsize in Mpc/h (comov)
     ! note: h's in Hz and boxsize cancel 
     z_old=zcurrent         !the 2 redshifts between which LOS is 1 boxsize
     zcurrent=zcurrent-Deltaz
     zmiddle=zcurrent+Deltaz/2. !halfway point
     !         print*,'z before loop check',z_old,zcurrent, Deltaz, zmiddle
     !     
     do nz=1,NumZred-1!between which output redshifts the halfway point lies?
        !            print*,'z before if check', zmiddle,zred_now(nz),zred_now(nz+1)
        if(zmiddle.lt.zred_now(nz) .and. zmiddle.ge.zred_now(nz+1))then 
           !               print*,'z check', zmiddle, zred_now(nz+1),zred_now(nz)
           call read_in(zred_now(nz),deltaT_old_real,tau_old_real)
           call read_in(zred_now(nz+1),deltaT_real,tau_real)
           deltaT_old=deltaT_old_real
           deltaT=deltaT_real
           tau_old=tau_old_real
           tau=tau_real
           weight=(zred_now(nz)-zmiddle)/(zred_now(nz)-zred_now(nz+1))   
           DeltaTadd = weight*DeltaT_old+(1.0d0-weight)*DeltaT
           tauadd = weight*tau_old+(1.0d0-weight)*tau
           mean1=sum(deltaTadd(1,:,:))/mesh(1)**2
           mean2=sum(deltaTadd(2,:,:))/mesh(1)**2
           mean3=sum(deltaTadd(3,:,:))/mesh(1)**2
           !               print*,'check means main',mean1,mean2,mean3

           rms1=sqrt(sum(deltaTadd(1,:,:)**2)/mesh(1)**2-mean1**2)
           rms2=sqrt(sum(deltaTadd(2,:,:)**2)/mesh(1)**2-mean2**2)
           rms3=sqrt(sum(deltaTadd(3,:,:)**2)/mesh(1)**2-mean3**2)

           !               print*,'check rms main',rms1,rms2,rms3


           !               print*,'check interp.',weight,sum(DeltaTadd)/mesh(1)**2
           !here we should switch directions (x->y->z) 
           !and translate box randomly
           direction=modulo(direction+1,3)
           !     print*,'check dir.',direction,mesh(1)
           call Random_int(ranx,0,mesh(1)-1)
           call Random_int(rany,0,mesh(1)-1)
           call Random_int(rant,0,100)
           call Random_int(rana,0,100)
           rantheta=real(rant)/100.*2.*pi
           ranampl=sqrt(-2.*log(real(rana)/100.))              
           print*,'check randoms',ranx,rany,rant,rana,rantheta,ranampl
           !               pause
           call vrms_singlez(zmiddle,boxsize,vrms,vrms_box)
           print*,'check random vel.',vrms,vrms_box,rantheta,ranampl
           !               pause
           !     vrandom=sqrt(vrms-vrms_box)*cos(rantheta)*1d5
           vrandom=sqrt(vrms-vrms_box)*ranampl*cos(rantheta)*1d5!Gaussian-distributed 
           !random number with zero mean and sigma=sqrt(vrms-vrms_box)
           !Box-Muller Transformation     
           !     ranx=0
           !     rany=0
           DeltaTadd1=DeltaTadd(direction,:,:)
           DeltaTadd2=cshift(DeltaTadd1,shift=ranx,dim=1)               
           DeltaTadd3=cshift(DeltaTadd2,shift=rany,dim=2)               
           tauadd1=tauadd(direction,:,:)
           tauadd2=cshift(tauadd1,shift=ranx,dim=1)               
           tauadd3=cshift(tauadd2,shift=rany,dim=2)            
           !     both velocities are in cm/s
           DeltaTtot=DeltaTtot + DeltaTadd3 + tauadd3*vrandom/c

           print*,'check vel.',vrandom/1e5,c/1e5
           !     save these partial integrals
           !     zmiddle is the redshift
           write(100,*) zmiddle
           write(zred_str,'(f6.3)') zmiddle
           file1='partial_integrals_w_random_vel/dTkSZ_'//trim(adjustl(zred_str))//'.bin'
           open(unit=54,file=file1,form='unformatted',status='unknown')
           write(54) mesh(1),mesh(1)
           write(54)((real(deltaTadd3(i,j)+tauadd3(i,j)*vrandom/c),i=1,mesh(1)),j=1,mesh(1))
           close(54)
           !               DeltaTtot=DeltaTtot + DeltaTadd1
           mean1=sum(deltaTadd(direction,:,:))/mesh(1)**2
           mean2=sum(deltaTadd1)/mesh(1)**2
           mean3=sum(deltaTtot)/mesh(1)**2
           rms1=sqrt(sum(deltaTadd(direction,:,:)**2)/mesh(1)**2-mean1**2)
           rms2=sqrt(sum(deltaTadd1**2)/mesh(1)**2-mean2**2)
           rms3=sqrt(sum(deltaTtot**2)/mesh(1)**2-mean3**2)
           !               print*,'check means mask',mean1,mean2,mean3
           !               print*,'check rms mask',rms1,rms2,rms3
        end if
     end do
     !         stop
     if(zcurrent.lt.zred_now(NumZred))go to 20
  end do

20 meandeltaT=sum(DeltaTtot)/mesh(1)**2
  rmsdeltaT=sqrt(sum(DeltaTtot**2)/mesh(1)**2-meandeltaT**2)
  print*,'means',meandeltaT,rmsdeltaT

  open(unit=3,file='result_w_random_vel.bin',form='unformatted',status='unknown')
  write(3) mesh(1),mesh(1)
  write(3)((real(DeltaTtot(i,j)),i=1,mesh(1)),j=1,mesh(1))

  stop
end program make_map

subroutine read_in(zr,deltaT,tau)

  use  sizes, only: mesh

  implicit none

!!$      include '../code/sizes.h'
!!$      include '../code/math.h'
!!$      include '../code/cgsatomic.h'
!!$      include '../code/astroconstants.h'
!!$      include '../code/cosmology.h'
!!$      include '../code/pmfast.h'
!!$      include '../code/abundances.h'

  integer :: i,j,dummy
  real*8 zr
  real*8 mean1,mean2,mean3,rms1,rms2,rms3
  character(len=6) :: zred_str
  character(len=180) :: file1, file2, file3
  real(kind=8) :: zred_now     
  real DeltaT(3,mesh(1),mesh(1))!assumes cubical mesh!!
  real tau(3,mesh(1),mesh(1))!assumes cubical mesh!!


  write(zred_str,'(f6.3)') zr
  file1='partial_integrals_w_random_vel/dTkSZ_x_'//trim(adjustl(zred_str))//'.bin'
  file2='partial_integrals_w_random_vel/dTkSZ_y_'//trim(adjustl(zred_str))//'.bin'
  file3='partial_integrals_w_random_vel/dTkSZ_z_'//trim(adjustl(zred_str))//'.bin'
  open(unit=54,file=file1,form='unformatted',status='unknown')
  open(unit=55,file=file2,form='unformatted',status='unknown')
  open(unit=56,file=file3,form='unformatted',status='unknown')
  !     x kSZ map
  read(54) dummy,dummy
  read(54)((deltaT(1,i,j),i=1,mesh(1)),j=1,mesh(1))
  !     y kSZ map
  read(55) dummy,dummy
  read(55)((deltaT(2,i,j),i=1,mesh(1)),j=1,mesh(1))
  !     y kSZ map
  read(56) dummy,dummy
  read(56)((deltaT(3,i,j),i=1,mesh(1)),j=1,mesh(1))
  close(54)
  close(55)
  close(56)

  file1='polarization/tau_x_'//trim(adjustl(zred_str))//'.bin'
  file2='polarization/tau_y_'//trim(adjustl(zred_str))//'.bin'
  file3='polarization/tau_z_'//trim(adjustl(zred_str))//'.bin'
  open(unit=54,file=file1,form='unformatted',status='unknown')
  open(unit=55,file=file2,form='unformatted',status='unknown')
  open(unit=56,file=file3,form='unformatted',status='unknown')
  !     x kSZ map
  read(54) dummy,dummy
  read(54)((tau(1,i,j),i=1,mesh(1)),j=1,mesh(1))
  !     y kSZ map
  read(55) dummy,dummy
  read(55)((tau(2,i,j),i=1,mesh(1)),j=1,mesh(1))
  !     y kSZ map
  read(56) dummy,dummy
  read(56)((tau(3,i,j),i=1,mesh(1)),j=1,mesh(1))
  close(54)
  close(55)
  close(56)

  mean1=sum(DeltaT(1,:,:))/mesh(1)**2
  mean2=sum(DeltaT(2,:,:))/mesh(1)**2
  mean3=sum(DeltaT(3,:,:))/mesh(1)**2
  !      print*,'check means',mean1,mean2,mean3

  rms1=sqrt(sum(DeltaT(1,:,:)**2)/mesh(1)**2-mean1**2)
  rms2=sqrt(sum(DeltaT(2,:,:)**2)/mesh(1)**2-mean2**2)
  rms3=sqrt(sum(DeltaT(3,:,:)**2)/mesh(1)**2-mean3**2)

  !      print*,'check rms',rms1,rms2,rms3

  return
end subroutine read_in

subroutine random_int (result, low, high)

  integer, intent (out) :: result
  integer, intent (in) :: low, high
  real uniform_random_value

  call random_number (uniform_random_value)
  result = int ((high - low + 1) * uniform_random_value + low)

end subroutine random_int


subroutine vrms_singlez(z,box,vrms,vrms_box)
  !     calculates the velocity power spectrum from the density one (CMBfast), 
  !     based on the linear theory. Note: returns vrms^2!!

  use cosmology_parameters
  use nbody

  implicit none

  integer :: k
  real*8 :: z, Ez, coeffz, growth, kr, v8, kmax, kmin,tophat 
  real*8 :: vrms, vrms_box,box

!  integer(4),   parameter :: nc     = n_box !       = 6144
!  integer(4),   parameter :: nc_rt  = mesh(1)!       = 256

  !! n_box is the number of cells per box length
  integer, parameter :: hc=n_box/2
  integer, parameter :: hc_rt=mesh(1)/2
  real, parameter    :: ncr=n_box
  real, parameter    :: hcr=hc

  !! Cosmo parameters
  !      real*8, parameter :: box= 100.0 !Mpc/h 
  real*8, parameter :: redshift= 300.0
  real*8, parameter :: scalefactor=1/(1+redshift)
!  real*8, parameter :: n_s=0.96 !WMAP3+
  !1.0 !0.95!!WMAP3 central values
  !      real*8, parameter :: s8=0.8 !0.9 !0.74
  real*8, parameter :: omegam=Omega0 !1.0-omegal !0.24
  real*8, parameter :: omegal=1-Omega0 !0.73 !0.76
  !      real*8, parameter :: h=0.7

  real*8 :: h00

!  real*8, parameter :: pi=3.14159
  real*8, parameter :: eps = 1.0d-03

  !     ! nk is the length of the initial power spectrum file
  integer, parameter      :: nk=437 !361 !408!204!407 !361 !406
  !      character(*), parameter :: fntf='cmbfast.z=0_WMAP3_fine'
  character(*), parameter :: fntf='CAMB_BAO_Tf_CAMB_BAO.dat'

  real*8, dimension(5,nk) :: tf
  !! tf(1,i) stores k
  !! tf(2,i) stores \Delta^2_m
  !! tf(3,i) stores \Delta^2_b
  !! tf(4,i) stores dk
  !! tf(5,i) will store \Delta^2_v (linear)  

  real*8 :: dummy

  !! Get transfer function from CAMB !CMBFAST
  write(*,*) 'Reading ',fntf
  open(11,file=fntf)
  do k=1,nk
     read(11,*) tf(1,k),tf(2,k),tf(3,k),tf(4,k),dummy,dummy,dummy
  end do
  !      read(11,*) tf
  close(11)

  !! Compute \Delta^2
  do k=1,nk
     kr     =tf(1,k)
     tf(2,k)=kr**(3+n_s)*tf(2,k)**2/(2*pi**2)
     tf(3,k)=kr**(3+n_s)*tf(3,k)**2/(2*pi**2)
  enddo

  !! Compute dk
  tf(4,1)=tf(1,2)/2
  do k=2,nk-1
     tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
  enddo
  tf(4,nk)=tf(1,nk)-tf(1,nk-1)

  !! Compute variance in 8 h^-1 Mpc spheres
  v8=0.
  kmax=2.*pi*sqrt(3.)*hc/box
  !      print*,'check',kmax,box,hc
  do k=1,nk
     if (tf(1,k) .gt. kmax) exit
     v8=v8+tf(2,k)*tophat(tf(1,k)*8.)**2*tf(4,k)/tf(1,k)
     !     print*,'v8',v8,tf(2,k),tf(1,k),tf(4,k),tophat(tf(1,k)*8.)
  enddo

  !! Normalize to \sigma_8
  !! Include growth factor
  tf(2:3,:)=tf(2:3,:)*(sigma8**2/v8)


  !      print*,'which redshift?'
  !      read(*,*) z

!  H0=100.*h                 !km/s/Mpc !no small h since k is in h/Mpc
  h00=100.*h
  Ez=sqrt(Omegam*(1.+z)**3+Omegal)
  coeffz=9.*h00**2*Omegam**2*(1.+z)**4/4./Ez**2*(1.-5./3./(1.+z)&
       /growth(omegam,omegal,z))**2

  do k=1,nk
     kr     =tf(1,k)/h
     !     tf(5,k)=tf(2,k)/kr**2*coeffz
     tf(2,k)=tf(2,k)/growth(omegam,omegal,z)**2
     tf(5,k)=tf(2,k)/kr**2*coeffz
     write(12,*) tf(1,k),sqrt(tf(5,k)),sqrt(tf(2,k))
     !     print*,'vel. ps', tf(1,k),sqrt(tf(5,k)),sqrt(tf(2,k))!Delta_v in km/s
  enddo

  !! Compute v_rms
  vrms=0.
  do k=1,nk
     vrms=vrms+tf(5,k)*tf(4,k)/tf(1,k)
  enddo

  !! Compute v_rms_box
  vrms_box=0.
  !     kmax=2.*pi*sqrt(3.)*hc/box
  kmax=2.*pi*hc_rt/box
  kmin=2.*pi/box
  !     print*,'check',kmax,kmin,box,tf(1,1),tf(1,nk)
  do k=nk,1,-1
     if (tf(1,k) .lt. kmax .and. tf(1,k) .gt. kmin)&
          vrms_box=vrms_box+tf(5,k)*tf(4,k)/tf(1,k)
  enddo

  print*,'v_rms=',vrms,vrms_box

  return
end subroutine vrms_singlez

function tophat(x)
  implicit none
  real*8 :: x,tophat

  if (x .ne. 0) then
     tophat=3.*(sin(x)-cos(x)*x)/x**3
  else
     tophat=1.
  endif

  return
end function tophat
