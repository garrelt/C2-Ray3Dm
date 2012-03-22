program make_maps_kSZ

  !     calculates the redshift-by-redshift integrated kSZ results
  !     from the simulation data

  use sizes, only: mesh
  use nbody  
  use grid, only: vol
  use material
  use cosmology_parameters
  use astroconstants
  use cgsconstants
  use velocity_mod
  use ion_frac_mod

  implicit none

  integer,parameter :: MaxNumZred=200
  integer ngrid1,ngrid2,ngrid3,i,j,k

  character(len=180) :: redshift_file ! name of file with list of redshifts
  character(len=180) :: results_file ! name of file with results
  real(kind=8) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter, convert
  integer :: nz, nz0 ! loop counter, starting slice

  real(kind=dp) :: zred_now !current redshift

  character(len=180) :: xfrac_file, dens_file, file1, file2, file3
  character(len=6) :: zred_str

  !operations done with doubles, otherwise sums of large arrays will be wrong
  real(kind=dp) :: ne(mesh(1),mesh(2),mesh(3))
  ! Data is double precision for newer runs!
  !  real(kind=si) :: xh1_real(mesh(1),mesh(2),mesh(3))

  !  real vel_real(3,mesh(1),mesh(2),mesh(3))
  !  real*8 vel(3,mesh(1),mesh(2),mesh(3))

  real DeltaT(3,mesh(1),mesh(1))!assumes cubical mesh!!

  real*8,parameter :: sigma_T=6.65e-25 !Thomson scattering cross-section
  real*8 :: dscale, x_cell,vnorm !, convert
!  real*8 :: lscale, tscale

  real*8 :: avg_mHII, avg_nHII, sum_xh1

  print*,'resolution', mesh

  call set_id

  allocate(ndens(mesh(1),mesh(2),mesh(3)))

  print*,id_str

  !just checks
  print*,c,H0, Omega0

!  lscale = boxsize/h/real(n_box)*Mpc ! size of a cell in cm, comoving
  dscale = Omega0*rho_crit_0 ! comoving gas density
!  tscale = 2./(3.*sqrt(omega0)*H0) ! time scale, when divided by (1+z)^2

  print*,'scales=',lscale,dscale,tscale,H0

  !box properties (n_box, boxsize) - now those are availale from cubep3m 
  ! module (but MAKE THEM PUBLIC!!!)
  rho_matter = dscale        ! mean matter density (g/cm^3)
 ! M_box      = rho_matter*(boxsize*Mpc/h)**3   ! mass in box (g, not M0) 
 ! M_particle = 8.0*M_box/(real(n_box)**3) ! mass per particle (g, not M0)
 ! M_grid = M_particle/8. 
  vol = (boxsize*Mpc/h)**3/(mesh(1)*mesh(2)*mesh(3))!volume per RT cell
  !     M_particle should be in g
  convert=M_grid*Omega_B/Omega0/(mu*m_p)/vol  ! density is CIC i.e. in M_grid
  x_cell=boxsize/h/mesh(1)  !rad. transfer cell size in comoving Mpc
  vnorm=lscale/tscale       !*(1.+z)FOR CONVERTING VELOCITIES FROM INTERNAL CODE UNITS TO CM/S

  print*,vnorm, x_cell, convert
  print*,n_box, M_box, M_grid

  !     Ask for redshift file
  write(*,'(A,$)') 'File with redshifts: '
  read(*,*) redshift_file
  write(*,'(A,$)') 'Initial snapshot number: '
  read(*,*) nz0
  write(*,'(A,$)') 'File with results: '
  read(*,*) results_file

  open(unit=2,file=results_file)

  !     Open and read redshift file
  open(unit=60,file=redshift_file,form='formatted',status='old')
  read(60,*) NumZred

  if (NumZred > MaxNumZred) write(*,*)&
       'Your program is about to die in a horrible way'

  do nz=nz0,NumZred !loop over redshifts
     read(60,*) zred_now

     print*,'doing z=',zred_now

     call read_density(zred_now)

     !     pause
     call read_velocity(zred_now)

     call read_xfrac(zred_now)

     ! Find the average ionized fraction: here just a check if 
     ! data is read correctly
     avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))

     avg_mHII=sum(ndens*xh1)/sum(ndens)

     print*,zred_now,avg_nHII,avg_mHII
     write(2,*) zred_now,avg_nHII,avg_mHII

     ne = xh1*ndens         !ionized gas number density (array)

     deltaT=0.0d0
     !NOTE: make sure indices correspond: vel. component 1 is along i, 2
     !     -along j, etc.

     do j=1,mesh(2)
        do k=1,mesh(3)
           do i=1,mesh(1)
              DeltaT(1,j,k)=DeltaT(1,j,k)+vel(1,i,j,k)*ne(i,j,k)
              !     &                 *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     do i=1,mesh(1)
        do k=1,mesh(3)
           do j=1,mesh(2)
              DeltaT(2,i,k)=DeltaT(2,i,k)+vel(2,i,j,k)*ne(i,j,k)
              !     &                *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     do i=1,mesh(1)
        do j=1,mesh(2)
           do k=1,mesh(3)
              DeltaT(3,i,j)=DeltaT(3,i,j)+vel(3,i,j,k)*ne(i,j,k)
              !     &                 *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     deltaT=deltaT*vnorm*x_cell*Mpc/c*sigma_T*convert*(1.+zred_now)**3 
     ! note that x_cell is comoving, so
     ! the 1+z normalization of the velocities drops out
     print*,'check TkSZ 1',minval(DeltaT),maxval(DeltaT),x_cell,&
          convert,maxval(ne*convert),minval(ne*convert)
     print*,'total norm.', vnorm*x_cell*Mpc/c*sigma_T*convert

     write(zred_str,'(f6.3)') zred_now
     file1='partial_integrals_w_random_vel/dTkSZ_x_'//trim(adjustl(zred_str))//'.bin'
     file2='partial_integrals_w_random_vel/dTkSZ_y_'//trim(adjustl(zred_str))//'.bin'
     file3='partial_integrals_w_random_vel/dTkSZ_z_'//trim(adjustl(zred_str))//'.bin'
     open(unit=54,file=file1,form='unformatted',status='unknown')
     open(unit=55,file=file2,form='unformatted',status='unknown')
     open(unit=56,file=file3,form='unformatted',status='unknown')
     !        x kSZ map
     write(54) mesh(1),mesh(1)
     write(54)((real(deltaT(1,i,j)),i=1,mesh(1)),j=1,mesh(1))
     !        y kSZ map
     write(55) mesh(1),mesh(1)
     write(55)((real(deltaT(2,i,j)),i=1,mesh(1)),j=1,mesh(1))
     !        z kSZ map
     write(56) mesh(1),mesh(1)
     write(56)((real(deltaT(3,i,j)),i=1,mesh(1)),j=1,mesh(1))
     close(54)
     close(55)
     close(56)

     print*,'done with ',zred_now

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     read(60,*) zred_now

     call read_xfrac(zred_now) 

     ! Find the average ionized fractions: here just a check if 
     ! data is read correctly
     avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))

     avg_mHII=sum(ndens*xh1)/sum(ndens)

     print*,zred_now,avg_nHII,avg_mHII
     write(2,*) zred_now,avg_nHII,avg_mHII

     ne = xh1*ndens         !ionized gas number density (array)

     deltaT=0.0d0
     !NOTE: make sure indices correspond: vel. component 1 is along i, 2
     !     -along j, etc.

     do j=1,mesh(2)
        do k=1,mesh(3)
           do i=1,mesh(1)
              DeltaT(1,j,k)=DeltaT(1,j,k)+vel(1,i,j,k)*ne(i,j,k)
              !     &                 *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     do i=1,mesh(1)
        do k=1,mesh(3)
           do j=1,mesh(2)
              DeltaT(2,i,k)=DeltaT(2,i,k)+vel(2,i,j,k)*ne(i,j,k)
              !     &                *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     do i=1,mesh(1)
        do j=1,mesh(2)
           do k=1,mesh(3)
              DeltaT(3,i,j)=DeltaT(3,i,j)+vel(3,i,j,k)*ne(i,j,k)
              !     &                 *exp(-sigma_T*x_cell*ne(i,j,k))*x_cell
           end do
        end do
     end do

     deltaT=deltaT*vnorm*x_cell*Mpc/c*sigma_T*convert*(1.+zred_now)**3 
     ! note that x_cell is comoving, so
     ! the 1+z normalization of the velocities drops out
     print*,'check TkSZ 1',minval(DeltaT),maxval(DeltaT),x_cell,&
          convert,maxval(ne*convert),minval(ne*convert)
     print*,'total norm.', vnorm*x_cell*Mpc/c*sigma_T*convert

     write(zred_str,'(f6.3)') zred_now
     file1='partial_integrals_w_random_vel/dTkSZ_x_'//trim(adjustl(zred_str))//'.bin'
     file2='partial_integrals_w_random_vel/dTkSZ_y_'//trim(adjustl(zred_str))//'.bin'
     file3='partial_integrals_w_random_vel/dTkSZ_z_'//trim(adjustl(zred_str))//'.bin'
     open(unit=54,file=file1,form='unformatted',status='unknown')
     open(unit=55,file=file2,form='unformatted',status='unknown')
     open(unit=56,file=file3,form='unformatted',status='unknown')
     !        x kSZ map
     write(54) mesh(1),mesh(1)
     write(54)((real(deltaT(1,i,j)),i=1,mesh(1)),j=1,mesh(1))
     !        y kSZ map
     write(55) mesh(1),mesh(1)
     write(55)((real(deltaT(2,i,j)),i=1,mesh(1)),j=1,mesh(1))
     !        z kSZ map
     write(56) mesh(1),mesh(1)
     write(56)((real(deltaT(3,i,j)),i=1,mesh(1)),j=1,mesh(1))
     close(54)
     close(55)
     close(56)

     print*,'done with ',zred_now
!     pause

  end do
  close(60)

  stop
end program make_maps_kSZ
