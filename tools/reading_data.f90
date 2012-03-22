module directories

 ! Where the data is and where the results will be written
  character(len=180),parameter :: &
       data_dir="../results/" 
  character(len=180),parameter :: &
       results_dir="./" 
  character(len=180),parameter :: &
       dir_dens_read="../" 

  character(len=3) :: id, id0 !output number string

end module directories
module velocity_mod

  use sizes, only: mesh
  use precision, only: dp,si

  implicit none

  real(kind=dp) :: vel(3,mesh(1),mesh(2),mesh(3))

end module velocity_mod

module ion_frac_mod

  use sizes, only: mesh
  use precision, only: dp,si

  implicit none

  real(kind=dp) :: xh1(mesh(1),mesh(2),mesh(3))

end module ion_frac_mod

module source_statistics

  use precision, only: dp,si

  integer :: NumSrc00,NumSrc01
  
  real(kind=dp) :: TotSrcFlux,TotSrcFlux0, LargeSrcFlux0,SmallSrcFlux0

end module source_statistics

subroutine set_id
  use sizes, only: mesh
  use nbody, only : id_str
  use material
 
  implicit none

  ! Set identifying string (resolution-dependent)
  ! Construct the file name
  if (mesh(1) == 784) id_str="coarse"
  if (mesh(1) == 504) id_str="coarser"
  if (mesh(1) == 252) id_str="coarsest"
  
  print*,'resolution:', id_str,mesh
  
  return
end subroutine set_id

subroutine read_density(z)
  use material
  use nbody, only : id_str
  use directories

  implicit none

  integer m1,m2,m3
  integer i,j,k
  character(len=6) :: zred_str
  real(kind=dp) :: z
  character(len=180) :: dens_file
!  real(kind=dp) :: avg_dens
  ! density in file is in 4B reals, read in via this array
  real(kind=si),dimension(:,:,:),allocatable :: ndens_real
  
  print*,'id, dens:', id_str, dir_dens_read

  write(zred_str,'(f6.3)') z
  
!  if (id_str /= "coarse") then
     write(zred_str,'(f6.3)') z
!     if (id_str /= "coarse") then
        ! Allocate array needed to read in data
        allocate(ndens_real(mesh(1),mesh(2),mesh(3)))

        write(30,*) 'Reading ',id_str,' input'
        dens_file=trim(adjustl(dir_dens_read))//"coarser_densities/"// &
        !dens_file="../../c2ray_114Mpc_f10_150S_256_Cz/coarser_densities/"// &
             trim(adjustl(zred_str))//"n_all.dat"

        print*,dens_file

        ! Open density file: note that it is in `unformatted' form
        open(unit=20,file=dens_file,form='binary',status='old')
        !          open(unit=20,file=dens_file,form='unformatted',status='old')
        
        read(20) m1, m2, m3 
        print*,m1,m2,m3
        if (m1.ne.mesh(1).or.m2.ne.mesh(2).or.m3.ne.mesh(3))then
           write(*,*) 'Warning: file with density has problems'
           print*,m1,m2,m3
        end if
        ! Read in data and store it in ndens
        read(20) ndens_real
        ndens(:,:,:)=ndens_real(:,:,:)
        
        ! close file
        close(20)

        ! Deallocate array needed for reading in the data.
        deallocate(ndens_real)

!!$     else
!!$        allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
!!$        dens_file="../"//&
!!$             trim(adjustl(zred_str))//"rho_c.dat"
!!$        write(30,*) dens_file
!!$        open(unit=20,file=dens_file,form='binary',status='old')
!!$        read(20) ndens_real
!!$        ndens(:,:,:)=ndens_real(:,:,:)
!!$        close(20)
!!$
!!$        deallocate(ndens_real)
!!$     endif
!!$  endif
  
  !     Find summed_and average density
  avg_dens=sum(1.0d0*ndens)/(real(mesh(1),8)*real(mesh(2),8)*real(mesh(3),8))
  
  !     Report on data: min, max, total
  write(*,*) 'minimum: ',minval(ndens)/8.
  write(*,*) 'maximum: ',maxval(ndens)/8.
  write(*,*) 'summed density: ',sum(1.0d0*ndens)/8.
  write(*,*) 'average dens.',avg_dens/8.

  return
end subroutine read_density

subroutine read_xfrac(z)
  use sizes, only: mesh  
  use material
  use nbody
  use ion_frac_mod

  implicit none

  integer i,j,k
  integer m1,m2,m3
  character(len=6) :: zred_str
  real(kind=dp) :: z
  character(len=180) :: xfrac_file
!  real(kind=dp) :: xh1(mesh(1),mesh(2),mesh(3))
  real(kind=dp),dimension(:,:,:),allocatable :: xh1_real
!  real(kind=si) :: xh1_r(mesh(1),mesh(2),mesh(3))

!  common /x/ xh1

  allocate(xh1_real(mesh(1),mesh(2),mesh(3))) 

  write(zred_str,'(f6.3)') z
!  xfrac_file= "../results/Ifront3_"//trim(adjustl(zred_str))//".bin"
  xfrac_file= "../results/xfrac3d_"//trim(adjustl(zred_str))//".bin"
  
  !     Open ionized fractions file
  open(unit=20,file=xfrac_file,form='unformatted',status='old')
  
  !     Read in data
  read(20) m1,m2,m3
  print*,m1,m2,m3
  if (m1.ne.mesh(1).or.m2.ne.mesh(2).or.m3.ne.mesh(3)) then
     write(*,*) 'Warning: file with ionization fractions unusable'
  else
     read(20) xh1_real
     xh1(:,:,:)=real(xh1_real(:,:,:),8)
     !            xh1=1.0d0*xh1_real
  endif
  !     close file
  close(20)

  deallocate(xh1_real)

  return
end subroutine read_xfrac

subroutine read_sources(z,supp)
  use sizes, only: mesh, Ndim  
  use material
  use nbody
  use source_statistics
  use c2ray_parameters
  use cosmology_parameters

  implicit none

  integer i,j,k, supp !supp=0 - sources without, supp=1 - with suppression
  integer :: ns0
  integer m1,m2,m3
  character(len=6) :: zred_str
  real(kind=dp) :: z
  character*180 sourcelistfile,sourcelistfilesuppress  
!  integer, parameter :: Nsrc=200000 !note that this is done for minimizing the storage
  real(kind=dp) :: rsrcpos,NormFlux,SrcMass !(Nsrc)
!  integer :: NumSrc0, NumSrc

  integer srcpos(Ndim)!note that this is done for minimizing the storage, change this if 
!  integer :: srcpos(Ndim,Nsrc) !all source positions are needed.

  real(kind=dp) :: SrcMass00, SrcMass01!massive (0) and less massive (01) sources
  
  real(kind=dp) :: dscale, rho_matter
  
  print*,'doing z=',z
  !     Sources are read from file
  !        
  !     Construct the file name

  dscale = Omega0*rho_crit_0 ! comoving gas density
  rho_matter = dscale        ! mean matter density (g/cm^3)
!  M_box      = rho_matter*(boxsize*Mpc/h)**3   ! mass in box (g, not M0) 
!  M_particle = 8.0*M_box/(real(n_box)**3) ! mass per particle (g, not M0)
!  M_grid = M_particle/8. 

  if(supp==0)then
     if (z.ge.10.0) then
        write(sourcelistfile,'(f6.3,3A)') z,&
             "-",trim(adjustl(id_str)),"_sources.dat"
     else
        write(sourcelistfile,'(f5.3,3A)') z,&
             "-",trim(adjustl(id_str)),"_sources.dat"
     endif
     open(unit=50,file='../sources/'//sourcelistfile,status='old')
     !        Number of sources without suppression
     read(50,*) NumSrc00
     !        Report
     write(*,*) 'Number of sources, no suppression: ',NumSrc00

     TotSrcFlux0=0.0
     LargeSrcFlux0=0.0
     SmallSrcFlux0=0.0

     do ns0=1,NumSrc00     
        
        read(50,*) srcpos(1),srcpos(2),srcpos(3),&
             SrcMass00,SrcMass01

        TotSrcFlux0=TotSrcFlux0+(SrcMass00*phot_per_atom(1)  & !massive sources
             +SrcMass01*phot_per_atom(2))* & !small sources   
             M_grid*Omega_B/(Omega0*m_p)   !normalization, note that the source lifetime 
                                           !is not included
        LargeSrcFlux0=LargeSrcFlux0+SrcMass00*phot_per_atom(1)*M_grid*Omega_B/(Omega0*m_p)
        SmallSrcFlux0=SmallSrcFlux0+SrcMass01*phot_per_atom(2)*M_grid*Omega_B/(Omega0*m_p)

     end do
  else
     if (z.ge.10.0) then
        write(sourcelistfilesuppress,'(f6.3,3A)') z,&
             "-",trim(adjustl(id_str)),"_sources_used_wfgamma.dat"
     else
        write(sourcelistfile,'(f5.3,3A)') z,&
             "-",trim(adjustl(id_str)),"_sources.dat"
        write(sourcelistfilesuppress,'(f5.3,3A)') z,&
             "-",trim(adjustl(id_str)),"_sources_used_wfgamma.dat"
     endif
     open(unit=49,file='../sources/'//sourcelistfilesuppress,status='old')
     !        Number of sources
     read(49,*) NumSrc01       
     !        Report
     write(*,*) 'Number of sources: ',NumSrc01 

     TotSrcFlux=0.

     do ns0=1,NumSrc01 
        read(49,*) srcpos(1),srcpos(2),srcpos(3),SrcMass
        !srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0),SrcMass(ns0)
        TotSrcFlux=TotSrcFlux+SrcMass* & !small sources   
             M_grid*Omega_B/(Omega0*m_p) !normalization, note that the source lifetime 
                                         !is not included, while f_gamma is in SrcMass here
     end do
  end if

  return
end subroutine read_sources

subroutine read_ion_rate(z)
  use sizes, only: mesh  
  use material
  use nbody

  implicit none

  integer i,j,k
  integer m1,m2,m3
  character(len=6) :: zred_str
  real(kind=dp) :: z
  character(len=180) :: ion_rate_file
  real(kind=dp) :: gamma(mesh(1),mesh(2),mesh(3))
  real(kind=si),dimension(:,:,:),allocatable :: gamma_real

  common /g/ gamma

  allocate(gamma_real(mesh(1),mesh(2),mesh(3))) 

  write(zred_str,'(f6.3)') z
  ion_rate_file= "../results/IonRates3_"//trim(adjustl(zred_str))//".bin"
  
  !     Open ionized fractions file
  open(unit=20,file=ion_rate_file,form='unformatted',status='old')
  
  !     Read in data
  read(20) m1,m2,m3
  print*,m1,m2,m3
  if (m1.ne.mesh(1).or.m2.ne.mesh(2).or.m3.ne.mesh(3)) then
     write(*,*) 'Warning: file with ionization rates unusable'
  else
     read(20) gamma_real
     gamma(:,:,:)=real(gamma_real(:,:,:),8)
  endif
  !     close file
  close(20)

  deallocate(gamma_real)

  return
end subroutine read_ion_rate

    
subroutine read_velocity(z)
  use material
  use nbody, only : id_str
  use velocity_mod

  implicit none

  integer m1,m2,m3
  integer i,j,k
  character(len=6) :: zred_str
  real(kind=dp),intent(in) :: z
  real(kind=dp) :: avg_vel(3) 
  character(len=180) :: vel_file

  ! velocity in file is in 4B reals, read in via this array
!  real(kind=si),dimension(:,:,:,:),allocatable :: vel_real
  real(kind=si),dimension(:,:,:,:), allocatable :: vel_real
!  real(kind=dp) :: vel(3,mesh(1),mesh(2),mesh(3))


  print*,'id, dens:', id_str

  write(zred_str,'(f6.3)') z
  
  
        ! Allocate array needed to read in data
  allocate(vel_real(3,mesh(1),mesh(2),mesh(3)))
  write(30,*) 'Reading ',id_str,' input'
  !          dens_file=trim(adjustl(dir_dens))//"coarser_densities/"// &
  !             "rho_"//trim(adjustl(id_str))//".dat"
  vel_file="../../../coarser_data_SPH_3072_114Mpc/halos_removed/"// &
       trim(adjustl(zred_str))//"v_all.dat"
  
  print*,vel_file
  
  ! Open density file: note that it is in `unformatted' form
  open(unit=20,file=vel_file,form='binary',status='old')
  !          open(unit=20,file=dens_file,form='unformatted',status='old')
  
  read(20) m1, m2, m3 
  print*,m1,m2,m3
  if (m1.ne.mesh(1).or.m2.ne.mesh(2).or.m3.ne.mesh(3))then
     write(*,*) 'Warning: file with density has problems'
     print*,m1,m2,m3
  end if
  ! Read in data and store it in vel
  read(20) vel_real
!!$  do k=1,m3
!!$     do j=1,m2
!!$        do i=1,m1
!!$           !read(20) vel_real(:,i,j,k)
!!$        enddo
!!$     enddo
!!$  enddo

  !this is because the velocity in SPH-smoothed form is not normalized per particle
  do i=1,3
     vel(i,:,:,:)=vel_real(i,:,:,:)*8/ndens(:,:,:)
  end do
!!$        
  ! close file
  close(20)
  
  ! Deallocate array needed for reading in the data.
  deallocate(vel_real)
  
!!$  endif
!!$  
!!$!  ndens=ndens/8. !to fix normalization (in particle masses, not grid) 
!!$
  !     Find summed_and average velocity (as a check)
  do i=1,3
     avg_vel(i)=sum(vel(i,:,:,:))/(real(mesh(1),8)*real(mesh(2),8)*real(mesh(3),8))
  end do

  !     Report on data: min, max, total
  write(*,*) 'minimum: ',minval(vel)
  write(*,*) 'maximum: ',maxval(vel)
  write(*,*) 'summed velocity: ',sum(vel)
  write(*,*) 'average vel.',avg_vel
  
  return
end subroutine read_velocity


