      program make_PDF_Gamma

!     calculates the PDF of the photoionization rates
!     from the simulation results. Also calculates and 
!     makes a histogram of the corresponding equilibrium 
!     ion. rates

!     Author: Ilian Iliev
!     Date: 20-Feb-2010

      use sizes, only: mesh
!      use mathconstants
!      use astroconstants
!      use cgsconstants
!      use cosmology
      use nbody
      use material
      use ion_frac_mod
      use c2ray_parameters, only: type_of_clumping,clumping_factor
      use radiative_cooling
      use grid, only: vol, grid_ini
      use cgsconstants, only: albpow,bh00,colh0,temph0

      implicit none
      
      integer,parameter :: MaxNumZred=200
      integer ngrid1,ngrid2,ngrid3,i,j,k,ii,jj,l, ll
      
      ! name of file with list of redshifts
      character(len=180) :: redshift_file 
      character(len=180) :: single_z_file, out_file
      character(len=180) :: results_file ! name of file with results
      real(kind=dp) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter
      character(len=6) :: z_str
      integer :: nz ! loop counter
      
!      real(kind=dp) ::  d(mesh(1)*mesh(2)*mesh(3)),x(mesh(1)*mesh(2)*mesh(3)),xd(mesh(1)*mesh(2)*mesh(3))
      real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)) :: d,x,gamma,lgamma
      real*8,dimension(mesh(1),mesh(2),mesh(3)) :: gamma_eq,lgamma_eq,gamma_eq_c1

      integer ans, size, sample,size_Mpc, arr_size
      character(len=2) :: size_str
      character(len=1) :: nfile_str
      
      ! PDF array 
      
      integer,parameter :: n_x=500 !how many bins to use
      real(kind=dp),dimension(0:n_x) :: distriblir, distriblir_eq
      real(kind=dp) :: maxlir, minlir, mean_lir, sigma_lir, dlir
      real(kind=dp) :: maxlir_eq, minlir_eq, mean_lir_eq, sigma_lir_eq, dlir_eq
      
      real*8 :: temp0, brech0, acolh0, sqrtt0, convert
!      real*8 :: colh0, bh00, albpow, temph0 !, clumping
      real*8 :: maxlir0, minlir0
      
      real(kind=dp) :: zred_now, z, z_next
      integer nz0

      character(len=180) :: xfrac_file, dens_file
      character(len=6) :: zred_str 
      
      real(kind=dp) :: avg_mHII, avg_nHII, sum_xh1, avg_gamma, avg_mgamma

      common /g/ gamma
      
      call grid_ini()

      allocate(ndens(mesh(1),mesh(2),mesh(3)))
      
      print*,'resolution:', mesh     
      call set_id

      out_file="PDF_means.dat"

      !     Open output file
      open(unit=23,file=out_file,status='unknown')


!     Ask for redshift file
      write(*,'(A,$)') 'File with redshifts: '
      read(*,*) redshift_file
      write(*,'(A,$)') 'Initial snapshot number: '
      read(*,*) nz0
!      write(*,'(A,$)') 'File with results: '
!      read(*,*) results_file

      write(*,'(A,$)') 'Size of regions (in cells) and'
      write(*,'(A,$)') ' sampling (every ... in x,y,x):'
      read(*,*) size,sample
      
      temp0=1.0d4
      minlir0=-21.5
      maxlir0=-8.
     
!     Open and read redshift file
      open(unit=60,file=redshift_file,form='formatted',status='old')
      read(60,*) NumZred

!      allocate(snap(NumZred))
!      allocate(zred_array(NumZred))

      if (NumZred > MaxNumZred) write(*,*)&
          'Your program is about to die in a horrible way'

!!$      do nz=1,nz0-1      
!!$         read(60,*) zred_array(nz)
!!$      end do
!!$
      do nz=nz0,NumZred

!!$         write(id0,'(i3)') nz 
!!$         if(nz<10) then	    
!!$            id="00"//trim(adjustl(id0))
!!$         elseif(10<=nz.and.nz<100) then
!!$            id="0"//trim(adjustl(id0))
!!$         else
!!$            id=trim(adjustl(id0))
!!$         end if

!         print*,'check ',id,id0,id_str,mesh

         read(60,*) zred_now

         print*,'doing z=',zred_now
         !
         call read_density(zred_now)
         
         call read_xfrac(zred_now)

         call read_ion_rate(zred_now)

!     Find the average ionized fraction
         avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         
         avg_mHII=sum(ndens*xh1)/sum(ndens)
         
         avg_gamma=sum(gamma)/real(mesh(1)*mesh(2)*mesh(3))
         avg_mgamma=sum(ndens*gamma)/sum(ndens)

         print*,zred_now,avg_nHII,avg_mHII,avg_gamma,avg_mgamma
!         write(2,*) zred_now,avg_nHII,avg_mHII
         
         !report averages for ion. rates
         write(*,*) 'minimum: ',minval(gamma)
         write(*,*) 'maximum: ',maxval(gamma)
         write(*,*) 'average ion. rate',avg_gamma

         write(102,*) zred_now, avg_gamma,avg_mgamma

         call set_clumping(zred_now) ! this comes from up-to-date versions 
                                     ! of mat_ini_cubep3m.F90; sets clumping
                                     ! based on type_of_clumping, set in 
                                     ! module c2ray_parameters in c2ray_parameters.f90

         print*,'clumping set to ',clumping

         print*,'check parameters',Omega0,mu,m_p,vol,temp0,  bh00, albpow, temph0, colh0

         !FOR VARIABLE CLUMPING WE NEED TO MAKE ALL THESE ARRAYS!

         brech0=clumping*bh00*(temp0/1e4)**albpow!recomb. rate
         sqrtt0=sqrt(temp0)
         acolh0=colh0*sqrtt0*exp(-temph0/temp0)!collisional ion. rate
         convert=M_grid*Omega_B/Omega0/(mu*m_p)/vol 
         ndens=ndens*convert*(1.0d0+zred_now)**3
         
         print*,'check equilib. ion. rates',sum(ndens)/(real(mesh(1))*real(mesh(2))*real(mesh(3))),brech0,acolh0

         do i=1,mesh(1)
            do j=1,mesh(2)
               do k=1,mesh(3)
                  gamma_eq(i,j,k)=xh1(i,j,k)*ndens(i,j,k)*(xh1(i,j,k)/(1.d0-xh1(i,j,k))*brech0-acolh0)
                  if(gamma_eq(i,j,k)<0.0d0)gamma_eq(i,j,k)=0.0d0
                  gamma_eq_c1(i,j,k)=xh1(i,j,k)*ndens(i,j,k)*(xh1(i,j,k)/(1.d0-xh1(i,j,k))*brech0/clumping-acolh0)
                  if(gamma_eq_c1(i,j,k)<0.0d0)gamma_eq_c1(i,j,k)=0.0d0
               end do
            end do
         end do

         write(*,*) 'minimum eq.: ',minval(gamma_eq)
         write(*,*) 'maximum eq.: ',maxval(gamma_eq)
         write(*,*) 'average eq. ion. rate',sum(gamma_eq)/real(mesh(1)*mesh(2)*mesh(3))
         
         print*,'check',clumping,sum(gamma_eq)/real(mesh(1)*mesh(2)*mesh(3)),sum(gamma_eq_c1)/real(mesh(1)*mesh(2)*mesh(3)),clumping*sum(gamma_eq_c1)/real(mesh(1)*mesh(2)*mesh(3))
         
         !start making the PDF
         
         !make histograms 
         
         lgamma=log10(gamma+1.d-50)
         lgamma_eq=log10(gamma_eq+1.d-50)
         
         maxlir_eq=maxval(lgamma_eq)
         minlir_eq=minval(lgamma_eq)
         mean_lir_eq=sum(lgamma_eq)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         sigma_lir_eq=sqrt(sum((lgamma_eq-mean_lir_eq)**2)/(real(mesh(1))*real(mesh(2))*real(mesh(3))-1))
         
         maxlir=maxval(lgamma)
         minlir=minval(lgamma)
         mean_lir=sum(lgamma)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         sigma_lir=sqrt(sum((lgamma-mean_lir)**2)/(real(mesh(1))*real(mesh(2))*real(mesh(3))-1))
         
!         minlir=minlir0
!         minlir_eq=minlir0
         !     bin sizes
         dlir_eq  = (maxlir_eq-minlir_eq)/real(n_x)
         dlir  = (maxlir-minlir)/real(n_x)
         print*,'check',maxlir,minlir,dlir
         print*,'mean, sigma',mean_lir, sigma_lir
         
         distriblir = 0.0d0
         distriblir_eq = 0.0d0

         !     print*,'check arrays0',maxlir,maxlir_eq,minlir,minlir_eq,dlir,dlir_eq
         
         do i=1,mesh(1)
            do j=1,mesh(2)
               do k=1,mesh(3)
                  !  print*,'check arrays',floor((lgamma(i,j,k)-minlir)/dlir),floor((lgamma_eq(i,j,k)-minlir_eq)/dlir_eq)
                  !  total number of cells in each density bin
                  
                  distriblir(max(floor((lgamma(i,j,k)-minlir)/dlir),0))=&
                       distriblir(max(floor((lgamma(i,j,k)-minlir)/dlir),0))+1. 
                  distriblir_eq(max(floor((lgamma_eq(i,j,k)-minlir_eq)/dlir_eq),0))=&
                       distriblir_eq(max(floor((lgamma_eq(i,j,k)-minlir_eq)/dlir_eq),0))+1. 
               end do
            end do
         end do
         !normalize
         distriblir=distriblir/sum(distriblir)
         distriblir_eq=distriblir_eq/sum(distriblir_eq)

         call openfile(zred_now)
         
         write(23,*) zred_now
         write(23,*) minval(gamma),maxval(gamma),sum(gamma)/real(mesh(1)*mesh(2)*mesh(3))
         write(23,*) minval(gamma_eq),maxval(gamma_eq),sum(gamma_eq)/real(mesh(1)*mesh(2)*mesh(3))
         
         write(20,*) mean_lir, sigma_lir,dlir
         do ll=0,n_x
            write(20,*)minlir+real(ll)*dlir,&
                 real(distriblir(ll))/dlir
         end do
         
         write(22,*) mean_lir_eq, sigma_lir_eq, dlir_eq
         do ll=0,n_x
            write(22,*)minlir_eq+real(ll)*dlir_eq,&
                 real(distriblir_eq(ll))/dlir_eq
         end do
         
         print*,'check normalization',sum(distriblir),sum(distriblir_eq)     
         
         call closefile
         
         ndens=ndens/(sum(ndens)/real(mesh(1)*mesh(2)*mesh(3)))-1.0d0 !overdensity

     call output_correl(zred_now)
     
  end do
  
  close(23)
  
contains
  
  
  subroutine openfile(z)

    real*8 :: z

    write(z_str,'(f6.3)') z    
    out_file="PDF_114Mpc_f10_150S_"//trim(adjustl(z_str))//".dat"
    
    !     Open output file
    open(unit=20,file=out_file,status='unknown')

    out_file="PDF_114Mpc_f10_150S_eq_"//trim(adjustl(z_str))//".dat"
    
    !     Open output file
    open(unit=22,file=out_file,status='unknown')
    
  end subroutine openfile

  subroutine closefile
    close(20)
    close(22)
  end subroutine closefile

  subroutine output_correl(z)
    
    real*8 :: z
    
    write(z_str,'(f6.3)') z
    
    out_file="ionrate_dens_correl_114Mpc_f10_150S_"//trim(adjustl(z_str))//".dat"
    
    !     Open output file
    open(unit=20,file=out_file,status='unknown')

     do i=1,mesh(1),2
        do j=1,mesh(2),2
           do k=1,mesh(3),2    
              write(20,*) ndens(i,j,k), lgamma(i,j,k)
           end do
        end do
     end do

     close(20)

   end subroutine output_correl

 end program make_PDF_Gamma
