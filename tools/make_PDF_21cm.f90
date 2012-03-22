      program make_PDF_21cm

!     calculates the volume- and mass-neutral fractions at a given z
!     from the simulation results 

      use sizes, only: mesh
!      use mathconstants
!      use astroconstants
!      use cgsconstants
!      use cosmology
      use nbody
      use material
      use ion_frac_mod

      implicit none

      integer,parameter :: MaxNumZred=200
      integer ngrid1,ngrid2,ngrid3,i,j,k,ii,jj,l

      character(len=180) :: redshift_file ! name of file with list of redshifts
      character(len=180) :: single_z_file
      character(len=180) :: results_file ! name of file with results
      real(kind=dp) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter
      integer :: nz ! loop counter
!      integer ::  NumZred

      ! density in file is in real(kind=4), read in via this array
!      real ndens_real(mesh(1),mesh(2),mesh(3))
!     real(kind=dp) ndens(mesh(1),mesh(2),mesh(3))

      real(kind=dp) ::  d(mesh(1)*mesh(2)*mesh(3)),x(mesh(1)*mesh(2)*mesh(3)),xd(mesh(1)*mesh(2)*mesh(3))

      real(kind=dp) :: mean_ion_rho_xm, mean_ion_m, mean_dens
      real(kind=dp) :: lmean_ion_rho_xm, lmean_ion_m, lmean_dens
      real(kind=dp) :: maxd,mind,maxx,minx,maxxd,minxd,dd,dx,dxd
      real(kind=dp) :: mean_d,mean_x,mean_xd,sigma_d,sigma_x,sigma_xd
      integer ans, size, sample,size_Mpc, arr_size
      character(len=2) :: size_str
      character(len=1) :: nfile_str
      integer,parameter :: n_x=100 !how many bins to use
      real(kind=dp),dimension(0:n_x) :: distribxm, distribrho, distribrhoxm

      real(kind=dp) :: zred_now
      integer nz0

      character(len=180) :: xfrac_file, dens_file
      character(len=6) :: zred_str


      !operations done with doubles, otherwise sums of large arrays will be wrong
!      real(kind=dp) :: xh1(mesh(1),mesh(2),mesh(3))
      ! Array needed to read in 4B reals
!      real(kind=si),dimension(:,:,:),allocatable :: xh1_real
!      real(kind=si) :: xh1_real(mesh(1),mesh(2),mesh(3))

      real(kind=dp) :: avg_mHII, avg_nHII, sum_xh1
      
!      common /x/ xh1

!      allocate(xh1_real(mesh(1),mesh(2),mesh(3))) 

      allocate(ndens(mesh(1),mesh(2),mesh(3)))

      print*,'resolution:', mesh     
      call set_id

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

!      open(unit=2,file=results_file)
     
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

!     Find the average ionized fraction
         avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         
         avg_mHII=sum(ndens*xh1)/sum(ndens)

         print*,zred_now,avg_nHII,avg_mHII
!         write(2,*) zred_now,avg_nHII,avg_mHII

!start making the PDF

         distribrho=0.0d0
         distribxm=0.0d0
         distribrhoxm=0.0d0
         x=0.0d0
         d=0.0d0
         xd=0.0d0

!     output files
         write(zred_str,'(f6.3)') zred_now
         size_Mpc=nint(size*boxsize/mesh(1))
         print*,'Doing regions of size (in Mpc)',size_Mpc
         if(size_Mpc .lt. 10)then
            write(size_str,'(i1)') size_Mpc
         elseif(size_Mpc .lt. 100)then
            write(size_str,'(i2)') size_Mpc
         else
            write(size_str,'(i3)') size_Mpc
         end if

         single_z_file="PDF_d_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=61,file='../subvolume_data/'//single_z_file)
         single_z_file="PDF_xm1_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=62,file='../subvolume_data/'//single_z_file)
         single_z_file="PDF_xm1d_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=63,file='../subvolume_data/'//single_z_file)

         ii=size
         jj=0
!     loop over the cubical regions of size ii along each direction
         do i=1,mesh(1)-ii+1,sample      
            do j=1,mesh(1)-ii+1,sample 
               do k=1,mesh(1)-ii+1,sample
                  jj=jj+1
                  d(jj)=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1))/&
                       real(ii)**3/avg_dens !normalized mean density per region  
                  x(jj)=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1)& !normalized mean mass ...
                       *xh1(i:i+ii-1,j:j+ii-1,k:k+ii-1))& !neutral fraction per region
                       /sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1))
                  x(jj)=(1.d0-x(jj))/(1.0d0-avg_mHII) !for doing the neutral mass  
               end do
            end do
         end do
         xd=x*d !normalized mean ionized/neutral density per region
               
         arr_size=jj
!     do not count 0's! Make arrays allocatable and allocate only part that is needed!
         maxd=maxval(d(1:arr_size))
         mind=minval(d(1:arr_size))
         mean_d=sum(d(1:arr_size))/real(arr_size)
         sigma_d=sqrt(sum((d(1:arr_size)-mean_d)**2)/real(arr_size-1))

         maxx=maxval(x(1:arr_size))
         minx=minval(x(1:arr_size))
         mean_x=sum(x(1:arr_size))/real(arr_size)
         sigma_x=sqrt(sum((x(1:arr_size)-mean_x)**2)/real(arr_size-1))

         maxxd=maxval(xd(1:arr_size))
         minxd=minval(xd(1:arr_size))
         mean_xd=sum(xd(1:arr_size))/real(arr_size)
         sigma_xd=sqrt(sum((xd(1:arr_size)-mean_xd)**2)/real(arr_size-1))

         
!     bin sizes
         dd  = (maxd-mind)/real(n_x)
         dx  = (maxx-minx)/real(n_x)
         dxd = (maxxd-minxd)/real(n_x)
         print*,'check',maxd,mind,maxxd,minxd,dd,dx,dxd
         print*,'means, sigmas',mean_d, sigma_d,mean_x, sigma_x,&
              mean_xd, sigma_xd
         
         do jj=1,arr_size
!     total number of cells in each density bin
            distribrho(floor((d(jj)-mind)/dd))=distribrho(floor((d(jj)-mind)/dd))+1 
            distribxm(floor((x(jj)-minx)/dx))=& !total number of cells in each ioniz. bin
                 distribxm(floor((x(jj)-minx)/dx))+1. !in log bins 
            distribrhoxm(floor((xd(jj)-minxd)/dxd))=& !total number of cells in each 
                 distribrhoxm(floor((xd(jj)-minxd)/dxd))+1. !ioniz. bin in log bins 
         end do
         
         distribrho=distribrho/sum(distribrho)
         distribxm=distribxm/sum(distribxm)
         distribrhoxm=distribrhoxm/sum(distribrhoxm)

         write(61,*) mean_d, sigma_d,dd
         write(62,*) mean_x, sigma_x,dx         
         write(63,*) mean_xd, sigma_xd,dxd

         do l=0,n_x
            write(61,*)mind+real(l)*dd,real(distribrho(l))/dd
            write(62,*)minx+real(l)*dx,real(distribxm(l))/dx
            write(63,*)minxd+real(l)*dxd,real(distribrhoxm(l))/dxd
         end do
         
         print*,'check normalization',sum(distribrho),sum(distribxm),sum(distribrhoxm)
         close(61)
         close(62)
         close(63)

         read(60,*) zred_now
         
         call read_xfrac(zred_now)

!     Find the average ionized fraction
         avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         
         avg_mHII=sum(ndens*xh1)/sum(ndens)

         print*,zred_now,avg_nHII,avg_mHII
!         write(2,*) zred_now,avg_nHII,avg_mHII

!start making the PDF

         distribrho=0.0d0
         distribxm=0.0d0
         distribrhoxm=0.0d0
         x=0.0d0
         d=0.0d0
         xd=0.0d0

!     output files
         write(zred_str,'(f6.3)') zred_now

         single_z_file="PDF_d_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=61,file='../subvolume_data/'//single_z_file)
         single_z_file="PDF_xm1_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=62,file='../subvolume_data/'//single_z_file)
         single_z_file="PDF_xm1d_"//trim(adjustl(size_str))//&
              "_z"//trim(adjustl(zred_str))//".dat"
         open(unit=63,file='../subvolume_data/'//single_z_file)

         ii=size
         jj=0
!     loop over the cubical regions of size ii along each direction
         do i=1,mesh(1)-ii+1,sample      
            do j=1,mesh(1)-ii+1,sample 
               do k=1,mesh(1)-ii+1,sample
                  jj=jj+1
                  d(jj)=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1))/&
                       real(ii)**3/avg_dens !normalized mean density per region  
                  x(jj)=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1)& !normalized mean mass ...
                       *xh1(i:i+ii-1,j:j+ii-1,k:k+ii-1))& !neutral fraction per region
                       /sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ii-1))
                  x(jj)=(1.d0-x(jj))/(1.0d0-avg_mHII) !for doing the neutral mass  
               end do
            end do
         end do
         xd=x*d !normalized mean ionized/neutral density per region
               
         arr_size=jj
!     do not count 0's! Make arrays allocatable and allocate only part that is needed!
         maxd=maxval(d(1:arr_size))
         mind=minval(d(1:arr_size))
         mean_d=sum(d(1:arr_size))/real(arr_size)
         sigma_d=sqrt(sum((d(1:arr_size)-mean_d)**2)/real(arr_size-1))

         maxx=maxval(x(1:arr_size))
         minx=minval(x(1:arr_size))
         mean_x=sum(x(1:arr_size))/real(arr_size)
         sigma_x=sqrt(sum((x(1:arr_size)-mean_x)**2)/real(arr_size-1))

         maxxd=maxval(xd(1:arr_size))
         minxd=minval(xd(1:arr_size))
         mean_xd=sum(xd(1:arr_size))/real(arr_size)
         sigma_xd=sqrt(sum((xd(1:arr_size)-mean_xd)**2)/real(arr_size-1))

         
!     bin sizes
         dd  = (maxd-mind)/real(n_x)
         dx  = (maxx-minx)/real(n_x)
         dxd = (maxxd-minxd)/real(n_x)
         print*,'check',maxd,mind,maxxd,minxd,dd,dx,dxd
         print*,'means, sigmas',mean_d, sigma_d,mean_x, sigma_x,&
              mean_xd, sigma_xd
         
         do jj=1,arr_size
!     total number of cells in each density bin
            distribrho(floor((d(jj)-mind)/dd))=distribrho(floor((d(jj)-mind)/dd))+1 
            distribxm(floor((x(jj)-minx)/dx))=& !total number of cells in each ioniz. bin
                 distribxm(floor((x(jj)-minx)/dx))+1. !in log bins 
            distribrhoxm(floor((xd(jj)-minxd)/dxd))=& !total number of cells in each 
                 distribrhoxm(floor((xd(jj)-minxd)/dxd))+1. !ioniz. bin in log bins 
         end do
         
         distribrho=distribrho/sum(distribrho)
         distribxm=distribxm/sum(distribxm)
         distribrhoxm=distribrhoxm/sum(distribrhoxm)

         write(61,*) mean_d, sigma_d,dd
         write(62,*) mean_x, sigma_x,dx         
         write(63,*) mean_xd, sigma_xd,dxd

         do l=0,n_x
            write(61,*)mind+real(l)*dd,real(distribrho(l))/dd
            write(62,*)minx+real(l)*dx,real(distribxm(l))/dx
            write(63,*)minxd+real(l)*dxd,real(distribrhoxm(l))/dxd
         end do
         
         print*,'check normalization',sum(distribrho),sum(distribxm),sum(distribrhoxm)
         close(61)
         close(62)
         close(63)

      enddo
      close(60)

      stop
    end program make_PDF_21cm


