      program fractions_calc

!     calculates the volume- and mass-neutral fractions at a given z
!     from the simulation results 

      use sizes, only: mesh
!      use mathconstants
!      use astroconstants
!      use cgsconstants
!      use cosmology
!      use lg
      use material
      use ion_frac_mod

      implicit none

      integer,parameter :: MaxNumZred=200
      integer ngrid1,ngrid2,ngrid3,i,j,k

      character(len=180) :: redshift_file ! name of file with list of redshifts
      character(len=180) :: results_file ! name of file with results
      real(kind=8) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter
      integer :: nz ! loop counter
!      integer ::  NumZred

      ! density in file is in real(kind=4), read in via this array
!      real ndens_real(mesh(1),mesh(2),mesh(3))
!     real*8 ndens(mesh(1),mesh(2),mesh(3))

      real(kind=dp) :: zred_now
      integer nz0

      character(len=180) :: xfrac_file, dens_file
      character(len=6) :: zred_str


      !operations done with doubles, otherwise sums of large arrays will be wrong
!      real(kind=dp) :: xh1(mesh(1),mesh(2),mesh(3))
      ! Array needed to read in 4B reals
!      real(kind=si),dimension(:,:,:),allocatable :: xh1_real
!      real(kind=si) :: xh1_real(mesh(1),mesh(2),mesh(3))

      real*8 :: avg_mHII, avg_nHII, sum_xh1
      
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
      write(*,'(A,$)') 'File with results: '
      read(*,*) results_file

      open(unit=2,file=results_file)
     
!     Open and read redshift file
      open(unit=60,file=redshift_file,form='formatted',status='old')
      read(60,*) NumZred

!      allocate(snap(NumZred))
!      allocate(zred_array(NumZred))

      if (NumZred > MaxNumZred) write(*,*)&
          'Your program is about to die in a horrible way'

      do nz=1,nz0-1     
!          read(60,*) zred_array(nz)
         read(60,*) zred_now
      end do

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
         write(2,*) zred_now,avg_nHII,avg_mHII

         read(60,*) zred_now
         
         call read_xfrac(zred_now)

!     Find the average ionized fraction
         avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
         
         avg_mHII=sum(ndens*xh1)/sum(ndens)

         print*,zred_now,avg_nHII,avg_mHII
         write(2,*) zred_now,avg_nHII,avg_mHII

      enddo
      close(60)

      stop
    end program fractions_calc


