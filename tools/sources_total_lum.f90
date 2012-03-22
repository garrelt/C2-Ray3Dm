      program sources_calc

!     calculates the number of sources (active and total)
!     from the simulation results 

      use sizes, only: mesh
      use sourceprops
      use source_statistics
!      use mathconstants
!      use astroconstants
!      use cgsconstants
!      use cosmology
!      use lg
      use material
      use c2ray_parameters

      implicit none

      integer,parameter :: MaxNumZred=200
      integer ngrid1,ngrid2,ngrid3,i,j,k

      character(len=180) :: redshift_file, redshift_file2 ! names of files with list of redshifts
      character(len=180) :: results_file ! name of file with results
      real(kind=8) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter
!      integer :: nz ! loop counter
      integer ::  NumZred0

      ! density in file is in real(kind=4), read in via this array
!      real ndens_real(mesh(1),mesh(2),mesh(3))
!     real*8 ndens(mesh(1),mesh(2),mesh(3))

      real(kind=dp) :: zred_now
      integer nz0, nz

      character(len=180) :: xfrac_file, dens_file
      character(len=6) :: zred_str

      print*,'resolution:', mesh     
      call set_id

!     Ask for redshift file
      write(*,'(A,$)') 'File with used redshifts: '
      read(*,*) redshift_file
      write(*,'(A,$)') 'File with all redshifts: '
      read(*,*) redshift_file2
      write(*,'(A,$)') 'Initial snapshot number: '
      read(*,*) nz0
      write(*,'(A,$)') 'File with results: '
      read(*,*) results_file

      open(unit=2,file=results_file)
     
!     Open and read redshift file
      open(unit=60,file=redshift_file,form='formatted',status='old')
      read(60,*) NumZred0

      open(unit=61,file=redshift_file2,form='formatted',status='old')
      read(61,*) NumZred

      print*,NumZred,' redshifts to read'


      if (NumZred > MaxNumZred) write(*,*)&
          'Your program is about to die in a horrible way'

      do nz=nz0,NumZred

         read(61,*) zred_now

         print*,'doing z=',zred_now
         !
         !sources with no suppression (i.e. total)
         call read_sources(zred_now,0) 
         
         !sources with suppression (i.e. active ones)
         if(nz<=NumZred0) then
            call read_sources(zred_now,1) 
         else
            TotSrcFlux=LargeSrcFlux0 !after overlap set total flux to the total coming from the large sources
         end if

         !report
         write(*,"(f10.3,2x,2i8,2x,4es12.3)") zred_now,NumSrc00,NumSrc01,TotSrcFlux0,TotSrcFlux,LargeSrcFlux0,SmallSrcFlux0
         write(2,"(f10.3,2x,2i8,2x,4es12.3)") zred_now,NumSrc00,NumSrc01,TotSrcFlux0,TotSrcFlux,LargeSrcFlux0,SmallSrcFlux0

      enddo
      close(60)

      stop
    end program sources_calc


