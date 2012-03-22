      program rms21cm

!     calculates the 21-cm mean background evolution and its fluctuations
!     from the simulation results 

      use sizes, only: mesh
      use cgsconstants
      use astroconstants
      use abundances
      use material
      use cosmology_parameters
      use ion_frac_mod

      implicit none

      integer,parameter :: MaxNumZred=200
      integer ngrid1,ngrid2,ngrid3,i,j,k

      character(len=180) :: redshift_file ! name of file with list of redshifts
      character(len=180) :: results_file1, results_file2 ! names of files with results
      real(kind=8) :: rho_crit_0_MpcM_sun,rho_bar,rho_matter

      real(kind=dp) :: zred_now
      integer nz0,nz

      character(len=180) :: xfrac_file, dens_file
      character(len=6) :: zred_str

      !operations done with doubles, otherwise sums of large arrays will be wrong
!      real(kind=dp) :: xh1(mesh(1),mesh(2),mesh(3))
      ! Array needed to read in 4B reals
!      real(kind=si) :: xh1_real(mesh(1),mesh(2),mesh(3))

!      common /x/ xh1

      integer :: nn, ii, ij, n_r, n_l, nfile, count ! loop counters

      real(kind=dp) :: r, l, dtheta, dtheta1, deltanu
      real(kind=dp) :: mean_ion_m, da, nu, H_z, n_bg_z, n_bar, n_bar_H
      real(kind=dp) :: delta_Tb,delta_Tb_mean, fluct

      real*8 :: kappa1
      real(kind=dp),parameter :: A10=2.85d-15 !g;s^-1
      real(kind=dp),parameter :: nu0=1420.40575d6,T_star=0.0681 !Hz;K

      real(kind=dp),parameter :: lambda0=1.0d0-Omega0
      real(kind=dp),parameter :: boxsize = 425. !Mpc/h
      
      real(kind=dp),parameter :: mu_H=1.0d0+4.0d0*abu_he/(1.0d0-abu_he), amass=m_p*mu
      !amass = mu*m_H
      !where mu=1.222 (neutral material), boltzk=k_B

      real(kind=dp) :: hz, omegaz
      real(kind=dp) :: z, z1, avg_mHII, avg_nHII !, avg_dens
      real(kind=dp) :: fluct_dens, mean_dens_i, mean_dens, delta_overdens

      real*8 :: dx

      allocate(ndens(mesh(1),mesh(2),mesh(3)))

!     set the scales at which to calculate the fluctuations
      dtheta=3.                 !arcmin
      deltanu=0.2e6             !Hz
      
      dtheta1=dtheta*2.90888e-04 !in rad
      dx =  real(boxsize/h/mesh(1))  !cell size in comoving Mpc

      n_bar=rho_crit_0*Omega_b/amass !mean gas number density at present
      n_bar_H = n_bar*mu/mu_H !hydrogen number density at present

!     Ask for redshift file
      write(*,'(A,$)') 'File with redshifts: '
      read(*,*) redshift_file
      write(*,'(A,$)') 'Initial snapshot number: '
      read(*,*) nz0
!!$      write(*,'(A,$)') 'Files with results: '
!!$      read(*,*) results_file1
!!$      read(*,*) results_file2  


      call set_id

!      open(unit=2,file=results_file)

!     Open and read redshift file
      open(unit=60,file=redshift_file,form='formatted',status='old')
      read(60,*) NumZred

      if (NumZred > MaxNumZred) write(*,*)&
           'Your program is about to die in a horrible way'

                                !21-cm kappa coefficient
      kappa1 = 3.0d0*c**3*A10*T_star/(32.*pi*nu0**3)
!     
      if(nz0>1)then
         do nz=1,nz0-1
            read(60,*) zred_now
         end do
      end if

      do nz=nz0,NumZred
         read(60,*) zred_now
         print*,'doing z=',zred_now
         write(zred_str,'(f6.3)') zred_now
!     
!     output files
         results_file1="Tb_fluct_z.dat"
         results_file2="Tb_mean_z.dat"
         
         open(unit=61,file=results_file1)
         open(unit=62,file=results_file2)

         call read_density(zred_now)
         
         !     Find summed_and average density
         avg_dens=sum(1.0d0*ndens)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))!particles/cell
         
         !     Report on data: min, max, total
         write(*,*) 'minimum: ',minval(ndens)/8.
         write(*,*) 'maximum: ',maxval(ndens)/8.
         write(*,*) 'summed density: ',sum(1.0d0*ndens)/8.
         write(*,*) 'average dens.',avg_dens/8.
         
         do nn=1,2 !since 2 ionized fraction outputs are done per each density file
            if(nn.eq.2)then
               read(60,*) zred_now
               print*,'doing z=',zred_now
            end if
            
            write(zred_str,'(f6.3)') zred_now
            
            call read_xfrac(zred_now)        
         
            !     Find the average ionized fraction                                      

            avg_nHII=sum(xh1)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))            
            avg_mHII=sum(ndens*xh1)/sum(ndens)
            
            print*,zred_now,avg_nHII,avg_mHII
            write(2,*) zred_now,avg_nHII,avg_mHII
!     calculate mean 21-cm bg
            H_z = 100.*h*sqrt(omega0*(1.+zred_now)**3+lambda0+ &
                 (1.-omega0-lambda0)*(1.+zred_now)**2)/3.086e19 !in s^-1 
            n_bg_z = n_bar_H*(1.+zred_now)**3 !physical number density of H at z
            delta_Tb_mean=kappa1*(1.0d0-avg_mHII)*n_bg_z/& 
                 (H_z*(1.+zred_now))*1e3 !in mK
            write(62,*) zred_now, delta_Tb_mean, avg_mHII            

            ! calculate the 21-cm (delta T_b) and density fluctuations
            
            z1=zred_now+1.
            z=zred_now
            omegaz=omega0*z1**3/((1.-lambda0+omega0*z)*z1**2+lambda0)
            hz=dsqrt(omega0*h**2*z1**3/omegaz)
            nu=nu0/(1.+z)
            l=z1*(c/1e7/hz)*deltanu/nu !comoving size along LOS in Mpc
            call ang_dist(z,da) !find the angular diameter distance
            r=dtheta1*z1*da !comoving size on the sky in Mpc
            
            n_l = anint(l/dx)     !same sizes in cell units
            n_r = anint(r/dx)
            print*,'check ',n_l,n_r, da
            print*,'check 2',l,r,dx
            ii=n_r
            ij=n_l
            count = 0
            fluct=0.0d0
            fluct_dens=0.0d0
            mean_dens = avg_dens*real(ii)**2*real(ij) !particles/volume considered (in cell volumes)

            do i=1,mesh(1)-ii+1,4  !loop over the regions of size ii x ii x ij along each direction
               do j=1,mesh(2)-ii+1,4   
                  do k=1,mesh(3)-ij+1,4
                     mean_dens_i=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ij-1))!mean density in the volume
                     mean_ion_m=sum(ndens(i:i+ii-1,j:j+ii-1,k:k+ij-1) &!mean ionization in the volume
                          *xh1(i:i+ii-1,j:j+ii-1,k:k+ij-1))/mean_dens_i
                     delta_overdens=mean_dens_i/mean_dens  !overdensity of the volume
                     delta_Tb=kappa1*(1.0d0-mean_ion_m)*n_bg_z* &
                          delta_overdens/(H_z*(1.+zred_now))*1e3 !in mK
                     fluct=fluct+(delta_Tb-delta_Tb_mean)**2
                     fluct_dens = fluct_dens+(mean_dens_i-mean_dens)**2
!                     print*,'check Tb',delta_Tb, delta_Tb_mean
                     count = count+1
                  end do
               end do
            end do
            fluct = sqrt(fluct/count)
            fluct_dens=sqrt(fluct_dens/count)/mean_dens
            print*,'result',z,fluct,delta_Tb_mean,fluct_dens
            write(61,*) z,fluct,fluct_dens
         end do
      end do
!      stop
      
      close(60)
      close(61)
      close(62)
         
      stop
      end
      
      subroutine ang_dist(z,dist)
                                ! calculates the angular diameter distance 
                                ! out to a given redshift z for LambdaCDM

!      use sizes, only: mesh
      use cgsconstants
      use astroconstants
!      use abundances
      use material
      use cosmology_parameters


      implicit none
      
      external integrant2
      real*8 :: z, H00
      real*8 :: area, integrant, dist
      integer :: i
      
      H00=100.                  !h^-1 km/s/Mpc
      call INTEG2(0.0d0,z,integrant2,area)
      dist = c/1e5/H00/(1.d0+z)*area !h^-1 Mpc, LambdaCDM

      return
      end subroutine ang_dist
!=========================================================================
      subroutine integ2(a,b,F,area)
                                !
                                ! This subroutine computes integrals using 
                                !  Simpson's rule.
                                !

      implicit none

      real*8 :: tol, a, b, c, area, areaold, sumend, sumeven, sumodd, h, F
      integer :: i, j, n


      data tol /1.e-08/

      if(a.eq.b) then
         area=0.
         return
      endif
      
      areaold=0.
      sumend=F(a)+F(b)
      sumeven=0.
      sumodd=0.
      
      do n=1,25
         i=2.**n
         h=(b-a)/dfloat(i)
         sumodd=0.
         do j=1,i-1,2
            c=a+j*h
            sumodd=sumodd+F(c)
         enddo
         area=(sumend+4.*sumodd+2.*sumeven)*h/3.
         if(dabs((area-areaold)/area).lt.tol) return
         sumeven=sumeven+sumodd
         areaold=area
      enddo
      
      write(6,1000)
 1000 format(/5x,'Error, no convergence.')
      stop
      
      end subroutine integ2
!=========================================================================
      function integrant2(x)
 
      use cosmology_parameters

      implicit none

      real*8 :: x, integrant2
            
      integrant2=Omega0**(-0.5)/((1.0d0+x)**3+1.0d0/Omega0-1.0d0)**.5
      
      return
      end function integrant2
      
