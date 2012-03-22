      program tau
                                !calculates the Thomson optical depth tau 
                                !based on the ionized fraction data
      implicit none

      integer :: i,n
      real*8,parameter :: omega0=0.27,lambda0=0.73, omegab=0.044, h=0.7
      real*8,parameter :: rho_crit_0=1.88e-29*h**2, m_p=1.67e-24
      real*8,parameter :: mu_H=1.32, sigma_T=6.65e-25, c=3e10
      real*8,parameter :: AHe=0.08, chi1=1.+AHe, chi2=1.+2.*AHe
      real*8,dimension(200) :: frac_m, frac_v, z, tau_es, tau_es_full
      real*8 :: tau0, coeff, zred

      open(1,file='frac_425Mpc_f2_10S_504.dat')
      open(2,file='tau_425Mpc_f2_10S_504.dat')

      tau_es=0.0d0

      do i=1,200
         read(1,*,end=10) z(i),frac_v(i),frac_m(i)  
      end do
 10   n=i-1
      print*,n,z(n)

!      pause
 
      coeff = 2.*c*sigma_T*Omegab*rho_crit_0/
     &     (3.*100.*h/3.086e19*m_p*Omega0)*chi1/mu_H
      do i=0,50
         zred=real(i)/50.*z(n)
         tau0 = coeff*(sqrt(omega0*(1.+zred)**3+lambda0) - 1.) 
         write(2,*) zred+1.,tau0,tau0
!         print*,zred.,tau0,tau0
      end do
      
      tau_es(n)=tau0

      print*,'start',z(n),tau_es(n)
      write(2,*) z(n)+1.,tau_es(n),tau_es(n)
      do i=n-1,1,-1
         tau_es(i) =tau_es(i+1)+1.5*coeff*Omega0*
     &        (frac_m(i)*(1.+z(i))**2/sqrt(omega0*(1.+z(i))**3+lambda0)
     &  + frac_m(i+1)*(1.+z(i+1))**2/sqrt(omega0*(1.+z(i+1))**3+lambda0)
     &        )*(-z(i+1)+z(i))/2.
         tau_es_full(i)= coeff*(sqrt(omega0*(1.+z(i))**3+lambda0) - 1.)
!         print*,z(i),tau_es(i)+0.011
!         write(2,*) z(i)+1.,tau_es(i)+0.011,tau_es_full(i)+0.011
         print*,z(i),tau_es(i)
         write(2,*) z(i)+1.,tau_es(i),tau_es_full(i)
      end do
      
      close(1)
      close(2)
      stop
      end
