      program contours

      implicit none
c     calculates the Gamma-delta correlation as a 2D image for making contour plots 
      
      integer,parameter :: size=100
      character(len=180) :: redshift_file ! name of file with list of redshifts 
      character(len=180) :: file1, file2, file3
      character(len=6) :: z_str
      character(len=1) :: nfile_str 

      real Gamma(2097152)
      real delta(2097152),ld(2097152)

      integer numzred,nz,i,j,k,xc,yc
      real zred_now, delmin,delmax, Gammamin, Gammamax
      real image(size,size)

c     Ask for redshift file
c      write(*,'(A,$)') 'File with redshifts: '
c      read(*,*) redshift_file
c      redshift_file='redshifts.dat'
c      redshift_file='red.dat'
      redshift_file='red_short.dat'
     
c     Open and read redshift file
      open(unit=60,file=redshift_file,form='formatted',status='old')
      read(60,*) NumZred
      do nz=1,NumZred
         read(60,*) zred_now
         print*,'doing z=',zred_now
c     
         write(z_str,'(f6.3)') zred_now
         
         image=0.0d0
         file1=
     &        'ionrate_dens_correl_114Mpc_f10_150S_'//
     &    trim(adjustl(z_str))//'.dat'
         open(unit=54,file=file1,
     &        status='old')
         file2=
     &        'ionrate_dens_correl_114Mpc_f10_150S_log_'//
     &        trim(adjustl(z_str))//'.dat'
         open(unit=55,file=file2,
     &        status='unknown')
         do i=1,2097152
            read(54,*) delta(i), Gamma(i)
            ld(i)=log10(delta(i)+1.)
            Gamma(i)=max(Gamma(i),-21.)
            write(55,*) ld(i), Gamma(i)
         end do
         
c         delmin=minval(log10(delta+1))
c         delmax=maxval(log10(delta+1))
c         Gammamin=minval(Gamma)
c         Gammamax=maxval(Gamma)
         delmin=-1.0
         delmax=2.0
         Gammamin=-21
         Gammamax=-8
         print*,'min,max',delmin,delmax,Gammamin,Gammamax
         print*,'actual min,max',minval(ld),maxval(ld),
     &        minval(Gamma),maxval(Gamma)
c         pause
         do i=1,2097152
            xc=int((ld(i)-delmin)/(delmax-delmin)*(size-1))+1
            yc=int((Gamma(i)-Gammamin)/(Gammamax-Gammamin)*(size-1))+1
            if(xc>size.or.yc>size.or.xc<=0.or.yc<=0)then
               print*,xc,yc,i
               xc=min(xc,size)
               yc=min(yc,size)
               xc=max(xc,1)
               yc=max(yc,1)
            end if
            image(xc,yc)=image(xc,yc)+1
!            print*,'check',xc,yc,Gamma(i),ld(i),i
!            pause
         end do
         
         close(54)
!         close(55)
         
         file1='imagec_'//trim(adjustl(z_str))//'.bin'         
         open(57,file=file1,FORM='UNFORMATTED')
         write(57) size,size
         write(57) image
         close(57)

      end do

      stop
      end
