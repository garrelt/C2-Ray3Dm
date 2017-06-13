program morph_dens

  use precision, only: dp, si
  use cosmosimple, only: zred2time, time2zred
  use vtk_2d, only: vtk2d

  implicit none
  integer, parameter             :: N_inter = 5
  real(kind=dp), dimension(N_inter) :: timeinter 
  integer,dimension(3),parameter :: mesh=(/250,250,250/)
  integer                        :: N_zred, m1, m2, m3, ii, i, j, iiii
  real(kind=dp)                  :: z_red, z_red_prev, z_red_inter, timenow, timeprev
  real(kind=dp)                  :: ndens_prev_av, ndens_av, frac_av

  character(len=512),parameter :: dir_dens="../../coarser_densities/halos_included/"
  character(len=512)           :: dens_prev_file, dens_file, dens_fine_file
  real(kind=si),dimension(mesh(1),mesh(2),mesh(3)) :: ndens_prev, ndens, frac, ndens_inter

  character(len=6)   :: z_prev_str, z_str, z_inter_str
  character(len=3)   :: ii_str
  
  open(unit=51, file="redshifts.dat", status="old")
  open(unit=52, file="redshifts_fine.dat", status="unknown")
  read (51,*) N_zred
  write(52,*) (N_zred-1)*(N_inter+1)

  read(51,*) z_red_prev
  iiii = 0
  do ii=2, N_zred
     write(52,'(f6.3)') z_red_prev
     read(51,*) z_red
     timeprev = zred2time(z_red_prev)
     timenow  = zred2time(z_red)
     print *, 'z_red_prev: ', z_red_prev, ' timeprev: ', timeprev
     print *, 'z_red:      ', z_red,      ' timenow:  ', timenow
     write(z_prev_str,'(f6.3)') z_red_prev
     write(z_str,     '(f6.3)') z_red
     


     ! do whatever analysis
     dens_prev_file=trim(adjustl(dir_dens))//trim(adjustl(z_prev_str))//"ntot_all.dat"
     dens_file     =trim(adjustl(dir_dens))//trim(adjustl(z_str))//"ntot_all.dat"
     
     open(unit=20,file=dens_prev_file,form="binary",status="old")
     read(20) m1,m2,m3
     print *, 'mesh number', m1, m2, m3
     if (m1 /= mesh(1)) then
        print *, 'mesh number not matched'
        close(20)
        goto 333
     endif
     read(20) ndens_prev
     close(20)

     ! density file spaced in finer time sequence: for z_prev
     dens_fine_file = trim(adjustl(z_prev_str))//"n_all.dat"
     open(unit=20,file=dens_fine_file,form="binary",status="unknown")
     write(20) m1,m2,m3
     write(20) ndens_prev
     close(20)

     ndens_prev_av = sum(dble(ndens_prev))/real(m1*m2*m3,dp)

     open(unit=20,file=dens_file,form="binary",status="old")
     read(20) m1,m2,m3
     print *, 'mesh number', m1, m2, m3
     if (m1 /= mesh(1)) then
        print *, 'mesh number not matched'
        close(20)
        goto 333
     endif
     read(20) ndens
     close(20)

     ndens_av = sum(dble(ndens))/real(m1*m2*m3,dp)

     frac     = ndens / ndens_prev
     frac_av  = sum(dble(frac))/real(m1*m2*m3,dp)

!!$     iiii = iiii + 1
!!$     if (100 <= iiii .and. iiii < 1000) then
!!$        write(ii_str,'(I3)') iiii
!!$     elseif (10 <= iiii .and. iiii < 100) then
!!$        write(ii_str,'(I2)') iiii
!!$        ii_str = "0"//trim(adjustl(ii_str))
!!$     else
!!$        write(ii_str,'(I1)') iiii
!!$        ii_str = "00"//trim(adjustl(ii_str))
!!$     endif
!!$     
!!$     open(unit=1,file="dens2d"//trim(adjustl(ii_str))//"-"//trim(adjustl(z_prev_str))//".vtk",status="unknown")
!!$     call vtk2d(real(ndens_prev(:,:,m3/2),dp), m1, m2, 1, "dens2d", 6)
!!$     close(1)

     ! Now morph for intermediate steps
     ! set intermediate times
     do i=1, N_inter
        timeinter(i)= timeprev+(timenow-timeprev)*real(i,dp)/real(N_inter+1,dp)
        ndens_inter =ndens_prev+(ndens-ndens_prev)*real(i,dp)/real(N_inter+1,dp)
        z_red_inter = time2zred(timeinter(i))
        write(52,'(f6.3)') z_red_inter
        write(z_inter_str, '(f6.3)') z_red_inter

        ! density file spaced in finer time sequence: for z_prev
        dens_fine_file = trim(adjustl(z_inter_str))//"n_all.dat"
        open(unit=20,file=dens_fine_file,form="binary",status="unknown")
        write(20) m1,m2,m3
        write(20) ndens_inter
        close(20)
!!$        iiii = iiii + 1
!!$        if (100 <= iiii .and. iiii < 1000) then
!!$           write(ii_str,'(I3)') iiii
!!$        elseif (10 <= iiii .and. iiii < 100) then
!!$           write(ii_str,'(I2)') iiii
!!$           ii_str = "0"//trim(adjustl(ii_str))
!!$        else
!!$           write(ii_str,'(I1)') iiii
!!$           ii_str = "00"//trim(adjustl(ii_str))
!!$        endif
!!$
!!$        open(unit=1,file="dens2d"//trim(adjustl(ii_str))//"-"//trim(adjustl(z_inter_str))//".vtk",status="unknown")
!!$        call vtk2d(real(ndens_inter(:,:,m3/2),dp), m1, m2, 1, "dens2d", 6)
!!$        close(1)
        

     enddo

     print *, 'maxfrac:  ', maxval(frac)
     print *, 'minfrac:  ', minval(frac) 
     print *, 'av frac:  ', frac_av
     print *, 'rms frac: ', sqrt(sum(dble((frac-frac_av)**2))/real(m1*m2*m3,dp))
     print *, 'growth fraction: ', (1d0+z_red_prev)/(1d0+z_red)


!!$     if (100 <= iiii .and. iiii < 1000) then
!!$        write(ii_str,'(I3)') iiii
!!$     elseif (10 <= iiii .and. iiii < 100) then
!!$        write(ii_str,'(I2)') iiii
!!$        ii_str = "0"//trim(adjustl(ii_str))
!!$     else
!!$        write(ii_str,'(I1)') iiii
!!$        ii_str = "00"//trim(adjustl(ii_str))
!!$     endif
!!$     
!!$     print *, 'iiii: ', iiii
!!$     print *, 'ii_str: ', ii_str, ' z_str: ', z_str
!!$     open(unit=1,file="densfrac"//trim(adjustl(ii_str))//"-"//trim(adjustl(z_str))//".vtk",status="unknown")
!!$
!!$     call vtk2d(real(frac(:,:,m3/2),dp), m1, m2, 1, "densfrac", 8)
!!$     
!!$     close(1)
     
     z_red_prev = z_red

  enddo
  close(51)
  close(52)

!!$  iiii = iiii + 1
!!$  if (100 <= iiii .and. iiii < 1000) then
!!$     write(ii_str,'(I3)') iiii
!!$  elseif (10 <= iiii .and. iiii < 100) then
!!$     write(ii_str,'(I2)') iiii
!!$     ii_str = "0"//trim(adjustl(ii_str))
!!$  else
!!$     write(ii_str,'(I1)') iiii
!!$     ii_str = "00"//trim(adjustl(ii_str))
!!$  endif
!!$  
!!$  open(unit=1,file="dens2d"//trim(adjustl(ii_str))//"-"//trim(adjustl(z_str))//".vtk",status="unknown")
!!$  call vtk2d(real(ndens(:,:,m3/2),dp), m1, m2, 1, "dens2d", 6)
!!$  close(1)
  

333 continue
end program morph_dens

